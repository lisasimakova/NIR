import os
import re

def parse_mutation_code(name):
    parts = name.split('_')
    if len(parts) < 4:
        raise ValueError(f"Не удалось разобрать имя {name}")
    chain_mut = parts[2]
    digits = re.findall(r'\d+', parts[-1])
    if not digits:
        raise ValueError(f"Не найдены цифры в коде мутации {parts[-1]} из {name}")
    pos_mut = int(digits[0])
    return chain_mut, pos_mut
def get_mutation_region_from_block(block, pos_mut, mutated_chain):
    pattern = fr'---\s*Chain\s+{mutated_chain}\s+---(.*?)(?:\n---|\n===|\Z)'
    m = re.search(pattern, block, flags=re.DOTALL)
    if not m:
        return "Chain block not found", ""
    chain_block = m.group(1)
    region_pattern = r'\s*(\S+):\s+RMSD\s*=\s*([0-9.]+)\s*Å\s*\(\s*(\d+)\s*-\s*(\d+)\s*residues\)'
    for line in chain_block.splitlines():
        rm = re.match(region_pattern, line)
        if rm:
            start, end = int(rm.group(3)), int(rm.group(4))
            if start <= pos_mut <= end:
                return rm.group(1), rm.group(2)
    return "Region not found", ""

def parse_logs(input_log):
    records = []
    blocks = re.split(r'\n\*{3}\s*', input_log)
    for block in blocks:
        if not block.strip():
            continue
        header = block.splitlines()[0].strip()
        hp = r'([^/]+)/([^*]+)\*{3}.*cur_h\s*=\s*(\d+),\s*cur_l\s*=\s*(\d+)'
        m = re.match(hp, header)
        if not m:
            print("Не удалось распарсить заголовок:", header)
            continue
        acid_type = m.group(1).strip()
        name_file = os.path.splitext(m.group(2).strip())[0]
        mutated_chain = 'H' if m.group(3)=="1" else ('L' if m.group(4)=="1" else None)
        total = re.search(r'===\s*Overall RMSD.*?:\s*([0-9.]+)', block)
        heavy = re.search(r'>>\s*Chain H global RMSD:\s*([0-9.]+)', block)
        light = re.search(r'>>\s*Chain L global RMSD:\s*([0-9.]+)', block)
        try:
            _, pos_mut = parse_mutation_code(name_file)
        except Exception as e:
            pos_mut = None
            print(f"Ошибка при разборе имени {name_file}: {e}")
        if pos_mut and mutated_chain:
            region, region_rmsd = get_mutation_region_from_block(block, pos_mut, mutated_chain)
        else:
            region, region_rmsd = "Undefined", ""
        records.append({
            'acid_type': acid_type,
            'name_file': name_file,
            'total_RMSD_python': total.group(1) if total else None,
            'heavy_RMSD_python': heavy.group(1) if heavy else None,
            'light_RMSD_python': light.group(1) if light else None,
            'REGION_python': region,
            'region_RMSD_python': region_rmsd,

            'total_RMSD_pml': None,
            'heavy_RMSD_pml': None,
            'light_RMSD_pml': None,
            'REGION_pml': None,

            'DDG': None
        })
    return records

def load_pml_data(pml_filepath):
    pml_data = {}
    with open(pml_filepath, "r", encoding="utf-8") as fin:
        fin.readline()  
        for line in fin:
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            key = parts[0]
            pml_data[key] = {
                'total_RMSD_pml': parts[1],
                'heavy_RMSD_pml': parts[2],
                'light_RMSD_pml': parts[3],
                'REGION_pml': parts[4]
            }
    return pml_data

def load_ddg_data(ddg_filepath):
    """
    Читает файл вида:
      name_file<TAB>DDG (с десятичной запятой)
    Возвращает dict: { name_file: float(DDG), … }
    """
    ddg_data = {}
    # print(ddg_filepath)
    with open(ddg_filepath, "r", encoding="utf-8") as fin:
        # print(fin)
        for line in fin:
            # print(line)
            parts = line.strip().split()
            # print(parts)
            if len(parts) != 2:
                continue
            key, val = parts
            key = key[:-3]

            try:
                ddg_data[key] = float(val.replace(',', '.'))
            except ValueError:
                ddg_data[key] = None
    # print(ddg_data)
    return ddg_data

def save_records_by_acid(records, out_dir):
    from collections import defaultdict
    groups = defaultdict(list)
    for rec in records:
        groups[rec['acid_type']].append(rec)
    os.makedirs(out_dir, exist_ok=True)
    header = (
        "name_file\ttotal_RMSD_python\theavy_RMSD_python\t"
        "light_RMSD_python\tREGION_python\tregion_RMSD_python\t"
        "total_RMSD_pml\theavy_RMSD_pml\tlight_RMSD_pml\tREGION_pml\tDDG"
    )
    for acid, recs in groups.items():
        chunks = [recs[i:i+10] for i in range(0, len(recs), 10)]
        for idx, chunk in enumerate(chunks, start=1):
            suffix = f"_{idx}" if len(chunks) > 1 else ""
            filename = f"{acid}{suffix}.tsv"
            path = os.path.join(out_dir, filename)
            with open(path, "w", encoding="utf-8") as fout:
                fout.write(header + "\n")
                for rec in chunk:
                    row = [
                        rec.get('name_file', ""),
                        rec.get('total_RMSD_python') or "",
                        rec.get('heavy_RMSD_python') or "",
                        rec.get('light_RMSD_python') or "",
                        rec.get('REGION_python') or "",
                        rec.get('region_RMSD_python') or "",
                        rec.get('total_RMSD_pml') or "",
                        rec.get('heavy_RMSD_pml') or "",
                        rec.get('light_RMSD_pml') or "",
                        rec.get('REGION_pml') or "",

                        f"{rec['DDG']:.2f}".replace('.', ',') if rec['DDG'] is not None else ""
                    ]
                    fout.write("\t".join(row) + "\n")
            print(f"Сохранено: {path}")

if __name__ == '__main__':
    # Пути
    input_log   = r"C:\Users\Redmi\Desktop\NIR\mutation\type_amin\all_rmsd.log"
    pml_filepath= r"C:\Users\Redmi\Desktop\NIR\mutation\tsv_tables\parsed_rmsd.tsv"
    ddg_filepath= r"C:\Users\Redmi\Desktop\NIR\mutation\DDG.txt"
    out_dir     = r"C:\Users\Redmi\Desktop\NIR\mutation\tsv_tables"

    with open(input_log, "r", encoding="utf-8") as fin:
        log_data = fin.read()
    records = parse_logs(log_data)

    pml_data = load_pml_data(pml_filepath)

    ddg_data = load_ddg_data(ddg_filepath)

    for rec in records:
        key = rec['name_file']

        rec.update(pml_data.get(key, {
            'total_RMSD_pml': "",
            'heavy_RMSD_pml': "",
            'light_RMSD_pml': "",
            'REGION_pml': ""
        }))

        rec['DDG'] = ddg_data.get(key)

    save_records_by_acid(records, out_dir)
