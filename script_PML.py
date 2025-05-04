import os
import re
from Bio.PDB import MMCIFParser
def parse_chothia_regions(summary_file):
    regions = {}
    current = None
    chain = None
    with open(summary_file, encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line.startswith('===') and line.endswith('==='):
                fn = line.strip('= ').split()[0]
                current = os.path.splitext(fn)[0]
                regions[current] = {'L': {}, 'H': {}}
            elif line.startswith('Chain '):
                chain = line.split()[1].rstrip(':')
            elif current and chain and ':' in line:
                name, rng = line.split(':', 1)
                rng = rng.strip().replace('–', '-')
                start, end = map(int, rng.split('-'))
                regions[current][chain][name.strip()] = (start, end)
    return regions

def parse_mutation_code(folder_name):
    parts = folder_name.split('_')
    if len(parts) < 4:
        raise ValueError(f"Неподдерживаемое имя папки: {folder_name}")
    mut = parts[3]  
    m = re.match(r"[A-Za-z](?P<chain>[A-Za-z])(?P<pos>\d+)[A-Za-z]", mut)
    if not m:
        raise ValueError(f"Не удалось разобрать код мутации '{mut}'")
    chain = m.group('chain').upper()
    pos = int(m.group('pos'))
    return chain, pos

def find_first_cif(root):
    """
    Ищет рекурсивно первый .cif в поддиректориях root и возвращает путь.
    """
    for dirpath, dirnames, filenames in os.walk(root):
        for fn in filenames:
            if fn.lower().endswith('.cif'):
                return os.path.join(dirpath, fn)
    return None

def parse_complex_name(name: str):
    """
    Извлекает два идентификатора цепей из имени комплекса вида XXX_AB, XXX_CD или XXX_HL.
    """
    m = re.search(r'(?<=_)([A-Z]{2})(?=_)', name) or \
        re.search(r'(?<=_)([A-Z]{2})$', name)
    if not m:
        raise ValueError(f"Не удалось из имени комплекса '{name}' извлечь пару идентификаторов цепей")
    return list(m.group(1))

def generate_pml(wt_cif, mt_cif, regions_for_id, struct_folder, out_dir, summary_file):
    pdb_id = struct_folder.split('_')[0].upper()
    chain_mut, pos_mut = parse_mutation_code(struct_folder)
    print(chain_mut, pos_mut)
    chain_pair = parse_complex_name(struct_folder)
    print(chain_pair)
    if pdb_id == '3NPS':
        chain_map = {'C': 'L', 'B': 'H'}
    else:
        if chain_mut in ['H', 'L']:
            chain_map =  {'H': 'H', 'L': 'L'}
        else:
            a, b = sorted(chain_pair)
            chain_map = {a: 'L', b: 'H'}

    # print(list(regions_for_id.keys()))
    # print(chain_map)
    chain_H = chain_map.get(chain_mut)
    # print(chain_H, "khkhchain_H")
    if chain_H is None:
        raise KeyError(f"Цепь мутации {chain_mut} не найдена в разметке для {pdb_id}")

    mut_region = None
    for region, (start, end) in regions_for_id[chain_H].items():
        if start <= pos_mut <= end:
            mut_region = (region, start, end)
            break
    if mut_region is None:
        raise ValueError(f"Позиция {pos_mut} не входит ни в один регион цепи {chain_H}")
    region_name, region_start, region_end = mut_region

    ranges = {}
    for ch, regs in regions_for_id.items():
        starts = [s for s, e in regs.values()]
        ends = [e for s, e in regs.values()]
        ranges[ch] = (min(starts), max(ends))
    start_H, end_H = ranges['H']
    start_L, end_L = ranges['L']

    wt_label = 'wt'
    mt_label = 'mt'
    pml_lines = []
    pml_lines.append(f'load "{wt_cif}", {wt_label}')
    pml_lines.append(f'load "{mt_cif}", {mt_label}')

    select_expr = (
        f"(({wt_label} and name CA and chain H and resi {start_H}-{end_H}) or "
        f"({wt_label} and name CA and chain L and resi {start_L}-{end_L}))"
    )
    pml_lines.append(f'select wt_sel, {select_expr}')

    select_expr_mt = (
        f"(({mt_label} and name CA and chain H and resi {start_H}-{end_H}) or "
        f"({mt_label} and name CA and chain L and resi {start_L}-{end_L}))"
    )
    pml_lines.append(f'select mt_sel, {select_expr_mt}')


    pml_lines.append('align mt_sel, wt_sel')

    pml_lines.append('hide everything')
    pml_lines.append(f'show cartoon, (wt_sel)')
    pml_lines.append(f'show cartoon, (mt_sel)')
    pml_lines.append('')

    pml_lines.append('color red, wt_sel and chain H')
    pml_lines.append('color orange, wt_sel and chain L')
    pml_lines.append('color blue, mt_sel and chain H')
    pml_lines.append('color cyan, mt_sel and chain L')
    pml_lines.append('')

    pml_lines.append(f'select wt_H, {wt_label} and name CA and chain H and resi {start_H}-{end_H}')
    pml_lines.append(f'select wt_L, {wt_label} and name CA and chain L and resi {start_L}-{end_L}')
    pml_lines.append(f'select mt_H, {mt_label} and name CA and chain H and resi {start_H}-{end_H}')
    pml_lines.append(f'select mt_L, {mt_label} and name CA and chain L and resi {start_L}-{end_L}')
    pml_lines.append('align mt_H, wt_H')
    pml_lines.append('align mt_L, wt_L')

    pml_lines.append(f'show sticks, {wt_label} and name CA and chain {chain_H} and resi {pos_mut}')
    pml_lines.append(f'show sticks, {mt_label} and name CA and chain {chain_H} and resi {pos_mut}')
    pml_lines.append(f'color yellow, {wt_label} and name CA and chain {chain_H} and resi {pos_mut}')
    pml_lines.append(f'color green, {mt_label} and name CA and chain {chain_H} and resi {pos_mut}')
    pml_lines.append('')

    pml_lines.append(f'select {region_name}_wt, {wt_label} and name CA and chain {chain_H} and resi {region_start}-{region_end}')
    pml_lines.append(f'select {region_name}_mt, {mt_label} and name CA and chain {chain_H} and resi {region_start}-{region_end}')
    pml_lines.append(f'align {region_name}_wt, {region_name}_mt')

    os.makedirs(out_dir, exist_ok=True)
    pml_path = os.path.join(out_dir, f"{struct_folder}.pml")
    with open(pml_path, 'w', encoding='utf-8') as f:
        f.write("\n".join(pml_lines))
    return pml_path

if __name__ == '__main__':
    regions_file = r"C:\Users\Redmi\Desktop\NIR\mutation\chotia_clean\regions_summary.txt"
    summary = parse_chothia_regions(regions_file)
    base_wt = r"C:\Users\Redmi\Desktop\NIR\mutation\str_WT\struct_WT"
    types_base = r"C:\Users\Redmi\Desktop\NIR\mutation\type_amin"
    out_dir = r"C:\Users\Redmi\Desktop\NIR\mutation\file_PML"
    for aa_type in os.listdir(types_base):
        aa_dir = os.path.join(types_base, aa_type)
        if not os.path.isdir(aa_dir):
            continue
        for struct_folder in os.listdir(aa_dir):
            struct_dir = os.path.join(aa_dir, struct_folder)
            if not os.path.isdir(struct_dir):
                continue
            pdb_id = struct_folder.split('_')[0].upper()
            wt_cif = None
            for wt_sub in os.listdir(base_wt):
                if wt_sub.upper().startswith(pdb_id + "_") and wt_sub.upper().endswith("_WT"):
                    wt_cif = find_first_cif(os.path.join(base_wt, wt_sub))
                    break
            mt_cif = None
            for sub in os.listdir(struct_dir):
                subp = os.path.join(struct_dir, sub)
                if os.path.isdir(subp) and "_MT" in sub.upper():
                    mt_cif = find_first_cif(subp)
                    break
            if wt_cif and mt_cif and pdb_id in summary:
                # print(wt_cif, mt_cif, summary[pdb_id], struct_folder, out_dir, regions_file)
                pml = generate_pml(wt_cif, mt_cif, summary[pdb_id], struct_folder, out_dir, regions_file)
                print(f"Generated {pml}")
