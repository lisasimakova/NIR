import os
import shutil
import json
import re
from typing import List, Dict


def parse_complex_name(name: str) -> List[str]:
    m = re.search(r'(?<=_)([A-Z]{2})(?=_)', name) or \
        re.search(r'(?<=_)([A-Z]{2})$', name)
    if not m:
        raise ValueError(f"Не удалось извлечь пару идентификаторов из имени комплекса '{name}'")
    return list(m.group(1))


def assign_H_L(chain_ids: List[str]) -> Dict[str, str]:
    if chain_ids == '3NPS_A_BC':
        return {'C': 'L', 'B': 'H'}
    if set(chain_ids) == {'H', 'L'}:
        return {'H': 'H', 'L': 'L'}
    a, b = sorted(chain_ids)
    return {a: 'L', b: 'H'}


def apply_mutation(json_path: str, complex_name: str, mutation: str) -> None:
    # парсим mutation
    m = re.match(r'^([A-Z])([A-Z])(\d+)([A-Z])$', mutation)
    if not m:
        raise ValueError(f"Неверный формат мутации: {mutation}")
    wt_aa, orig_chain, pos_str, new_aa = m.groups()
    pos = int(pos_str)

    with open(json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    new_name = os.path.splitext(os.path.basename(json_path))[0]
    data['name'] = new_name

    chain_pair = parse_complex_name(complex_name)
    if complex_name == '3NPS_A_BC':
        mapping = {'C': 'L', 'B': 'H'}
    else:
        mapping = assign_H_L(chain_pair)
    role = mapping.get(orig_chain)
    if not role:
        raise KeyError(f"Chain '{orig_chain}' не найден в комплексe {complex_name}")

    # ищем запись с нужным role
    seq_entry = None
    for entry in data.get('sequences', []):
        if entry.get('protein', {}).get('id') == role:
            seq_entry = entry
            break
    if not seq_entry:
        raise KeyError(f"В JSON нет цепи роли '{role}' для {complex_name}")

    seq = seq_entry['protein']['sequence']
    if pos < 1 or pos > len(seq):
        raise IndexError(f"Позиция {pos} выходит за длину последовательности ({len(seq)})")
    if seq[pos-1] != wt_aa:
        raise ValueError(f"На позиции {pos} ожидалась {wt_aa}, но в последовательности '{seq[pos-1]}'")

    mutated = seq[:pos-1] + new_aa + seq[pos:]
    seq_entry['protein']['sequence'] = mutated

    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


def main():
    doc_file_list = [r"C:\Users\Redmi\Desktop\НИР\Мутации\str_WT\мутации_ароматические.txt",
r"C:\Users\Redmi\Desktop\НИР\Мутации\str_WT\мутации_неполярные.txt",
r"C:\Users\Redmi\Desktop\НИР\Мутации\str_WT\мутации_особенные.txt",
r"C:\Users\Redmi\Desktop\НИР\Мутации\str_WT\мутации_отрицательно_заряженные.txt",
r"C:\Users\Redmi\Desktop\НИР\Мутации\str_WT\мутации_положительно_заряженные.txt",
r"C:\Users\Redmi\Desktop\НИР\Мутации\str_WT\мутации_полярные.txt"]
    
    output_root_list = [r"C:\Users\Redmi\Desktop\НИР\Мутации\типы_аминокислот\Ароматические\сгенерированные",
r"C:\Users\Redmi\Desktop\НИР\Мутации\типы_аминокислот\Неполярные\сгенерированные",
r"C:\Users\Redmi\Desktop\НИР\Мутации\типы_аминокислот\Особенные",
r"C:\Users\Redmi\Desktop\НИР\Мутации\типы_аминокислот\Отрицательно_заряженные",
r"C:\Users\Redmi\Desktop\НИР\Мутации\типы_аминокислот\Положительно_заряженные\сгенерированные",
r"C:\Users\Redmi\Desktop\НИР\Мутации\типы_аминокислот\Полярные_незаряженные"]
    
    json_dir = r"C:\Users\Redmi\Desktop\НИР\Мутации\str_WT\file_script_PDB-json_clean"
    for i in range(6):
        doc_file = doc_file_list[i]

        output_root = output_root_list[i]
        os.makedirs(output_root, exist_ok=True)

        with open(doc_file, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                complex_name, mutation, mutation2 = line.split()
                # if complex_name == "2NYY_DC_A" and  mutation2 == "YD30D":
                #     mutation2 = mutation 
                mutation = mutation2
                folder = f"{complex_name}_{mutation}"
                out_dir = os.path.join(output_root, folder)
                os.makedirs(out_dir, exist_ok=True)

                src_json = os.path.join(json_dir, f"{complex_name}_WT.json")
                if not os.path.isfile(src_json):
                    print(f"WARNING: JSON для {complex_name} не найден: {src_json}")
                    continue

                # dst_json = os.path.join(out_dir, os.path.basename(src_json))
                # shutil.copy(src_json, dst_json)

                new_filename = f"{complex_name}_{mutation}_MT.json"
                dst_json = os.path.join(out_dir, new_filename)
                shutil.copy(src_json, dst_json)
                try:
                    apply_mutation(dst_json, complex_name, mutation)
                    print(f"Processed {complex_name} {mutation}")
                except Exception as e:
                    print(f"ERROR processing {complex_name} {mutation}: {e}")


if __name__ == '__main__':
    main()



# # a = "DVQLQESGPSLVKPSQTLSLTCSVTGDSITSDYWSWIRKFPGNRLEYMGYVSYSGSTYYNPSLKSRISITRDTSKNQYYLDLNSVTTEDTATYYCANWDGDYWGQGTLVTVSAAKTTPPSVYPLAPGSAAQTNSMVTLGCLVKGYFPEPVTVTWNSGSLSSGVHTFPAVLQSDLYTLSSSVTVPSSPRPSETVTCNVAHPASSTKVDKKIVPRDC"
# b = "DIVLTQSPATLSVTPGNSVSLSCRASQSIGENLHWYQQKSHESPRLLIKYASQSISGIPSRFSGSGSGTDFTLSINSVETEDFGMYFCQQSNSWPYTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFYPKDINVKWKIDGSERQNGVLNSWTDQDSKDSTYSMSSTLTLTKDEYERHNSYTCEATHKTSTSPIVKSFNRNEC"
# print(b[30])