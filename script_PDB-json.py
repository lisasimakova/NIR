import os
import re
import json
from typing import List, Dict, Any
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import PPBuilder


def parse_complex_name(name: str) -> List[str]:
    """
    Извлекаем идентификаторы цепей
    """
    m = re.search(r'(?<=_)([A-Z]{2})(?=_)', name) or \
        re.search(r'(?<=_)([A-Z]{2})$', name)
    if not m:
        raise ValueError(f"Не удалось из имени комплекса '{name}' извлечь пару идентификаторов цепей")
    return list(m.group(1))


def assign_H_L(chain_ids: List[str]) -> Dict[str, str]:
    """
    Вход — два chain ID (например ['C','F']).
    Если это ['H','L'] или ['L','H'] — оставляем напрямую.
    Иначе сортируем и даём первой L, второй H.
    """


    if set(chain_ids) == {'H', 'L'}:
        return {'H':'H', 'L':'L'}
    a, b = sorted(chain_ids)
    return {a: 'L', b: 'H'}


def parse_pdb_sequences(pdb_path: str) -> List[Dict[str, Any]]:
    """
    Читает PDB и возвращает список послеодвательностей
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(os.path.basename(pdb_path), pdb_path)
    ppb = PPBuilder()
    sequences: Dict[str, str] = {}

    # Работаем с первой моделью
    model = next(structure.get_models())
    for chain in model:
        seq = ""
        for pp in ppb.build_peptides(chain):
            seq += str(pp.get_sequence())
        if seq:
            sequences[chain.id] = seq

    return [{'ids': [chain_id], 'seq': seq} for chain_id, seq in sequences.items()]


def process_pdb(pdb_path: str, out_dir: str, complex_name: str) -> None:
    # пара chain IDs из имени комплекса
    chain_pair = parse_complex_name(complex_name)

    # парсим PDB 
    records = parse_pdb_sequences(pdb_path)
    print(f"Пары цепей в структуре: {[r['ids'] for r in records]}")

    # ищем порядок, чтобы каждая буква была в разной записи
    found_order = None
    for candidate in (chain_pair, chain_pair[::-1]):
        idxs = []
        for c in candidate:
            idx = next((i for i, r in enumerate(records) if c in r['ids']), None)
            if idx is None:
                break
            idxs.append(idx)
        if len(idxs) == 2 and idxs[0] != idxs[1]:
            found_order = candidate
            break

    if not found_order:
        raise KeyError(f"Не удалось сопоставить цепи {chain_pair} структуре {os.path.basename(pdb_path)}")

    # Этот комплекс отличается от других
    if complex_name == '3NPS_A_BC':
        mapping = {'C': 'L', 'B': 'H'}
    else:
        mapping = assign_H_L(found_order)

    # собираем записи для JSON
    entries = []
    for role in ['H', 'L']:
        orig_chain = next(k for k, v in mapping.items() if v == role)
        rec = next(r for r in records if orig_chain in r['ids'])
        entries.append({
            "protein": {
                "id": role,
                "sequence": rec['seq']
            }
        })

    # формируем и сохраняем JSON
    obj = {
        "name": f"{complex_name}_WT",
        "sequences": entries,
        "modelSeeds": [1],
        "dialect": "alphafold3",
        "version": 1
    }
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{complex_name}_WT.json")
    with open(out_path, 'w') as fw:
        json.dump(obj, fw, indent=2)
    print(f"Сохранён {out_path}")


def main():
    pdb_dir = r"C:\Users\Redmi\Desktop\НИР\Мутации\SKEMPI2_PDBs\PDBs\PDB"
    output_dir = r"C:\Users\Redmi\Desktop\НИР\Мутации\str_WT\file_script_PDB-json_clean"
    complexes_file = r"C:\Users\Redmi\Desktop\НИР\Мутации\str_WT\Обработка_файлов.txt"

    # Читаем список комплексов, которые нужно обработать 
    with open(complexes_file, 'r') as f:
        complexes = [line.strip() for line in f if line.strip()]

    # Список всех PDB-файлов
    pdb_files = [os.path.join(pdb_dir, fn)
                 for fn in os.listdir(pdb_dir)
                 if fn.lower().endswith('.pdb')]

    for complex_name in complexes:
        prefix = complex_name.split('_')[0]
        # Ищем соответствующий PDB
        for pdb_path in pdb_files:
            base = os.path.splitext(os.path.basename(pdb_path))[0]
            if prefix == base or prefix in base:
                process_pdb(pdb_path, output_dir, complex_name)

if __name__ == '__main__':
    main()
