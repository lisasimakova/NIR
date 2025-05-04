import json
import os
from pathlib import Path

json_folder = r"C:\Users\Redmi\Desktop\НИР\Мутации\str_WT\file_script_PDB-json_clean"
output_folder = r"C:\Users\Redmi\Desktop\НИР\Мутации\str_WT\fasta_clean"

Path(output_folder).mkdir(parents=True, exist_ok=True)

for json_file in Path(json_folder).glob('*.json'):
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    base_name = data['name']
    fasta_filename = f"{base_name}.fasta"
    fasta_path = os.path.join(output_folder, fasta_filename)
    
    # Собираем все последовательности из файла
    sequences = []
    for seq_data in data['sequences']:
        chain_id = seq_data['protein']['id']
        sequence = seq_data['protein']['sequence']
        sequences.append((chain_id, sequence))
    
    # Записываем в FASTA-файл
    with open(fasta_path, 'w') as fasta_file:
        for chain_id, sequence in sequences:
            fasta_file.write(f">{base_name}_{chain_id}\n{sequence}\n\n")

print(f"Все файлы преобразованы и сохранены в папку {output_folder}")