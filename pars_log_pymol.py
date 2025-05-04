import re
import os
def pars_pml_file(pml_file):

    name_file = os.path.splitext(os.path.basename(pml_file))[0]
    return name_file
# Путь к входному файлу
input_file = r"C:\Users\Redmi\Desktop\NIR\mutation\output_PML.txt"
# Путь к выходному файлу
output_file =r"C:\Users\Redmi\Desktop\NIR\mutation\tsv_tables\parsed_rmsd.tsv"

# Открываем лог
with open(input_file, 'r', encoding='utf-8') as f:
    log_data = f.read()

# Регулярка для каждого случая
block_pattern = r"Running .*?\.pml.*?(?=Running |\Z)"  # делим по pml-скриптам
blocks = re.findall(block_pattern, log_data, flags=re.DOTALL)

# Ищем RMSD
rmsd_pattern = r"Executive: RMSD =\s+([0-9.]+) \(\d+ to \d+ atoms\)"
results = []

for block in blocks:
    script_match = re.search(r"Running\s+(.*\.pml)", block)
    script_name = script_match.group(1).strip() if script_match else "unknown"

    rmsds = re.findall(rmsd_pattern, block)
    
    if len(rmsds) >= 4:
        name_file = (pars_pml_file(script_name))
        total, heavy, light, region = rmsds[:4]
        results.append((name_file, total, heavy, light, region))
    else:
        print(f"Warning: not enough RMSDs found in {script_name}, found {len(rmsds)}")


                               

# Сохраняем в файл
with open(output_file, 'w', encoding='utf-8') as f:

    
    f.write("name_file\ttotal_RMSD_pml\theavy_RMSD_pml\tlight_RMSD_pml\tREGION_pml\n")
    for row in results:
        f.write("\t".join(row) + "\n")

print(f"Done. Parsed RMSDs saved to {output_file}")
