import os
import json
import pytest

# Путь к папке с «правильными» JSON
# EXPECTED_DIR = r"C:\Users\Redmi\Desktop\НИР\Мутации\str_WT\json_проверенные"
EXPECTED_DIR = r"C:\Users\Redmi\Desktop\НИР\Мутации\str_WT\файлы_script_no_clean"

# Путь к папке, где лежат только что сгенерированные JSON
OUTPUT_DIR   = r"C:\Users\Redmi\Desktop\НИР\Мутации\str_WT\file_script_PDB-json_clean"

def normalize(obj: dict) -> dict:
    if "sequences" in obj and isinstance(obj["sequences"], list):
        obj["sequences"] = sorted(
            obj["sequences"],
            key=lambda entry: entry["protein"]["id"]
        )
    return obj

@pytest.fixture(params=[
    fn for fn in os.listdir(EXPECTED_DIR)
    if fn.lower().endswith(".json")
])
def pair(request):
    fname = request.param
    exp_path = os.path.join(EXPECTED_DIR, fname)
    out_path = os.path.join(OUTPUT_DIR, fname)
    if not os.path.exists(out_path):
        pytest.skip(f"Нет сгенерированного файла {fname}")
    with open(exp_path, "r", encoding="utf-8") as f:
        expected = json.load(f)
    with open(out_path, "r", encoding="utf-8") as f:
        output = json.load(f)
    return fname, normalize(expected), normalize(output)

def test_jsons_match(pair):
    fname, expected, output = pair
    assert expected == output, (
        f"\n\nФайл {fname} НЕ совпадает!\n"
        f"Ожидалось:\n{json.dumps(expected, indent=2, ensure_ascii=False)}\n\n"
        f"Получено:\n{json.dumps(output, indent=2, ensure_ascii=False)}"
    )
if __name__ == "__main__":
    import subprocess
    from datetime import datetime

    log_dir = "pytest_logs"
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, f"log_{datetime.now():%Y-%m-%d_%H-%M-%S}.txt")

    with open(log_file, "w", encoding="utf-8") as f:
        result = subprocess.run(
            ["pytest", __file__, "--tb=long", "--capture=tee-sys"],
            stdout=f,
            stderr=subprocess.STDOUT
        )
        exit(result.returncode)
