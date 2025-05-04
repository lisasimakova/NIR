import os
import json
import pytest

EXPECTED_DIR = r"C:\Users\Redmi\Desktop\НИР\Мутации\str_WT\Структуры_WT"

OUTPUT_DIR = r"C:\Users\Redmi\Desktop\НИР\Мутации\типы_аминокислот\Положительно_заряженные\готовые\проверка"

def normalize(obj: dict) -> dict:

    if "sequences" in obj and isinstance(obj["sequences"], list):
        obj["sequences"] = sorted(
            obj["sequences"],
            key=lambda entry: entry["protein"]["id"]
        )
    return obj

def get_all_json_files(directory: str):
    for root, _, files in os.walk(directory):
        for fname in files:
            if fname.lower().endswith(".json"):
                yield os.path.join(root, fname)

def find_expected_file(structure_id: str, directory: str) -> str:
    for root, _, files in os.walk(directory):
        for fname in files:
            if fname.lower().endswith(".json") and structure_id in fname and "WT" in fname:
                return os.path.join(root, fname)
    return None

def extract_structure_id(output_fname: str) -> str:

    parts = output_fname.split("_")
    if len(parts) >= 3:
        return "_".join(parts[:3])  # первые три части
    else:
        raise ValueError(f"Неожиданное имя файла: {output_fname}")

@pytest.fixture(params=list(get_all_json_files(OUTPUT_DIR)))
def pair(request):
    output_path = request.param
    output_fname = os.path.basename(output_path)

    structure_id = extract_structure_id(output_fname)

    expected_path = find_expected_file(structure_id, EXPECTED_DIR)
    if expected_path is None:
        pytest.skip(f"Референсный файл с ID '{structure_id}' не найден для файла {output_fname}")

    with open(expected_path, "r", encoding="utf-8") as f:
        expected = json.load(f)
    with open(output_path, "r", encoding="utf-8") as f:
        output = json.load(f)

    return output_fname, normalize(expected), normalize(output)

def test_chains_match(pair):
    fname, expected, output = pair
    if expected.get("sequences") != output.get("sequences"):
        diff = {
            "expected": expected.get("sequences"),
            "output": output.get("sequences")
        }
        pytest.fail(
            f"\n\nФайл {fname} НЕ совпадает по цепям!\nРазличия:\n{json.dumps(diff, indent=2, ensure_ascii=False)}"
        )
    assert expected.get("sequences") == output.get("sequences")

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
