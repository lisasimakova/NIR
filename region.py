import os
import re
from collections import defaultdict

FOLDER = r"C:\Users\Redmi\Desktop\НИР\Мутации\chotia_clean"

CHOTHIA_REGIONS = {
    'L': {
        'FR1':    (1,  23),
        'CDR-L1': (24, 34),
        'FR2':    (35, 49),
        'CDR-L2': (50, 56),
        'FR3':    (57, 88),
        'CDR-L3': (89, 97),
        'FR4':    (98,107),
    },
    'H': {
        'FR1':    (1,  25),
        'CDR-H1': (26, 35),
        'FR2':    (36, 51),
        'CDR-H2': (52, 56),
        'FR3':    (57, 94),
        'CDR-H3': (95,102),
        'FR4':    (103,113),
    }
}

def parse_chothia_file(path):
    chains = defaultdict(list)
    seq_counter = defaultdict(int)

    with open(path, encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('//'):
                continue
            parts = line.split()
            chain = parts[0]
            if chain not in ('H', 'L'):
                continue
            print(parts)
            pos_label = parts[1]
            m = re.match(r"^(\d+)", pos_label)
            if not m:
                continue
            chothia_num = int(m.group(1))

            if parts[2].isalpha():
                seq_counter[chain] += 1
            chains[chain].append((chothia_num, seq_counter[chain]))

    # print(chains)
    return chains

def write_summary(folder):
    summary_path = os.path.join(folder, "regions_summary.txt")
    with open(summary_path, 'w', encoding='utf-8') as out:
        for fn in sorted(os.listdir(folder)):
            # if "1JRH" not in fn:
            #     continue
            if not fn.lower().endswith(".txt"):
                continue

            inp = os.path.join(folder, fn)
            chains = parse_chothia_file(inp)

            out.write(f"=== {fn} ===\n")
            for chain_label in ('L', 'H'):
                if chain_label not in chains:
                    out.write(f"Chain {chain_label}: (not found)\n\n")
                    continue

                out.write(f"Chain {chain_label}:\n")
                positions = chains[chain_label]
                for region, (start, end) in CHOTHIA_REGIONS[chain_label].items():
                    # из всех (chothia_num, seq_idx) выбираем по диапазону chothia_num
                    idxs = [seq_idx for ch_num, seq_idx in positions if start <= ch_num <= end]
                    if idxs:
                        # такое условие нужно при случае, что первый остаток = "-"
                        if min(idxs) == 0: 
                            out.write(f"  {region}: {1}–{max(idxs)}\n")
                        else:
                            out.write(f"  {region}: {min(idxs)}–{max(idxs)}\n")
                    else:
                        out.write(f"  {region}: (no residues)\n")
                out.write("\n")
            out.write("\n")

    print(f"Сводный отчёт записан в: {summary_path}")

if __name__ == "__main__":
    write_summary(FOLDER)
