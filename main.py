import os
import sys
from Bio.PDB import MMCIFParser, Superimposer

from skrip_pars_cif import extract_sequence_from_cif

class Tee:
    def __init__(self, *files):
        self.files = files
    def write(self, msg):
        for f in self.files:
            f.write(msg)
    def flush(self):
        for f in self.files:
            f.flush()


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
                name, rng = line.split(':',1)
                rng = rng.strip().replace('–','-')
                start, end = map(int, rng.split('-'))
                regions[current][chain][name.strip()] = (start, end)
    return regions


def find_first_cif(root):

    for dirpath, dirnames, filenames in os.walk(root):
        for fn in filenames:
            if fn.lower().endswith('.cif'):
                return os.path.join(dirpath, fn)
    return None

def find_first_json(root):

    for dirpath, dirnames, filenames in os.walk(root):
        for fn in filenames:
            if fn.lower().endswith('.json'):
                return os.path.join(dirpath, fn)
    return None

def compute_rmsd(wt_cif, mt_cif, regions_for_id, out):
    parser = MMCIFParser(QUIET=True)
    struct_wt = parser.get_structure('WT', wt_cif)
    struct_mt = parser.get_structure('MT', mt_cif)

    global_wt_atoms = []
    global_mt_atoms = []

    for chain_label in ('L','H'):
        out.write(f"--- Chain {chain_label} ---\n")
        ranges = regions_for_id.get(chain_label, {})
        try:
            ch_wt = struct_wt[0][chain_label]
            ch_mt = struct_mt[0][chain_label]
        except KeyError:
            out.write(f"Chain {chain_label} not found in one of the structures\n\n")
            continue

        seq_wt = [res['CA'] for res in ch_wt if 'CA' in res]
        seq_mt = [res['CA'] for res in ch_mt if 'CA' in res]
        # print(seq_wt)
        chain_wt_atoms = []
        chain_mt_atoms = []
        cur = 1
        for region, (start, end) in ranges.items():
            pos = [i for i in range(start, end+1)
                   if 1 <= i <= len(seq_wt) and 1 <= i <= len(seq_mt)]
            if not pos:
                out.write(f"  {region}: (no residues)\n")
                continue

            wt_sel = [seq_wt[i-1] for i in pos]
            mt_sel = [seq_mt[i-1] for i in pos]

            sup = Superimposer()
            sup.set_atoms(mt_sel, wt_sel)
            out.write(f"  {region}: RMSD = {sup.rms:.3f} Å ({cur} - {cur + len(pos)-1} residues)\n")
            cur += len(pos)
            chain_wt_atoms += wt_sel
            chain_mt_atoms += mt_sel
            global_wt_atoms += wt_sel
            global_mt_atoms += mt_sel

        if chain_wt_atoms:
            sup_chain = Superimposer()
            sup_chain.set_atoms(chain_mt_atoms, chain_wt_atoms)
            out.write(f"  >> Chain {chain_label} global RMSD: {sup_chain.rms:.3f} Å\n")
        out.write("\n")

    if global_wt_atoms:
        sup_all = Superimposer()
        sup_all.set_atoms(global_mt_atoms, global_wt_atoms)
        out.write(f"=== Overall RMSD (all chains): {sup_all.rms:.3f} Å ===\n\n")

def compare(wt_cif, mt_cif):
    dict_wt = extract_sequence_from_cif(wt_cif)
    dict_mt = extract_sequence_from_cif(mt_cif)
    # print(dict_wt)
    cur_h = 0
    for i in range(len(dict_wt["H"])):
        if dict_mt["H"][i]  != dict_wt["H"][i]:
            cur_h += 1
    cur_l = 0
    for i in range(len(dict_wt["L"])):
        if dict_mt["L"][i]  != dict_wt["L"][i]:
            cur_l += 1 
    try:
        if cur_h + cur_l == 1:
            return cur_h, cur_l
        else:
            raise ValueError(f"Сумма значений ({cur_h} + {cur_l} = {cur_h + cur_l}) не равна 1")
    except TypeError as e:
        raise TypeError(f"Один из аргументов не является числом: cur_h={type(cur_h)}, cur_l={type(cur_l)}") from e
    except Exception as e:
        raise Exception(f"Произошла непредвиденная ошибка: {str(e)}") from e





def main():
    types_base = r"C:\Users\Redmi\Desktop\НИР\Мутации\типы_аминокислот"
    wt_base    = r"C:\Users\Redmi\Desktop\НИР\Мутации\str_WT\Структуры_WT"
    regions_file = r"C:\Users\Redmi\Desktop\НИР\Мутации\chotia_clean\regions_summary.txt"

    # парсим все разметки
    regions = parse_chothia_regions(regions_file)

    # общий лог
    all_log_path = os.path.join(types_base, "all_rmsd.log")
    all_log = open(all_log_path, "w", encoding="utf-8")

    for aa_type in os.listdir(types_base):
        aa_dir = os.path.join(types_base, aa_type)
        if not os.path.isdir(aa_dir):
            continue

        for struct_folder in os.listdir(aa_dir):
            struct_dir = os.path.join(aa_dir, struct_folder)
            if not os.path.isdir(struct_dir):
                continue
            print(struct_dir)

            pdb_id = struct_folder.split('_')[0].upper()
            wt_cif = None
            for wt_sub in os.listdir(wt_base):
                if wt_sub.upper().startswith(pdb_id + "_") and wt_sub.upper().endswith("_WT"):
                    wt_cif = find_first_cif(os.path.join(wt_base, wt_sub))
                    break
            print(wt_cif)

            # wt_json = None
            # for wt_sub in os.listdir(wt_base):
            #     if wt_sub.upper().startswith(pdb_id + "_") and wt_sub.upper().endswith("_WT"):
            #         wt_json = find_first_json(os.path.join(wt_base, wt_sub))
            #         break
            # print(wt_json)

            mt_cif = None
            for sub in os.listdir(struct_dir):
                subp = os.path.join(struct_dir, sub)
                if os.path.isdir(subp) and "_MT" in sub.upper():
                    mt_cif = find_first_cif(subp)
                    break
            print(mt_cif)
            # for sub in os.listdir(struct_dir):
            #     subp = os.path.join(struct_dir, sub)
            #     if os.path.isdir(subp) and "_MT" in sub.upper():
            #         mt_json = find_first_json(subp)
            #         break
            # print(mt_json)

            
            if not (wt_cif and mt_cif):
                all_log.write(f"Skip {struct_folder!r}: WT or MT .cif not found\n")
                continue
            cur_h, cur_l = compare(wt_cif, mt_cif)
            if pdb_id not in regions:
                all_log.write(f"Warning: regions for {pdb_id} not found, skipping regions\n")
                continue
            print(pdb_id)
            mutant_log_path = os.path.join(struct_dir, "rmsd.log")
            mlog = open(mutant_log_path, "w", encoding="utf-8")

            tee = Tee(sys.stdout, mlog, all_log)
            sys.stdout = tee

            print(f"\n*** {aa_type}/{struct_folder} ***, cur_h = {cur_h}, cur_l = {cur_l}")
            compute_rmsd(wt_cif, mt_cif, regions[pdb_id], out=tee)

            sys.stdout = sys.__stdout__
            mlog.close()
        #     break
        # break
    

    all_log.close()

if __name__ == "__main__":
    main()
