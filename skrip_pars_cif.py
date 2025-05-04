from pathlib import Path
from Bio.PDB import MMCIFParser, Polypeptide

def extract_sequence_from_cif(cif_file: str) -> dict:

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('protein', cif_file)
    
    sequences = {}
    
    for model in structure:
        for chain in model:
            ppb = Polypeptide.PPBuilder()
            peptides = ppb.build_peptides(chain)
            sequence = "".join([str(p.get_sequence()) for p in peptides])
            sequences[chain.id] = sequence
    
    return sequences

# def print_sequences(cif_file: str):
#     """
#     Печатает последовательности цепей из CIF файла в удобном формате.
#     """
#     sequences = extract_sequence_from_cif(cif_file)
    
#     print(f"\nSequences from {Path(cif_file).name}:")
#     for chain_id, seq in sequences.items():
#         print(f"Chain {chain_id}: {seq}")
#     print(f"\nTotal chains: {len(sequences)}")
