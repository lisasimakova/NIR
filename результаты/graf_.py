import pandas as pd
import matplotlib.pyplot as plt
import os
import sys

excel_path = r"C:\Users\Redmi\Desktop\NIR\mutation\tsv_tables\result.xlsx"
path = r"C:\Users\Redmi\Desktop\NIR\mutation\tsv_tables"
log_path = os.path.join(path, "log.txt")
sys.stdout = open(log_path, 'w', encoding='utf-8') 
sheets = pd.read_excel(excel_path, sheet_name=None)

sheet_names = list(sheets.keys())[:6]

def plot_rmsd_vs_ddg_from_sheets(sheet_dict, sheet_list, nrows, ncols, figure_title, rmsd_col_number):
    print(f'============={figure_title}==============')
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols,
                            figsize=(4*ncols, 3*nrows), squeeze=False)
    fig.suptitle(figure_title, fontsize=16)

    for idx, name in enumerate(sheet_list):
        df = sheet_dict[name]

        rmsd = pd.to_numeric(df.iloc[:, rmsd_col_number], errors='coerce')
        ddg = pd.to_numeric(df.iloc[:, 11], errors='coerce')

        valid_idx = rmsd.dropna().index.intersection(ddg.dropna().index)
        rmsd_valid = rmsd.loc[valid_idx]
        ddg_valid = ddg.loc[valid_idx]

        i = idx // ncols
        j = idx % ncols
        ax = axes[i][j]

        ax.scatter(ddg_valid, rmsd_valid, alpha=0.5, s=10)
        ax.set_title(name)
        ax.set_xlabel('ddG')
        ax.set_ylabel('RMSD')

        ax.grid(True)
    for idx in range(len(sheet_list), nrows*ncols):
        i = idx // ncols
        j = idx % ncols
        axes[i][j].axis('off')

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    safe_title = figure_title.replace(" ", "_")
    file_path = os.path.join(path, f"{safe_title}_ddg_vs_rmsd.png")
    plt.savefig(file_path, dpi=300)
    plt.show()

a = ['RMSD по всей структуре', 'RMSD мутировавшей цепи', 'RMSD мутировавшего региона']
for i, z in enumerate([1, 4, 6]):
    plot_rmsd_vs_ddg_from_sheets(sheets, sheet_names, nrows=2, ncols=3,
                            figure_title=a[i], rmsd_col_number=z)

sys.stdout.close()
sys.stdout = sys.__stdout__
