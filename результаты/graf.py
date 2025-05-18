import pandas as pd
import matplotlib.pyplot as plt
import os

excel_path = r"C:\Users\Redmi\Desktop\NIR\mutation\tsv_tables\result.xlsx"
path = r"C:\Users\Redmi\Desktop\NIR\mutation\tsv_tables"

sheets = pd.read_excel(excel_path, sheet_name=None)

sheet_names = list(sheets.keys())[:6]

def plot_boxplots_from_sheets(sheet_dict, sheet_list, nrows, ncols, figure_title, number):

        data = []
        titles = []
        for name in sheet_list:
            df = sheet_dict[name]

            col = pd.to_numeric(df.iloc[:, number], errors='coerce').dropna()
            data.append(col)
            titles.append(name)

        fig, axes = plt.subplots(nrows=nrows, ncols=ncols,
                                figsize=(4*ncols, 3*nrows), squeeze=False)
        fig.suptitle(figure_title, fontsize=16)

        for idx, series in enumerate(data):
            i = idx // ncols
            j = idx % ncols
            axes[i][j].boxplot(series)
            axes[i][j].set_title(titles[idx])
            axes[i][j].set_ylabel("RMSD")

        for idx in range(len(data), nrows*ncols):
            i = idx // ncols
            j = idx % ncols
            axes[i][j].axis('off')

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        
        safe_title = figure_title.replace(" ", "_")
        file_path = os.path.join(path, f"{safe_title}.png")
        plt.savefig(file_path, dpi=300)
        plt.show()

a = ['RMSD по всей структуре',	'RMSD мутировавшей цепи','RMSD мутировавшего региона']
for i, z in enumerate([1, 4, 6]):
    plot_boxplots_from_sheets(sheets, sheet_names, nrows=2, ncols=3,
                            figure_title=a[i], number=z)
