import seaborn as sns
import pandas as pd
import string
import matplotlib.pyplot as plt



def mk_heatmap(df, title, outfile):

    values = df.median(axis=0)
    heatmap_df = pd.DataFrame(index=list(string.ascii_uppercase)[0:16])
    columns = [str.zfill(str(x), 2) for x in range(1, 25)]

    for col in columns:
        curr_column = values[[x[-2:] == col for x in values.index]]
        curr_column.index = [y[-3] for y in curr_column.index]
        heatmap_df[col] = curr_column

    sns.heatmap(heatmap_df, linewidths=.1, cmap="coolwarm")
    plt.yticks(rotation=1)
    plt.title(title)
    plt.savefig(outfile)
    plt.clf()