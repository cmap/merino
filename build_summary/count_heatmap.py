import sys
import pandas as pd
import seaborn as sns
import string
from broadinstitute_cmap.io.pandasGEXpress import parse as pe
import matplotlib.pyplot as plt
import os


def make_count_heatmap(count_gct, title, output)
    concat_list = []
    for x in count_gct.col_metadata_df['prism_replicate'].unique():
        index = count_gct.col_metadata_df[count_gct.col_metadata_df['prism_replicate'] == x].index
        temp_df = count_gct.data_df[index]
        temp_df.columns = [x[-3:] for x in temp_df.columns]
        concat_list.append(temp_df)

    data_stack = pd.concat(concat_list)

    values = data_stack.median(axis=0)
    heatmap_df = pd.DataFrame(index = list(string.ascii_uppercase)[0:16])
    columns = [str.zfill(str(x) ,2) for x in range(1 ,25)]
    for col in columns:
        curr_column = values[[x[-2:] == col for x in values.index]]
        curr_column.index = [y[0] for y in curr_column.index]
        heatmap_df[col] = curr_column
    sns.heatmap(heatmap_df, linewidths=.1, cmap="Reds", vmin=20, vmax=100)
    plt.yticks(rotation=1)
    plt.title(title)
    plt.savefig(output)
    plt.show()
    plt.clf()