import sys
import pandas as pd
import seaborn as sns
import string
from broadinstitute_cmap.io.pandasGEXpress import parse as pe
import matplotlib.pyplot as plt
import os

'''
Make heatmap of median beadcount per well - used for cohort level analysis. Takes a build gct or gctx file and
associated metadata (most likely inst_info.txt)
'''

def make_count_heatmap(count_gct_path, metadata_path, title, output):

    print 'Reading GCT: {}'.format(count_gct_path)

    count_gct = pe.parse(count_gct_path)
    col_metadata_df = pd.read_table(metadata_path)

    concat_list = []
    print 'Rearranging by well'
    for x in col_metadata_df['prism_replicate'].unique():
        print x
        index = col_metadata_df[col_metadata_df['prism_replicate'] == x].index
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
    plt.clf()

    return