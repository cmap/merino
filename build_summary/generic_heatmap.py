import seaborn as sns
import pandas as pd
import string
import matplotlib.pyplot as plt


def mk_heatmap(df, title, outfile, lims=[], colormap='coolwarm'):
    values = df.median(axis=0)
    #values = df
    heatmap_df = pd.DataFrame(index=list(string.ascii_uppercase)[0:16])
    columns = [str.zfill(str(x), 2) for x in range(1, 25)]

    for col in columns:
        curr_column = values[[x[-2:] == col for x in values.index]]
        curr_column.index = [y[-3] for y in curr_column.index]
        try:
            heatmap_df[col] = curr_column
        except:
            print 'skipping {}'.format(title)
            return
    if len(lims) == 0:
        sns.heatmap(heatmap_df, linewidths=.1, cmap=colormap)
    elif len(lims) == 2:
        sns.heatmap(heatmap_df, linewidths=.1, cmap=colormap, vmin=lims[0], vmax=lims[1])
    else:
        raise Exception('Must pass exactly 2 color scale limits or none at all')
    plt.yticks(rotation=1)
    plt.title(title)
    plt.savefig(outfile)
    plt.clf()