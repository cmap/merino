import sys
sys.path.append('/Users/elemire/Workspace/l1ktools')
import python.broadinstitute_cmap.io.pandasGEXpress.parse as pe
import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns
import pandas as pd
import string

invariants = ['661', '662', '663', '664', '665', '666', '667', '668', '669', '670', '671', '672', '673', '674', '675', '676', '677', '678', '679', '680']



def data_distribution(data_df, xlabel, title, outfile, xlim):

    print title
    bins = np.linspace(xlim[0], xlim[1], 50)
    ax = sns.distplot(data_df.unstack().dropna(), bins=bins, color='g')
    ax.set(xlabel=xlabel, ylabel='Frequency')
    plt.title(title)
    plt.savefig(outfile)
    plt.clf()



def stacked_heatmap(data_df, column_metadata, title, outfile, lims):

    my_palette = sns.cubehelix_palette(8, start=0.5, rot=-.8, dark=0.2, light=.90)


    data_df.columns = column_metadata['pert_well']
    values = data_df.unstack().reset_index().pivot_table(columns='pert_well', values=0, aggfunc=np.median)

    values_unstacked = values.unstack().reset_index().drop('level_1', axis=1)
    values_unstacked['well_row'] = [x[1:] for x in values_unstacked['pert_well']]
    values_unstacked['well_col'] = [x[0] for x in values_unstacked['pert_well']]
    heatmap_df = values_unstacked.pivot_table(values=0, index='well_row', columns='well_col')

    sns.heatmap(heatmap_df, linewidths=.1, cmap='coolwarm', vmin=lims[0], vmax=lims[1])
    plt.yticks(rotation=1)
    plt.title(title)
    plt.savefig(outfile)
    plt.clf()

    return

