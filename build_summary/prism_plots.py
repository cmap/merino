import sys
import matplotlib
matplotlib.use("Agg")
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



def stacked_heatmap(df, column_metadata, title, outfile, lims, reduce_upper_limit=False):

    my_palette = sns.cubehelix_palette(8, start=0.5, rot=-.8, dark=0.2, light=.90)
    data_df = df

    if reduce_upper_limit is True:
        lims[1] = lims[1] / 2

    data_df.columns = column_metadata['pert_well']
    values = data_df.unstack().reset_index().pivot_table(columns='pert_well', values=0, aggfunc=np.median)

    values_unstacked = values.unstack().reset_index().drop('level_1', axis=1)
    values_unstacked['well_row'] = [x[1:] for x in values_unstacked['pert_well']]
    values_unstacked['well_col'] = [x[0] for x in values_unstacked['pert_well']]
    heatmap_df = values_unstacked.pivot_table(values=0, index='well_row', columns='well_col').T

    sns.heatmap(heatmap_df, linewidths=.1, cmap='coolwarm', vmin=lims[0], vmax=lims[1])
    plt.yticks(rotation=1)
    plt.title(title)
    plt.savefig(outfile)
    plt.clf()

    data_df.columns = column_metadata.index

    return


def sc_plot(signature_info, outfile):
    _types = signature_info['pert_type'].unique().tolist()
    colors = ['b', 'g', 'r', 'orange', 'y']
    data_df = signature_info[(signature_info.pert_type.isin(_types))]
    ax = sns.relplot(data=data_df, x='cc_q75', y='ss_ltn2', size='pert_idose', hue='pert_type', col='pert_type', palette=colors[0:len(_types)])
    ax.set(xlabel='Replicate Correlation (CC_Q75)', ylabel='Num. Sens. Cell Lines (SS_ltn2)')
    plt.savefig(outfile)
    plt.cla()
    plt.clf()
    ax.fig.set_size_inches(10, 10)

def modz_distribution(df, col_metadata_df, outfile):


    pos_dex = col_metadata_df[col_metadata_df['pert_type'] == 'trt_poscon'].index.tolist()
    neg_dex = col_metadata_df[col_metadata_df['pert_type'] == 'ctl_vehicle'].index.tolist()

    treatment = df.data_df[df.data_df.columns[~df.data_df.columns.isin(neg_dex + pos_dex)]]

    neg_df = df.data_df[neg_dex]
    pos_df = df.data_df[pos_dex]


    full_dist = treatment.unstack().dropna()
    pos_dist = pos_df.unstack().dropna()
    neg_dist = neg_df.unstack().dropna()
    full_dist.replace([np.inf, -np.inf], np.nan, inplace=True)
    full_dist.dropna(inplace=True)

    bines = np.linspace(min(full_dist), 4, 50)

    sns.distplot(neg_dist, bines, color='blue', norm_hist=True, label='DMSO')
    sns.distplot(full_dist, bines, color='orange', norm_hist=True, label='Treatment')
    sns.distplot(pos_dist, bines, color='red', norm_hist=True, label='Bortezomib')

    plt.xlim(-10, 10)
    plt.xlabel('MODZSPC')
    plt.ylabel('Frequency')
    plt.title('Distribution of MODZSPC by Compound')
    axes = plt.gca()
    axes.legend(bbox_to_anchor=(0., 0.8, 0.8, .102), loc=3, borderaxespad=0.)

    plt.savefig(outfile)

    plt.clf()

