import sys
import bisect
import python.broadinstitute_cmap.io.pandasGEXpress.parse as pe
import matplotlib.pyplot as plt
from statsmodels.distributions.empirical_distribution import ECDF
import numpy as np
import os
import pandas as pd
import seaborn as sns

def median_ZSPC_histogram(filepath, outfile, det ='', norm_data=False):

    # Make histogram of median ZSPC by compound and mark the location of bortezomib
    df = pe.parse(filepath)

    pos_dex = df.col_metadata_df[df.col_metadata_df['pert_type'] == 'trt_poscon'].index.tolist()



    full_dist = df.data_df.median(axis=0).dropna()
    pos_dist = df.data_df[pos_dex].median(axis=0).dropna()
    full_dist.replace([np.inf, -np.inf], np.nan, inplace=True)
    full_dist.dropna(inplace=True)
    bines = np.linspace(min(full_dist), max(full_dist), 30)


    n, bins, patches = plt.hist(full_dist, bines, facecolor='green', alpha=1)
    n0, bins0, patches0 = plt.hist(pos_dist, bines, facecolor='red', alpha=1, label='Bortezomib')

    plt.xlabel('ZSPC')
    plt.ylabel('Frequency')
    plt.title('{} Distribution of Median ZSPC by Compound'.format(det))
    axes = plt.gca()
    axes.legend(bbox_to_anchor=(0., 0.9, 0.8, .102), loc=3, borderaxespad=0.)

    plt.savefig(os.path.join(outfile, '{}median_zspc_hist.png'.format(det)))

    plt.clf()

def median_ZSPC_ecdf(filepath, outfile, det='', norm_data=False):
    # Make ECDF of median ZSPC by compound and mark the location of bortezomib
    df = pe.parse(filepath)


    pos_dex = df.col_metadata_df[df.col_metadata_df['pert_type'] == 'trt_poscon'].index.tolist()


    ecdf = ECDF(df.data_df.median(axis=0))
    poscons_x = []
    poscons_y = []

    for x in df.data_df[pos_dex].median(axis=0):
        y = bisect.bisect(ecdf.x, x)
        poscons_x.append(x)
        if y > 384:
            y = 384
        poscons_y.append(ecdf.y[y])

    plt.plot(ecdf.x, ecdf.y)
    plt.scatter(poscons_x, poscons_y, marker='o', color='r', label='Bortezomib')
    plt.xlabel('Median ZSPC')
    plt.ylabel('Portion of Compounds')
    plt.title('{} ECDF of Median ZSPC by Compound'.format(det))

    axes = plt.gca()
    axes.set_ylim([0, 1])
    axes.legend(bbox_to_anchor=(0., 0.9, 0.8, .102), loc=3, borderaxespad=0.)
    plt.savefig(os.path.join(outfile, '{}median_zspc_ecdf.png'.format(det)))

    plt.clf()


def kill_count(filepath, outfile, det=''):
    # Calculate number of cell lines killed per compound based on a ZSPC threshold of -1.5 and make ECDF
    narrow_kills = ['BRD-K66175015', 'BRD-K13154216', 'BRD-K19687926', 'BRD-K02113016', 'BRD-K51313569',
                    'BRD-K56343971', 'BRD-K05804044', 'BRD-A12230535', 'BRD-K59369769', 'BRD-K29905972',
                    'BRD-K09951645',
                    'BRD-K70401845', 'BRD-K64052750']

    df = pe.parse(filepath)

    kill_counts = df.data_df[df.data_df < -1].count()
    ecdf = ECDF(kill_counts)

    poscons_x = []
    poscons_y = []

    broad_x = []
    broad_y = []

    narrow_x = []
    narrow_y = []

    pos_dex = df.col_metadata_df[df.col_metadata_df['pert_type'] == 'trt_poscon'].index.tolist()

    broad_dex = df.col_metadata_df[(df.col_metadata_df['pert_id'] == 'BRD-K12343256') & (df.col_metadata_df['pert_dose'] == 10)].index.tolist()

    narrow_dex = df.col_metadata_df[df.col_metadata_df['pert_id'].isin(narrow_kills) & (df.col_metadata_df['pert_dose'] == df.col_metadata_df[df.col_metadata_df['pert_id'].isin(narrow_kills)]['pert_dose'].max())].index.tolist()

    bortez = df.data_df[pos_dex]
    trametinib_high = df.data_df[broad_dex]
    narrow_df = df.data_df[narrow_dex]

    for pos in bortez[bortez < -1].count():
        y = bisect.bisect_left(ecdf.x, pos)
        if y >= len(ecdf.y):
            y = 384
        poscons_x.append(pos)
        poscons_y.append(ecdf.y[y])

    for x in trametinib_high[trametinib_high < -1].count():
        y = bisect.bisect_left(ecdf.x, x)
        if y >= len(ecdf.y):
            y = 384
        broad_x.append(x)
        broad_y.append(ecdf.y[y])

    for x in narrow_df[narrow_df < -1].count():
        y = bisect.bisect_left(ecdf.x, x)
        if y >= len(ecdf.y):
            y = 384
        narrow_x.append(x)
        narrow_y.append(ecdf.y[y])

    unexpected_broad_killers = df.col_metadata_df.loc[narrow_df[narrow_df < -1].count()[narrow_df[narrow_df < -1].count() > 200].index][
        'pert_iname']
    kill_amount = narrow_df[narrow_df < -1].count()[narrow_df[narrow_df < -1].count() > 200]

    unexpected_broad_killers_file = pd.concat([unexpected_broad_killers, kill_amount], axis=1)

    unexpected_broad_killers_file.to_csv(os.path.join(outfile, 'unexpected_broad_killers.txt'), sep='\t')

    plt.plot(ecdf.x, ecdf.y, color='b')
    plt.scatter(narrow_x, narrow_y, marker='o', color='g', label='Narrow Killers Max Dose')
    plt.scatter(poscons_x, poscons_y, marker='o', color='r', label='Bortezomib')
    plt.scatter(broad_x, broad_y, marker='o', color='y', label='Trametinib High Dose (broad killer)')

    plt.xlabel('Compound Kill Count Based on ZSPC Threshold of -1')
    plt.ylabel('Portion of Compounds \n n={}'.format(len(kill_counts)))
    plt.title('{} ECDF of Cell Lines Killed by Compound'.format(det))

    axes = plt.gca()
    axes.set_xlim([0, max(ecdf.x)])
    axes.set_ylim([0, 1])
    axes.legend(bbox_to_anchor=(0.7, 0.65, 0.8, .102), loc=3, borderaxespad=0.)
    plt.savefig(os.path.join(outfile, '{}compound_strength_ecdf.png'.format(det)))

    plt.clf()

    return bortez[bortez < -1].count(), narrow_df[narrow_df < -1].count()