import pandas as pd
import python.broadinstitute_cmap.io.pandasGEXpress.parse as pe
from statsmodels.distributions.empirical_distribution import ECDF
import matplotlib.pyplot as plt
import os


def make_plots(gct_filepath, sensitivity_map, outfile):

    gct = pe.parse(gct_filepath)
    zscores = []
    for key in sensitivity_map:
        for y in sensitivity_map[key]:
            zscores.append(gct.data_df[gct.col_metadata_df[gct.col_metadata_df['pert_id'] == key].index].ix[y].tolist())

    flat_list = [item for sublist in zscores for item in sublist]


    ecdf = ECDF(flat_list)

    plt.plot(ecdf.x, ecdf.y)
    plt.xlabel('Zscore')
    plt.ylabel('Portion of Compounds')
    plt.title('ECDF of Zscores for Expected Sensitivities')

    axes = plt.gca()
    axes.legend(bbox_to_anchor=(0., 0.9, 0.8, .102), loc=3, borderaxespad=0.)
    plt.savefig(os.path.join(outfile, 'sensor_zscore_ecdf.png'))

    ranks = []
    for key in sensitivity_map:
        for y in sensitivity_map[key]:
            ranks.append(
                gct.data_df.rank(axis=1)[gct.col_metadata_df[gct.col_metadata_df['pert_id'] == key].index].ix[y].tolist())

    flat_rank_list = [item for sublist in ranks for item in sublist]

    ecdf = ECDF(flat_rank_list)

    plt.plot(ecdf.x, ecdf.y)
    plt.xlabel('Zscore Rank (Vs. All other cell lines exposed to the compound)')
    plt.ylabel('Portion of Compounds')
    plt.title('ECDF of Zscore Ranks for Expected Sensitivities')

    axes = plt.gca()
    axes.legend(bbox_to_anchor=(0., 0.9, 0.8, .102), loc=3, borderaxespad=0.)
    plt.savefig(os.path.join(outfile, 'sensor_zscore_ranks_ecdf.png'))