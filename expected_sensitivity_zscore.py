import pandas as pd
import python.broadinstitute_cmap.io.pandasGEXpress.parse as pe
from statsmodels.distributions.empirical_distribution import ECDF
import matplotlib.pyplot as plt
import os
import numpy as np
import bisect
import matlab.engine
eng1 = matlab.engine.start_matlab()


invariants = ['661', '662', '663', '664', '665', '666', '667', '668', '669', '670', '671', '672', '673', '674', '675', '676', '677', '678', '679', '680']



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

    plt.clf()

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

    plt.clf()

def wtks(gct_filepath, sensitivity_map, outfile):
    gct = pe.parse(gct_filepath)

    invariants = ['661', '662', '663', '664', '665', '666', '667', '668', '669', '670']

    s_qi_map = {}

    for key in sensitivity_map:
        s = gct.data_df[
            gct.col_metadata_df[(gct.col_metadata_df['pert_id'] == key) & (gct.col_metadata_df['pert_dose'] == 2.5)].index]
        s = s[~s.index.isin(invariants)]
        qi = [thing + 1 for thing in np.where(s.index.isin(sensitivity_map[key]))[0].tolist()]
        s_qi_map[key] = [s.unstack().tolist(), qi]

    for key in s_qi_map:
        os.system(
            "/Applications/MATLAB_R2016a.app/bin/matlab -nodesktop -nosplash -nodisplay -r \"[wtks, maxi, rs, r] = compute_wtks({}, {}'); plot_es(wtks, maxi, rs, r); saveas(gcf, \'{}.png\'); clf; exit;\"".format(
                s_qi_map[key][0], s_qi_map[key][1], os.path.join(outfile, key)))
        sensitivity_score = eng1.compute_wtks(matlab.double(s_qi_map[key][0]), matlab.double(s_qi_map[key][1]))
        all_scores = []
        for column in gct.data_df:
            y = eng1.compute_wtks(matlab.double(x.data_df[column].tolist()), matlab.double(s_qi_map[key][1]))
            all_scores.append(abs(y))

        ecdf = ECDF(all_scores)
        mark = bisect.bisect(ecdf.x, sensitivity_score)
        plt.figure()
        plt.plot(ecdf.x, ecdf.y)
        plt.scatter(sensitivity_score, ecdf.y[mark], marker='o', color='r', label='Sensitivity')
        plt.xlabel('Enrichment Score for {} Sensitivities'.format(key))
        plt.ylabel('Fraction of Compounds')
        plt.title('ECDF of Enrichment Score by Compound')

        axes = plt.gca()
        axes.set_xlim([0, max(ecdf.x)])
        axes.set_ylim([0, 1])
        axes.legend(bbox_to_anchor=(0., 0.9, 0.8, .102), loc=3, borderaxespad=0.)
        plt.savefig(os.path.join(outfile, '{}_ECDF_of_enrichment_for_all_other_compounds.png'.format(key)))
