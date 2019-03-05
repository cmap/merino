import argparse
import bisect
import logging
import os
import sys

import compute_wtcs as compute_wtcs
import plot_enrichment_score as plot_enrichment_score
import cmapPy.pandasGEXpress.parse as pe
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import merino.setup_logger as setup_logger
import numpy as np
import pandas as pd
import pestle.cmap.io.gmt as gmt
from statsmodels.distributions.empirical_distribution import ECDF

logger = logging.getLogger(setup_logger.LOGGER_NAME)

invariants = ['661', '662', '663', '664', '665', '666', '667', '668', '669', '670', '671', '672', '673', '674', '675', '676', '677', '678', '679', '680']

def build_parser():

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # The following arguments are required. These are files that are necessary for assembly and which change
    # frequently between cohorts, replicates, etc.
    parser.add_argument("-gct",  help="path to data file",
                        type=str, required=True)
    parser.add_argument("-meta", "-m", help="path to metadata file",
                        type=str, required=True)
    parser.add_argument("-out", "-o",
                        help="Our folder for results",
                        type=str, required=True)
    parser.add_argument("-sense", "-s",
                        help="Path to expected sensitivity gmt if required",
                        type=str, default=None, required=False)
    parser.add_argument("-prefix_name", "-pref",
                        help="Column name for weave_prefix/prism_replicate",
                        type=str, default='weave_prefix',required=False)
    parser.add_argument("-index_col", "-ind",
                        help="Column name for index - sig_id/profile_id",
                        type=str, default='sig_id', required=False)
    parser.add_argument("-verbose", '-v', help="Whether to print a bunch of output", action="store_true", default=False)


    return parser

def reformat_gmt(gmt):

    sense = {}
    for x in gmt:

        sense[x['desc']] = x['sig']

    return sense


def make_sqi_map(sensitivity_map, col_meta, data):
    s_qi_map = {}

    for key in sensitivity_map:
        key = key.replace('_UP', '')

        if key in col_meta['pert_id'].tolist():

            for dose in col_meta[col_meta['pert_id'] == key]['pert_idose'].unique():
                for well in col_meta[(col_meta['pert_id'] == key) & (col_meta['pert_idose'] == dose)]['pert_well'].unique():

                    if type(dose) is float:
                        dose = 'ctl_vehicle'
                    s_value = data[col_meta[(col_meta['pert_id'] == key) & (col_meta['pert_idose'] == dose) & (col_meta['pert_well'] == well)].index]
                    s_value = s_value[~s_value.index.isin(invariants)]

                    if s_value.shape[1] != 1:
                        raise Exception('S Value has more than one column, should correspond to a single well')

                    s_value = s_value.iloc[:,0]

                    qi = [thing for thing in
                          np.where(s_value.dropna().index.astype(str).isin(sensitivity_map[key]))[0].tolist()]


                    s_qi_map[key + '_' + str(dose) + '_' + well] = [s_value.dropna().tolist(), qi]

    return s_qi_map


def wtks(gct, metadata, outfolder, gmt_path='/Users/elemire/Workspace/merino/full_sensitivities.gmt', group_col='weave_prefix'):
    gmt_file = gmt.read(gmt_path)

    sensitivity_map = reformat_gmt(gmt_file)

    metadata = metadata.loc[gct.data_df.columns]

    expected_sensitivity_ranks = []


    for rep in metadata[group_col].dropna().unique():


        if not os.path.exists(os.path.join(outfolder, rep)):
            os.mkdir(os.path.join(outfolder, rep))
            os.mkdir(os.path.join(outfolder, rep, 'enrichment_ecdf'))
            os.mkdir(os.path.join(outfolder, rep, 'mountain_plots'))

        rep_folder = os.path.join(outfolder, rep)

        ids = metadata[metadata[group_col] == rep].index

        data = gct.data_df[ids]

        data.index = [str(x) for x in data.index]

        col_meta = metadata.loc[ids]

        marks = {}

        s_qi_map = make_sqi_map(sensitivity_map, data=data, col_meta=col_meta)

        for key in s_qi_map:
            print key.split('_')[0]
            if key.split('_')[0] in col_meta['pert_id'].tolist():


                if len(s_qi_map[key][0]) > 0 and len(s_qi_map[key][1]) > 0:
                    #calculate enrichment score
                    sensitivity_score, cumsum = compute_wtcs.compute_wtcs(pd.Series(s_qi_map[key][0]), pd.Series(s_qi_map[key][1]))
                    #make mountain plot of enrichment score
                    plot_enrichment_score.plot_enrichment_score(sensitivity_score, cumsum,
                                                                title='Enrichment Score of {}, Set Size = {}'.format
                                                                    (col_meta[col_meta['pert_id'] == key.split('_')[0]]
                                                                        ['pert_iname'][0] ,len(s_qi_map[key][1])),
                                                                outfile =  os.path.join(rep_folder, 'mountain_plots', '{}_{}_dose'.format(col_meta
                                                                    [col_meta['pert_id'] == key.split('_')[0]]['pert_iname'][0], key.split('_')[1]) + '.png'))



                    bortezomib_scores = []
                    dmso_scores = []
                    pos_dex = col_meta[col_meta['pert_type'] == 'trt_poscon'].index.tolist()
                    neg_dex = col_meta[col_meta['pert_type'] == 'ctl_vehicle'].index.tolist()
                    bortez = data[pos_dex]
                    dmso = data[neg_dex]

                    brd = key.split('_')[0]
                    rids = sensitivity_map[brd]

                    ################################

                    for column in bortez:
                        temp_s = bortez[column].dropna()
                        temp_qi = [thing for thing in np.where(temp_s.dropna().index.isin(rids))[0].tolist()]

                        if (len(temp_s) - len(temp_qi)) < 1:
                            continue

                        y, x = compute_wtcs.compute_wtcs(temp_s, pd.Series(temp_qi))
                        bortezomib_scores.append(y)

                    for column in dmso:
                        temp_s = dmso[~dmso.index.isin(invariants)][column].dropna()
                        temp_qi = [thing for thing in np.where(temp_s.dropna().index.isin(rids))[0].tolist()]

                        if (len(temp_s) - len(temp_qi)) < 1:
                            continue

                        y, x = compute_wtcs.compute_wtcs(temp_s, pd.Series(temp_qi))
                        dmso_scores.append(y)

                    all_scores = []
                    for column in data:
                        temp_s = data[~data.index.isin(invariants)][column].dropna()
                        temp_qi = [thing for thing in np.where(temp_s.dropna().index.isin(rids))[0].tolist()]

                        if (len(temp_s) - len(temp_qi)) < 1:
                            continue

                        y, x = compute_wtcs.compute_wtcs(temp_s, pd.Series(temp_qi))
                        all_scores.append(y)

                    ecdf = ECDF(all_scores)
                    mark = bisect.bisect_left(ecdf.x, sensitivity_score)
                    if mark == len(ecdf.y):
                        mark -= 1
                    marks[col_meta[col_meta['pert_id'] == key.split('_')[0]]['pert_iname'][0] + '_' + str(key.split('_')[1])] = mark

                    poscons_x = []
                    poscons_y = []

                    for pos in bortezomib_scores:
                        y = bisect.bisect_left(ecdf.x, pos)
                        poscons_x.append(pos)
                        poscons_y.append(ecdf.y[y - 1])

                    neg_x = []
                    neg_y = []

                    for neg in dmso_scores:
                        y = bisect.bisect_left(ecdf.x, neg)
                        neg_x.append(neg)
                        neg_y.append(ecdf.y[y - 1])

                    plt.figure()
                    plt.plot(ecdf.x, ecdf.y)
                    plt.scatter(poscons_x, poscons_y, marker='o', color='g', label='Bortezomib')
                    # plt.scatter(neg_x, neg_y, marker='o', color='b', label='DMSO')
                    plt.scatter(sensitivity_score, ecdf.y[mark], marker='o', color='r', label='Sensitivity {}_dose'.format(key.split('_')[1]))
                    plt.xlabel('Enrichment Score for {} Sensitivities Compound Rank = {}'.format
                        (col_meta[col_meta['pert_id'] == key.split('_')[0]]['pert_iname'][0], mark))
                    plt.ylabel('Fraction of Compounds')
                    plt.title('ECDF of Enrichment Score by Compound, Set Size = {}'.format(len(s_qi_map[key][1])))

                    axes = plt.gca()
                    axes.set_xlim([np.nanmin(ecdf.x[ecdf.x != -np.inf]), max(ecdf.x)])
                    axes.set_ylim([0, 1])
                    axes.legend(bbox_to_anchor=(0., 0.8, 0.8, .102), loc=3, borderaxespad=0.)
                    plt.savefig(os.path.join(rep_folder, 'enrichment_ecdf', '{}_Competitive_Enrichment_ECDF_{}_dose.png'.format
                                                 (col_meta[col_meta['pert_id'] == key.split('_')[0]]['pert_iname'][0], key.split('_')[1])))
                    plt.clf()



        plate_summary = pd.Series(marks)
        # plate_summary.set_index('det_plate', inplace=True)
        plate_summary.rename(rep, inplace=True)
        pd.DataFrame(plate_summary).to_csv(os.path.join(rep_folder, '{}_expected_sensitivity_ranks.txt'.format(rep)), sep='\t')
        expected_sensitivity_ranks.append(marks)

        marks['det_plate'] = rep


    summary = pd.DataFrame(expected_sensitivity_ranks)
    summary.set_index('det_plate', inplace=True)

    summary.to_csv(os.path.join(outfolder, 'expected_sensitivity_ranks.txt'), sep='\t')


def main(args):
    data = pe.parse(args.gct)
    meta = pd.read_table(args.meta, index_col=args.index_col)
    if args.sense is not None:
        wtks(data, meta, args.out, args.sense, group_col=args.prefix_name)
    else:
        wtks(data, meta, args.out, group_col=args.prefix_name)

if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    logger.debug("args:  {}".format(args))

    main(args)

