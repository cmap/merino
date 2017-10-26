import prism_pipeline.compute_wtcs as compute_wtcs
import prism_pipeline.plot_enrichment_score as plot_enrichment_score
import pandas as pd
import cmapPy.pandasGEXpress.parse as pe
import os
import matplotlib.pyplot as plt
import numpy as np
import bisect
from statsmodels.distributions.empirical_distribution import ECDF

invariants = ['661', '662', '663', '664', '665', '666', '667', '668', '669', '670', '671', '672', '673', '674', '675', '676', '677', '678', '679', '680']
standard_sensitivity_map = {'BRD-K66175015': ['12', '64', '108', '134', '142', '182', '185', '223', '236', '294'], 'BRD-K19687926': ['142', '185', '294'],
                            'BRD-K51313569':  ['28', '88', '237', '286', '294'], 'BRD-K12343256': ['53', '65', '121', '135', '161', '264', '270', '318', '411', '459'],
                            'BRD-K56343971': ['53' ,'59', '65', '121' ,'127' ,'135' ,'146' ,'160' ,'161' ,'180' ,'216'
                                              ,'285' ,'434' ,'484'],
                            'BRD-K05804044': ['59', '65', '121', '135', '180', '320'], 'BRD-A12230535': ['108', '121', '244', '328', '345', '374', '511'],
                            'BRD-K59369769': ['187', '207', '212', '299', '303', '379', '428'], 'BRD-K29905972': ['52', '57', '88', '110', '166', '339'],
                            'BRD-K09951645': ['53', '59', '65', '102', '207', '316', '484', '512'], 'BRD-K70401845': ['12', '58', '142', '182', '339', '408', '448', '464', '485'],
                            'BRD-K64052750': ['12', '64', '68', '142', '185', '298', '359'], 'BRD-K49328571': ['6', '22', '193', '255', '368'],
                            'BRD-K88560311': ['52', '410'], 'BRD-K51544265': ['3', '57', '149', '166', '374'], 'BRD-K99023089': ['86', '187', '294'],
                            'BRD-K77625799': ['142', '377', '464'], 'BRD-K44227013': ['53', '57', '110', '166', '187', '212', '223', '339', '512']}


def make_sqi_map(sensitivity_map, col_meta, data):

    s_qi_map = {}

    for key in sensitivity_map:
        if key in col_meta['pert_id'].tolist():
            for dose in col_meta[col_meta['pert_id'] == key]['pert_dose'].unique():
                for well in col_meta[(col_meta['pert_id'] == key) & (col_meta['pert_dose'] == dose)]['pert_well'].unique():

                    if type(dose) is float:
                        dose = 'ctl_vehicle'
                    s_value = data[col_meta[(col_meta['pert_id'] == key) & (col_meta['pert_dose'] == dose) & (col_meta['pert_well'] == well)].index]
                    s_value = s_value[~s_value.index.isin(invariants)]

                    if s_value.shape[1] != 1:
                        raise Exception('S Value has more than one column, should correspond to a single well')

                    s_value = s_value.iloc[:,0]

                    qi = [thing for thing in
                          np.where(s_value.dropna().index.astype(str).isin(sensitivity_map[key]))[0].tolist()]


                    s_qi_map[key + '_' + str(dose) + '_' + well] = [s_value.dropna().tolist(), qi]

    return s_qi_map


def wtks(gct, metadata, outfolder, sensitivity_map=standard_sensitivity_map):

    expected_sensitivity_ranks = []

    print 'meta'
    print 'data'

    for rep in metadata['prism_replicate'].unique():

        print rep

        if not os.path.exists(os.path.join(outfolder, rep)):
            os.mkdir(os.path.join(outfolder, rep))
            os.mkdir(os.path.join(outfolder, rep, 'enrichment_ecdf'))
            os.mkdir(os.path.join(outfolder, rep, 'mountain_plots'))

        rep_folder = os.path.join(outfolder, rep)

        ids = metadata[metadata['prism_replicate'] == rep].index

        data = gct.data_df[ids]

        col_meta = metadata.loc[ids]

        marks = {}

        s_qi_map = make_sqi_map(sensitivity_map, data=data, col_meta=col_meta)


        for key in s_qi_map:
            if key.split('_')[0] in col_meta['pert_id'].tolist():
                if len(s_qi_map[key][0]) > 0 and len(s_qi_map[key][1]) > 0:
                    print key
                    #calculate enrichment score
                    sensitivity_score, cumsum = compute_wtcs.compute_wtcs(pd.Series(s_qi_map[key][0]), pd.Series(s_qi_map[key][1]))
                    print sensitivity_score
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
                    print rids

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

                    marks[col_meta[col_meta['pert_id'] == key.split('_')[0]]['pert_iname'][0] + '_' + key.split('_')[1] + '_' + key.split('_')[2]] = mark

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

                    print mark
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