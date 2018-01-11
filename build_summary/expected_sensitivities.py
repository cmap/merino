import merino.compute_wtcs as compute_wtcs
import merino.plot_enrichment_score as plot_enrichment_score
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

expanded_sense = {'BRD-A12230535': ['108', '345', '328', '121', '244', '374', '511'],
 'BRD-A60245366': ['240', '414', '501', '328', '309', '33'],'BRD-A62182663': ['40','156', '296', '144', '410', '3', '262', '298', '88', '26', '24', '77', '29','149','458'],
 'BRD-K05104363': ['160','127','161','31','121','96', '411', '135','318','68','53','100','14','284','168','216','233','392','473','65','529','153','101'],
 'BRD-K05804044': ['59', '320', '121', '76', '180', '135'],
 'BRD-K09951645': ['53', '65', '59', '102', '207', '484', '512', '316', '460', '160', '127', '379', '161', '180', '291', '121', '96', '81', '135', '255',
  '216', '309', '478', '534', '473', '459'],'BRD-K12343256': [],'BRD-K12787259': ['85','15','293','151','234','157','77','349','353','362'],
 'BRD-K14991967': ['296', '240', '88', '339', '150'],'BRD-K16478699': ['127','201','121','135','61','487','53','407','216','473','484','434','285','160','161','59','65'],
 'BRD-K19687926': ['185', '142', '294', '550', '31', '268', '301', '12', '414', '199', '553', '320', '395', '20', '189', '33', '368', '209'],
 'BRD-K20696416': ['90','491','218','135','171','438','328','339','362','353','256','173','244'],
 'BRD-K23984367': ['374', '162', '118', '107', '406', '352'], 'BRD-K24132293': ['317', '88', '53', '309'],'BRD-K24496482': ['121', '59', '135', '53', '216', '65'],
 'BRD-K29395450': ['52', '238', '433', '213'],
 'BRD-K29905972': ['166', '57', '339', '52', '88', '110', '190'],
 'BRD-K34581968': ['127', '369', '145', '359', '399', '39', '65', '256'],
 'BRD-K36740062': ['7','127','85', '169', '531', '13', '269', '373', '51', '133', '330','52','61','183','433','339','277','353'],
 'BRD-K41312087': ['161', '52', '57', '400'],
 'BRD-K41895714': ['127','237','232','284','29','439','163','220','65'],
 'BRD-K42828737': [],
 'BRD-K44227013': ['57','110','512','53','223','339','212','187','166'],
 'BRD-K48112880': ['95', '240', '100', '281', '477', '395', '216', '484'],
 'BRD-K48488978': ['52', '433', '157'],
 'BRD-K49328571': ['255', '6', '22', '368', '193'],
 'BRD-K49669041': ['374', '85', '52'],
 'BRD-K50140147': ['90','201','491','12','380','218','272','87','312','171','438','52','454','406','339','352','362','353','190','277','256','101','285'],
 'BRD-K51313569': ['28', '286', '294', '237', '88'],
 'BRD-K51544265': ['57', '166', '3', '374', '149', '371', '355'],
 'BRD-K51575138': ['238', '352', '274'],
 'BRD-K54606188': ['258', '411', '242', '102'],
 'BRD-K56343971': ['285', '484', '121', '53', '180', '434', '59', '160','65','146','216','127','206','161','135'],
 'BRD-K57282030': ['476','252','85','174', '335', '95', '42', '68', '52', '57', '433', '157', '264', '138', '213'],
 'BRD-K58772419': ['156','161','6', '99', '262', '13', '289', '411', '358', '293', '120', '238', '359', '155', '363', '539', '118', '105',
  '307','393','190','372','480','215','213'],
 'BRD-K59369769': ['379','412','212','207','187','299','303','428','428'],
 'BRD-K61860715': ['410','301','88','52','223','303','413','299','433'],
 'BRD-K64052750': ['12', '359', '298', '142', '68', '64', '185'],
 'BRD-K64800655': ['85','18','269','293','22','76','51','52','250','234','29','264'],
 'BRD-K66175015': ['108','142','294','134','223','12','182', '185', '236', '64', '316', '31', '335', '311', '346', '68', '553', '239', '492', '413',
  '20', '395', '359', '172', '283'],'BRD-K66538826': ['174','228','128','259','373','76','84','295','449','149'],
 'BRD-K68065987': ['394','86','474','161','99','108','240', '289', '253', '293', '236', '88', '68', '307', '171', '186', '57','25','303','295','97','157',
  '518', '148', '208', '213'],'BRD-K70401845': ['485', '12', '464', '448', '142', '408', '182', '58', '339', '90', '31', '414', '16', '185',
  '413', '350', '113', '372', '33', '42', '171', '320', '393', '209'],'BRD-K72636697': ['63', '476', '252', '40','550', '530', '144','95','13','68',
  '330', '52', '449', '234', '433', '246', '172', '214', '264', '309', '277', '242', '39', '213', '104'],
 'BRD-K74514084': ['126', '374', '110', '539', '166', '115', '26', '302'],
 'BRD-K76894938': ['476', '336', '114', '52', '255'],
 'BRD-K77625799': ['142','464','377', '90','31','182','12','162','257','448',  '510', '185', '413', '350'],
 'BRD-K78431006': ['491', '371', '84', '438', '355', '300'],
 'BRD-K79404599': ['476','170', '293', '227', '355', '250', '449', '234', '249', '189'],
 'BRD-K81418486': ['86', '161', '31', '96', '240', '135', '397', '68', '100', '509', '365', '185', '137', '518', '148', '39', '473', '529', '150', '101',
  '283'],'BRD-K81528515': ['491', '380', '482', '413'],'BRD-K82746043': ['122','86', '180', '291', '371', '289', '318', '457', '356', '426', '409', '260', '365',
  '238', '185', '433', '518', '208', '353'],'BRD-K85402309': ['126','40','374','87','482','312','171','339','190','378','173'],
 'BRD-K86118762': ['127','369',  '296',  '170','108','272','373','356','347','537','501','232','234','107','399','328','246','555','39','138','534','353','33'],
 'BRD-K88560311': [],
 'BRD-K90497590': ['252','267','317','262','15','346','271','88','254','77', '138', '153', '275', '155'],'BRD-K92991072': ['252', '127', '121', '108', '178', '239', '483'],
 'BRD-K93725829': ['127','85','161','201','411', '135', '257', '269', '270', '61', '171', '359', '77', '29', '233', '264','274', '220', '65'],
 'BRD-K93747373': ['369', '410', '88', '24', '303', '256'],
 'BRD-K93918653': ['374', '530', '166', '373', '186', '57'],
 'BRD-K95435023': ['371', '257', '355', '33'],
 'BRD-K97764662': ['32', '114', '289', '356'],
 'BRD-K99023089': ['294', '86', '187'],
 'BRD-K99545815': ['114', '240', '302', '208', '277'],
 'BRD-K99749624': ['374', '110', '530', '166', '57', '256'],
 'BRD-K99964838': ['127','316','161','31','301','81','129','411','253','12','76','51','414','58','100','340','413','137','359','208','350','207','153','159'],
 'BRD-U07805514': ['22','182','12','380','162','253','448','510','16','20','185','350','352','122','320','275']}


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


def wtks(gct, metadata, outfolder, sensitivity_map=expanded_sense):

    expected_sensitivity_ranks = []

    print 'meta'
    print 'data'


    for rep in metadata['weave_prefix'].unique():

        print rep

        if not os.path.exists(os.path.join(outfolder, rep)):
            os.mkdir(os.path.join(outfolder, rep))
            os.mkdir(os.path.join(outfolder, rep, 'enrichment_ecdf'))
            os.mkdir(os.path.join(outfolder, rep, 'mountain_plots'))

        rep_folder = os.path.join(outfolder, rep)

        ids = metadata[metadata['weave_prefix'] == rep].index

        data = gct.data_df[ids]

        data.index = [str(x) for x in data.index]

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


