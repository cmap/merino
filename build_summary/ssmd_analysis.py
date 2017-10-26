import sys
import python.broadinstitute_cmap.io.pandasGEXpress.parse as pe
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
from statsmodels.distributions.empirical_distribution import ECDF as ECDF
import glob

invariants = ['661', '662', '663', '664', '665', '666', '667', '668', '669', '670', '671', '672', '673', '674', '675', '676', '677', '678', '679', '680']


def ssmd_matrix(search_pattern, exclude_wells=[], outpath = ''):

    ssmd_report = pd.DataFrame()
    for x in glob.glob(search_pattern):
        scores = get_ssmd(x, exclude_wells=exclude_wells)
        ssmd_report[os.path.basename(x)] = scores

    ssmd_report.to_csv(outpath, sep='\t', index=False)


def get_ssmd(gctoo, metadata, exclude_wells=[]):


    metadata.set_index('cid', inplace=True)

    neg_dex = metadata[metadata['pert_type'] == 'ctl_vehicle'].index.tolist()
    pos_dex = metadata[(metadata['pert_type'] == 'trt_poscon')].index.tolist()

    neg_df = gctoo.data_df[neg_dex]
    pos_df = gctoo.data_df[pos_dex]

    remove = metadata[metadata['pert_well'].isin(exclude_wells)].index

    neg_remove = [x for x in remove if x in neg_df.columns]
    pos_remove = [x for x in remove if x in pos_df.columns]

    neg_df.drop(neg_remove, axis=1, inplace=True)
    pos_df.drop(pos_remove, axis=1, inplace=True)

    print neg_df.columns

    dmso_medians = neg_df[~neg_df.index.isin(invariants)].median(axis=1, skipna=True)
    poscon_medians = pos_df[~neg_df.index.isin(invariants)].median(axis=1, skipna=True)

    dmso_mad_af = neg_df[~neg_df.index.isin(invariants)].mad(axis=1, skipna=True)
    poscon_mad_af = pos_df[~neg_df.index.isin(invariants)].mad(axis=1, skipna=True)

    ssmds = (dmso_medians - poscon_medians) / np.sqrt((dmso_mad_af * dmso_mad_af) + (poscon_mad_af * poscon_mad_af))

    return ssmds


def ssmd_overview(data_path, pert_outfile, pert, threshold=2):

    ssmds = pd.DataFrame()
    for path in norm_paths:
        x = get_ssmd(path)
        ssmds[os.path.basename(path)] = x

    get_names_gct = pe.parse(norm_paths[0])
    names = get_names_gct.row_metadata_df['name'][~get_names_gct.row_metadata_df.index.isin(invariants)]

    all_ssmds = ssmds.unstack().tolist()

    failures = ssmds.loc[ssmds.median(axis=1) < threshold]

    failure_names = names.loc[ssmds.median(axis=1) < threshold]

    failures.set_index(failure_names, inplace=True)

    data = [all_ssmds]
    labels = ['All SSMDs \n n={}'.format(len(all_ssmds))]

    meds = failures.median(axis=1)
    meds.sort(ascending=True, inplace=True)
    failures2 = failures.loc[meds.index]

    iterrator = failures2.iterrows()

    for index, row in iterrator:
        data.append(row.values)
        labels.append(index)

    plt.boxplot(data, labels=labels)
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.xticks(rotation=70)
    plt.ylim(0, 7)
    plt.ylabel('SSMD Values')

    if len(data) > 1:

        plt.title('{} SSMD Distribution (n={}) with {} Failed Cell Lines (n={})'.format(pert, len(all_ssmds), len(meds),
                                                                                    len(data[1])))
    else:
        plt.title('{} SSMD Distribution (n={}) with {} Failed Cell Lines'.format(pert, len(all_ssmds), len(meds)))

    plt.savefig(os.path.join(pert_outfile, '{}_SSMD_overview.png'.format(pert)))

    plt.clf()

    return len(data) - 1