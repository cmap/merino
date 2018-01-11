import sys
import cmapPy.pandasGEXpress.parse as pe
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import cmapPy.pandasGEXpress.GCToo as GCToo
import glob
import seaborn as sns

from statsmodels.distributions.empirical_distribution import ECDF

invariants = ['661', '662', '663', '664', '665', '666', '667', '668', '669', '670', '671', '672', '673', '674', '675', '676', '677', '678', '679', '680']

def ssmd_matrix(norm_paths, exclude_wells=[], log=True, outpath = ''):

    ssmd_report = pd.DataFrame()
    for path in norm_paths:
        print path
        x = pe(path)
        if log==True:
            unlog = np.power(2, x.data_df)
            x = GCToo.GCToo(data_df=unlog, row_metadata_df=x.row_metadata_df, col_metadata_df=x.col_metadata_df)
        scores = get_ssmd(x, exclude_wells=exclude_wells)
        ssmd_report[os.path.basename(os.path.dirname(path))] = scores

    return ssmd_report
    #ssmd_report.to_csv(outpath, sep='\t', index=False)


def get_ssmd(df, exclude_wells = []):

    neg_dex = df.col_metadata_df[df.col_metadata_df['pert_type'] == 'ctl_vehicle'].index.tolist()
    pos_dex = df.col_metadata_df[(df.col_metadata_df['pert_type'] == 'trt_poscon')].index.tolist()

    neg_df = df.data_df[neg_dex]
    pos_df = df.data_df[pos_dex]

    remove = df.col_metadata_df[df.col_metadata_df['pert_well'].isin(exclude_wells)].index

    neg_remove = [x for x in remove if x in neg_df.columns]
    pos_remove = [x for x in remove if x in pos_df.columns]

    neg_df.drop(neg_remove, axis=1, inplace=True)
    pos_df.drop(pos_remove, axis=1, inplace=True)


    dmso_medians = neg_df[~neg_df.index.isin(invariants)].median(axis=1, skipna=True)
    poscon_medians = pos_df[~neg_df.index.isin(invariants)].median(axis=1, skipna=True)

    dmso_mad = neg_df[~neg_df.index.isin(invariants)].mad(axis=1, skipna=True)
    poscon_mad = pos_df[~neg_df.index.isin(invariants)].mad(axis=1, skipna=True)

    ssmds = (dmso_medians - poscon_medians) / np.sqrt((dmso_mad * dmso_mad) + (poscon_mad * poscon_mad))

    return ssmds

def get_failed_cell_lines(filepaths, exclude_wells=[]):
    failed_lines_2_map = {}
    failed_lines_1_5_map = {}

    for x in filepaths:
        gct = pe(x)
        ssmds = get_ssmd(gct, exclude_wells=exclude_wells)
        failed_lines2 = gct.row_metadata_df.loc[ssmds[ssmds < 2].index]['name']
        failed_lines1_5 = gct.row_metadata_df.loc[ssmds[ssmds < 1.5].index]['name']
        failed_lines_2_map[os.path.basename(x)] = failed_lines2
        failed_lines_1_5_map[os.path.basename(x)] = failed_lines1_5

    return failed_lines_2_map, failed_lines_1_5_map


def norm_v_mfi_ssmd(norm_gct, median_gct, outfile):

    norm_ssmd = get_ssmd(norm_gct)
    median_ssmd = get_ssmd(median_gct)


    data = [norm_ssmd.dropna(), median_ssmd.dropna()]

    if len(norm_ssmd.dropna()) > 0:

        labels = ['Norm SSMDs, n={}'.format(len(norm_ssmd)), 'MFI SSMDs, n={}'.format(len(norm_ssmd))]
        plt.boxplot(data, labels=labels)

        plt.ylabel('SSMD Values')
        plt.title('Comparison of SSMDs Between Norm and MFI Data')

        plt.savefig(os.path.join(outfile, 'NORMvMFI_SSMD_Boxplot.png'))

        plt.clf()


        bines = np.linspace(min(norm_ssmd.dropna().append(median_ssmd.dropna())), max(norm_ssmd.dropna().append(median_ssmd.dropna())), 50)
        n, bins, patches = plt.hist(norm_ssmd.dropna(), bines, facecolor='blue', alpha=.4, label='NORM, n={}'.format(len(norm_ssmd)))
        n1, bins1, patches1 = plt.hist(median_ssmd.dropna(), bines, facecolor='orange', alpha=.4,
                                   label='MFI, n={}'.format(len(median_ssmd)))

        plt.xlabel('SSMD Values')
        plt.ylabel('Frequency')
        plt.title('Distribution of SSMDs for Norm and MFI Data')

        axes = plt.gca()
        axes.legend(bbox_to_anchor=(.615, 0.81, 0.8, .6), loc=3, borderaxespad=0.)

        plt.savefig(os.path.join(outfile, 'SSMD_Distributions.png'))
        plt.clf()


def ssmd_overview(norm_paths, pert_outfile, pert, threshold=2, exclude_wells = []):

    ssmds = pd.DataFrame()
    for path in norm_paths:
        gct = pe(path)
        x = get_ssmd(gct, exclude_wells=exclude_wells)
        ssmds[os.path.basename(path)] = x
    get_names_gct = pe(norm_paths[0])
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
    plt.ylim(0, 10)
    plt.ylabel('SSMD Values')

    if len(data) > 1:

        plt.title('{} SSMD Distribution (n={}) with {} Failed Cell Lines (n={})'.format(pert, len(all_ssmds), len(meds),
                                                                                    len(data[1])))
    else:
        plt.title('{} SSMD Distribution (n={}) with {} Failed Cell Lines'.format(pert, len(all_ssmds), len(meds)))

    plt.savefig(os.path.join(pert_outfile, '{}_SSMD_overview.png'.format(pert)))

    plt.clf()

    return len(data) - 1

def failures_ecdf(filepaths, outfile, exclude_wells = []):

    fail_2, fail_1_5 = get_failed_cell_lines(filepaths, exclude_wells=exclude_wells)

    cell_line_map = {}
    for jey in fail_2:
        for cell in fail_2[jey]:
            if cell not in cell_line_map.keys():
                cell_line_map[cell] = []
                cell_line_map[cell].append(jey)
            else:
                cell_line_map[cell].append(jey)

    cell_line_map15 = {}
    for jey in fail_1_5:
        for cell in fail_1_5[jey]:
            if cell not in cell_line_map15.keys():
                cell_line_map15[cell] = []
                cell_line_map15[cell].append(jey)
            else:
                cell_line_map15[cell].append(jey)


    cell_failure_numbers = []
    for key in cell_line_map:
        cell_failure_numbers.append(pd.Series(cell_line_map[key]).count())

    cell_failure_numbers15 = []
    for key in cell_line_map15:
        cell_failure_numbers15.append(pd.Series(cell_line_map15[key]).count())

    cell_ecdf2 = ECDF(cell_failure_numbers)

    cell_ecdf15 = ECDF(cell_failure_numbers15)

    plt.plot(cell_ecdf2.x, cell_ecdf2.y, label='Failures Per Cell Line (Threshold 2)')
    plt.plot(cell_ecdf15.x, cell_ecdf15.y, label='Failures Per Cell Line (Threshold 1.5)')
    plt.xlabel('Number Of Failures Across Cohort By Cell Line')
    plt.ylabel('Portion of Cell Lines')
    plt.xlim(0,70)
    plt.title('ECDF of Cell Line SSMD Failures Across Cohort For PCAL')

    axes = plt.gca()
    axes.set_ylim([0, 1])
    axes.legend(bbox_to_anchor=(0.55, 0.5, 0.8, .102), loc=3, borderaxespad=0.)
    plt.savefig(os.path.join(outfile, 'CellSSMDFailureECDF.png'))
    plt.show()
    plt.clf()

    failure_numbers1_5 = []
    import seaborn as sns

    for key in fail_1_5:
        failure_numbers1_5.append(fail_1_5[key].count())

    ecdf15 = ECDF(failure_numbers1_5)

    failure_numbers2 = []
    for key in fail_2:
        failure_numbers2.append(fail_2[key].count())

    ecdf2 = ECDF(failure_numbers2)

    plt.plot(ecdf2.x, ecdf2.y, label='Failures Per Plate(Threshold 2)')
    plt.plot(ecdf15.x, ecdf15.y, label='Failures Per Plate (Threshold 1.5)')
    plt.xlabel('Cell Lines Failed Per Plate')
    plt.xlim(0,300)
    plt.ylabel('Portion of Plates')
    plt.title('ECDF of Cell Line SSMD Failures by Plate for PCAL')

    axes = plt.gca()
    axes.set_ylim([0, 1])
    axes.legend(bbox_to_anchor=(0.55, 0.5, 0.8, .102), loc=3, borderaxespad=0.)
    plt.savefig(os.path.join(outfile, 'PlateSSMDFailureECDF.png'))
    plt.show()
    plt.clf()