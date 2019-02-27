import sys
import cmapPy.pandasGEXpress.parse as pe
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import cmapPy.pandasGEXpress.GCToo as GCToo
import glob
import seaborn as sns

from statsmodels.distributions.empirical_distribution import ECDF

invariants = ['661', '662', '663', '664', '665', '666', '667', '668', '669', '670', '671', '672', '673', '674', '675', '676', '677', '678', '679', '680']

def ssmd_matrix(norm_paths, exclude_wells=[], outpath = ''):

    ssmd_report = pd.DataFrame()

    for path in norm_paths:
        print path
        x = pe.parse(path)
        scores = get_ssmd(x, exclude_wells=exclude_wells, unlog=True)
        ssmd_report[os.path.basename(os.path.dirname(path))] = scores
    if len(outpath) > 0:
        ssmd_report.to_csv(outpath, sep='\t')
    return ssmd_report



def get_ssmd(df, exclude_wells = [], unlog=False, pos_field = 'pert_type', pos_val='trt_poscon', neg_val='ctl_vehicle'):


    if unlog==True:
        new_data = np.power(2, df.data_df)
        df = GCToo.GCToo(data_df=new_data, col_metadata_df=df.col_metadata_df, row_metadata_df=df.row_metadata_df)

    neg_dex = df.col_metadata_df[df.col_metadata_df['pert_type'] == neg_val].index.tolist()
    pos_dex = df.col_metadata_df[(df.col_metadata_df[pos_field] == pos_val)].index.tolist()

    neg_df = df.data_df[neg_dex]
    pos_df = df.data_df[pos_dex]


    remove = df.col_metadata_df[df.col_metadata_df['pert_well'].isin(exclude_wells)].index

    neg_remove = [x for x in remove if x in neg_df.columns]
    pos_remove = [x for x in remove if x in pos_df.columns]

    neg_df.drop(neg_remove, axis=1, inplace=True)
    pos_df.drop(pos_remove, axis=1, inplace=True)

    dmso_medians = neg_df[~neg_df.index.isin(invariants)].median(axis=1, skipna=True)
    #print 'dmso median {}'.format(dmso_medians.values[0])
    poscon_medians = pos_df[~neg_df.index.isin(invariants)].median(axis=1, skipna=True)
    #print 'pos median {}'.format(poscon_medians.values[0])


    pos_median_devs = abs(pos_df.subtract(poscon_medians, axis=0))
    poscon_mad = pos_median_devs.median(axis=1) * 1.4826
    #print 'pos mad {}'.format(poscon_mad.values[0])

    dmso_median_devs = abs(neg_df.subtract(dmso_medians, axis=0))
    dmso_mad = dmso_median_devs.median(axis=1) * 1.4826
    #print 'dmso mad {}'.format(dmso_mad.values[0])
    #print 'diff {}'.format(dmso_medians.values[0] - poscon_medians.values[0])
    #print 'mads {}'.format(dmso_mad.values[0] + poscon_mad.values[0])

    #dmso_mad = neg_df[~neg_df.index.isin(invariants)].mad(axis=1, skipna=True)
    #poscon_mad = pos_df[~neg_df.index.isin(invariants)].mad(axis=1, skipna=True)

    ssmds = (dmso_medians - poscon_medians) / np.sqrt((dmso_mad * dmso_mad) + (poscon_mad * poscon_mad))
    #print 'ssmd {}'.format(ssmds.values[0])


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


def ssmd_ecdf(norm_gct, median_gct, title, outfile):

    norm_ssmd = pd.Series()
    median_ssmd = pd.Series()


    for rep in norm_gct.col_metadata_df['prism_replicate'].unique():
        norm_temp = norm_gct.data_df[norm_gct.col_metadata_df[norm_gct.col_metadata_df['prism_replicate'] == rep].index]
        norm_temp_col = norm_gct.col_metadata_df.loc[norm_gct.col_metadata_df['prism_replicate'] == rep]
        temp_norm_ssmd = get_ssmd(GCToo.GCToo(data_df=norm_temp, col_metadata_df=norm_temp_col, row_metadata_df=norm_gct.row_metadata_df),
                             unlog=True)
        med_temp = median_gct.data_df[median_gct.col_metadata_df[median_gct.col_metadata_df['prism_replicate'] == rep].index]
        med_temp_col = median_gct.col_metadata_df.loc[median_gct.col_metadata_df['prism_replicate'] == rep]
        temp_median_ssmd = get_ssmd(GCToo.GCToo(data_df=med_temp, col_metadata_df=med_temp_col, row_metadata_df=median_gct.row_metadata_df))
        norm_ssmd = norm_ssmd.append(temp_norm_ssmd)
        median_ssmd = median_ssmd.append(temp_median_ssmd)

    norm_ecdf = ECDF(norm_ssmd)
    med_ecdf = ECDF(median_ssmd)

    plt.plot(norm_ecdf.x, norm_ecdf.y, label='NORM')
    plt.plot(med_ecdf.x, med_ecdf.y, label='MFI')
    plt.xlim(-1,10)
    plt.xlabel('SSMD Values')
    plt.title(title)
    axes = plt.gca()
    axes.legend(bbox_to_anchor=(.615, 0.81, 0.8, .6), loc=3, borderaxespad=0.)
    plt.savefig(os.path.join(outfile, 'SSMD_ECDF.png'))


def norm_v_mfi_ssmd(norm_gct, median_gct, outfile):

    norm_ssmd = get_ssmd(norm_gct, unlog=True)
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
    plt.clf()

def ssmd_by_pool(ssmd, cell_info, outfile):
    ssmd_count = pd.DataFrame()
    ssmd_median = pd.DataFrame()
    cell_info.index = [str(x) for x in cell_info.index]
    cell_info = cell_info.loc[ssmd.index.values]
    for pool in cell_info['pool_id'].unique():
        print pool
        pool_dex = cell_info[cell_info['pool_id'] == pool].index
        pool_ssmd = ssmd.loc[pool_dex]
        ssmd_count[pool] = pool_ssmd[pool_ssmd < 2].count() / float(pool_ssmd.shape[0])

    for pool in cell_info['pool_id'].unique():
        print pool
        pool_dex = cell_info[cell_info['pool_id'] == pool].index
        pool_ssmd = ssmd.loc[pool_dex]
        ssmd_median[pool] = pool_ssmd.median()

    ssmd_count.to_csv(os.path.join(outfile, 'ssmd_failure%_by_pool.txt'), sep='\t')
    ssmd_median.to_csv(os.path.join(outfile, 'ssmd_median_by_pool.txt'), sep='\t')
