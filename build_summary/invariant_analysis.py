import sys
import cmapPy.pandasGEXpress.parse as pe
import numpy
import os
from scipy import stats
import matplotlib.pyplot as plt
from collections import OrderedDict
import seaborn as sns
import pandas as pd
import numpy as np
import cmapPy.pandasGEXpress.GCToo as GCToo


def spearmanwrapper(sequence):
    # TODO: Make sure sequence is of equal length to sequence passed?
    if len(sequence.dropna()) > 2:
        mono_value = stats.spearmanr(sequence, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        return mono_value[0]
    else:
        return

def extract_invariants(gctoo):

    invariants = range(661,671)

    invs = gctoo.data_df.loc[invariants]

    row_meta = gctoo.row_metadata_df.loc[invariants]

    new_gctoo = GCToo.GCToo(data_df=invs, row_metadata_df=row_meta, col_metadata_df=gctoo.col_metadata_df)

    return new_gctoo



def invariant_table(gct, count_gct):

    #TODO Include median count in the well and median count of the invariants, median of the rest of the analytes
    inv_gct = extract_invariants(gct)
    invs = inv_gct.data_df
    vals = invs.T
    table = pd.concat([invs.median(), np.log2(invs.median()), invs.max() - invs.min(),
                       invs.apply(spearmanwrapper, axis=0)], axis=1)
    table.columns = ['Median', 'Log Median', 'Range', 'Monotonicity']
    final = table.join(vals)

    return final

def invariant_heatmap(gct):
    """
    Make heatmap of invariant medians here
    :return:
    """
    inv_gct = extract_invariants(gct)

def invariant_monotonicity(filepath, metadata_path, outfile):
    x = pe.parse(filepath)
    col_metadata = pd.read_table(metadata_path)

    # Find your controls
    neg_dex = col_metadata[col_metadata['pert_type'] == 'ctl_vehicle'].index.tolist()
    pos_dex = col_metadata[col_metadata['pert_type'] == 'trt_poscon'].index.tolist()

    # Extract them
    neg_invariant_df = x.data_df[neg_dex]
    pos_invariant_df = x.data_df[pos_dex]

    # Calculate monotonicity of invariants for positive and negative controls
    pos_mono_values = pos_invariant_df.apply(spearmanwrapper, axis=0)

    neg_mono_values = neg_invariant_df.apply(spearmanwrapper, axis=0)

    treatment = x.data_df[x.data_df.columns[~x.data_df.columns.isin(neg_dex + pos_dex)]]

    trt_mono_values = treatment.apply(spearmanwrapper, axis=0)

    high_dose_dex = []
    for compound in col_metadata[~x.data_df.columns.isin(neg_dex + pos_dex)]['pert_id'].unique():
        pert_df = col_metadata[
            (col_metadata['pert_id'] == compound) & (col_metadata['pert_type'] == 'trt_cp')]
        high_dose_dex = high_dose_dex + pert_df[(pert_df['pert_dose'] == pert_df['pert_dose'].max())].index.tolist()

    trt_high_dose = treatment[high_dose_dex]

    trt_high_mono_values = trt_high_dose.apply(spearmanwrapper, axis=0)


    mono_data = [trt_mono_values.dropna(), trt_high_mono_values.dropna(), pos_mono_values.dropna(),
                 neg_mono_values.dropna()]
    labels = ['All Treatment, n={}'.format(len(trt_mono_values.dropna())), 'Max Dose, n={}'.format(len(trt_high_mono_values.dropna())),
              'Poscons, n={}'.format(len(pos_mono_values.dropna())), 'DMSO, n={}'.format(len(neg_mono_values.dropna()))]
    plt.boxplot(mono_data, labels=labels)

    plt.xlabel('Compounds Type')
    plt.ylabel('Monotonicity of Invariants')
    plt.title('Monotonicity of Invariants')

    plt.savefig(os.path.join(outfile, 'invariant_mono.png'))
    plt.clf()

    return {'all treatment': trt_mono_values.dropna(), 'high_dose': trt_high_mono_values.dropna(), 'poscons': pos_mono_values.dropna(),
    'DMSO': neg_mono_values.dropna()}


def invariant_curves_plot(df, full_col_meta, outfile):

    master_inv_df = extract_invariants(df)

    for x in full_col_meta['prism_replicate'].unique():

        col_metadata_df = full_col_meta[full_col_meta['prism_replicate'] == x]

        inv_df = master_inv_df.data_df[col_metadata_df.index]

        neg_dex = col_metadata_df[col_metadata_df['pert_type'] == 'ctl_vehicle'].index.tolist()
        pos_dex = col_metadata_df[(col_metadata_df['pert_iname'] == 'Bortezomib') | (col_metadata_df['pert_iname'] == 'MG-132')].index.tolist()


        dmso_invariants = inv_df[neg_dex]
        poscon_invariants = inv_df[pos_dex]


        remainder_index = neg_dex + pos_dex
        remainder_invariants = inv_df.drop(inv_df[remainder_index], axis=1)


        plt.plot(remainder_invariants, 'y', label='Treatment')
        plt.plot(dmso_invariants, 'b', label='DMSO')
        plt.plot(poscon_invariants, 'r', label='poscons')

        plt.xlabel('Invariant ID')
        plt.ylabel('MFI')
        plt.title('PR500 Invariant Curves')

        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = OrderedDict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys())

        if not os.path.exists(os.path.join(outfile, x)):
            os.mkdir(os.path.join(outfile, x))

        plt.savefig(os.path.join(outfile, x, 'invariant_curves.png'))

        plt.clf()

def invariant_range_distributions(gctoo, full_col_meta, outfile):
    # We are making IQR boxplots and distributions for POSCON, DMSO, trt, and trt highest dose

    master_inv_df = extract_invariants(gctoo)

    for y in full_col_meta['prism_replicate'].unique():

        col_metadata_df = full_col_meta[full_col_meta['prism_replicate'] == y]

        inv_df = master_inv_df.data_df[col_metadata_df.index]


        neg_dex = col_metadata_df[col_metadata_df['pert_type'] == 'ctl_vehicle'].index.tolist()
        pos_dex = col_metadata_df[col_metadata_df['pert_type'] == 'trt_poscon'].index.tolist()

        neg_df = inv_df[neg_dex]
        pos_df = inv_df[pos_dex]

        # Return if we are using an all DMSO plate
        if len(pos_df.unstack()) == 0:
            return

        treatment = inv_df[inv_df.columns[~inv_df.columns.isin(neg_dex + pos_dex)]]

        all_trt_IQRs = treatment.quantile(0.9, axis=0) - treatment.quantile(0.1, axis=0)
        all_trt_q9_ratios = treatment.quantile(0.1, axis=0) / treatment.quantile(0.9, axis=0)

        ####################################################################################################
        # Find columns where treatment compounds are used at the maximum dose

        high_dose_dex = []
        for compound in col_metadata_df[~inv_df.columns.isin(neg_dex + pos_dex)]['pert_id'].unique():
            pert_df = col_metadata_df[(col_metadata_df['pert_id'] == compound) & (col_metadata_df['pert_type'] == 'trt_cp')]
            high_dose_dex = high_dose_dex + pert_df[(pert_df['pert_dose'] == pert_df['pert_dose'].max())].index.tolist()

        trt_high_dose = treatment[high_dose_dex]
        all_high_trt_IQRs = trt_high_dose.quantile(0.9, axis=0) - trt_high_dose.quantile(0.1, axis=0)
        all_high_trt_q9_ratios = trt_high_dose.quantile(0.1, axis=0) / trt_high_dose.quantile(0.9,axis=0)

        ########################################################################################################

        trt_lower_dose = treatment.drop(high_dose_dex, axis=1)
        all_lower_trt_IQRs = trt_lower_dose.quantile(0.9, axis=0) - trt_lower_dose.quantile(0.1,axis=0)
        all_lower_trt_q9_ratios = trt_lower_dose.quantile(0.1, axis=0) / trt_lower_dose.quantile(0.9,axis=0)
        #########################################################

        all_pos_invariant_IQRs = pos_df.quantile(0.9, axis=0) - pos_df.quantile(0.1, axis=0)
        all_pos_invariant_q9_ratios = pos_df.quantile(0.1, axis=0) / pos_df.quantile(0.9, axis=0)

        all_neg_invariant_IQRs = neg_df.quantile(0.9, axis=0) - neg_df.quantile(0.1, axis=0)
        all_neg_invariant_q9_ratios = neg_df.quantile(0.1, axis=0) / neg_df.quantile(0.9, axis=0)

        all_IQRs = all_trt_IQRs + all_high_trt_IQRs + all_pos_invariant_IQRs + all_neg_invariant_IQRs
        bines = numpy.linspace(0, max(all_IQRs), 50)
        sns.distplot(all_trt_IQRs.dropna().tolist(), bines, hist=True, kde=False, norm_hist=True, color='green',  label='Trt, n={}'.format(len(all_trt_IQRs)))
        sns.distplot(all_high_trt_IQRs.dropna().tolist(), bines, hist=True, kde=False, norm_hist=True, color='yellow', label='Trt_high_dose, n={}'.format(len(all_high_trt_IQRs)))
        sns.distplot(all_pos_invariant_IQRs.dropna().tolist(), bines, hist=True, kde=False, norm_hist=True, color='red', label='Poscon, n={}'.format(len(all_pos_invariant_IQRs)))
        sns.distplot(all_neg_invariant_IQRs.dropna().tolist(), bines, hist=True, kde=False, norm_hist=True,color='blue', label='DMSO, n={}'.format(len(all_neg_invariant_IQRs)))

        plt.xlabel('Interquartile Range of Invariants')
        plt.ylabel('Frequency')
        plt.title('Distribution of IQRs for Invariants')

        axes = plt.gca()
        axes.legend(bbox_to_anchor=(.7, 0.7, 0.6, .6), loc=3, borderaxespad=0.)


        if not os.path.exists(os.path.join(outfile, y)):
            os.mkdir(os.path.join(outfile, y))

        plt.savefig(os.path.join(outfile, y, 'Invariant_IQR_Distributions.png'))
        plt.clf()

        all_q9 = all_trt_q9_ratios + all_high_trt_q9_ratios + all_pos_invariant_q9_ratios + all_neg_invariant_q9_ratios
        bines = numpy.linspace(0, max(all_q9), 50)
        sns.distplot(all_trt_q9_ratios.dropna().tolist(), bines, hist=True, kde=False, norm_hist=True, color='green', label='Trt, n={}'.format(len(all_trt_q9_ratios)))
        sns.distplot(all_high_trt_q9_ratios.dropna().tolist(), bines, hist=True, kde=False, norm_hist=True,color='yellow', label='Trt_high_dose, n={}'.format(len(all_high_trt_q9_ratios)))
        sns.distplot(all_pos_invariant_q9_ratios.dropna().tolist(), bines, hist=True, kde=False, norm_hist=True, color='red', label='Poscon, n={}'.format(len(all_pos_invariant_q9_ratios)))
        sns.distplot(all_neg_invariant_q9_ratios.dropna().tolist(), bines, hist=True, kde=False, norm_hist=True,color='blue', label='DMSO, n={}'.format(len(all_neg_invariant_q9_ratios)))

        plt.xlabel('Ratio of q1 to q9 of Invariants')
        plt.ylabel('Frequency')
        plt.title('Distribution of q1/q9 for Invariants')

        axes = plt.gca()
        axes.legend(bbox_to_anchor=(0.7, 0.7, 0.6, .6), loc=3, borderaxespad=0.)

        plt.savefig(os.path.join(outfile, y, 'Invariant_Fold_Change_Distributions.png'))
        plt.clf()

        IQR_data = [all_trt_IQRs, all_high_trt_IQRs, all_pos_invariant_IQRs, all_neg_invariant_IQRs]
        labels = ['All Treatment, n={}'.format(len(all_trt_IQRs)),
                  'Max Dose, n={}'.format(len(all_high_trt_IQRs)),
                  'Poscons, n={}'.format(len(all_pos_invariant_IQRs)),
                  'DMSO, n={}'.format(len(all_neg_invariant_IQRs))]
        plt.boxplot(IQR_data, labels=labels)

        plt.xlabel('Compounds Type')
        plt.ylabel('Interquartile Range of Invariants')
        plt.title('Distribution of IQRs for Invariants')

        plt.savefig(os.path.join(outfile, y, 'Invariant_IQR_Boxplots.png'))

        plt.clf()

        q9_data = [all_trt_q9_ratios, all_high_trt_q9_ratios, all_pos_invariant_q9_ratios, all_neg_invariant_q9_ratios]
        labels = ['All Treatment, n={}'.format(len(all_trt_q9_ratios)),
                  'Max Dose, n={}'.format(len(all_high_trt_q9_ratios)),
                  'Poscons, n={}'.format(len(all_pos_invariant_q9_ratios)),
                  'DMSO, n={}'.format(len(all_neg_invariant_q9_ratios))]
        plt.boxplot(q9_data, labels=labels)

        plt.xlabel('Compounds Type')
        plt.ylabel('Ratio of q1 to q9 of Invariants')
        plt.title('Distribution of q1/q9 for Invariants')

        plt.savefig(os.path.join(outfile, y, 'Invariant_Fold_Change_Boxplots.png'))

        plt.clf()

        all_high_low_q9 = all_high_trt_q9_ratios + all_lower_trt_q9_ratios
        bines = numpy.linspace(0, max(all_high_low_q9), 50)
        sns.distplot(all_lower_trt_q9_ratios.dropna(), bines, hist=True, kde=False, norm_hist=True, color='green',
                     label='Trt_lower_doses, n={}'.format(len(all_lower_trt_q9_ratios)))
        sns.distplot(all_high_trt_q9_ratios.dropna(), bines, hist=True, kde=False, norm_hist=True, color='yellow',
                     label='Trt_high_dose, n={}'.format(len(all_high_trt_q9_ratios)))

        plt.xlabel('Ratio of q1 to q9 of Invariants')
        plt.ylabel('Frequency')
        plt.title('Distribution of q1/q9 for Invariants')

        axes = plt.gca()
        axes.legend(bbox_to_anchor=(0.7, 0.7, 0.6, .6), loc=3, borderaxespad=0.)

        plt.savefig(os.path.join(outfile, 'Low_vs_High_Fold_Change_Distributions.png'))
        plt.clf()

        all_high_low_IQRs = all_high_trt_IQRs + all_lower_trt_IQRs
        bines = numpy.linspace(0, max(all_high_low_IQRs), 50)
        sns.distplot(all_lower_trt_IQRs.dropna().tolist(), bines, hist=True, kde=False, norm_hist=True, color='green',
                     label='Trt_lower_doses, n={}'.format(len(all_lower_trt_IQRs)))
        sns.distplot(all_high_trt_IQRs.dropna().tolist(), bines, hist=True, kde=False, norm_hist=True, color='yellow',
                     label='Trt_high_dose, n={}'.format(len(all_high_trt_IQRs)))


        plt.xlabel('Interquartile Range of Invariants')
        plt.ylabel('Frequency')
        plt.title('Distribution of IQRs for Invariants High Dose vs. Low Dose')

        axes = plt.gca()
        axes.legend(bbox_to_anchor=(.7, 0.7, 0.6, .6), loc=3, borderaxespad=0.)

        plt.savefig(os.path.join(outfile, y, 'Low_vs_High_Invariant_IQR_Distributions.png'))
        plt.clf()
        ##############################################################################################
        all_lower_invariants = trt_lower_dose.unstack().dropna()
        all_high_invariants = trt_high_dose.unstack().dropna()

        all_invariants_low_v_high = all_lower_invariants.tolist() + all_high_invariants.tolist()
        bines = numpy.linspace(0, max(all_invariants_low_v_high), 50)

        sns.distplot(all_lower_invariants.dropna(), bines, hist=True, kde=False, norm_hist=True, color='green',
                     label='Trt_lower_doses, n={}'.format(len(all_lower_invariants)))
        sns.distplot(all_high_invariants.dropna(), bines, hist=True, kde=False, norm_hist=True, color='yellow',
                     label='Trt_high_dose, n={}'.format(len(all_high_invariants)))

        plt.xlabel('Invariant MFI Values')
        plt.ylabel('Frequency')
        plt.title('Distribution of All Invariant Values High Dose vs. Low Dose')

        axes = plt.gca()
        axes.legend(bbox_to_anchor=(.7, 0.7, 0.6, .6), loc=3, borderaxespad=0.)

        plt.savefig(os.path.join(outfile, y, 'Low_vs_High_Invariant_Values_Distributions.png'))
        plt.clf()

    return all_lower_invariants, all_high_invariants, all_high_trt_IQRs, all_lower_trt_IQRs, q9_data, IQR_data
