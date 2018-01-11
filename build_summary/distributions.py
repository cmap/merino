import sys
sys.path.append('/Users/elemire/Workspace/l1ktools')
import python.broadinstitute_cmap.io.pandasGEXpress.parse as pe
import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns

invariants = ['661', '662', '663', '664', '665', '666', '667', '668', '669', '670', '671', '672', '673', '674', '675', '676', '677', '678', '679', '680']

def distributions(norm_gct,median_gct,count_gct, zscore_gct,viability_gct, col_metadata_df, outfile):

    ###################################################################################
    # Making distribution of normalised MFI values

    # The replacement of inf with NaN was only necessary because of issues with DP11/12 Data. May just remove this in the future
    norm_df = norm_gct.data_df

    norm_data = norm_df.unstack()
    norm_data.replace([np.inf, -np.inf], np.nan, inplace=True)

    norm_data.dropna(inplace=True)
    norm_bins = np.linspace(-10, max(norm_data), 100)

    # Matplotlib Histogram
    n, bins, patches = plt.hist(norm_data, norm_bins, facecolor='green', alpha=1)


    plt.xlabel('Normalized Data, n={}'.format(len(norm_data)))
    plt.ylabel('Frequency')
    plt.title('Distribution of Norm Data')
    plt.grid(True)
    axes = plt.gca()
    axes.set_xlim(xmin=-10, xmax=max(norm_data))

    plt.savefig(os.path.join(outfile, 'norm_dist.png'))

    plt.clf()

    ###########################################################################################
    # Making distribution of normalized DMSO values

    neg_dex = col_metadata_df[col_metadata_df['pert_type'] == 'ctl_vehicle'].index.tolist()
    pos_dex = col_metadata_df[(col_metadata_df['pert_iname'] == 'Bortezomib') | (col_metadata_df['pert_iname'] == 'MG-132')].index.tolist()

    neg_dex = [y for y in neg_dex if y in norm_df.columns]
    neg_df = norm_df[neg_dex]
    neg_data = neg_df.unstack()
    neg_data.dropna(inplace=True)

    neg_bins = np.linspace(min(neg_data), max(neg_data), 100)
    n, bins, patches = plt.hist(neg_data, neg_bins, facecolor='blue', alpha=1)

    plt.xlabel('Normalized DMSO Values, n={}'.format(len(neg_df.unstack().dropna())))
    plt.ylabel('Frequency')
    plt.title('Distribution of Norm DMSO Data')
    plt.grid(True)

    plt.savefig(os.path.join(outfile, 'dmso_dist.png'))

    plt.clf()

    #############################################################################################
    # Making distribution of normalized POSCON values

    pos_dex = [y for y in pos_dex if y in norm_df.columns]
    pos_df = norm_df[pos_dex]
    pos_data = pos_df.unstack()
    pos_data.dropna(inplace=True)

    if len(pos_data) > 0:
        pos_bins = np.linspace(min(pos_data), max(pos_data), 100)
        n, bins, patches = plt.hist(pos_data, pos_bins, facecolor='red', alpha=1)

        plt.xlabel('Normalized Poscon Values, n={}'.format(len(pos_data)))
        plt.ylabel('Frequency')
        plt.title('Distribution of Norm Poscon Data')
        plt.grid(True)
        axes = plt.gca()
        axes.set_xlim(xmin=min(pos_data), xmax=max(pos_data))

        plt.savefig(os.path.join(outfile, 'poscon_dist.png'))

        plt.clf()

    ##################################################################################################
    # Plotting normalized POSCON and DMSO distributions together

    if len(pos_data) > 0:
        control_bins = np.linspace(min(pos_data), max(neg_data), 100)

    else:
        control_bins = np.linspace(min(neg_data), max(neg_data), 100)
    sns.distplot(neg_data, control_bins, color='blue', hist=True, kde=False, norm_hist=True, label='DMSO, n={}'.format(len(neg_df.unstack().dropna())))
    sns.distplot(pos_data, control_bins, color='red', hist=True, kde=False, norm_hist=True, label='POSCON, n={}'.format(len(pos_df.unstack().dropna())))

    if len(pos_data) > 0:
        plt.xlim((min(pos_data), max(neg_data)))

    else:
        plt.xlim((min(neg_data), max(neg_data)))

    plt.xlabel('Normalized Control Values')
    plt.ylabel('Frequency')
    plt.title('Distribution of Normalized Control Data')
    plt.grid(True)
    axes = plt.gca()
    axes.legend(bbox_to_anchor=(0.57, 0.85, 0.8, .102), loc=3, borderaxespad=0.)

    plt.savefig(os.path.join(outfile, 'norm_control_dist.png'))

    plt.clf()

    #####################################################################################################
    # Making distribution of RAW MFI values.

    median_df = median_gct.data_df[~median_gct.data_df.index.isin(invariants)]
    median_data = median_df.unstack()
    median_data.dropna(inplace=True)
    mfi_bins = np.linspace(0, max(median_data), 100)
    n, bins, patches = plt.hist(median_data, mfi_bins, facecolor='blue', alpha=1)

    plt.xlabel('MFI Data, n={}'.format(len(median_data)))
    plt.ylabel('Frequency')
    plt.title('Distribution of Median Data')
    plt.grid(True)
    axes = plt.gca()
    axes.set_xlim(xmin=0, xmax=max(median_data))

    plt.savefig(os.path.join(outfile, 'mfi_dist.png'))

    plt.clf()
    #################################################################################################
    # Plotting MFI control density together

    mfi_neg_df = median_df[neg_dex]
    mfi_neg_data = mfi_neg_df.unstack()
    mfi_neg_data.dropna(inplace=True)

    mfi_pos_df = median_df[pos_dex]
    mfi_pos_data = mfi_pos_df.unstack()
    mfi_pos_data.dropna(inplace=True)


    mfi_control_bins = np.linspace(0, max(mfi_neg_data), 100)
    sns.distplot(mfi_neg_data, mfi_control_bins, hist=True, kde=False, norm_hist=True, color='blue', label='DMSO, n={}'.format(len(mfi_neg_df.unstack().dropna())))
    sns.distplot(mfi_pos_data, mfi_control_bins, hist=True, kde=False, norm_hist=True, color='red',  label='POSCON, n={}'.format(len(mfi_pos_df.unstack().dropna())))

    plt.xlim((0, max(mfi_neg_data)))
    plt.xlabel('MFI Control Values')
    plt.ylabel('Frequency')
    plt.title('Distribution of Control Data')

    axes = plt.gca()
    axes.legend(bbox_to_anchor=(0.57, 0.85, 0.8, .102), loc=3, borderaxespad=0.)

    plt.savefig(os.path.join(outfile, 'mfi_control_dist.png'))
    plt.clf()

    ###############################################################################################
    # Making distribution of ALL count values
    count_df = count_gct.data_df
    count_data = count_df.unstack()
    count_data.dropna(inplace=True)
    ber = np.linspace(0, 150, 150)

    n, bins, patches = plt.hist(count_data, ber, facecolor='yellow', alpha=1)

    plt.xlabel('Bead Count Data, n={}'.format(len(count_data)))
    plt.ylabel('Frequency')
    plt.title('Distribution of Bead Count Data')
    plt.xlim(0,150)

    plt.savefig(os.path.join(outfile, 'bead_count_dist.png'))

    plt.clf()

    ####################################################################################################
    # Making Dist of Median Count Values per well
    median_well_count = count_df.median()
    median_well_count.dropna(inplace=True)
    bins = np.linspace(0,100,50)

    n, b, patches = plt.hist(median_well_count, bins, facecolor='orange', alpha=1)

    plt.xlabel('Median Bead Count by Well, n={}'.format(len(median_well_count)))
    plt.ylabel('Frequency')
    plt.title('Distribution of Median Count By Well')
    plt.xlim(0,150)

    plt.savefig(os.path.join(outfile, 'median_count_dist.png'))

    plt.clf()

    #############################################################################


    viability_df = viability_gct.data_df[~viability_gct.data_df.index.isin(invariants)]
    viability_data = viability_df.unstack()
    viability_data.replace([np.inf, -np.inf], np.nan, inplace=True)
    viability_data.dropna(inplace=True)
    vib_bins = np.linspace(min(viability_data.dropna()), 4, 100)
    n, bins, patches = plt.hist(viability_data, vib_bins, facecolor='pink', alpha=1)

    plt.xlabel('Viability Data, n={}'.format(len(viability_data)))
    plt.ylabel('Frequency')
    plt.title('Distribution of Viability Data')
    plt.grid(True)
    axes = plt.gca()
    axes.set_xlim(xmin=min(viability_data.dropna()), xmax=max(viability_data.dropna()))

    plt.savefig(os.path.join(outfile, 'viability_dist.png'))

    plt.clf()

    #############################################################################

    zscore_df = zscore_gct.data_df[~zscore_gct.data_df.index.isin(invariants)]
    zscore_data = zscore_df.unstack()
    zscore_data.replace([np.inf, -np.inf], np.nan, inplace=True)
    zscore_data.dropna(inplace=True)
    bins = np.linspace(-10, 10, 100)

    n, bins, patches = plt.hist(zscore_data, bins, facecolor='brown', alpha=1)

    plt.xlabel('ZSCORE Data, n={}, median={}'.format(len(zscore_data), zscore_data.median()))
    plt.ylabel('Frequency')
    plt.title('Distribution of ZSPC Data')
    plt.grid(True)
    axes = plt.gca()

    axes.set_xlim(xmin=-10, xmax=10)

    plt.savefig(os.path.join(outfile, 'zscore_dist.png'))

    plt.clf()


    ##############################################################################

    #modzscore_df = modzscore_gct.data_df[~modzscore_gct.data_df.index.isin(invariants)]
    #modzscore_data = modzscore_df.unstack()
    #modzscore_data.replace([np.inf, -np.inf], np.nan, inplace=True)
    #modzscore_data.dropna(inplace=True)

    #bins = np.linspace(-10, 10, 50)

    #n, bins, patches = plt.hist(modzscore_data, bins, facecolor='brown', alpha=1)

    #plt.xlabel('MODZSCORE Data, n={}, median={}'.format(len(zscore_data), zscore_data.median()))
    #plt.ylabel('Frequency')
    #plt.title('Distribution of MODZSPC Data')
    #plt.grid(True)
    #axes = plt.gca()

    #axes.set_xlim(xmin=-10, xmax=10)

    #plt.savefig(os.path.join(outfile, 'modz_dist.png'))

    #plt.clf()



