import seaborn as sns
import matplotlib.pyplot as plt

def modz_dist(df, col_metadata_df, exclude_wells, outfile):

    #Make histogram of median ZSPC by compound and mark the location of bortezomib
    import numpy as np
    #df = pe.parse('/Users/elemire/Workspace/PCAL_proj_dir/modz_master.gct')

    pos_dex = col_metadata_df[col_metadata_df['pert_type'] == 'trt_poscon'].index.tolist()

    neg_dex = col_metadata_df[col_metadata_df['pert_type'] == 'ctl_vehicle'].index.tolist()

    treatment = df.data_df[df.data_df.columns[~df.data_df.columns.isin(neg_dex + pos_dex)]]

    neg_df = df.data_df[neg_dex]
    pos_df = df.data_df[pos_dex]

    remove = col_metadata_df[col_metadata_df['pert_well'].isin(exclude_wells)].index



    neg_remove = [x for x in remove if x in neg_df.columns]
    pos_remove = [x for x in remove if x in pos_df.columns]

    neg_df.drop(neg_remove, axis=1, inplace=True)
    pos_df.drop(pos_remove, axis=1, inplace=True)


    full_dist = treatment.unstack().dropna()
    pos_dist = pos_df.unstack().dropna()
    neg_dist = neg_df.unstack().dropna()
    full_dist.replace([np.inf, -np.inf], np.nan, inplace=True)
    full_dist.dropna(inplace=True)

    bines = np.linspace(min(full_dist), 4, 50)

    sns.distplot(neg_dist, bines, color='blue', norm_hist=True, label='DMSO')
    sns.distplot(full_dist, bines, color='orange', norm_hist=True, label='Treatment')
    sns.distplot(pos_dist, bines, color='red', norm_hist=True, label='Bortezomib')


    plt.xlim(-10, 10)
    plt.xlabel('MODZSPC')
    plt.ylabel('Frequency')
    plt.title('Distribution of MODZSPC by Compound')
    axes = plt.gca()
    axes.legend(bbox_to_anchor=(0., 0.8, 0.8, .102), loc=3, borderaxespad=0.)

    plt.savefig(outfile)

    plt.clf()