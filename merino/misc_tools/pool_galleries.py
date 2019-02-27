import glob
import cmapPy.pandasGEXpress.parse as pe
import os


for path in glob.glob('/Volumes/cmap_obelix/pod/custom/PCAL/elwork/PCAL_T3A/*ZSPC*/*/*.gct'):
    plate_name = os.path.basename(os.path.dirname(path))
    gct = pe(path)
    for pool in gct.row_metadata_df['pool_id'].unique():

        if not os.path.exists('/Volumes/cmap_obelix/pod/custom/PCAL/elwork/PCAL_T3A/QC/pool_heatmaps/{}'.format(pool)):
            os.mkdir('/Volumes/cmap_obelix/pod/custom/PCAL/elwork/PCAL_T3A/QC/pool_heatmaps/{}'.format(pool))

        pool_df = gct.data_df.loc[gct.row_metadata_df[gct.row_metadata_df['pool_id'] == pool].index]

        pool_df.columns = [x.split(':')[0] + ':' + x.split(':')[1] for x in pool_df.columns]

        map.mk_heatmap(pool_df, 'Heatmap of Median ZSCORE Values for {}'.format(plate_name),
                       '/Volumes/cmap_obelix/pod/custom/PCAL/elwork/PCAL_T3A/QC/pool_heatmaps/{}/{}.png'.format(pool,
                                                                                                                plate_name),
                       lims=[-10, 10], colormap='coolwarm')