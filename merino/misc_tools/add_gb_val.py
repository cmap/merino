import glob
import cmapPy.pandasGEXpress.parse as pe
import cmapPy.pandasGEXpress.write_gct as wgct


for path in glob.glob('/Volumes/cmap_obelix/pod/custom/PCAL/elwork/PCAL_T3A/assemble/*/*.gct'):
    print path
    ztest = pe(path)
    ztest.col_metadata_df['gb_id'] = ztest.col_metadata_df['pert_id']
    ztest.col_metadata_df.loc[ztest.col_metadata_df['pert_iname'] == 'Bortezomib', 'pert_type'] = 'trt_poscon'
    ztest.col_metadata_df.loc[ztest.col_metadata_df['pert_type'] == 'trt_poscon', 'gb_id'] = ztest.col_metadata_df.loc[ztest.col_metadata_df['pert_type'] == 'trt_poscon', 'gb_id'] + [':' + str(x) for x in range(1, len(ztest.col_metadata_df.loc[ztest.col_metadata_df['pert_type'] == 'trt_poscon', 'gb_id']) + 1)]
    ztest.col_metadata_df.loc[ztest.col_metadata_df['pert_type'] == 'ctl_vehicle', 'gb_id'] = ztest.col_metadata_df.loc[ztest.col_metadata_df['pert_type'] == 'ctl_vehicle', 'gb_id'] + [':' + str(x) for x in range(1, len(ztest.col_metadata_df.loc[ztest.col_metadata_df['pert_type'] == 'ctl_vehicle', 'gb_id']) + 1)]
    wgct.write(ztest, path)