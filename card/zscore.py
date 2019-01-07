import cmapPy.pandasGEXpress.parse as pe
import cmapPy.pandasGEXpress.GCToo as GCToo


def zscore(mat, ctrl_mat=None):

    if ctrl_mat is not None:
        medians = ctrl_mat.median(axis=1)
        median_devs = abs(ctrl_mat.subtract(medians, axis=0))

    else:
        medians = mat.median(axis=1)
        median_devs = abs(mat.subtract(medians, axis=0))

    sub = mat.subtract(medians, axis='index')
    mads = median_devs.median(axis=1)
    mads[mads < .1] = .1
    zscore_data = sub.divide(mads * 1.4826, axis='index')

    return zscore_data



def calculate_zscore(df, plate_control=False):

    # Calculate level 4 data from level 3

    if plate_control == False:
        neg_dex = df.col_metadata_df[df.col_metadata_df['pert_type'] == 'ctl_vehicle'].index.tolist()
        neg_df = df.data_df[neg_dex]
        zscore_data = zscore(df.data_df, neg_df)
        df.col_metadata_df['data_level'] = 'ZSVC'
        df.col_metadata_df['provenance'] = [x + ' | ZSVC' for x in df.col_metadata_df['provenance']]

    elif plate_control == True:
        zscore_data = zscore(df.data_df)
        df.col_metadata_df['data_level'] = 'ZSPC'
        df.col_metadata_df['provenance'] = [x + ' | ZSPC' for x in df.col_metadata_df['provenance']]

    row_metadata_df = df.row_metadata_df

    zscore_data[zscore_data < -10] = -10

    zscore_data[zscore_data > 10] = 10

    zscore_data.sort_index(inplace=True)
    row_metadata_df.sort_index(inplace=True)
    zscore_gctoo = GCToo.GCToo(data_df=zscore_data, row_metadata_df=row_metadata_df, col_metadata_df=df.col_metadata_df)

    return zscore_gctoo

