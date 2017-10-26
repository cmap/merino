import cmapPy.pandasGEXpress.parse as pe
import cmapPy.pandasGEXpress.GCToo as GCToo


def calculate_zscore(df, plate_control=False):

    # Calculate level 4 data from level 3

    if plate_control == False:

        neg_dex = df.col_metadata_df[df.col_metadata_df['pert_type'] == 'ctl_vehicle'].index.tolist()
        neg_df = df.data_df[neg_dex]

        sub = df.data_df.subtract(neg_df.median(axis=1), axis='index')

        #sub = df.data_df.subtract(
        #    neg_df.median(axis=1), axis='index')

        #zscore_data = sub.divide(neg_df.mad(axis=1) * 1.4826, axis='index')

        dmso_medians = neg_df.median(axis=1)
        dmso_median_devs = abs(neg_df.subtract(dmso_medians, axis=0))
        dmso_mads = dmso_median_devs.median(axis=1) * 1.4826
        dmso_mads[dmso_mads < .1] = .1

        zscore_data = sub.divide(dmso_mads, axis='index')

        df.col_metadata_df['data_level'] = 'ZSVC'
        df.col_metadata_df['provenance'] = [x + ' | ZSVC' for x in df.col_metadata_df['provenance']]

    elif plate_control == True:

        medians = df.data_df.median(axis=1)

        sub = df.data_df.subtract(medians, axis='index')

        median_devs = abs(df.data_df.subtract(medians, axis=0))
        mads = median_devs.median(axis=1) * 1.4826
        mads[mads < .1] = .1

        zscore_data = sub.divide(mads, axis='index')

        #sub = df.data_df.subtract(
        #    df.data_df.median(axis=1), axis='index')

        #zscore_data = sub.divide(df.data_df.mad(axis=1) * 1.4826, axis='index')

        df.col_metadata_df['data_level'] = 'ZSPC'
        df.col_metadata_df['provenance'] = [x + ' | ZSPC' for x in df.col_metadata_df['provenance']]

    row_metadata_df = df.row_metadata_df

    zscore_data[zscore_data<-10] = -10

    zscore_data[zscore_data > 10] = 10

    zscore_data.sort_index(inplace=True)
    row_metadata_df.sort_index(inplace=True)

    zscore_gctoo = GCToo.GCToo(data_df=zscore_data, row_metadata_df=row_metadata_df, col_metadata_df=df.col_metadata_df)

    #if write is True:
     #   wg.write(Zscore_GCT, outfile + '_ZSCORE.gct')

    return zscore_gctoo

