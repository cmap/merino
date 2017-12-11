import cmapPy.pandasGEXpress.GCToo as GCToo

def calculate_viability(df, plate_control=False):

    # Calculation of level 3 to level 4 data
    # TODO make this use a metadata file


    neg_dex = df.col_metadata_df[df.col_metadata_df['pert_type'] == 'ctl_vehicle'].index.tolist()
    neg_df = df.data_df[neg_dex]
    neg_df.median(axis=1)

    if plate_control is False:
        DNorm_Data = df.data_df.divide(neg_df.median(axis=1), axis='index')

    else:
        DNorm_Data = df.data_df.divide(df.data_df.median(axis=1), axis='index')

    row_metadata_df = df.row_metadata_df

    DNorm_Data.sort_index(inplace=True)

    row_metadata_df.sort_index(inplace=True)

    df.col_metadata_df['data_level'] = 'viability'
    df.col_metadata_df['provenance'] = [x + ' | viability' for x in df.col_metadata_df['provenance']]

    DNorm_GCT = GCToo.GCToo(data_df=DNorm_Data, row_metadata_df=row_metadata_df, col_metadata_df=df.col_metadata_df, make_multiindex=True)

    #if write is True:
     #   wg.write(DNorm_GCT, outfile + '_VIABILITY.gct')

    return DNorm_GCT

def log_viability(df, plate_control=False):
    # Calculation of level 3 to level 4 data

    neg_dex = df.col_metadata_df[df.col_metadata_df['pert_type'] == 'ctl_vehicle'].index.tolist()
    neg_df = df.data_df[neg_dex]
    neg_df.median(axis=1)

    if plate_control is False:
        DNorm_Data = df.data_df.subtract(neg_df.median(axis=1), axis='index')

    else:
        DNorm_Data = df.data_df.subtract(df.data_df.median(axis=1), axis='index')

    row_metadata_df = df.row_metadata_df

    DNorm_Data.sort_index(inplace=True)

    row_metadata_df.sort_index(inplace=True)

    df.col_metadata_df['data_level'] = 'viability'
    df.col_metadata_df['provenance'] = [x + ' | viability' for x in df.col_metadata_df['provenance']]

    DNorm_GCT = GCToo.GCToo(data_df=DNorm_Data, row_metadata_df=row_metadata_df, col_metadata_df=df.col_metadata_df,
                            make_multiindex=True)

    # if write is True:
    #   wg.write(DNorm_GCT, outfile + '_VIABILITY.gct')

    return DNorm_GCT