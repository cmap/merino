import cmapPy.pandasGEXpress.parse as parse
import cmapPy.pandasGEXpress.GCToo as GCToo
import numpy as np
import pandas as pd

# TODO add to metadata
invariant_rids = ['c-661', 'c-662', 'c-663', 'c-664', 'c-665', 'c-666', 'c-667', 'c-668', 'c-669', 'c-670']

def normalize(mfi_gctoo, log=True, inv=True, inv_threshold = 600):

    # Level 2-3 Normalization based on prism invariant

    if inv is True:
        mfi_gctoo = remove_outlier_invariants(mfi_gctoo, inv_threshold)


    mfi_gctoo.data_df[mfi_gctoo.data_df < 1] = 1
    data_df = mfi_gctoo.data_df

    data_df = data_df.round(2)

    if inv is True:
        invs = data_df.loc[invariant_rids].median(axis=0)
        data_df = data_df.divide(invs, axis='columns')

    if log is True:
        data_df = np.log2(data_df)
        mfi_gctoo.col_metadata_df['provenance'] = 'assembled | log2 | median normalized'

    else:
        data_df = data_df
        mfi_gctoo.col_metadata_df['provenance'] = 'assembled | median normalized'

    mfi_gctoo.col_metadata_df['data_level'] = 'normalized'

    data_df = data_df.loc[~data_df.index.isin(invariant_rids)]
    row_metadata_df = mfi_gctoo.row_metadata_df.loc[~mfi_gctoo.row_metadata_df.index.isin(invariant_rids)]

    new_gctoo = GCToo.GCToo(data_df=data_df, row_metadata_df=row_metadata_df, col_metadata_df=mfi_gctoo.col_metadata_df)

    return new_gctoo


def remove_low_bead_wells(mfi_gct, count_gct):

    medians = count_gct.data_df.median(axis=0)

    bad_wells = medians[medians < 20].index

    bad_wells = [x for x in bad_wells if x in mfi_gct.data_df.columns]

    data = mfi_gct.data_df.drop(bad_wells, axis=1)
    col_data = mfi_gct.col_metadata_df.drop(bad_wells, axis=0)

    new_gctoo = GCToo.GCToo(data_df=data, col_metadata_df=col_data, row_metadata_df=mfi_gct.row_metadata_df)

    return new_gctoo


def remove_outlier_invariants(gctoo, inv_threshold):

    invdata = gctoo.data_df.loc[invariant_rids]

    bad_wells = invdata.median()[invdata.median() < inv_threshold].index

    data = gctoo.data_df.drop(bad_wells, axis=1).dropna(axis=1, how='all')
    col_data = gctoo.col_metadata_df.loc[data.columns]

    new_gctoo = GCToo.GCToo(data_df=data, col_metadata_df=col_data, row_metadata_df=gctoo.row_metadata_df)

    return new_gctoo


def dp_normalize(filepath, outfile):

    # For use with DP11/12 protocol (deprecated)

    df = parse.parse(filepath)
    dp11_invariant_rids = ['661', '662', '663', '664', '665', '666', '667', '668', '669', '670']
    dp12_invariant_rids = ['671', '672', '673', '674', '675', '676', '677', '678', '679', '680']

    dp11_dex = df.row_metadata_df[df.row_metadata_df['davepool_id'] == 'DP7'].index.tolist()
    dp12_dex = df.row_metadata_df[df.row_metadata_df['davepool_id'] == 'DP8.1'].index.tolist()

    dp11_df = df.data_df.loc[dp11_dex]

    for column in dp11_df:
        median = dp11_df[column][dp11_invariant_rids].median()

        for index in dp11_df[column].index:
            if index not in dp11_invariant_rids:
                dp11_df[column][index] = dp11_df[column][index] / median

    dp12_df = df.data_df.loc[dp12_dex]

    for column in dp12_df:
        median = dp12_df[column][dp12_invariant_rids].median()

        for index in dp12_df[column].index:
            if index not in dp12_invariant_rids:
                dp12_df[column][index] = dp12_df[column][index] / median

    recombine = dp11_df.append(dp12_df)
    recombine.sort_index(inplace=True)
    df.row_metadata_df.sort_index(inplace=True)
    df.col_metadata_df['provenance'] = 'assembled | median normalized'
    df.col_metadata_df['data_level'] = 'normalized'

    my_gctoo = GCToo.GCToo(data_df=recombine, row_metadata_df=df.row_metadata_df, col_metadata_df=df.col_metadata_df)

    #write_gct.write(my_gctoo, outfile)


def no_inv_norm(mfi_gctoo, log=True, inv=True):
    mfi_gctoo.data_df[mfi_gctoo.data_df < 1] = 1
    log_data = np.log2(mfi_gctoo.data_df)

    dp_data = []

    for x in mfi_gctoo.row_metadata_df['davepool_id'].unique():

        temp = log_data.loc[mfi_gctoo.row_metadata_df[mfi_gctoo.row_metadata_df['davepool_id'] == x].index]

        medians = temp.median()
        group1_dex = medians[medians < np.log2(1500)].index
        group2_dex = medians[medians >= np.log2(1500)].index

        group1_data = temp.loc[:, group1_dex]
        group2_data = temp.loc[:, group2_dex]


        group1_norm = group1_data - 14.2551
        group2_norm = group2_data.subtract(group2_data.median(), axis='columns')

        norm_data = pd.concat([group1_norm, group2_norm], axis=1)

        dp_data.append(norm_data)

    concat_data = pd.concat(dp_data, axis=0)

    concat_data = concat_data.loc[mfi_gctoo.data_df.index,mfi_gctoo.data_df.columns]

    norm_gct = GCToo.GCToo(data_df=concat_data, row_metadata_df=mfi_gctoo.row_metadata_df,
                           col_metadata_df=mfi_gctoo.col_metadata_df)

    return norm_gct
