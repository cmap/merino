import broadinstitute_cmap.io.GCToo.parse as parse
import broadinstitute_cmap.io.GCToo.write_gctoo as write_gctoo
import sys
sys.path.append('/Users/elemire/Workspace/l1ktools')
from python.broadinstitute_cmap.io.pandasGEXpress import GCToo

def normalize(filepath, outfile=None):

    # Level 2-3 Normalization based on prism invariant

    my_gctoo = parse.parse(filepath)

    invariant_rids = ['661', '662', '663', '664', '665', '666', '667', '668', '669', '670']

    for column in my_gctoo.data_df:
        median = my_gctoo.data_df[column][invariant_rids].median()

        for index in my_gctoo.data_df[column].index:
            if index not in invariant_rids:
                my_gctoo.data_df[column][index] = my_gctoo.data_df[column][index] / median

    my_gctoo.col_metadata_df['provenance'] = 'assembled | median normalized'
    my_gctoo.col_metadata_df['data_level'] = 'normalized'

    if outfile is None:
        outfile = filepath[:-10] + 'NORM.gct'


    write_gctoo.write(my_gctoo, outfile)

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

    write_gctoo.write(my_gctoo, outfile)