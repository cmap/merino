import sys
sys.path.append('/Users/elemire/Workspace/l1ktools')
import python.broadinstitute_cmap.io.pandasGEXpress.parse as pe
import python.broadinstitute_cmap.io.pandasGEXpress.write_gct as wg
from python.broadinstitute_cmap.io.pandasGEXpress import GCToo

def calculate_zscorevc(filepath, write=False, outfile=''):

    # Calculate level 4 data from level 3 using vehicle control

    invariants = ['661', '662', '663', '664', '665', '666', '667', '668', '669', '670']
    df = pe.parse(filepath)
    neg_dex = df.col_metadata_df[df.col_metadata_df['pert_type'] == 'ctl_vehicle'].index.tolist()
    neg_df = df.data_df[neg_dex]
    neg_df.median(axis=1)

    sub = df.data_df[~df.data_df.index.isin(invariants)].subtract(
        neg_df[~df.data_df.index.isin(invariants)].median(axis=1), axis='index')

    zscore_data = sub.divide(neg_df[~df.data_df.index.isin(invariants)].mad(axis=1) * 1.4826, axis='index')

    zscore_data = zscore_data.append(df.data_df[df.data_df.index.isin(invariants)])
    
    zscore_data.sort_index(inplace=True)
    df.row_metadata_df.sort_index(inplace=True)

    df.col_metadata_df['data_level'] = 'ZSVC'
    df.col_metadata_df['provenance'] = [x + ' | ZSVC' for x in df.col_metadata_df['provenance']]

    Zscore_GCT = GCToo.GCToo(data_df=zscore_data, row_metadata_df=df.row_metadata_df, col_metadata_df=df.col_metadata_df)

    if write is True:
        wg.write(Zscore_GCT, outfile)

    return Zscore_GCT

def calculate_zscorepc(filepath, write=False, outfile=''):

    # Calculate level 4 data from level 3 using plate control

    invariants = ['661', '662', '663', '664', '665', '666', '667', '668', '669', '670']
    df = pe.parse(filepath)

    sub = df.data_df[~df.data_df.index.isin(invariants)].subtract(
        df.data_df[~df.data_df.index.isin(invariants)].median(axis=1), axis='index')

    zscore_data = sub.divide(df.data_df[~df.data_df.index.isin(invariants)].mad(axis=1) * 1.4826, axis='index')

    row_metadata_df = df.row_metadata_df[~df.row_metadata_df.index.isin(invariants)]

    row_metadata_df.sort_index(inplace=True)

    zscore_data.sort_index(inplace=True)

    df.col_metadata_df['data_level'] = 'ZSPC'
    df.col_metadata_df['provenance'] = [x + ' | ZSPC' for x in df.col_metadata_df['provenance']]


    Zscore_GCT = GCToo.GCToo(data_df=zscore_data, row_metadata_df=row_metadata_df, col_metadata_df=df.col_metadata_df)

    if write is True:
        wg.write(Zscore_GCT, outfile)

    return Zscore_GCT


