import python.broadinstitute_cmap.io.pandasGEXpress.parse as pe
import python.broadinstitute_cmap.io.pandasGEXpress.write_gct as wg
from python.broadinstitute_cmap.io.pandasGEXpress import GCToo

def calculate_davenorm(filepath, write=False, outfile=''):
    invariants = ['661', '662', '663', '664', '665', '666', '667', '668', '669', '670']
    df = pe.parse(filepath)
    neg_dex = df.col_metadata_df[df.col_metadata_df['pert_type'] == 'ctl_vehicle'].index.tolist()
    neg_df = df.data_df[neg_dex]
    neg_df.median(axis=1)

    DNorm_Data = df.data_df[~df.data_df.index.isin(invariants)].divide(
        neg_df[~df.data_df.index.isin(invariants)].median(axis=1), axis='index')
    DNorm_Data = DNorm_Data.append(df.data_df[df.data_df.index.isin(invariants)])

    DNorm_GCT = GCToo.GCToo(data_df=DNorm_Data, row_metadata_df=df.row_metadata_df, col_metadata_df=df.col_metadata_df, make_multiindex=True)

    if write is True:
        wg.write(DNorm_GCT, outfile)

    return DNorm_GCT