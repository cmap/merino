import sys
sys.path.append('/Users/elemire/Workspace/l1ktools')
import python.broadinstitute_cmap.io.pandasGEXpress.parse as pe
import python.broadinstitute_cmap.io.pandasGEXpress.write_gct as wg
from python.broadinstitute_cmap.io.pandasGEXpress import GCToo

def calculate_viability(filepath, write=False, outfile=''):

    # Calculation of level 3 to level 4 data

    invariants = ['661', '662', '663', '664', '665', '666', '667', '668', '669', '670', '671', '672', '673', '674',
                  '675', '676', '677', '678', '679', '680']
    df = pe.parse(filepath)
    neg_dex = df.col_metadata_df[df.col_metadata_df['pert_type'] == 'ctl_vehicle'].index.tolist()
    neg_df = df.data_df[neg_dex]
    neg_df.median(axis=1)

    DNorm_Data = df.data_df[~df.data_df.index.isin(invariants)].divide(
        neg_df[~df.data_df.index.isin(invariants)].median(axis=1), axis='index')

    row_metadata_df = df.row_metadata_df[~df.row_metadata_df.index.isin(invariants)]

    DNorm_Data.sort_index(inplace=True)

    row_metadata_df.sort_index(inplace=True)

    df.col_metadata_df['data_level'] = 'viability'
    df.col_metadata_df['provenance'] = [x + ' | viability' for x in df.col_metadata_df['provenance']]

    DNorm_GCT = GCToo.GCToo(data_df=DNorm_Data, row_metadata_df=row_metadata_df, col_metadata_df=df.col_metadata_df, make_multiindex=True)

    if write is True:
        wg.write(DNorm_GCT, outfile)

    return DNorm_GCT