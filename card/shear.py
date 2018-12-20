import cmapPy.pandasGEXpress.GCToo as GCToo
import pandas as pd

def shear(gctoo, bad_wells):
    remove = gctoo.col_metadata_df[gctoo.col_metadata_df['pert_well'].isin(bad_wells)].index
    new_data_df = gctoo.data_df.drop(remove, axis=1)
    new_col_df = gctoo.col_metadata_df.drop(remove)
    new_gctoo = GCToo.GCToo(data_df=new_data_df, col_metadata_df=new_col_df, row_metadata_df=gctoo.row_metadata_df)

    return new_gctoo