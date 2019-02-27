import glob
from broadinstitute_cmap.io.pandasGEXpress import concat_gctoo
from broadinstitute_cmap.io.pandasGEXpress import parse
import pandas as pd
import os
from scipy import stats
import plotly.plotly as py
import plotly.graph_objs as go
import plotly


def calc_mono(pert_plates, project_folder):
    mono_list = []
    for plate in pert_plates:
        gct_list = []
        for gct in glob.glob(os.path.join(project_folder, plate + '.gct')):
            curr_gct = parse.parse(gct, make_multiindex=True)
            gct_list.append(curr_gct)
        import pdb
        pdb.set_trace()
        if len(gct_list) == 3:
            merged_gcts = concat_gctoo.hstack(gct_list,
                                        fields_to_remove=['det_plate', 'det_plate_scan_time', 'assay_plate_barcode'])

        # take unique sets of pert ids and bring them out into their own
        pert_info = merged_gcts.col_metadata_df[['pert_id', 'pert_dose']]
        transpose_df = merged_gcts.data_df.transpose()
        merged_df = pert_info.merge(transpose_df, left_index=True, right_index=True)

        melted_df = pd.melt(merged_df, id_vars=['pert_id', 'pert_dose'])
        no_nan = melted_df[pd.notnull(melted_df["pert_dose"])]

        no_nan['Dose Level'] = no_nan.groupby(['variable', 'pert_id'])['pert_dose'] \
            .apply(lambda x: (x.rank(method='dense')))

        replicate_collapsed = no_nan.groupby(["variable", "pert_id", "Dose Level"], axis=0).median()


        replicate_collapsed = replicate_collapsed.drop('pert_dose', axis=1)
        cornrelations_df = replicate_collapsed.unstack(level='Dose Level')

        def spearmanwrapper(sequence):
            x = stats.spearmanr(sequence.dropna(), range(1, len(sequence.dropna()) + 1))
            return x[0]

        mono_df = cornrelations_df.apply(spearmanwrapper, axis=1).unstack(level='pert_id')

        poscon_ids = merged_gcts.col_metadata_df[merged_gcts.col_metadata_df['pert_type'] == 'trt_poscon'][
            'pert_id'].unique()
        dmso_ids = merged_gcts.col_metadata_df[merged_gcts.col_metadata_df['pert_type'] == 'trt_ctl'][
            'pert_id'].unique()

        for poscon in poscon_ids:
            mono_df = mono_df.drop(poscon, axis=1)
        for dmso in dmso_ids:
            mono_df.drop(dmso, axis=1)

        mono_list.append(mono_df)

    return mono_list



