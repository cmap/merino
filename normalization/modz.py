import sys
import python.broadinstitute_cmap.io.pandasGEXpress.parse as pe
import python.broadinstitute_cmap.io.pandasGEXpress.write_gct as wg
from python.broadinstitute_cmap.io.pandasGEXpress import GCToo
from broadinstitute_cmap.io.pandasGEXpress import concat_gctoo as cg
import pandas as pd
import numpy as np
import os

def calculate_modz(zscore_paths, project_folder):

    gct_list = []
    for path in zscore_paths:
        gct = pe.parse(path)
        gct_list.append(gct)

    master_gct = cg.hstack(gct_list, fields_to_remove = ['det_plate', 'det_plate_scan_time', 'assay_plate_barcode'])

    modZ_mat = pd.DataFrame(index = master_gct.data_df.index)
    all_weights = pd.Series()
    raw_weights = pd.Series()
    raw_wells = set([x[-3:] for x in master_gct.data_df.columns])
    wells = sorted(list(raw_wells))
    for well in wells:
        mat = master_gct.data_df[[x for x in master_gct.data_df.columns if x[-3:] == well]]
        corr_mat = mat.corr(method='spearman')
        corr_mat[corr_mat < 0] = 0
        np.fill_diagonal(corr_mat.values, 0)
        raw_weights = 0.5 * corr_mat.sum(axis=1)
        raw_weights[raw_weights < .01] = .01
        weights = raw_weights / sum(raw_weights.abs())
        all_weights = all_weights.append(weights)
        raw_weights = raw_weights.append(raw_weights)
        weighted_values = mat * weights
        modZ_mat[mat.columns[0][:11] + mat.columns[0][-4:]] = weighted_values.sum(axis=1)

    gct.col_metadata_df.index = [x[:11] + x[-4:] for x in gct.col_metadata_df.index]
    gct.col_metadata_df['prism_replicate'] = [x[:11] for x in gct.col_metadata_df['prism_replicate']]
    gct.col_metadata_df['data_level'] = 'modZ'
    gct.col_metadata_df['provenance'] = [x + ' | modZ' for x in gct.col_metadata_df['provenance']]

    modZ_GCT = GCToo.GCToo(data_df=modZ_mat, row_metadata_df=master_gct.row_metadata_df, col_metadata_df=gct.col_metadata_df)

    brew_prefix = gct.col_metadata_df['prism_replicate'].iloc[0]

    outfile = os.path.join(project_folder, 'modZ', brew_prefix)

    all_weights.to_csv(os.path.join(outfile, brew_prefix + '_norm_weights.txt'), sep='\t')
    raw_weights.to_csv(os.path.join(outfile, brew_prefix + '_raw_weights.txt'), sep='\t')

    if not os.path.exists(outfile):
        os.mkdir(outfile)
    wg.write(modZ_GCT, os.path.join(outfile, brew_prefix + '_MODZ.gct'))