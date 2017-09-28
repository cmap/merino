import sys
sys.path.append('/Users/elemire/Workspace/l1ktools')
import python.broadinstitute_cmap.io.pandasGEXpress.parse as pe
import python.broadinstitute_cmap.io.pandasGEXpress.write_gct as wg
from python.broadinstitute_cmap.io.pandasGEXpress import GCToo
from broadinstitute_cmap.io.pandasGEXpress import concat_gctoo as cg
import pandas as pd
import numpy as np
import os
import math

def upper_triangle(correlation_matrix):
    upper_triangle = correlation_matrix.where(np.triu(np.ones(correlation_matrix.shape), k=1).astype(np.bool))

    # convert matrix into long form description
    upper_tri_series = upper_triangle.stack().reset_index(level=1)

    upper_tri_series.columns = ['rid', 'spearman_corr']

    # Index at this point is CID, it now becomes a column
    upper_tri_series.reset_index(level=0, inplace=True)

    return upper_tri_series


def calculate_weights(correlation_matrix):

    np.fill_diagonal(correlation_matrix.values, 0)
    correlation_matrix[correlation_matrix < 0] = 0
    raw_weights = 0.5 * correlation_matrix.sum(axis=1)
    raw_weights[raw_weights < .01] = .01
    weights = raw_weights / sum(raw_weights.abs())

    return raw_weights, weights

def calculate_sig_strength(modz_values, n_reps):

    modz_values_adjusted = modz_values * math.sqrt(n_reps)

    ss1_5 = modz_values_adjusted[modz_values_adjusted < -3].count()
    ss1 = modz_values_adjusted[modz_values_adjusted < -2].count()
    ss_5 = modz_values_adjusted[modz_values_adjusted < -1].count()

    return ss1_5, ss1, ss_5


def calculate_cis(ss1_5, ss1, ss_5, q75, num_cells):
    q75_div = max(q75, .001) / num_cells
    cis1_5 = np.sqrt(max(ss1_5, 1) * q75_div)
    cis1 = np.sqrt(max(ss1, 1) * q75_div)
    cis_5 = np.sqrt(max(ss_5, 1) * q75_div)

    return cis1_5, cis1, cis_5


def calculate_q75(spearman):

    return spearman.quantile(0.75)


def calculate_modz(zscorepc_paths, project_folder, group_by='pert_well'):

    gct_list = []
    for path in zscorepc_paths:
        gct = pe.parse(path)
        print gct.data_df.shape
        gct_list.append(gct)

    fields_to_remove = [x for x in gct_list[0].row_metadata_df.columns if x in ['det_plate', 'det_plate_scan_time', 'assay_plate_barcode']]
    master_gct = cg.hstack(gct_list, fields_to_remove=fields_to_remove)

    #TODO Change to replicate set ID when we have it in assemble
    #TODO change prism_replicate to replicate_id in assemble
    replicate_set_id = gct.col_metadata_df.index[0].split('_')[0] + '_' + gct.col_metadata_df.index[0].split('_')[1]

    cc_q75_df = pd.DataFrame(
        columns=['weave_prefix', 'det_well', 'profile_ids', 'cc_ut', 'cc_q75', 'nprofile', 'ss_ltn3', 'ss_ltn2',
                 'ss_ltn1', 'cis_ltn3', 'cis_ltn12', 'cis_ltn1'])

    cc_q75_df.index.name = 'sig_id'

    modZ_mat = pd.DataFrame(index=master_gct.data_df.index)
    all_weights = pd.Series()
    all_raw_weights = pd.Series()

    all_corr_values = pd.DataFrame(columns=['weave_prefix', 'cid', 'rid', 'spearman_corr'])

    raw_groupby_vals = set(master_gct.col_metadata_df[group_by])
    groupby_vals = sorted(list(raw_groupby_vals))

    #TODO Abstract out a method which just takes the GCT objects. Find a way to get rid of any field that doesn't match up in cg.hstack.

    for gv_val in groupby_vals:
        dex = master_gct.col_metadata_df[master_gct.col_metadata_df[group_by] == gv_val].index.tolist()
        mat = master_gct.data_df[dex]

        corr_mat = mat.corr(method='spearman')

        upper_tri_series = upper_triangle(corr_mat)

        upper_tri_series['weave_prefix'] = replicate_set_id

        all_corr_values = all_corr_values.append(upper_tri_series)
        upper_tri_series['spearman_corr'][upper_tri_series['spearman_corr'] < 0] = 0

        raw_weights, weights = calculate_weights(corr_mat)

        all_weights = all_weights.append(weights)
        all_raw_weights = all_raw_weights.append(raw_weights)

        weighted_values = mat * weights

        modz_values = weighted_values.sum(axis=1)

        ss1_5, ss1, ss_5 = calculate_sig_strength(modz_values, n_reps=len(corr_mat.index))

        q75 = calculate_q75(upper_tri_series['spearman_corr'])

        cis1_5, cis1, cis_5 = calculate_cis(ss1_5, ss1, ss_5, q75, len(weighted_values))

        cc_q75_df.loc[replicate_set_id + ':' + gv_val] = [replicate_set_id,
                                                   gv_val,
                                                   ','.join(corr_mat.columns.values.tolist()),
                                                   ','.join([str(x) for x in
                                                             upper_tri_series['spearman_corr'].values.tolist()]),
                                                   q75,
                                                   len(corr_mat.index),
                                                   ss1_5,
                                                   ss1,
                                                   ss_5,
                                                   cis1_5,
                                                   cis1,
                                                   cis_5
                                                   ]
        #TODO just use sig_id here
        modz_values[modz_values < -10] = -10
        modz_values[modz_values > 10] = 10
        modZ_mat[mat.columns[0][:11] + mat.columns[0][-4:]] = modz_values


    col_meta = master_gct.col_metadata_df.drop_duplicates(subset=group_by, keep="last")
    col_meta.sort('pert_well', inplace=True)
    col_meta.index = modZ_mat.columns
    col_meta['data_level'] = 'modZ'
    col_meta['provenance'] = gct.col_metadata_df['provenance'] + ' | modZ'

    modZ_mat.index = modZ_mat.index.astype(str)
    master_gct.row_metadata_df.index = master_gct.row_metadata_df.index.astype(str)

    modZ_GCT = GCToo.GCToo(data_df=modZ_mat, row_metadata_df=master_gct.row_metadata_df,
                           col_metadata_df=col_meta)

    outfile = os.path.join(project_folder, 'modZ', replicate_set_id)

    if not os.path.exists(outfile):
        os.mkdir(outfile)

    all_weights.to_csv(os.path.join(outfile, replicate_set_id + '_norm_weights.txt'), sep='\t')
    all_raw_weights.to_csv(os.path.join(outfile, replicate_set_id + '_raw_weights.txt'), sep='\t')
    all_corr_values.set_index(all_corr_values['weave_prefix'], inplace=True)
    del all_corr_values['weave_prefix']
    all_corr_values.to_csv(os.path.join(outfile, replicate_set_id + '_corr_values.txt'), sep='\t')

    wg.write(modZ_GCT, os.path.join(outfile, replicate_set_id + '_MODZ.gct'))
    cc_q75_df.to_csv(os.path.join(outfile, replicate_set_id + 'cc_q75.txt'), sep='\t')