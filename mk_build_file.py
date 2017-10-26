import cut_to_l2
import glob
import sys
sys.path.append('/Users/elemire/Workspace/cmapPy/cmapPy')
sys.path.append('/Users/elemire/Workspace/l1ktools')
from broadinstitute_cmap.io.pandasGEXpress import concat_gctoo as cg
#import pandasGEXpress.concat_gctoo as cg
#import pandasGEXpress.GCToo as GCToo
from broadinstitute_cmap.io.pandasGEXpress import GCToo as GCToo
#import pandasGEXpress.parse as pe
from broadinstitute_cmap.io.pandasGEXpress import parse as pe
import broadinstitute_cmap.io.pandasGEXpress.write_gctx as wg
#import pandasGEXpress.write_gct as wg
import pandas as pd
import os

def build(search_pattern, outfile, file_suffix):
    gct_list = glob.glob(search_pattern)
    old_len = len(gct_list)
    #gct_list = cut_to_l2.cut_l1(gct_list)
    new_len = len(gct_list)

    print 'Number of old lysate plates removed = {}'.format(old_len - new_len)

    gcts = []
    for gct in gct_list:
        temp = pe.parse(gct)
        gcts.append(temp)



    concat_gct = cg.hstack(gcts, fields_to_remove=['assay_plate_barcode', 'det_plate', 'det_plate_scan_time'])

    concat_gct_wo_meta = GCToo.GCToo(data_df = concat_gct.data_df, row_metadata_df = pd.DataFrame(index=concat_gct.data_df.index),
                                     col_metadata_df=pd.DataFrame(index=concat_gct.col_metadata_df.index))

    print concat_gct_wo_meta.data_df.shape

    wg.write(concat_gct_wo_meta, outfile + 'n{}x{}'.format(concat_gct.data_df.shape[1], concat_gct.data_df.shape[0]) + file_suffix)

    return concat_gct

def complete_build(pod_dir, search_pattern, cohort_name, build_folder):

    modz_path = os.path.join(pod_dir, 'modz', search_pattern, '*.gct')
    print modz_path
    modz_out_path = os.path.join(build_folder, cohort_name + '_LEVEL5_MODZ_')
    sig_data = build(modz_path,modz_out_path, '.gctx')

    zscorepc_path = os.path.join(pod_dir, 'zscorepc', search_pattern, '*ZSPC.gct')
    zscorepc_out_path = os.path.join(build_folder, cohort_name + '_LEVEL4_ZSPC_')
    build(zscorepc_path, zscorepc_out_path, '.gctx')

    zscorevc_path = os.path.join(pod_dir, 'zscorevc', search_pattern, '*ZSVC.gct')
    zscorevc_out_path = os.path.join(build_folder, cohort_name + '_LEVEL4_ZSVC_')
    build(zscorevc_path, zscorevc_out_path, '.gctx')

    viabilitypc_path = os.path.join(pod_dir, 'viabilitypc', search_pattern, '*.gct')
    viabilitypc_out_path = os.path.join(build_folder, cohort_name + '_LEVEL4_FCPC_')
    build(viabilitypc_path, viabilitypc_out_path, '.gctx')

    viabilityvc_path = os.path.join(pod_dir, 'viabilityvc', search_pattern, '*.gct')
    viabilityvc_out_path = os.path.join(build_folder, cohort_name + '_LEVEL4_FCVC_')
    build(viabilityvc_path, viabilityvc_out_path, '.gctx')

    norm_path = os.path.join(pod_dir, 'normalize', search_pattern, '*NORM.gct')
    norm_out_path = os.path.join(build_folder, cohort_name + '_LEVEL3_NORM_')
    inst_data = build(norm_path, norm_out_path, '.gctx')

    mfi_path = os.path.join(pod_dir, 'assemble', search_pattern, '*MEDIAN.gct')
    mfi_out_path = os.path.join(build_folder, cohort_name + '_LEVEL2_MFI_')
    build(mfi_path, mfi_out_path, '.gctx')

    count_path = os.path.join(pod_dir, 'assemble', search_pattern, '*COUNT.gct')
    count_out_path = os.path.join(build_folder, cohort_name + '_LEVEL2_COUNT_')
    build(count_path, count_out_path, '.gctx')

    inst_info = inst_data.col_metadata_df

    for x in ['data_level', 'provenance']:
        del inst_info[x]

    inst_info.to_csv(os.path.join(build_folder, 'inst_info.txt'), sep='\t')


    jon = pd.DataFrame()

    for y in glob.glob(os.path.join(pod_dir, 'modz', search_pattern, '*cc_q75.txt')):

        temp = pd.read_table(y)

        jon = jon.append(temp)

    jon.set_index('sig_id', inplace=True)

    sig_info = sig_data.col_metadata_df.join(jon)

    for x in ['data_level', 'prism_replicate', 'provenance', 'det_well']:
        del sig_info[x]

    sig_info.to_csv(os.path.join(build_folder, 'sig_info.txt'), sep='\t')





