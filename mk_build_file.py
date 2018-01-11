import cut_to_l2
import glob
import sys
import cmapPy.pandasGEXpress.concat_gctoo as cg
import cmapPy.pandasGEXpress.GCToo as GCToo
import cmapPy.pandasGEXpress.parse as pe
import cmapPy.pandasGEXpress.write_gctx as wg
import pandas as pd
import glob
import merino.cut_to_l2
import cmapPy.pandasGEXpress.write_gct as wgx
import cmapPy.pandasGEXpress.parse as pe
import functools
import merino.setup_logger as setup_logger
import logging
import argparse
import sys
import ConfigParser
import os

logger = logging.getLogger(setup_logger.LOGGER_NAME)


def build_parser():

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # The following arguments are required. These are files that are necessary for assembly and which change
    # frequently between cohorts, replicates, etc.
    parser.add_argument("-proj_dir", "-pd", help="path to the pod directory you want to run card on",
                        type=str, required=True)
    parser.add_argument("-cohort_name", "-cn", help="string designating the prefix to each build file eg. PCAL075-126_T2B",
                        type=str, required=True)
    parser.add_argument("-build_folder", "-bf", help="outfolder for build files",
                        type=str, required=True)
    parser.add_argument("-search_pattern", "-sp",
                        help="Search for this string in the directory, only run plates which contain it. "
                             "Default is wildcard",
                        type=str, default='*', required=False)
    parser.add_argument("-verbose", '-v', help="Whether to print a bunch of output", action="store_true", default=False)
    parser.add_argument("-bad_wells", "-wells", help="List of wells to be excluded from processing", type=list,
                        default=[])
    parser.add_argument("-log_tf", "-log", help="True or false, if true log transform the data",
                        action="store_true", default=True)

    return parser

def build(search_pattern, outfile, file_suffix):
    gct_list = glob.glob(search_pattern)
    old_len = len(gct_list)
    reload(cut_to_l2)
    gct_list = cut_to_l2.cut_l1(gct_list)
    new_len = len(gct_list)

    print 'Number of old lysate plates removed = {}'.format(old_len - new_len)

    gcts = []
    for gct in gct_list:
        if gct.endswith('X4'):
            continue
        temp = pe(gct)
        gcts.append(temp)

    for ct in gcts:
        ct.row_metadata_df = gcts[0].row_metadata_df

    fields_to_remove = [x for x in gcts[0].row_metadata_df.columns if
                        x in ['det_plate', 'det_plate_scan_time', 'assay_plate_barcode']]
    concat_gct = cg.hstack(gcts, fields_to_remove=fields_to_remove)

    concat_gct_wo_meta = GCToo.GCToo(data_df = concat_gct.data_df, row_metadata_df = pd.DataFrame(index=concat_gct.data_df.index),
                                     col_metadata_df=pd.DataFrame(index=concat_gct.col_metadata_df.index))

    print concat_gct_wo_meta.data_df.shape

    wg.write(concat_gct_wo_meta, outfile + 'n{}x{}'.format(concat_gct.data_df.shape[1], concat_gct.data_df.shape[0]) + file_suffix)

    return concat_gct

def main(args):

    modz_path = os.path.join(args.proj_dir, 'modz.ZSPC', args.search_pattern, '*.gct')
    print modz_path
    modz_out_path = os.path.join(args.build_folder, args.cohort_name + '_LEVEL5_MODZ.ZSPC_')
    print 'MODZ'
    zspc_sig_data = build(modz_path,modz_out_path, '.gctx')

    modz_path = os.path.join(args.proj_dir, 'modz.LFCPC', args.search_pattern, '*.gct')
    print modz_path
    modz_out_path = os.path.join(args.build_folder, args.cohort_name + '_LEVEL5_MODZ.LFCPC_')
    print 'MODZ'
    fcpc_sig_data = build(modz_path,modz_out_path, '.gctx')

    zscorepc_path = os.path.join(args.proj_dir, 'ZSPC', args.search_pattern, '*ZSPC.gct')
    zscorepc_out_path = os.path.join(args.build_folder, args.cohort_name + '_LEVEL4_ZSPC_')
    print 'ZSPC'
    build(zscorepc_path, zscorepc_out_path, '.gctx')

    zscorevc_path = os.path.join(args.proj_dir, 'ZSVC', args.search_pattern, '*ZSVC.gct')
    zscorevc_out_path = os.path.join(args.build_folder, args.cohort_name + '_LEVEL4_ZSVC_')
    print 'ZSVC'
    build(zscorevc_path, zscorevc_out_path, '.gctx')

    viabilitypc_path = os.path.join(args.proj_dir, 'LFCPC', args.search_pattern, '*.gct')
    viabilitypc_out_path = os.path.join(args.build_folder, args.cohort_name + '_LEVEL4_LFCPC_')
    print 'LFCPC'
    build(viabilitypc_path, viabilitypc_out_path, '.gctx')

    viabilityvc_path = os.path.join(args.proj_dir, 'LFCVC', args.search_pattern, '*.gct')
    viabilityvc_out_path = os.path.join(args.build_folder, args.cohort_name + '_LEVEL4_LFCVC_')
    print 'LFCVC'
    build(viabilityvc_path, viabilityvc_out_path, '.gctx')

    norm_path = os.path.join(args.proj_dir, 'normalize', args.search_pattern, '*NORM.gct')
    norm_out_path = os.path.join(args.build_folder, args.cohort_name + '_LEVEL3_NORM_')
    print 'NORM'
    build(norm_path, norm_out_path, '.gctx')

    mfi_path = os.path.join(args.proj_dir, 'assemble', args.search_pattern, '*MEDIAN.gct')
    mfi_out_path = os.path.join(args.build_folder, args.cohort_name + '_LEVEL2_MFI_')
    print 'MFI'
    inst_data = build(mfi_path, mfi_out_path, '.gctx')

    count_path = os.path.join(args.proj_dir, 'assemble', args.search_pattern, '*COUNT.gct')
    count_out_path = os.path.join(args.build_folder, args.cohort_name + '_LEVEL2_COUNT_')
    print 'COUNT'
    build(count_path, count_out_path, '.gctx')

    inst_info = inst_data.col_metadata_df

    for x in ['data_level', 'provenance']:
        del inst_info[x]

    inst_info.to_csv(os.path.join(args.build_folder, 'inst_info.txt'), sep='\t')

    jon = pd.DataFrame()

    for y in glob.glob(os.path.join(args.proj_dir, 'modz.ZSPC', args.search_pattern, '*cc_q75.txt')):

        temp = pd.read_table(y)

        jon = jon.append(temp)

    jon.set_index('sig_id', inplace=True)

    zspc_sig_info = zspc_sig_data.col_metadata_df.join(jon)

    for x in ['data_level', 'prism_replicate', 'det_well']:
        del zspc_sig_info[x]

    zspc_sig_info.to_csv(os.path.join(args.build_folder, args.cohort_name + '_MODZ.ZSPC_sig_metrics.txt'), sep='\t')

    jon = pd.DataFrame()

    for y in glob.glob(os.path.join(args.proj_dir, 'modz.LFCPC', args.search_pattern, '*cc_q75.txt')):
        temp = pd.read_table(y)

        jon = jon.append(temp)

    jon.set_index('sig_id', inplace=True)

    fcpc_sig_info = fcpc_sig_data.col_metadata_df.join(jon)

    for x in ['data_level', 'prism_replicate', 'det_well']:
        del fcpc_sig_info[x]

    fcpc_sig_info.to_csv(os.path.join(args.build_folder, args.cohort_name + '_MODZ.LFCPC_sig_metrics.txt'), sep='\t')



if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    logger.debug("args:  {}".format(args))

    main(args)
