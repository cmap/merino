import os
import sys
import glob
import logging
import argparse
import functools
import ConfigParser
import pandas as pd

import cmapPy.pandasGEXpress.concat as cg
import cmapPy.pandasGEXpress.GCToo as GCToo
import cmapPy.pandasGEXpress.parse as pe
import cmapPy.pandasGEXpress.write_gct as wg
import cmapPy.pandasGEXpress.write_gctx as wgx

import merino.misc_tools.cut_to_l2 as cut_to_l2
import merino.setup_logger as setup_logger
import merino.build_summary.ssmd_analysis as ssmd

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
    parser.add_argument("-aggregate_out", "-agg", help="whether weave used aggregate_out flag", action="store_true", default=False)
    parser.add_argument("-verbose", '-v', help="Whether to print a bunch of output", action="store_true", default=False)
    parser.add_argument("-bad_wells", "-wells", help="List of wells to be excluded from processing", type=list,
                        default=[])
    parser.add_argument("-log_tf", "-log", help="True or false, if true log transform the data",
                        action="store_true", default=True)

    return parser

def build(search_pattern, outfile, file_suffix, cut=True):
    gct_list = glob.glob(search_pattern)
    old_len = len(gct_list)

    if cut==True:
        gct_list = cut_to_l2.cut_l1(gct_list)

    new_len = len(gct_list)

    print 'Number of old lysate plates removed = {}'.format(old_len - new_len)

    if new_len == 0:
        return

    gcts = []
    for gct in gct_list:
        if gct.endswith('X4'):
            continue
        temp = pe.parse(gct)
        gcts.append(temp)

    for ct in gcts:
        ct.row_metadata_df = gcts[0].row_metadata_df

    fields_to_remove = [x for x in gcts[0].row_metadata_df.columns if
                        x in ['det_plate', 'det_plate_scan_time', 'assay_plate_barcode']]


    concat_gct = cg.hstack(gcts, False, None, fields_to_remove=fields_to_remove)

    concat_gct_wo_meta = GCToo.GCToo(data_df = concat_gct.data_df, row_metadata_df = pd.DataFrame(index=concat_gct.data_df.index),
                                     col_metadata_df=pd.DataFrame(index=concat_gct.col_metadata_df.index))

    logger.debug("gct shape without metadata: {}".format(concat_gct_wo_meta.data_df.shape))

    wgx.write(concat_gct_wo_meta, outfile + 'n{}x{}'.format(concat_gct.data_df.shape[1], concat_gct.data_df.shape[0]) + file_suffix)

    return concat_gct


def mk_cell_metadata(args):
    if args.aggregate_out:
        paths = glob.glob(os.path.join((args.proj_dir, args.search_pattern, 'card', '*', '*NORM.gct')))
    else:
        paths = glob.glob(os.path.join(args.proj_dir, 'card', args.search_pattern, '*NORM.gct'))

    cell_temp = pe.parse(paths[0])
    cell_temp.row_metadata_df.to_csv(os.path.join(args.build_folder, args.cohort_name + '_cell_info.txt'), sep='\t')

    # Calculate SSMD matrix using paths that were just grabbed and write out
    ssmd_mat = ssmd.ssmd_matrix(paths)

    ssmd_gct = GCToo.GCToo(data_df=ssmd_mat, row_metadata_df=pd.DataFrame(index=ssmd_mat.index),
                           col_metadata_df=pd.DataFrame(index=ssmd_mat.columns))
    wg.write(ssmd_gct, os.path.join(args.build_folder, args.cohort_name + '_ssmd_matrix.gct'))

    ssmd_gct = GCToo.GCToo(data_df=ssmd_mat, col_metadata_df=pd.DataFrame(index=ssmd_mat.columns),
                           row_metadata_df=pd.DataFrame(index=ssmd_mat.index))
    wg.write(ssmd_gct, os.path.join(args.build_folder, args.cohort_name + '_ssmd_matrix.gct'))

def mk_inst_info(inst_data, norm_data, args):

    inst_info = inst_data.col_metadata_df
    inst_info['profile_id'] = inst_info.index

    for x in ['data_level', 'provenance']:
        del inst_info[x]

    inst_info.set_index('profile_id', inplace=True)
    inst_info['is_well_failure'] = False
    inst_info.loc[[x for x in inst_info.index if x not in norm_data.data_df.columns], 'is_well_failure'] = True
    inst_info.to_csv(os.path.join(args.build_folder, args.cohort_name + '_inst_info.txt'), sep='\t')


def mk_sig_info(search_pattern_dict, data_dict, args):

    for key in search_pattern_dict:

        if 'MODZ' in key:

            jon = pd.DataFrame()

            data_id = key.split('.g')[0].replace('*', '')

            if args.aggregate_out:
                meta_paths = glob.glob(os.path.join(args.proj_dir, args.search_pattern, 'weave', '*{}_cc_q75.txt'.format(data_id)))
            else:
                meta_paths = glob.glob(os.path.join(args.proj_dir, 'weave', args.search_pattern, '*{}_cc_q75.txt'.format(data_id)))


            if len(meta_paths) == 0:
                logger.info("No metadata found for {}".format(data_id))
                continue

            for y in meta_paths:
                temp = pd.read_table(y)

                jon = jon.append(temp)

            jon.set_index('sig_id', inplace=True)

            sig_data = data_dict[key]

            fcpc_sig_info = jon.join(sig_data.col_metadata_df)

            for x in ['data_level', 'prism_replicate', 'det_well']:
                del fcpc_sig_info[x]

            fcpc_sig_info.to_csv(os.path.join(args.build_folder, args.cohort_name + '_sig_metrics_{}.txt'.format(data_id)), sep='\t')


def main(args):
    search_pattern_dict = {'*MEDIAN.gct': ['assemble', '_LEVEL2_MFI_'],
                           '*COUNT.gct': ['assemble', '_LEVEL2_COUNT_'],
                           '*NORM.gct': ['card', '_LEVEL3_NORM_'],
                           '*ZSPC.gct': ['card', '_LEVEL4_ZSPC_'],
                           '*ZSPC.COMBAT.gct': ['card', '_LEVEL4_ZSPC.COMBAT_'],
                           '*ZSVC.gct': ['card', '_LEVEL4_ZSVC_'],
                           '*LFCPC.gct': ['card', '_LEVEL4_LFCPC_'],
                           '*LFCPC.COMBAT.gct': ['card', '_LEVEL4_LFCPC.COMBAT_'],
                           '*LFCVC.gct': ['card', '_LEVEL4_LFCVC_'],
                           '*MODZ.ZSPC.gct':['weave', '_LEVEL5_MODZ.ZSPC_'],
                           '*MODZ.LFCPC.gct':['weave', '_LEVEL5_MODZ.LFCPC_'],
                           '*MODZ.ZSPC.COMBAT.gct': ['weave', '_LEVEL5_MODZ.ZSPC.COMBAT_'],
                           '*MODZ.LFCPC.COMBAT.gct': ['weave', '_LEVEL5_MODZ.LFCPC.COMBAT_']}

    data_dict = {}

    for key in search_pattern_dict:
        if args.aggregate_out:
            if search_pattern_dict[key][0] == "weave":
                path = os.path.join(args.proj_dir, args.search_pattern, search_pattern_dict[key][0], key)
            else:
                path = os.path.join(args.proj_dir, args.search_pattern, search_pattern_dict[key][0], '*', key)
        else:
            path = os.path.join(args.proj_dir, search_pattern_dict[key][0],
                            args.search_pattern,key)


        out_path = os.path.join(args.build_folder, args.cohort_name + search_pattern_dict[key][1])

        logger.info("working on {}".format(path))
        if 'MODZ' in key:
            data = build(path, out_path, '.gctx', cut=False)
        else:
            data = build(path, out_path, '.gctx', cut=True)
        data_dict[key] = data

    mk_inst_info(data_dict['*MEDIAN.gct'], data_dict['*NORM.gct'], args)

    mk_sig_info(search_pattern_dict, data_dict, args)

    mk_cell_metadata(args)

def dep(args):

    modz_path = os.path.join(args.proj_dir, 'weave', args.search_pattern,'*MODZ.ZSPC.gct')
    print modz_path
    modz_out_path = os.path.join(args.build_folder, args.cohort_name + '_LEVEL5_MODZ.ZSPC_')
    print 'MODZ'
    zspc_sig_data = build(modz_path,modz_out_path, '.gctx', cut=False)

    cb_modz_path = os.path.join(args.proj_dir, 'weave', args.search_pattern, '*MODZ.ZSPC.COMBAT.gct')
    print cb_modz_path
    cb_modz_out_path = os.path.join(args.build_folder, args.cohort_name + '_LEVEL5_MODZ.ZSPC.COMBAT_')
    print 'MODZ'
    cb_zspc_sig_data = build(cb_modz_path, cb_modz_out_path, '.gctx', cut=False)

    modz_path = os.path.join(args.proj_dir, 'weave', args.search_pattern, '*MODZ.LFCPC.gct')
    print modz_path
    modz_out_path = os.path.join(args.build_folder, args.cohort_name + '_LEVEL5_MODZ.LFCPC_')
    print 'MODZ'
    fcpc_sig_data = build(modz_path,modz_out_path, '.gctx', cut=False)

    modz_path = os.path.join(args.proj_dir, 'weave', args.search_pattern, '*MODZ.LFCPC.COMBAT.gct')
    print modz_path
    modz_out_path = os.path.join(args.build_folder, args.cohort_name + '_LEVEL5_MODZ.LFCPC.COMBAT_')
    print 'MODZ'
    cb_fcpc_sig_data = build(modz_path, modz_out_path, '.gctx', cut=False)

    zscorepc_path = os.path.join(args.proj_dir, 'card', args.search_pattern, '*_ZSPC.gct')
    zscorepc_out_path = os.path.join(args.build_folder, args.cohort_name + '_LEVEL4_ZSPC_')
    print 'ZSPC'
    build(zscorepc_path, zscorepc_out_path, '.gctx')

    zscorepc_path = os.path.join(args.proj_dir, 'card', args.search_pattern, '_ZSPC.COMBAT.gct')
    zscorepc_out_path = os.path.join(args.build_folder, args.cohort_name + '_LEVEL4_ZSPC.COMBAT_')
    print 'ZSPC.CB'
    build(zscorepc_path, zscorepc_out_path, '.gctx')

    zscorevc_path = os.path.join(args.proj_dir, 'card', args.search_pattern, '*_ZSVC.gct')
    zscorevc_out_path = os.path.join(args.build_folder, args.cohort_name + '_LEVEL4_ZSVC_')
    print 'ZSVC'
    build(zscorevc_path, zscorevc_out_path, '.gctx')

    viabilitypc_path = os.path.join(args.proj_dir, 'card', args.search_pattern, '*_LFCPC.gct')
    viabilitypc_out_path = os.path.join(args.build_folder, args.cohort_name + '_LEVEL4_LFCPC_')
    print 'LFCPC'
    build(viabilitypc_path, viabilitypc_out_path, '.gctx')

    viabilitypc_path = os.path.join(args.proj_dir, 'card', args.search_pattern, '*_LFCPC.COMBAT.gct')
    viabilitypc_out_path = os.path.join(args.build_folder, args.cohort_name + '_LEVEL4_LFCPC.COMBAT_')
    print 'LFCPC'
    build(viabilitypc_path, viabilitypc_out_path, '.gctx')

    viabilityvc_path = os.path.join(args.proj_dir, 'card', args.search_pattern, '*_LFCVC.gct')
    viabilityvc_out_path = os.path.join(args.build_folder, args.cohort_name + '_LEVEL4_LFCVC_')
    print 'LFCVC'
    build(viabilityvc_path, viabilityvc_out_path, '.gctx')

    norm_path = os.path.join(args.proj_dir, 'card', args.search_pattern, '*_NORM.gct')
    norm_out_path = os.path.join(args.build_folder, args.cohort_name + '_LEVEL3_NORM_')
    print 'NORM'
    norm_data = build(norm_path, norm_out_path, '.gctx')


    mfi_path = os.path.join(args.proj_dir, 'assemble', args.search_pattern, '*MEDIAN.gct')
    mfi_out_path = os.path.join(args.build_folder, args.cohort_name + '_LEVEL2_MFI_')
    print 'MFI'
    inst_data = build(mfi_path, mfi_out_path, '.gctx')

    count_path = os.path.join(args.proj_dir, 'assemble', args.search_pattern, '*COUNT.gct')
    count_out_path = os.path.join(args.build_folder, args.cohort_name + '_LEVEL2_COUNT_')
    print 'COUNT'
    build(count_path, count_out_path, '.gctx')

    inst_info = inst_data.col_metadata_df
    inst_info['profile_id'] = inst_info.index

    for x in ['data_level', 'provenance']:
        del inst_info[x]

    inst_info.set_index('profile_id', inplace=True)
    inst_info['is_well_failure'] = False
    inst_info.loc[[x for x in inst_info.index if x not in norm_data.data_df.columns], 'is_well_failure'] = True
    inst_info.to_csv(os.path.join(args.build_folder, args.cohort_name + '_inst_info.txt'), sep='\t')

    ################################################################################

    paths = glob.glob(os.path.join(args.proj_dir, 'normalize', args.search_pattern, '*.gct'))
    cell_temp = pe.parse(paths[0])
    cell_temp.row_metadata_df.to_csv(os.path.join(args.build_folder, args.cohort_name + '_cell_info.txt'), sep='\t')

    #Calculate SSMD matrix using paths that were just grabbed and write out
    ssmd_mat = ssmd.ssmd_matrix(paths)

    ssmd_gct = GCToo.GCToo(data_df=ssmd_mat, row_metadata_df=pd.DataFrame(index=ssmd_mat.index),
                           col_metadata_df=pd.DataFrame(index=ssmd_mat.columns))
    wg.write(ssmd_gct, os.path.join(args.build_folder, args.cohort_name + '_ssmd_matrix.gct'))


    ssmd_gct = GCToo.GCToo(data_df=ssmd_mat, col_metadata_df=pd.DataFrame(index=ssmd_mat.columns),
                row_metadata_df=pd.DataFrame(index=ssmd_mat.index))
    wg.write(ssmd_gct, os.path.join(args.build_folder, args.cohort_name + '_ssmd_matrix.gct'))
    ###################################################################################
    jon = pd.DataFrame()
    for y in glob.glob(os.path.join(args.proj_dir, 'modz.ZSPC', args.search_pattern, '*cc_q75.txt')):

        temp = pd.read_table(y)

        jon = jon.append(temp)
    jon.set_index('sig_id', inplace=True)

    zspc_sig_info = jon.join(zspc_sig_data.col_metadata_df)

    for x in ['data_level', 'prism_replicate', 'det_well']:
        del zspc_sig_info[x]

    zspc_sig_info.to_csv(os.path.join(args.build_folder, args.cohort_name + '_sig_metrics_MODZ.ZSPC.txt'), sep='\t')
    ################################################################################
    jon = pd.DataFrame()

    for y in glob.glob(os.path.join(args.proj_dir, 'modz.LFCPC', args.search_pattern, '*cc_q75.txt')):
        temp = pd.read_table(y)

        jon = jon.append(temp)

    jon.set_index('sig_id', inplace=True)

    fcpc_sig_info = jon.join(fcpc_sig_data.col_metadata_df)

    for x in ['data_level', 'prism_replicate', 'det_well']:
        del fcpc_sig_info[x]

    fcpc_sig_info.to_csv(os.path.join(args.build_folder, args.cohort_name + '_sig_metrics_MODZ.LFCPC.txt'), sep='\t')

    #####################################################
    jon = pd.DataFrame()

    for y in glob.glob(os.path.join(args.proj_dir, 'modz.LFCPC.COMBAT', args.search_pattern, '*cc_q75.txt')):
        temp = pd.read_table(y)

        jon = jon.append(temp)

    jon.set_index('sig_id', inplace=True)

    fcpc_sig_info_cb = jon.join(cb_fcpc_sig_data.col_metadata_df)

    for x in ['data_level', 'prism_replicate', 'det_well']:
        del fcpc_sig_info_cb[x]

    #fcpc_sig_info_cb.drop(fcpc_sig_info_cb[fcpc_sig_info_cb['nprofile'] < 2].index, inplace=True)


    fcpc_sig_info_cb.to_csv(os.path.join(args.build_folder, args.cohort_name + '_sig_metrics_MODZ.LFCPC.COMBAT.txt'), sep='\t')

    ####################################################################

    ###################################################################################
    jon = pd.DataFrame()
    for y in glob.glob(os.path.join(args.proj_dir, 'modz.ZSPC.COMBAT', args.search_pattern, '*cc_q75.txt')):

        temp = pd.read_table(y)

        jon = jon.append(temp)
    jon.set_index('sig_id', inplace=True)

    zspc_sig_info_cb = jon.join(cb_zspc_sig_data.col_metadata_df)

    for x in ['data_level', 'prism_replicate', 'det_well']:
        del zspc_sig_info_cb[x]

    #zspc_sig_info_cb.drop(zspc_sig_info_cb[zspc_sig_info_cb['nprofile'] < 2].index, inplace=True)

    zspc_sig_info_cb.to_csv(os.path.join(args.build_folder, args.cohort_name + '_sig_metrics_MODZ.ZSPC.COMBAT.txt'), sep='\t')

    ###################################################################################


    pert_info = zspc_sig_info[['pert_id','pert_id_vendor','pert_iname','pert_mfc_desc','pert_mfc_id','pert_type','x_avg_mol_weight','x_purity','x_smiles']]
    pert_info.set_index('pert_id', inplace=True)
    pert_info.drop_duplicates(inplace=True)
    pert_info.to_csv(os.path.join(args.build_folder, args.cohort_name + '_sig_metrics_pert_info.txt'), sep='\t')

if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    logger.debug("args:  {}".format(args))

    main(args)
