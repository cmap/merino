import os
import glob
import sys
import json
import logging
import argparse

import cmapPy.pandasGEXpress.parse as pe
import cmapPy.pandasGEXpress.GCToo as GCToo
import cmapPy.pandasGEXpress.write_gct as wg

import merino.setup_logger as setup_logger
import merino.normalization.distil as distil
import merino.normalization.batch_adjust as batch_adjust
import merino.misc_tools.cut_to_l2 as cut
import merino.utils.exceptions as merino_exception

logger = logging.getLogger(setup_logger.LOGGER_NAME)

def build_parser():

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # The following arguments are required. These are files that are necessary for assembly and which change
    # frequently between cohorts, replicates, etc.
    parser.add_argument("-proj_dir", "-pd", help="path to the pod directory you want to run card on",
                        type=str, required=True)

    input = parser.add_mutually_exclusive_group(required=True)
    input.add_argument("-input_folder", "-if", help="Search for this string in the directory, only run plates which contain it.",
                        type=str, default='ZSPC')
    input.add_argument("-all_inputs", "-ai", help="Run on every card output", action="store_true", default=False)

    replicate_sets = parser.add_mutually_exclusive_group(required=True)
    replicate_sets.add_argument("-all_plates", "-ap", help="Flag to run on all plates",action="store_true", default=False)
    replicate_sets.add_argument("-replicate_set_name", "-rsn", help="Run on specified replicate_set", type=str, default=None)
    replicate_sets.add_argument("-search_pattern", "-sp", help="Run on plates which contain this pattern", type=str, default=None)


    parser.add_argument("-verbose", '-v', help="Whether to print a bunch of output", action="store_true", default=False)
    parser.add_argument("-bad_wells", "-wells", help="List of wells to be excluded from processing", type=list,
                        default=[])
    parser.add_argument("-group_by", "-gb", help="Field(s) to group by for MODZ. If you are using more than one field separate columns names with a comma"
                        ,type=str,default='pert_well')
    parser.add_argument("-skip", "-sk", help="Dictionary indicating which columns to exclude from the MODZ calculation "
                                             "eg. {'pert_type': ['ctl_vehicle', 'trt_poscon']}, to exclude controls",
                        type=str,default=None)
    #todo: changed this action to store_false, check that this is intended use
    parser.add_argument("-log_tf", "-log", help="True or false, if true log transform the data",
                        action="store_false", default=True)
    parser.add_argument("-nprofile_drop", "-nd",
                        help="Drop sigs from MODZ with less than two profiles",
                        action="store_false")
    parser.add_argument("-davepool_combat", "-dc",
                        help="Perform combat on the two detection plates - pertains to older data format",
                        action="store_true")

    return parser

def main(args):
    # NB required to choose one of replicate set options: all_plates, replicate_set_name, search_pattern
    if args.search_pattern is not None:
        glob_value = args.search_pattern
    if args.replicate_set_name is not None:
        glob_value = args.replicate_set_name + "*"
    if args.all_plates:
        glob_value = "*"

    # NB required choose one of input_folder options: all_inputs, input_folder
    if args.all_inputs:
        for input in ["ZSPC", "ZSVC", "LFCVC", "LFCPC"]:
            replicate_sets_search_results = glob.glob(os.path.join(args.proj_dir, input, glob_value))
            if not replicate_sets_search_results:
                logger.info("Insufficient inputs for input folder {} to run weave. Skipping... ".format(input))
                pass
            # Setup directories after passing checks for work to do
            setup_directories(args.proj_dir, input)

            # Naming convention: pertPlate_assayType_pertTime_replicateNum_beadBatch
            # weave group defined by first three tokens
            n_tokens = len(os.path.basename(replicate_sets_search_results[0]).split("_"))
            all_replicate_sets = set([ os.path.basename(x).rsplit("_",n_tokens-3)[0] for x in replicate_sets_search_results])
            if len(all_replicate_sets) > 1:
                for replicate_set_name in all_replicate_sets:
                    logger.info("Weaving together {} {} files".format(replicate_set_name, input))
                    weave(args.proj_dir, replicate_set_name, input_folder=input, nprofile_drop=args.nprofile_drop, args=args)
            else:
                replicate_set_name = all_replicate_sets.pop()
                logger.info("Weaving together {} {} files".format(replicate_set_name, input))
                weave(args.proj_dir, replicate_set_name, input_folder=input, nprofile_drop=args.nprofile_drop, args=args)

    else:
        setup_directories(args.proj_dir, args.input_folder)
        replicate_sets_search_results = glob.glob(os.path.join(args.proj_dir, args.input_folder, glob_value))
        if not replicate_sets_search_results:
            msg = "Unable to find sufficient inputs for input folder {} to run weave using '{}' to search".format(
                args.input_folder, glob_value)
            raise merino_exception.ReplicateSetSearchFailure(msg)

        n_tokens = len(os.path.basename(replicate_sets_search_results[0]).split("_"))
        all_replicate_sets = set([ os.path.basename(x).rsplit("_",n_tokens-3)[0] for x in replicate_sets_search_results])
        if len(all_replicate_sets) > 1:
            for replicate_set_name in all_replicate_sets:
                logger.info("Weaving together {} {} files".format(replicate_set_name, args.input_folder))
                weave(args.proj_dir, replicate_set_name, input_folder=args.input_folder, nprofile_drop=args.nprofile_drop, args=args)
        else:
            replicate_set_name = all_replicate_sets.pop()
            logger.info("Weaving together {} {} files".format(replicate_set_name, args.input_folder))
            weave(args.proj_dir, replicate_set_name, input_folder=args.input_folder, nprofile_drop=args.nprofile_drop, args=args)


def setup_directories(project_dir, input_folder):
    if not os.path.exists(os.path.join(project_dir, 'MODZ.{}'.format(input_folder))):
        os.mkdir(os.path.join(project_dir, 'MODZ.{}'.format(input_folder)))
    if not os.path.exists(os.path.join(project_dir, 'MODZ.{}.COMBAT'.format(input_folder))):
        os.mkdir(os.path.join(project_dir, 'MODZ.{}.COMBAT'.format(input_folder)))


def weave(proj_dir, replicate_set_name, args, input_folder='ZSPC', nprofile_drop=True):

    gct_list = define_replicate_set_files_and_parse(proj_dir, input_folder, replicate_set_name)

    if not os.path.exists(os.path.join(proj_dir, 'MODZ.{}'.format(input_folder), replicate_set_name)):
        os.mkdir(os.path.join(proj_dir, 'MODZ.{}'.format(input_folder), replicate_set_name))

    reload(distil)

    group_by_list = [x for x in args.group_by.split(',')]

    if args.davepool_combat == True:
        all_ds, pre_list = batch_adjust.combat_by_group(gct_list, col_group=group_by_list, batch_field='davepool_id')
        all_ds, combat_adjusted_gct_list = batch_adjust.combat_by_group(pre_list, col_group=group_by_list, batch_field='pool_id')

    else:
        all_ds, combat_adjusted_gct_list = batch_adjust.combat_by_group(gct_list, col_group=group_by_list, batch_field='pool_id')
        logger.debug("sample combat adjusted gct shape {}".format(combat_adjusted_gct_list[0].data_df.shape))


    for combat_adjusted_gct in combat_adjusted_gct_list:
        replicate_name = combat_adjusted_gct.col_metadata_df['prism_replicate'].unique()[0]

        if not os.path.exists(os.path.join(proj_dir, '{}.COMBAT'.format(input_folder))):
            os.mkdir(os.path.join(proj_dir, '{}.COMBAT'.format(input_folder)))

        if not os.path.exists(os.path.join(proj_dir, '{}.COMBAT'.format(input_folder), replicate_name)):
            os.mkdir(os.path.join(proj_dir, '{}.COMBAT'.format(input_folder), replicate_name))

        wg.write(combat_adjusted_gct, os.path.join(proj_dir, input_folder + '.COMBAT', replicate_name, replicate_name + '_' + input_folder + '.COMBAT.gct'))


    if args.skip is not None:
        modZ_GCT, cc_q75_df, weights = distil.calculate_modz(gct_list, group_by=group_by_list, skip=json.loads(args.skip))
        cb_modZ_GCT, cb_cc_q75_df, cb_weights = distil.calculate_modz(combat_adjusted_gct_list, group_by=group_by_list, skip=json.loads(args.skip))

    else:
        modZ_GCT, cc_q75_df, weights = distil.calculate_modz(gct_list, group_by=group_by_list)
        cb_modZ_GCT, cb_cc_q75_df, cb_weights = distil.calculate_modz(combat_adjusted_gct_list, group_by=group_by_list)

    # Drop sigs where nprofile =1
    if nprofile_drop==True:
        (modZ_GCT, cc_q75_df, cb_modZ_GCT, cb_cc_q75_df) = drop_less_than_2_replicates(modZ_GCT, cc_q75_df, cb_modZ_GCT, cb_cc_q75_df)

    outfile = os.path.join(proj_dir, 'MODZ.{}'.format(input_folder), replicate_set_name)
    cb_outfile = os.path.join(proj_dir, 'MODZ.{}.COMBAT'.format(input_folder), replicate_set_name)

    write_outputs(weights, cb_weights, modZ_GCT, cb_modZ_GCT, cc_q75_df, cb_cc_q75_df, outfile, replicate_set_name, input_folder, cb_outfile, replicate_set_name)


def define_replicate_set_files_and_parse(proj_dir, input_folder, replicate_set_name):
    logger.info("defining replicate set files for {}".format(replicate_set_name))
    plate_directories = glob.glob(os.path.join(proj_dir, input_folder, replicate_set_name + '*', '*'))
    keep_files = cut.cut_l1(plate_directories)
    if len(keep_files) <= 1:
        msg = "Insufficient number of replicates for replicate set {} in {} input folder".format(replicate_set_name, input_folder)
        raise merino_exception.ReplicateSetSearchFailure(msg)

    gct_list = []
    for path in keep_files:
        gct = pe.parse(path)
        logger.info("parsed GCT {} with data_df.shape {}".format(os.path.basename(path), gct.data_df.shape))
        gct_list.append(gct)


    return gct_list


def drop_less_than_2_replicates(modZ_GCT, cc_q75_df, cb_modZ_GCT, cb_cc_q75_df):
    sub_mat = modZ_GCT.data_df.drop(cc_q75_df[cc_q75_df['nprofile'] < 2].index, axis=1)
    sub_col_mat = modZ_GCT.col_metadata_df.drop(cc_q75_df[cc_q75_df['nprofile'] < 2].index)
    modZ_GCT = GCToo.GCToo(data_df=sub_mat, col_metadata_df=sub_col_mat, row_metadata_df=modZ_GCT.row_metadata_df)
    cc_q75_df.drop(cc_q75_df[cc_q75_df['nprofile'] < 2].index, inplace=True)

    sub_mat = cb_modZ_GCT.data_df.drop(cb_cc_q75_df[cb_cc_q75_df['nprofile'] < 2].index, axis=1)
    sub_col_mat = cb_modZ_GCT.col_metadata_df.drop(cb_cc_q75_df[cb_cc_q75_df['nprofile'] < 2].index)
    cb_modZ_GCT = GCToo.GCToo(data_df=sub_mat, col_metadata_df=sub_col_mat,
                              row_metadata_df=cb_modZ_GCT.row_metadata_df)
    cb_cc_q75_df.drop(cb_cc_q75_df[cb_cc_q75_df['nprofile'] < 2].index, inplace=True)

    return (modZ_GCT, cc_q75_df, cb_modZ_GCT, cb_cc_q75_df)

def write_outputs(weights, cb_weights, modZ_GCT, cb_modZ_GCT, cc_q75_df, cb_cc_q75_df, outfile, rep_set, input_folder, cb_outfile, pert):

    if not os.path.exists(outfile):
        os.mkdir(outfile)
    if not os.path.exists(cb_outfile):
        os.mkdir(cb_outfile)

    weights[0].to_csv(os.path.join(outfile, rep_set + '_norm_weights.txt'), sep='\t')
    cb_weights[0].to_csv(os.path.join(cb_outfile, rep_set + '_norm_weights.txt'), sep='\t')
    weights[1].to_csv(os.path.join(outfile, rep_set + '_raw_weights.txt'), sep='\t')
    cb_weights[1].to_csv(os.path.join(cb_outfile, rep_set + '_raw_weights.txt'), sep='\t')
    wg.write(modZ_GCT, os.path.join(outfile, pert + '_MODZ.{}'.format(input_folder)))
    wg.write(cb_modZ_GCT, os.path.join(cb_outfile, pert + '_MODZ.{}.COMBAT'.format(input_folder)))
    cc_q75_df.to_csv(os.path.join(outfile, rep_set + '_cc_q75.txt'), sep='\t')
    cb_cc_q75_df.to_csv(os.path.join(cb_outfile, rep_set + '_cc_q75.txt'), sep='\t')

if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    logger.debug("args:  {}".format(args))

    main(args)