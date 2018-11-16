import os
import glob
import distil
import functools
import shear
import pandas as pd
import cmapPy.pandasGEXpress.parse as pe
import cmapPy.pandasGEXpress.GCToo as GCToo
import merino.setup_logger as setup_logger
import logging
import argparse
import sys
import cmapPy.pandasGEXpress.write_gct as wg
import json
import ConfigParser
import merino.normalization.batch_adjust as batch_adjust


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
    parser.add_argument("-group_by", "-gb", help="Field(s) to group by for modZ. If you are using more than one field separate columns names with a comma"
                        ,type=str,default='pert_well')
    parser.add_argument("-skip", "-sk", help="Dictionary indicating which columns to exclude from the modZ calculation "
                                             "eg. {'pert_type': ['ctl_vehicle', 'trt_poscon']}, to exclude controls",
                        type=str,default=None)
    parser.add_argument("-log_tf", "-log", help="True or false, if true log transform the data",
                        action="store_true", default=True)
    parser.add_argument("-nprofile_drop", "-nd",
                        help="Drop sigs from modZ with less than two profiles",
                        action="store_false")
    parser.add_argument("-davepool_combat", "-dc",
                        help="Perform combat on the two detection plates - pertains to older data format",
                        action="store_true")

    return parser

def define_pert(proj_dir, input_folder, rep_set):
    print rep_set
    plate_directories = glob.glob(os.path.join(proj_dir, input_folder, rep_set + '_*'))
    plate_names = [os.path.basename(x) for x in plate_directories]
    replicate_ids = [x.split("_")[2] for x in plate_names]
    short_reps = [x.split('.')[0] for x in replicate_ids]
    keep = []

    for r in set(short_reps):
        temp = [y for y in replicate_ids if y.startswith(r)]
        print temp

        if len(temp) == 1:
            keep.append(temp[0])
        else:

            temp2 = [z for z in temp if "." in z]

            max_l = max([int(x[-1]) for x in temp2])
            temp3 = [b for b in temp2 if b.endswith(str(max_l))]
            keep.append(temp3[0])

    keep_perts = [x for x in plate_names if x.split('_')[2] in keep]
    keep_files = [glob.glob(path + '/*')[0] for path in plate_directories if os.path.basename(path) in keep_perts]

    pert = plate_names[0].split('_')[0] + '_' + plate_names[0].split('_')[1]

    gct_list = []
    for path in keep_files:
        gct = pe.parse(path)
        print gct.data_df.shape
        gct_list.append(gct)

    return (pert, gct_list)

def weave(proj_dir, rep_set, args, input_folder='ZSPC', nprofile_drop=True):

    (pert, gct_list) = define_pert(proj_dir, input_folder, rep_set)

    if not os.path.exists(os.path.join(proj_dir, 'modz.{}'.format(input_folder), pert)):
        os.mkdir(os.path.join(proj_dir, 'modz.{}'.format(input_folder), pert))

    reload(distil)


    gb = [x for x in args.group_by.split(',')]

    if args.davepool_combat == True:
        all_ds, pre_list = batch_adjust.combat_by_group(gct_list, col_group=gb, batch_field='davepool_id')
        all_ds, adj_list = batch_adjust.combat_by_group(pre_list, col_group=gb, batch_field='pool_id')

    else:
        print 'here'
        all_ds, adj_list = batch_adjust.combat_by_group(gct_list, col_group=gb, batch_field='pool_id')
        print adj_list[0].data_df.shape

    new_list = []
    for g in adj_list:
        g.data_df = g.data_df.astype(float)
        new_list.append(g)

    for thing in new_list:
        replicate_name = thing.col_metadata_df['prism_replicate'].unique()[0]

        if not os.path.exists(os.path.join(args.proj_dir, '{}.COMBAT'.format(args.input_folder))):
            os.mkdir(os.path.join(args.proj_dir, '{}.COMBAT'.format(args.input_folder)))

        if not os.path.exists(os.path.join(args.proj_dir, '{}.COMBAT'.format(args.input_folder), replicate_name)):
            os.mkdir(os.path.join(args.proj_dir, '{}.COMBAT'.format(args.input_folder), replicate_name))

        wg.write(thing, os.path.join(proj_dir, input_folder + '.COMBAT', replicate_name, replicate_name + '_' + input_folder + '.CB.gct'))

    #if len(keep_files) < 2:
    #    return

    if args.skip is not None:
        modZ_GCT, cc_q75_df, weights = distil.calculate_modz(gct_list, group_by=gb, skip=json.loads(args.skip))
        cb_modZ_GCT, cb_cc_q75_df, cb_weights = distil.calculate_modz(new_list, group_by=gb, skip=json.loads(args.skip))
    else:
        modZ_GCT, cc_q75_df, weights = distil.calculate_modz(gct_list, group_by=gb)
        cb_modZ_GCT, cb_cc_q75_df, cb_weights = distil.calculate_modz(new_list, group_by=gb)

    # Drop sigs where nprofile =1
    if nprofile_drop==True:
        (modZ_GCT, cc_q75_df, cb_modZ_GCT, cb_cc_q75_df) = drop_less_than_2_replicates(modZ_GCT, cc_q75_df, cb_modZ_GCT, cb_cc_q75_df)

    outfile = os.path.join(proj_dir, 'modz.{}'.format(input_folder), pert)
    cb_outfile = os.path.join(proj_dir, 'modz.{}.COMBAT'.format(input_folder), pert)

    if not os.path.exists(outfile):
        os.mkdir(outfile)
    if not os.path.exists(cb_outfile):
        os.mkdir(cb_outfile)


    write_outputs(weights, cb_weights, modZ_GCT, cb_modZ_GCT, cc_q75_df, cb_cc_q75_df, outfile, rep_set, input_folder, cb_outfile, pert)


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

    weights[0].to_csv(os.path.join(outfile, rep_set + '_norm_weights.txt'), sep='\t')
    cb_weights[0].to_csv(os.path.join(cb_outfile, rep_set + '_norm_weights.txt'), sep='\t')
    weights[1].to_csv(os.path.join(outfile, rep_set + '_raw_weights.txt'), sep='\t')
    cb_weights[1].to_csv(os.path.join(cb_outfile, rep_set + '_raw_weights.txt'), sep='\t')
    wg.write(modZ_GCT, os.path.join(outfile, pert + '_MODZ.{}'.format(input_folder)))
    wg.write(cb_modZ_GCT, os.path.join(cb_outfile, pert + '_MODZ.{}.COMBAT'.format(input_folder)))
    cc_q75_df.to_csv(os.path.join(outfile, rep_set + '_cc_q75.txt'), sep='\t')
    cb_cc_q75_df.to_csv(os.path.join(cb_outfile, rep_set + '_cc_q75.txt'), sep='\t')

def setup_directories(project_dir, input_folder):
    if not os.path.exists(os.path.join(project_dir, 'modz.{}'.format(input_folder))):
        os.mkdir(os.path.join(project_dir, 'modz.{}'.format(input_folder)))
    if not os.path.exists(os.path.join(project_dir, 'modz.{}.COMBAT'.format(input_folder))):
        os.mkdir(os.path.join(project_dir, 'modz.{}.COMBAT'.format(input_folder)))

def main(args):
    # replicate set options: all_plates, replicate_set_name, search_pattern; NB required to chose one
    if args.search_pattern is not None:
        glob_value = args.search_pattern
    if args.replicate_set_name is not None:
        glob_value = args.replicate_set_name + "*"
    if args.all_plates:
        glob_value = "*"

    # input_folder options: all_inputs, input_folder; NB required to chose one
    if args.all_inputs:
        for input in ["ZSPC", "ZSVC", "LFCVC", "LFCPC"]:
            replicate_sets_search_results = glob.glob(os.path.join(args.proj_dir, input, glob_value))
            if not replicate_sets_search_results:
                continue

            # Setup directories after passing checks for work to do
            setup_directories(args.proj_dir, input)

            n_tokens = len(os.path.basename(replicate_sets_search_results[0]).split("_"))
            replicate_sets = set([ os.path.basename(x).rsplit("_",n_tokens-2)[0] for x in replicate_sets_search_results])
            if len(replicate_sets) > 1:
                for replicate_set in replicate_sets:
                    print replicate_set, input
                    weave(args.proj_dir, replicate_set, input_folder=input, nprofile_drop=args.nprofile_drop, args=args)
            else:
                replicate_set = replicate_sets.pop()
                print replicate_set, input
                weave(args.proj_dir, replicate_set, input_folder=input, nprofile_drop=args.nprofile_drop, args=args)

    else:
        setup_directories(args.proj_dir, args.input_folder)
        replicate_sets_search_results = glob.glob(os.path.join(args.proj_dir, args.input_folder, glob_value))
        n_tokens = len(os.path.basename(replicate_sets_search_results[0]).split("_"))
        replicate_sets = set([ os.path.basename(x).rsplit("_",n_tokens-2)[0] for x in replicate_sets_search_results])
        if len(replicate_sets) > 1:
            for replicate_set in replicate_sets:
                print replicate_set, args.input_folder
                weave(args.proj_dir, replicate_set, input_folder=args.input_folder, nprofile_drop=args.nprofile_drop, args=args)
        else:
            replicate_set = replicate_sets.pop()
            print replicate_set, args.input_folder
            weave(args.proj_dir, replicate_set, input_folder=args.input_folder, nprofile_drop=args.nprofile_drop, args=args)

    # for x in set([os.path.basename(y).split('_')[0]  + '_' + os.path.basename(y).split('_')[1] for y in glob.glob(os.path.join(args.proj_dir, '{}/*'.format(args.input_folder)))]):
    #     weave(args.proj_dir, x, input_folder=args.input_folder, nprofile_drop=args.nprofile_drop, args=args)

def run_all_weaves():
    pass

if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    logger.debug("args:  {}".format(args))

    main(args)