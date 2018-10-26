import normalize_with_prism_invariant as norm
import viability_normalization as viability
import zscore
import os
import glob
import sys
#todo : um what
sys.path.append('/Users/elemire/Workspace/cmapPy')
import cmapPy.pandasGEXpress.write_gct as wgx
import functools
import shear
import pandas as pd
import merino.setup_logger as setup_logger
import logging
import argparse
import sys

import cmapPy.pandasGEXpress.parse as pe


logger = logging.getLogger(setup_logger.LOGGER_NAME)

def build_parser():

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # The following arguments are required. These are files that are necessary for assembly and which change
    # frequently between cohorts, replicates, etc.
    parser.add_argument("-proj_dir", "-pd", help="path to the pod directory you want to run card on",
                        type=str, required=True)
    plates_group = parser.add_mutually_exclusive_group(required=True)
    plates_group.add_argument("-search_pattern", "-sp",
                        help="Search for this string in the directory, only run plates which contain it. "
                             "Default is wildcard",
                        type=str, default='*')
    plates_group.add_argument("-plate_name", "-pn", help="name of individual plate to run on", type=str)


    parser.add_argument("-verbose", '-v', help="Whether to print a bunch of output", action="store_true", default=False)
    parser.add_argument("-bad_wells", "-wells", help="List of wells to be excluded from processing", type=list,
                        default=[])
    parser.add_argument("-log_tf", "-log", help="True or false, if true log transform the data",
                        action="store_false")

    return parser


def reader_writer(input_file, output_file, function, check_size=False):
    plate_failure = False
    # Read in input file
    gctoo = pe.parse(input_file)
    # Call normalizing function on gctoo
    new_gctoo = function(gctoo)

    # If told to, check size of new_gctoo and flag if too small
    if new_gctoo.data_df.shape[1] <= 349 and check_size==True:
        print '{} Plate Failure With {} Failed Wells'.format(os.path.basename(os.path.dirname(input_file)),
                                                             384 - new_gctoo.data_df.shape[1])
        plate_failure = True

    # write out new gctoo
    wgx.write(new_gctoo, out_fname=output_file)
    print output_file

    return plate_failure


def card(proj_dir, plate_name, log_tf=True, inv_tf=True, bad_wells=[], dp=False):

    # Make Level 3 data folder
    if not os.path.exists(os.path.join(proj_dir, 'normalize')):
        os.mkdir(os.path.join(proj_dir, 'normalize'))

    # Get path to raw mfi
    assemble_path = os.path.join(proj_dir, 'assemble', plate_name, plate_name + '_MEDIAN.gct')
    # Get path to beadcount values
    count_path = os.path.join(proj_dir, 'assemble', plate_name, plate_name + '_COUNT.gct')
    # Get path to LEVEL3 norm values
    norm_path = os.path.join(proj_dir, 'normalize', plate_name, plate_name + '_NORM.gct')
    # Set plate_failure variable to false
    plate_failure = False

    if not os.path.exists(os.path.join(proj_dir, 'normalize', plate_name)):
        # Create norm folder for plate if it doesn't exist already
        os.mkdir(os.path.join(proj_dir, 'normalize', plate_name))

        # Create norm file
        if dp == True:
            reader_writer(assemble_path, norm_path, norm.no_inv_norm)
        else:
            reader_writer(assemble_path, norm_path, functools.partial(norm.normalize, log=log_tf, inv=inv_tf))

        # Read in count file
        count_gctoo = pe.parse(count_path)
        # Remove low bead count wells and check GCT size, if too many wells have been stripped it will qualify as a failure
        plate_failure = reader_writer(norm_path, norm_path, functools.partial(norm.remove_low_bead_wells, count_gct=count_gctoo), check_size=True)
        # Shear predetermined bad wells (if any exist)
        reader_writer(norm_path, norm_path, functools.partial(shear.shear, bad_wells=bad_wells))

    reload(viability)

    if log_tf==True:
    # Map denoting each type of LEVEL4 data, its folder name, the function to create it, and the file ending.
        lvl4_card_map = {'ZSVC': [zscore.calculate_zscore, '_ZSVC.gct'],
                     'ZSPC': [functools.partial(zscore.calculate_zscore, plate_control=True), '_ZSPC.gct'],
                     'LFCPC': [functools.partial(viability.log_viability, plate_control=True, log=True), '_FCPC.gct'],
                     'LFCVC': [viability.log_viability, '_FCVC.gct']}

    else:
        lvl4_card_map = {'ZSVC': [zscore.calculate_zscore, '_ZSVC.gct'],
                         'ZSPC': [functools.partial(zscore.calculate_zscore, plate_control=True), '_ZSPC.gct'],
                         'LFCPC': [functools.partial(viability.log_viability, plate_control=True, log=False), '_FCPC.gct'],
                         'LFCVC': [functools.partial(viability.log_viability, plate_control=False, log=False), '_FCVC.gct']}

    # Loop through this map to output all level 4 data
    for x in lvl4_card_map.keys():
        if not os.path.exists(os.path.join(proj_dir, x)):
            os.mkdir(os.path.join(proj_dir, x))
        if not os.path.exists(os.path.join(proj_dir, x, plate_name)):
            os.mkdir(os.path.join(proj_dir, x, plate_name))

            output_path = os.path.join(proj_dir, x, plate_name, plate_name + lvl4_card_map[x][1])

            reader_writer(input_file=norm_path, output_file=output_path, function=lvl4_card_map[x][0])

    # Return status of plate failure

    return plate_failure


def main(args):
    if args.search_pattern:
        failure_list = []
        for folder in glob.glob(os.path.join(args.proj_dir, 'assemble', args.search_pattern)):
            name = os.path.basename(folder)
            print name
            plate_failure = card(args.proj_dir, name, log_tf=args.log_tf, inv_tf=args.inv_tf, bad_wells=args.bad_wells, dp=args.no_invariants)
            if plate_failure == True:
                failure_list.append(name)

        print failure_list
        pd.Series(failure_list).to_csv(os.path.join(args.proj_dir, 'failed_plates.txt'), sep='\t')

    else:
        failure = card(args.proj_dir, args.plate_name, log_tf=args.log_tf, inv_tf=args.inv_tf, bad_wells=args.bad_wells, dp=args.no_invariants)
        if failure:
            plate_failure_path = os.path.join(args.proj_dir, "assemble", args.plate_name, "failure.txt")
            with open(plate_failure_path, "w") as file:
                file.write("{} failed size checks")


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    logger.debug("args:  {}".format(args))

    main(args)


def oldmain(proj_dir, search_pattern='*', log_tf=True, bad_wells=[]):
    failure_list = []
    for folder in glob.glob(os.path.join(proj_dir, 'assemble', search_pattern)):
        name = os.path.basename(folder)
        plate_failure = card(proj_dir, name, log_tf = log_tf, bad_wells=bad_wells)
        if plate_failure == True:
            failure_list.append(name)

    print failure_list
    pd.Series(failure_list).to_csv(os.path.join(proj_dir, 'failed_plates.txt'), sep='\t')

