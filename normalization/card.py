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
                        type=str, default=None)
    plates_group.add_argument("-plate_name", "-pn", help="name of individual plate to run on", type=str, default=None)


    parser.add_argument("-verbose", '-v', help="Whether to print a bunch of output", action="store_true", default=False)
    parser.add_argument("-bad_wells", "-wells", help="List of wells to be excluded from processing", type=list,
                        default=[])
    parser.add_argument("-log_tf", "-log", help="True or false, if true log transform the data",
                        action="store_false")
    parser.add_argument("-inv_tf", "-inv", help="True or false, if true normalize to invariants",
                        action="store_false")
    parser.add_argument("-no_invariants", "-ni", help="True or false, if true log transform the data",
                        action="store_true")
    parser.add_argument("-flattened", help="whether directory structure is flattened with platenames at top level", action="store_true", default=False)

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


def setup_paths(proj_dir, plate_name, flattened):
    if flattened:
        top_level_dir = plate_name
        assemble_path = os.path.join(top_level_dir, "assemble", plate_name + '_MEDIAN.gct')
        count_path = os.path.join(top_level_dir, "assemble", plate_name + '_COUNT.gct')
        norm_dir = os.path.join(top_level_dir, "normalize")

    else :
    # Make Level 3 data folder
        top_level_dir = proj_dir
        # Get path to raw mfi
        assemble_path = os.path.join(top_level_dir, 'assemble', plate_name, plate_name + '_MEDIAN.gct')
        # Get path to beadcount values
        count_path = os.path.join(top_level_dir, 'assemble', plate_name, plate_name + '_COUNT.gct')
        # Get path to LEVEL3 norm values
        norm_dir = os.path.join(top_level_dir, 'normalize', plate_name)

    norm_path = os.path.join(norm_dir, plate_name + '_NORM.gct')

    if not os.path.exists(os.path.join(top_level_dir, 'normalize')):
        os.mkdir(os.path.join(top_level_dir, 'normalize'))


    return (top_level_dir, assemble_path, count_path, norm_dir, norm_path)


def card(proj_dir, plate_name, log_tf=True, inv_tf=True, bad_wells=[], dp=False, flattened=False):
    # Set plate_failure variable to false
    plate_failure = False

    (top_level_dir, assemble_path, count_path, norm_dir, norm_path) = setup_paths(proj_dir, plate_name, flattened)

    if not os.path.exists(norm_dir):
        # Create norm folder for plate if it doesn't exist already
        os.mkdir(norm_dir)
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
    for dir_name in lvl4_card_map.keys():
        if not os.path.exists(os.path.join(top_level_dir, dir_name)):
            os.mkdir(os.path.join(top_level_dir, dir_name))
        # flattened = False
        if top_level_dir == proj_dir and not os.path.exists(os.path.join(top_level_dir, dir_name, plate_name)):
            os.mkdir(os.path.join(top_level_dir, dir_name, plate_name))
            output_path = os.path.join(top_level_dir, dir_name, plate_name, plate_name + lvl4_card_map[dir_name][1])

        else : # flattened = True
            output_path = os.path.join(top_level_dir, dir_name, plate_name + lvl4_card_map[dir_name][1])

        reader_writer(input_file=norm_path, output_file=output_path, function=lvl4_card_map[dir_name][0])

    # Return status of plate failure

    return plate_failure


def main(args):
    # N.B. search pattern is not setup to handle flattened directory structure
    if args.search_pattern:
        failure_list = []
        for folder in glob.glob(os.path.join(args.proj_dir, 'assemble', args.search_pattern)):
            plate_name = os.path.basename(folder)
            print plate_name
            plate_failure = card(args.proj_dir, plate_name, log_tf=args.log_tf, inv_tf=args.inv_tf, bad_wells=args.bad_wells, dp=args.no_invariants, flattened=args.flattened)
            if plate_failure == True:
                failure_list.append(plate_name)

        print failure_list
        if len(failure_list) > 0:
            pd.Series(failure_list).to_csv(os.path.join(args.proj_dir, 'failed_plates.txt'), sep='\t')
            return

    else:
        failure = card(args.proj_dir, args.plate_name, log_tf=args.log_tf, inv_tf=args.inv_tf, bad_wells=args.bad_wells, dp=args.no_invariants, flattened=args.flattened)
        if failure:
            plate_failure_path = os.path.join(args.proj_dir, "normalize", args.plate_name, "failure.txt")
            with open(plate_failure_path, "w") as file:
                file.write("{} failed size checks".format(args.plate_name))
                return

    if args.flattened and args.plate_name:
        plate_success_path = os.path.join(args.plate_name, "normalize", "success.txt")
    elif args.flattened is False and args.plate_name:
        plate_success_path = os.path.join(args.proj_dir, "normalize", args.plate_name, "success.txt")
    else:
        plate_success_path = os.path.join(args.proj_dir, "assemble_success.txt")
    with open(plate_success_path, "w") as file:
        if args.plate_name:
            file.write("{} successfully processed".format(args.plate_name))
        else: file.write("all plates in project {} successfully processed".format(args.proj_dir))
if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    logger.debug("args:  {}".format(args))

    main(args)