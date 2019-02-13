import os
import glob
import sys
import functools
import logging

import argparse
import pandas as pd

import cmapPy.pandasGEXpress.parse as pe
import cmapPy.pandasGEXpress.write_gct as wgx

import merino.setup_logger as setup_logger
import merino.card.normalize_with_prism_invariant as norm
import merino.card.viability_normalization as viability
import merino.card.shear as shear
import merino.card.zscore as zscore

logger = logging.getLogger(setup_logger.LOGGER_NAME)

def build_parser():

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # The following arguments are required. These are files that are necessary for assembly and which change
    # frequently between cohorts, replicates, etc.
    parser.add_argument("-proj_dir", "-pd", help="path to the pod directory you want to run card on",
                        type=str, required=True)
    plates_group = parser.add_mutually_exclusive_group(required=True)
    plates_group.add_argument("-search_pattern", "-sp",
                        help="Search for this string in the directory, only run plates which contain it. ",
                        type=str, default=None)
    plates_group.add_argument("-plate_name", "-pn", help="name of individual plate to run on", type=str, default=None)
    parser.add_argument("-verbose", '-v', help="Whether to print a bunch of output", action="store_true", default=False)
    parser.add_argument("-bad_wells", "-wells", help="List of wells to be excluded from processing", type=list,
                        default=[])
    parser.add_argument("-log_tf", "-log", help="True or false, if true log transform the data",
                        action="store_false", default=True)
    parser.add_argument("-inv_tf", "-inv", help="True or false, if true normalize to invariants",
                        action="store_false")
    parser.add_argument("-no_invariants", "-ni", help="True or false, if true log transform the data",
                        action="store_true")

    return parser


def reader_writer(input_file, output_file, function, check_size=False):
    plate_failure = False
    # Read in input file
    gctoo = pe.parse(input_file)
    # Call normalizing function on gctoo
    new_gctoo = function(gctoo)

    # If told to, check size of new_gctoo and flag if too small
    if new_gctoo.data_df.shape[1] <= 349 and check_size==True:
        logger.debug('{} Plate Failure With {} Failed Wells'.format(os.path.basename(os.path.dirname(input_file)),
                                                             384 - new_gctoo.data_df.shape[1]))
        plate_failure = True

    # write out new gctoo
    wgx.write(new_gctoo, out_fname=output_file)
    logger.debug("{} file written.".format(output_file))

    return plate_failure


def card(proj_dir, plate_name, log_tf=True, inv_tf=True, bad_wells=[], dp=False):

    # Make Level 3 data folder
    if not os.path.exists(os.path.join(proj_dir, 'card')):
        os.mkdir(os.path.join(proj_dir, 'card'))

    # Get path to raw mfi
    assemble_path = os.path.join(proj_dir, 'assemble', plate_name, plate_name + '_MEDIAN.gct')
    # Get path to beadcount values
    count_path = os.path.join(proj_dir, 'assemble', plate_name, plate_name + '_COUNT.gct')
    # Get path to LEVEL3 norm values
    norm_path = os.path.join(proj_dir, 'card', plate_name, plate_name + '_NORM.gct')
    # Set plate_failure variable to false
    plate_failure = False

    if not os.path.exists(os.path.join(proj_dir, 'card', plate_name)):
        # Create norm folder for plate if it doesn't exist already
        os.mkdir(os.path.join(proj_dir, 'card', plate_name))

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
        lvl4_card_jobs = {'ZSVC': [zscore.calculate_zscore, '_ZSVC.gct'],
                     'ZSPC': [functools.partial(zscore.calculate_zscore, plate_control=True), '_ZSPC.gct'],
                     'LFCPC': [functools.partial(viability.log_viability, plate_control=True, log=True), '_LFCPC.gct'],
                     'LFCVC': [functools.partial(viability.log_viability, plate_control=False,log=True), '_LFCVC.gct']}

    else:
        lvl4_card_jobs = {'ZSVC': [zscore.calculate_zscore, '_ZSVC.gct'],
                         'ZSPC': [functools.partial(zscore.calculate_zscore, plate_control=True), '_ZSPC.gct'],
                         'LFCPC': [functools.partial(viability.log_viability, plate_control=True, log=False), '_LFCPC.gct'],
                         'LFCVC': [functools.partial(viability.log_viability, plate_control=False, log=False), '_LFCVC.gct']}

    # Loop through this map to output all level 4 data
    for job in lvl4_card_jobs.keys():

        output_path = os.path.join(proj_dir, "card", plate_name, plate_name + lvl4_card_jobs[job][1])

        reader_writer(input_file=norm_path, output_file=output_path, function=lvl4_card_jobs[job][0])

    # Return status of plate failure
    return plate_failure


def main(args):
    # NB: automation sets project_dir to project_dir/prism_replicate_set_name to set up fs for s3

    if args.search_pattern:
        failure_list = []
        for folder in glob.glob(os.path.join(args.proj_dir, 'assemble', args.search_pattern)):
            name = os.path.basename(folder)
            logger.info("Carding {}".format(name))
            plate_failure = card(args.proj_dir, name, log_tf=args.log_tf, inv_tf=args.inv_tf, bad_wells=args.bad_wells, dp=args.no_invariants)
            if plate_failure == True:
                logger.debug("Carding failed for {}".format(name))
                failure_list.append(name)

        if len(failure_list) > 0:
            logger.info("Carding failed for following plates: {}".format(failure_list))
            pd.Series(failure_list).to_csv(os.path.join(args.proj_dir, 'failed_plates.txt'), sep='\t')
        else:
            logger.info("Carding succesfully completed on all plates")
            with open(os.path.join(args.proj_dir, "success.txt"), "w") as file: file.write("successfully processed all plates")

    else:
        failure = card(args.proj_dir, args.plate_name, log_tf=args.log_tf, inv_tf=args.inv_tf, bad_wells=args.bad_wells, dp=args.no_invariants)
        if failure:
            plate_failure_path = os.path.join(args.proj_dir, "card", args.plate_name, "failure.txt")
            with open(plate_failure_path, "w") as file:
                file.write("{} failed size checks".format(args.plate_name))
        else:
            success_path = os.path.join(args.proj_dir, "card", args.plate_name, "success.txt")
            with open(success_path, "w") as file:
                file.write("succesfully processed {}".format(args.plate_name))

if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    logger.debug("args:  {}".format(args))

    main(args)
