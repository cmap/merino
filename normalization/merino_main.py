import card
import weave
import merino.setup_logger as setup_logger
import logging
import argparse
import sys
import merino.mk_build_file as mk
import ConfigParser

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
    parser.add_argument("-input_folder", "-if",
                        help="The directory from which to produce level 5 data.",
                        type=str, default='ZSPC', required=False)

    return parser


def main(args):
    card.main(args)
    weave.main(args)
    args.input_folder = 'LFCPC'
    weave.main(args)
    mk.main(args)


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    logger.debug("args:  {}".format(args))

    main(args)


