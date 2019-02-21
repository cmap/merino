import glob
import os
import cmapPy.pandasGEXpress.parse as pe
import cmapPy.pandasGEXpress.GCToo as GCToo
import prism_plots
import invariant_analysis as inv
import pandas as pd
import generic_heatmap as map
import ssmd_analysis as ssmd
import compound_strength as cp
import comp_strength_overview as comp
import merino.setup_logger as setup_logger
import logging
import expected_sensitivities as expected_sense
import count_heatmap as c_map
import argparse
import matplotlib.pyplot as plt
import numpy as np
import sc_plot
import sys
import seaborn as sns
import make_gallery as galleries
import merino.build_summary.build_summary as build_summary

logger = logging.getLogger(setup_logger.LOGGER_NAME)

def build_parser():

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # The following arguments are required. These are files that are necessary for assembly and which change
    # frequently between cohorts, replicates, etc.
    parser.add_argument("-project_folder", "-pf", help="path to the project directory for running QC on",
                        type=str, required=True)
    plates_group = parser.add_mutually_exclusive_group(required=True)
    plates_group.add_argument("-search_pattern", "-sp",
                              help="Search for this string in the directory, only run plates which contain it. ",
                              type=str, default=None)
    plates_group.add_argument("-plate_name", "-pn", help="name of individual plate to run on", type=str, default=None)
    parser.add_argument("-qc_folder", "-qc", help="string designating the prefix to each build file eg. PCAL075-126_T2B",
                        type=str, required=True)
    parser.add_argument("-invar", "-inv",
                        help="Flag to turn off invariant QC",
                        action="store_false")
    parser.add_argument("-verbose", '-v', help="Whether to print a bunch of output", action="store_true", default=False)
    parser.add_argument("-aggregate_out", "-agg", help="whether weave used aggregate_out flag", action="store_true",
                        default=False)


    return parser


def read_build_data(assemble_path, card_path):
    """

    Args:
        proj_dir: Folder containing all build files

    Returns: data_df: Dictionary linking an identifier of each data level to its corresponding GCToo object
             metadata_df: Dictionary of metadata with identifiers linked to pandas dataframes.

    """
    print "Reading Data"

    # Standard search patterns for a merino build
    assemble_data_pattern_list = ['*MEDIAN*.gct', '*COUNT*.gct']
    card_data_pattern_list =  ['*NORM.gct', '*_ZSPC.gct', '*_ZSPC.COMBAT.gct',
                    '*_ZSVC.gct', '*_LFCPC.gct', '*_LFCVC.gct']

    # Link identifiers to path to build folder
    data_search_patterns = [os.path.join(assemble_path, x) for x in assemble_data_pattern_list] + [os.path.join(card_path, y) for y in card_data_pattern_list]

    # Search for and reaad in each file in the build
    gcts = [build_summary.globandparse(x) for x in data_search_patterns]

    # Make dictionary and link each GCToo object to a corresponding key
    data_map = dict(zip(['mfi', 'count', 'norm', 'zspc', 'combat_zspc',
                'zsvc', 'lfcpc', 'lfcvc'], gcts))


    return data_map

def mk_distributions(data_map, project_name, out_dir):

    plot_map = {data_map['mfi']: ['MFI', 0, 50000], data_map['count']: ['BEAD', 0, 300],
              data_map['norm']: ['NORM', -10, 10],
              data_map['zspc']: ['ZSPC', -10, 10], data_map['zsvc']: ['ZSVC', -10, 10],
              data_map['lfcpc']: ['LFCPC', -10, 10], data_map['lfcvc']: ['LFCVC', -10, 10]}

    for df in plot_map.keys():


        prism_plots.data_distribution(data_df=df.data_df, xlabel=plot_map[df][0],
                                      title='Distribution for {} {}'.format(project_name, plot_map[df][0]),
                                      outfile=os.path.join(out_dir, 'distributions', 'histogram_{}.png'.format(plot_map[df][0])),
                                      xlim=[plot_map[df][1], plot_map[df][2]])

        prism_plots.stacked_heatmap(df=df.data_df, column_metadata=df.col_metadata_df,
                                    title='Median {} Across {}'.format(plot_map[df][0],project_name),
                                    outfile=os.path.join(out_dir, 'heatmaps', 'heatmap_{}.png'.format(plot_map[df][0])),
                                                         lims=[plot_map[df][1], plot_map[df][2]], reduce_upper_limit=True)


def plate_qc(proj_dir, out_dir, plate_name, invar=True, agg=True):

    if agg==True:
        card_path = os.path.join(proj_dir, plate_name, 'card')
        assemble_path = os.path.join(proj_dir, plate_name, 'assemble')

    else:
        card_path = os.path.join(proj_dir, 'card', plate_name)
        assemble_path = os.path.join(proj_dir, 'assemble', plate_name)


    plate_data_map = read_build_data(assemble_path, card_path)

    build_summary.mk_folders(out_dir, [plate_name])
    build_summary.mk_folders(os.path.join(out_dir, plate_name), ['invariants', 'distributions', 'heatmaps', 'ssmd', 'cp_strength'])

    if invar is True:
        for func in [inv.invariant_monotonicity, inv.invariant_curves_plot, inv.invariant_range_distributions]:
            func(plate_data_map['mfi'], plate_data_map['mfi'].col_metadata_df,
                 os.path.join(out_dir, plate_name, 'invariants'))

        inv.invariant_heatmap(plate_data_map['mfi'], os.path.join(out_dir, plate_name, 'invariants', 'inv_heatmap.png'),
                              lims=[0, 25000])

    mk_distributions(data_map=plate_data_map, project_name=plate_name,
                     out_dir=os.path.join(out_dir, plate_name))

    ssmd.norm_v_mfi_ssmd(plate_data_map['norm'], plate_data_map['mfi'], os.path.join(out_dir, plate_name, 'ssmd'))

    ssmd.ssmd_ecdf(plate_data_map['norm'], plate_data_map['mfi'], 'SSMD ECDF for {}'.format(plate_name),
                   os.path.join(out_dir, plate_name, 'ssmd'))

    cp.median_ZSPC_histogram(plate_data_map['zspc'], plate_data_map['zspc'].col_metadata_df,
                             os.path.join(out_dir, plate_name, 'cp_strength'), det=plate_name)

    cp.median_ZSPC_ecdf(plate_data_map['zspc'], plate_data_map['zspc'].col_metadata_df,
                        os.path.join(out_dir, plate_name, 'cp_strength'), det=plate_name)

    build_summary.qc_galleries(out_dir, plate_name)


def main(args):
    # NB: automation sets project_dir to project_dir/prism_replicate_set_name to set up fs for s3

    if args.search_pattern:
        for folder in glob.glob(os.path.join(args.project_folder, 'card', args.search_pattern)):
            name = os.path.basename(folder)
            logger.info("QCing {}".format(name))
            plate_qc(args.project_folder, args.qc_folder, name, invar=args.invar, agg = args.aggregate_out)
    else:

        plate_qc(args.proj_dir, args.qc_folder, args.plate_name, invar=args.invar, agg = args.aggregate_out)


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    logger.debug("args:  {}".format(args))

    main(args)