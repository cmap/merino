import os
import sys
import glob
import logging
import argparse

import merino.setup_logger as setup_logger

import build_summary as build_summary
import invariant_analysis as inv
import ssmd_analysis as ssmd
import compound_strength as cp
import make_gallery as galleries
import prism_plots
import pandas as pd

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
    if len(glob.glob(card_path + '/*_ZSPC.COMBAT.gct')) > 0:
        card_data_pattern_list =  ['*NORM.gct', '*_ZSPC.gct',
                        '*_ZSVC.gct', '*_LFCPC.gct', '*_LFCVC.gct', '*_ZSPC.COMBAT.gct']

    else:
        card_data_pattern_list = ['*NORM.gct', '*_ZSPC.gct',
                                  '*_ZSVC.gct', '*_LFCPC.gct', '*_LFCVC.gct']

    # Link identifiers to path to build folder
    data_search_patterns = [os.path.join(assemble_path, x) for x in assemble_data_pattern_list] + [os.path.join(card_path, y) for y in card_data_pattern_list]

    # Search for and reaad in each file in the build
    gcts = [build_summary.globandparse(x) for x in data_search_patterns]

    # Make dictionary and link each GCToo object to a corresponding key
    if '*_ZSPC.COMBAT.gct' in card_data_pattern_list:
        data_map = dict(zip(['mfi', 'count', 'norm', 'zspc',
                'zsvc', 'lfcpc', 'lfcvc', 'zspc.combat'], gcts))
    else:
        data_map = dict(zip(['mfi', 'count', 'norm', 'zspc',
                             'zsvc', 'lfcpc', 'lfcvc'], gcts))


    return data_map

def mk_distributions(data_map, project_name, out_dir):

    if 'zspc.combat' in  data_map.keys():
        plot_map = {data_map['mfi']: ['MFI', 0, 50000], data_map['count']: ['BEAD', 0, 300],
              data_map['norm']: ['NORM', -10, 10],
              data_map['zspc']: ['ZSPC', -10, 10], data_map['zsvc']: ['ZSVC', -10, 10],
              data_map['lfcpc']: ['LFCPC', -10, 10], data_map['lfcvc']: ['LFCVC', -10, 10],
                    data_map['zspc.combat']: ['ZSPC.COMBAT', -10, 10]}

    else:
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


def get_plate_qc_data_map(proj_dir, plate_name):
    card_path = os.path.join(proj_dir, 'card', plate_name)
    assemble_path = os.path.join(proj_dir, 'assemble', plate_name)

    plate_data_map = read_build_data(assemble_path, card_path)
    return plate_data_map

def plate_qc(out_dir, plate_name, plate_data_map, invar=True):

    build_summary.mk_folders(out_dir, [plate_name])
    build_summary.mk_folders(os.path.join(out_dir, plate_name), ['invariants', 'distributions', 'heatmaps', 'ssmd', 'cp_strength'])

    mk_distributions(data_map=plate_data_map, project_name=plate_name,
                     out_dir=os.path.join(out_dir, plate_name))

    if invar is True:
        for func in [inv.invariant_monotonicity, inv.invariant_curves_plot]:
            func(plate_data_map['mfi'], plate_data_map['mfi'].col_metadata_df,os.path.join(out_dir, plate_name, 'invariants'))

        inv.invariant_heatmap(plate_data_map['mfi'], os.path.join(out_dir, plate_name, 'invariants', 'inv_heatmap.png'),
                              lims=[0, 15000])



    ssmd.norm_v_mfi_ssmd(plate_data_map['norm'], plate_data_map['mfi'], os.path.join(out_dir, plate_name, 'ssmd'))

    ssmd.ssmd_ecdf(plate_data_map['norm'], plate_data_map['mfi'], 'SSMD ECDF for {}'.format(plate_name),
                   os.path.join(out_dir, plate_name, 'ssmd'))

    cp.median_ZSPC_histogram(plate_data_map['zspc'], plate_data_map['zspc'].col_metadata_df,
                             os.path.join(out_dir, plate_name, 'cp_strength'), det=plate_name)

    cp.median_ZSPC_ecdf(plate_data_map['zspc'], plate_data_map['zspc'].col_metadata_df,
                        os.path.join(out_dir, plate_name, 'cp_strength'), det=plate_name)

    make_gallery(out_dir, plate_name)

    mk_report(out_dir, plate_name, plate_data_map)

def make_gallery(qc_dir, plate_name):
    plate_qc_dir = os.path.join(qc_dir, plate_name)

    images = ['invariants/inv_heatmap.png', 'invariants/invariant_curves.png', 'invariants/invariant_mono.png',
              'ssmd/SSMD_ECDF.png', 'ssmd/NORMvMFI_SSMD_Boxplot.png', 'distributions/histogram_BEAD.png',
              'heatmaps/heatmap_BEAD.png', 'heatmaps/heatmap_MFI.png', 'heatmaps/heatmap_NORM.png',
              'heatmaps/heatmap_ZSCORE.png']
    outfile = os.path.join(plate_qc_dir, 'gallery.html')
    galleries.mk_gal(images, outfile)

def mk_report(qc_dir, plate_name, plate_data_map):
    ssmds = ssmd.get_ssmd(plate_data_map['norm'])
    ssmd_median = ssmds.median()
    ssmd_failures = ssmds[ssmds < 2].count()
    ssmd_pct_failure = (float(ssmds[ssmds < 2].count()) / len(ssmds)) * 100
    plate_shapes = plate_data_map['norm'].data_df.shape[1]
    well_dropouts = 384 - plate_data_map['norm'].data_df.shape[1]
    ss_ltn2 = plate_data_map['zspc'].data_df[plate_data_map['zspc'].data_df < -2].count()
    n_active = ss_ltn2[ss_ltn2 > (len(ss_ltn2) / 20)].count()
    unique_perts = len(plate_data_map['norm'].col_metadata_df.loc[plate_data_map['norm'].col_metadata_df['pert_type'] == 'trt_cp', 'pert_id'].unique())
    invariants = ['c-' + str(x) for x in range(661,671)]
    median_invariant = plate_data_map['mfi'].data_df.loc[[x for x in invariants if x in plate_data_map['mfi'].data_df.index]].median(axis=1).median()

    headers = ['plate','median SSMD', 'n SSMD failures', 'pct SSMD failures', 'median_invariant','n wells', 'n dropouts',  'n active wells','n unique perts']
    report = pd.DataFrame(
        [plate_name, ssmd_median, ssmd_failures, ssmd_pct_failure, median_invariant,plate_shapes, well_dropouts, n_active,
         unique_perts]).T
    report.columns = headers

    report.to_csv(os.path.join(qc_dir, plate_name, 'report.txt'), sep='\t', index=False)


def main(args):
    # NB: automation sets project_dir to project_dir/prism_replicate_set_name to set up fs for s3

    if args.search_pattern:
        for folder in glob.glob(os.path.join(args.project_folder, 'card', args.search_pattern)):
            name = os.path.basename(folder)
            logger.info("QCing {}".format(name))
            plate_data = get_plate_qc_data_map(args.project_folder, name)
            plate_qc(args.qc_folder, name, plate_data, invar=args.invar)
    else:
        plate_data = get_plate_qc_data_map(args.project_folder, args.plate_name)
        plate_qc(args.qc_folder, args.plate_name, plate_data, invar=args.invar)


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    logger.debug("args:  {}".format(args))

    main(args)