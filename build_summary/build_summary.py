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
plt.rcParams['figure.figsize'] = (10.0,8.0)

logger = logging.getLogger(setup_logger.LOGGER_NAME)

def build_parser():

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # The following arguments are required. These are files that are necessary for assembly and which change
    # frequently between cohorts, replicates, etc.
    parser.add_argument("-build_folder", "-bf", help="path to the pod directory you want to run card on",
                        type=str, required=True)
    parser.add_argument("-qc_folder", "-qc", help="string designating the prefix to each build file eg. PCAL075-126_T2B",
                        type=str, required=True)
    parser.add_argument("-project_name", "-pn",
                        help="Code for project eg. PCAL",
                        type=str, required=False)
    parser.add_argument("-invar", "-inv",
                        help="Flag to turn off invariant QC",
                        action="store_false")
    parser.add_argument("-sensitivities", "-sense",
                        help="Perform expected sensitivity analysis (time consuming)",
                        action="store_true")
    parser.add_argument("-verbose", '-v', help="Whether to print a bunch of output", action="store_true", default=False)


    return parser

def globandparse(search_pattern):
    """
    Search for inputs based on a search pattern and read in file
    Args:
        search_pattern: Identifier to find a given build file using glob search

    Returns: GCToo object found through glob search and read in

    """
    path = glob.glob(search_pattern)[0]
    gct = pe.parse(path)
    return gct


def read_build_data(proj_dir):
    """

    Args:
        proj_dir: Folder containing all build files

    Returns: data_df: Dictionary linking an identifier of each data level to its corresponding GCToo object
             metadata_df: Dictionary of metadata with identifiers linked to pandas dataframes.

    """
    print "Reading Data"

    # Standard search patterns for a merino build
    data_pattern_list = ['*MFI*.gctx', '*COUNT*.gctx', '*NORM*.gctx', '*_ZSPC_*.gctx', '*_ZSPC.COMBAT_*.gctx',
                    '*_ZSVC_*.gctx', '*4_LFCPC_*.gctx', '*4_LFCVC_*.gctx', '*MODZ.ZSPC_*.gctx', '*MODZ.ZSPC.COMBAT_*.gctx']

    # Link identifiers to path to build folder
    data_search_patterns = [os.path.join(proj_dir, x) for x in data_pattern_list]

    # Search for and reaad in each file in the build
    gcts = [globandparse(x) for x in data_search_patterns]

    # Make dictionary and link each GCToo object to a corresponding key
    data_map = dict(zip(['mfi', 'count', 'norm', 'zspc', 'combat_zspc',
                'zsvc', 'lfcpc', 'lfcvc', 'modz', 'combat_modz'], gcts))

    print "Reading Metadata"
    # Read in each metadata file
    inst_info = pd.read_table(glob.glob(os.path.join(proj_dir, '*inst_info*.txt'))[0], index_col='profile_id')
    sig_info = pd.read_table(glob.glob(os.path.join(proj_dir, '*sig_metrics_MODZ.ZSPC.txt'))[0], index_col='sig_id')
    cb_sig_info = pd.read_table(glob.glob(os.path.join(proj_dir, '*sig_metrics_MODZ.ZSPC.COMBAT.txt'))[0], index_col='sig_id')

    cell_info = pd.read_table(glob.glob(os.path.join(proj_dir, '*cell_info*.txt'))[0], index_col='rid')
    ssmd_info = pe.parse(glob.glob(os.path.join(proj_dir, '*ssmd*.gct'))[0]).data_df

    # Make metadata map
    metadata_map = {'inst': inst_info, 'sig': sig_info, 'cell': cell_info, 'ssmd': ssmd_info, 'cb_sig': cb_sig_info}

    return data_map, metadata_map


def mk_folders(out_dir, folders):

    for folder in folders:
        if not os.path.exists(os.path.join(out_dir, folder)):
            os.mkdir(os.path.join(out_dir, folder))


def mk_distributions(data_map, metadata_map,project_name, out_dir):

    plot_map = {data_map['mfi']: ['MFI', 0, 50000], data_map['count']: ['BEAD', 0, 300],
              data_map['norm']: ['NORM', -10, 10],
              data_map['zspc']: ['ZSPC', -10, 10], data_map['zsvc']: ['ZSVC', -10, 10],
              data_map['lfcpc']: ['LFCPC', -10, 10], data_map['lfcvc']: ['LFCVC', -10, 10],
              data_map['modz']: ['MODZ', -10, 10], data_map['combat_modz']: ['COMBAT MODZ', -10, 10],
              data_map['combat_zspc']: ['COMBAT ZSPC', -10, 10]}

    for df in plot_map.keys():

        if 'MODZ' in plot_map[df][0]:
            meta = metadata_map['sig']
        else:
            meta = metadata_map['inst']

        prism_plots.data_distribution(data_df=df.data_df, xlabel=plot_map[df][0],
                                      title='Distribution for {} {}'.format(project_name, plot_map[df][0]),
                                      outfile=os.path.join(out_dir, 'distributions', 'histogram_{}.png'.format(plot_map[df][0])),
                                      xlim=[plot_map[df][1], plot_map[df][2]])

        prism_plots.stacked_heatmap(df=df.data_df, column_metadata=meta.loc[df.data_df.columns],
                                    title='Median {} Across {}'.format(plot_map[df][0],project_name),
                                    outfile=os.path.join(out_dir, 'heatmaps', 'heatmap_{}.png'.format(plot_map[df][0])),
                                                         lims=[plot_map[df][1], plot_map[df][2]], reduce_upper_limit=True)


def plate_qc(data_map, metadata_map, norm_cell_metadata, project_name, out_dir, invar):

    for plate in metadata_map['inst']['prism_replicate'].unique():
        print plate
        plate_data_map = {}
        for key in data_map:
            if key == 'count' or key =='mfi':
                plate_data_map[key] = GCToo.GCToo(data_df=data_map[key].data_df[metadata_map['inst'][metadata_map['inst']['prism_replicate'] == plate].index],
                                col_metadata_df=metadata_map['inst'][metadata_map['inst']['prism_replicate'] == plate],
                                row_metadata_df=metadata_map['cell'])
            else:
                plate_data_map[key]= GCToo.GCToo(data_df=data_map[key].data_df.loc[:, [x for x in metadata_map['inst'][
                    metadata_map['inst']['prism_replicate'] == plate].index if x in data_map[key].data_df.columns]],
                                         col_metadata_df=metadata_map['inst'].loc[[x for x in metadata_map['inst'][
                                             metadata_map['inst']['prism_replicate'] == plate].index if x in data_map[key].data_df.columns]],
                                         row_metadata_df=norm_cell_metadata)

        print plate_data_map['norm'].data_df.shape

        mk_folders(out_dir, [plate])
        mk_folders(os.path.join(out_dir, plate), ['invariants', 'distributions', 'heatmaps', 'ssmd', 'cp_strength'])

        if invar is True:
            for func in [inv.invariant_monotonicity, inv.invariant_curves_plot, inv.invariant_range_distributions]:
                func(plate_data_map['mfi'], plate_data_map['mfi'].col_metadata_df, os.path.join(out_dir, plate, 'invariants'))

            inv.invariant_heatmap(plate_data_map['mfi'], os.path.join(out_dir, plate, 'invariants', 'inv_heatmap.png'), lims=[0,25000])

        mk_distributions(data_map=plate_data_map,  metadata_map=metadata_map, project_name=project_name,
                         out_dir=os.path.join(out_dir, plate))

        ssmd.norm_v_mfi_ssmd(plate_data_map['norm'], plate_data_map['mfi'], os.path.join(out_dir, plate, 'ssmd'))

        ssmd.ssmd_ecdf(plate_data_map['norm'], plate_data_map['mfi'], 'SSMD ECDF for {}'.format(plate),os.path.join(out_dir, plate, 'ssmd'))

        cp.median_ZSPC_histogram(plate_data_map['zspc'], plate_data_map['zspc'].col_metadata_df,
                                 os.path.join(out_dir, plate, 'cp_strength'), det=plate)

        cp.median_ZSPC_ecdf(plate_data_map['zspc'], plate_data_map['zspc'].col_metadata_df,
                                 os.path.join(out_dir, plate, 'cp_strength'), det=plate)


def qc_galleries(proj_dir, proj_name):
    for x in glob.glob(os.path.join(proj_dir, proj_name + '*')):
        print x
        images = ['invariants/inv_heatmap.png', 'invariants/invariant_curves.png', 'invariants/invariant_mono.png',
                  'ssmd/SSMD_ECDF.png', 'ssmd/NORMvMFI_SSMD_Boxplot.png', 'distributions/histogram_BEAD.png',
                  'heatmaps/heatmap_BEAD.png', 'heatmaps/heatmap_MFI.png', 'heatmaps/heatmap_NORM.png',
                  'heatmaps/heatmap_ZSCORE.png']
        outfolder = os.path.join(x, 'gallery.html')
        galleries.mk_gal(images, outfolder)


def main(proj_dir, out_dir, project_name, invar=True, sense=False):
    # Read in the data
    data_map, metadata_map = read_build_data(proj_dir=proj_dir)

    # Make folders for different outputs
    mk_folders(out_dir=out_dir, folders=['sensitivities', 'distributions', 'heatmaps'])

    # Check if project name arg is filled, if not use base folder name
    if project_name is None:
        project_name = os.path.basename(os.path.dirname(proj_dir))

    # Make distributions and heatmaps of all data at each data level
    mk_distributions(data_map, metadata_map, project_name, out_dir)

    # If expected sensitivities arg is set to true, run expected sensitivities analysis
    # TODO add argument for defining sensitivity cell set
    if sense is True:
        expected_sense.wtks(data_map['combat_modz'], metadata_map['sig'], os.path.join(out_dir, 'sensitivities'))

    # Make standard SC plot for whole dataset, signal strength vs correlation
    prism_plots.sc_plot(metadata_map['sig'], os.path.join(out_dir,'sc_modz.zspc.png'))

    # Make modz distribuions split by pert type
    comp.modz_dist(data_map['combat_modz'], metadata_map['cb_sig'], [], os.path.join(out_dir, 'modz_dist.png'))

    # If running on data with control barcodes, plot monotonicity of curves
    if invar is True:
        inv.invariant_monotonicity(data_map['mfi'], metadata_map['inst'], out_dir)

    # Calculate median SSMD by pool and output in table
    ssmd.ssmd_by_pool(metadata_map['ssmd'], metadata_map['cell'], out_dir)

    # Get cell metadata without control barcodes for later use
    norm_cell_metadata = metadata_map['cell'].loc[[x for x in metadata_map['cell'].index if x in data_map['norm'].row_metadata_df.index]]

    # ECDF of SSMD Scores in Norm Data and MFI Data
    ssmd.ssmd_ecdf(GCToo.GCToo(data_df=data_map['norm'].data_df, col_metadata_df=metadata_map['inst'].loc[data_map['norm'].data_df.columns],
                               row_metadata_df=norm_cell_metadata),
                   GCToo.GCToo(data_df=data_map['mfi'].data_df,
                               col_metadata_df=metadata_map['inst'].loc[data_map['mfi'].data_df.columns],
                               row_metadata_df=metadata_map['cell']), 'SSMD ECDF for {}'.format(os.path.dirname(proj_dir))
                   ,os.path.join(out_dir))

    # Make a bunch of plots at the plate level for each plate in cohort
    plate_qc(data_map, metadata_map, norm_cell_metadata, project_name, out_dir, invar)

    # Put plate level plots into html galleries

    qc_galleries(out_dir, project_name)



if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    logger.debug("args:  {}".format(args))

    main(proj_dir=args.build_folder, out_dir=args.qc_folder, project_name=args.project_name, invar=args.invar, sense=args.sensitivities)





