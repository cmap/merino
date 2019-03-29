import os
import sys
import glob
import logging
import argparse
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import cmapPy.pandasGEXpress.parse as pe
import cmapPy.pandasGEXpress.GCToo as GCToo

import merino.setup_logger as setup_logger

import prism_plots
import invariant_analysis as inv
import ssmd_analysis as ssmd
import comp_strength_overview as comp
import make_gallery as galleries
import plate_summary

plt.rcParams['figure.figsize'] = (10.0,8.0)

logger = logging.getLogger(setup_logger.LOGGER_NAME)


data_mapper = dict(zip(['mfi', 'count', 'norm', 'zspc', 'combat_zspc',
                'zsvc', 'lfcpc', 'lfcvc', 'modz', 'combat_modz'], ['*MFI*.gctx', '*COUNT*.gctx', '*NORM*.gctx',
                                                                   '*_ZSPC_*.gctx', '*_ZSPC.COMBAT_*.gctx',
                    '*_ZSVC_*.gctx', '*4_LFCPC_*.gctx', '*4_LFCVC_*.gctx', '*MODZ.ZSPC_*.gctx', '*MODZ.ZSPC.COMBAT_*.gctx']))


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
    parser.add_argument("-plate_qc", "-pq",
                        help="Perform QC at a plate level",
                        action="store_true")
    parser.add_argument("-data_types", "-dt",
                        help="Comma separated list of data types to use. Valid options are any combination of"
                             ": mfi,count,norm,zspc,combat_zspc,zsvc,lfcpc,lfcvc,modz,combat_modz", type=str,
                        default='mfi,count,norm,zspc,combat_zspc,zsvc,lfcpc,lfcvc,modz,combat_modz')
    parser.add_argument("-verbose", '-v', help="Whether to print a bunch of output", action="store_true", default=False)


    return parser

def globandparse(search_pattern):
    """
    Search for inputs based on a search pattern and read in file
    Args:
        search_pattern: Identifier to find a given build file using glob search

    Returns: GCToo object found through glob search and read in

    """
    logger.info("globbing for : {}".format(search_pattern))
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

    keep_data_type = args.data_types.split(',')
    data_mapper_subset = dict((k, data_mapper[k]) for k in keep_data_type)

    # Standard search patterns for a merino build
    data_pattern_list = [data_mapper_subset[x] for x in data_mapper_subset.keys()]

    # Link identifiers to path to build folder
    data_search_patterns = [os.path.join(proj_dir, x) for x in data_pattern_list]

    # Search for and reaad in each file in the build
    gcts = [globandparse(x) for x in data_search_patterns]

    # Make dictionary and link each GCToo object to a corresponding key
    data_map = dict(zip(data_mapper_subset.keys(), gcts))

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
    plot_map = {'mfi': ['MFI', 0, 50000], 'count': ['BEAD', 0, 300],
              'norm': ['NORM', -10, 10],
              'zspc': ['ZSPC', -10, 10], 'zsvc': ['ZSVC', -10, 10],
              'lfcpc': ['LFCPC', -10, 10], 'lfcvc': ['LFCVC', -10, 10],
              'modz': ['MODZ', -10, 10], 'combat_modz': ['COMBAT MODZ', -10, 10],
              'combat_zspc': ['COMBAT ZSPC', -10, 10]}

    plot_map_subset = dict((data_map[k], plot_map[k]) for k in data_map.keys())

    for df in plot_map_subset.keys():

        if 'MODZ' in plot_map_subset[df][0]:
            meta = metadata_map['sig']
        else:
            meta = metadata_map['inst']

        prism_plots.data_distribution(data_df=df.data_df, xlabel=plot_map_subset[df][0],
                                      title='Distribution for {} {}'.format(project_name, plot_map_subset[df][0]),
                                      outfile=os.path.join(out_dir, 'distributions', 'histogram_{}.png'.format(plot_map_subset[df][0])),
                                      xlim=[plot_map_subset[df][1], plot_map_subset[df][2]])

        prism_plots.stacked_heatmap(df=df.data_df, column_metadata=meta.loc[df.data_df.columns],
                                    title='Median {} Across {}'.format(plot_map_subset[df][0],project_name),
                                    outfile=os.path.join(out_dir, 'heatmaps', 'heatmap_{}.png'.format(plot_map_subset[df][0])),
                                                         lims=[plot_map_subset[df][1], plot_map_subset[df][2]], reduce_upper_limit=True)



def get_plate_qc_data_map_and_run(data_map, metadata_map, norm_cell_metadata, project_name, out_dir, invar):

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
        plate_summary.plate_qc(out_dir, plate, plate_data_map, invar=invar)


def qc_galleries(out_dir, proj_name, metadata_map, data_map):
    local_paths = glob.glob(os.path.join(out_dir, proj_name + '*', '*.html'))
    dex = [os.path.basename(os.path.dirname(x)) for x in local_paths]

    ssmd_medians = metadata_map['ssmd'].median().loc[dex]
    ssmd_failures = metadata_map['ssmd'][metadata_map['ssmd'] < 2].count().loc[dex]
    ssmd_pct_failure = (metadata_map['ssmd'][metadata_map['ssmd'] < 2].count().loc[dex] / metadata_map['ssmd'].shape[0]) * 100
    plate_shapes = []
    well_dropouts = []
    signal_strengths = []
    sig_stength_75 = []
    correlations = []
    unique_perts = []
    n_kills = []
    for plate in dex:
        temp_sig = metadata_map['sig'].loc[
            [x for x in metadata_map['sig'].index if x.startswith(plate.rsplit('_', 2)[0])]]
        temp_inst = metadata_map['inst'].loc[[x for x in metadata_map['inst'].index if x.startswith(plate)]]
        temp_norm = data_map['norm'].data_df.loc[:,[x for x in data_map['norm'].data_df.columns if x.startswith(plate)]]
        plate_shapes.append(temp_norm.shape[1])
        well_dropouts.append(384 - temp_norm.shape[1])
        signal_strengths.append(temp_sig['ss_ltn2'].median())
        sig_stength_75.append(temp_sig['ss_ltn2'].quantile(.75))
        n_kills.append(temp_sig['ss_ltn2'][temp_sig['ss_ltn2'] > (len(temp_sig['ss_ltn2']) / 20)].count())
        correlations.append(temp_sig['cc_q75'].median())
        unique_perts.append(len(temp_inst.loc[temp_inst['pert_type'] == 'trt_cp', 'pert_id'].unique()))


    def make_url(ref, name):
        ref = os.path.relpath(ref, out_dir)
        return '<a target="_blank" href="{}">{}</a>'.format(ref, name)

    premadelinks = [make_url(x, dex[i]) for i, x in enumerate(local_paths)]

    headers = ['plate','median SSMD', 'n SSMD failures', 'pct SSMD failures', 'n wells',
               'n dropouts', 'median signal strength', 'q75 signal strength', 'n_active_wells','median ccQ75', 'n unique perts']

    index_info = zip(premadelinks, ssmd_medians,
                     ssmd_failures, ssmd_pct_failure, plate_shapes,
                     well_dropouts, signal_strengths, sig_stength_75, n_kills, correlations,
                     unique_perts)
    #print index_info
    made_gallery = galleries.mk_index(table_headers=headers, table_tuples=index_info, outfolder=out_dir, project_name=proj_name)
    if made_gallery:
        logger.info("successfully made QC gallery")

def main(args, proj_dir, out_dir, project_name,invar=True):

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

    # Make standard SC plot for whole dataset, signal strength vs correlation
    prism_plots.sc_plot(metadata_map['sig'], os.path.join(out_dir,'sc_modz.zspc.png'))

    # Make modz distribuions split by pert type
    if 'combat_modz' in data_map.keys():
        comp.modz_dist(data_map['combat_modz'], metadata_map['cb_sig'], [], os.path.join(out_dir, 'modz_dist.png'))

    # If running on data with control barcodes, plot monotonicity of curves
    #if invar is True:
    #    inv.invariant_monotonicity(data_map['mfi'], metadata_map['inst'], out_dir)

    # Calculate median SSMD by pool and output in table
    ssmd.ssmd_by_pool(metadata_map['ssmd'], metadata_map['cell'], out_dir)

    # Get cell metadata without control barcodes for later use
    norm_cell_metadata = metadata_map['cell'].loc[[x for x in metadata_map['cell'].index if x in data_map['norm'].row_metadata_df.index]]

    # ECDF of SSMD Scores in Norm Data and MFI Data
    if 'norm' and 'mfi' in data_map.keys():
        ssmd.ssmd_ecdf(GCToo.GCToo(data_df=data_map['norm'].data_df, col_metadata_df=metadata_map['inst'].loc[data_map['norm'].data_df.columns],
                               row_metadata_df=norm_cell_metadata),
                   GCToo.GCToo(data_df=data_map['mfi'].data_df,
                               col_metadata_df=metadata_map['inst'].loc[data_map['mfi'].data_df.columns],
                               row_metadata_df=metadata_map['cell']), 'SSMD ECDF for {}'.format(os.path.dirname(proj_dir))
                   ,os.path.join(out_dir))

    # Make a bunch of plots at the plate level for each plate in cohort
    if args.plate_qc:
        get_plate_qc_data_map_and_run(data_map, metadata_map, norm_cell_metadata, project_name, out_dir, invar)
    #todo: add a check for get_plate_qc_data_map_and_run before running qc_galleries --> dependent
        qc_galleries(out_dir, project_name, metadata_map, data_map)


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    logger.debug("args:  {}".format(args))

    main(args, proj_dir=args.build_folder, out_dir=args.qc_folder, project_name=args.project_name, invar=args.invar)





