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

logger = logging.getLogger(setup_logger.LOGGER_NAME)

def build_parser():

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # The following arguments are required. These are files that are necessary for assembly and which change
    # frequently between cohorts, replicates, etc.
    parser.add_argument("-build_folder", "-bf", help="path to the pod directory you want to run card on",
                        type=str, required=True)
    parser.add_argument("-qc_folder", "-qc", help="string designating the prefix to each build file eg. PCAL075-126_T2B",
                        type=str, required=True)
    parser.add_argument("-invar", "-inv",
                        help="Drop sigs from modZ with less than one profile",
                        action="store_false")
    parser.add_argument("-sensitivities", "-sense",
                        help="Perform expected sensitivity analysis (time consuming)",
                        action="store_true")
    parser.add_argument("-verbose", '-v', help="Whether to print a bunch of output", action="store_true", default=False)


    return parser

def globandparse(search_pattern):
    path = glob.glob(search_pattern)[0]
    gct = pe.parse(path)
    return gct


def read_build_data(proj_dir):
    print "Reading LEVEL2 Data"

    data_pattern_list = ['*MFI*.gctx', '*COUNT*.gctx', '*NORM*.gctx', '*_ZSPC_*.gctx', '*_ZSPC.COMBAT_*.gctx',
                    '*_ZSVC_*.gctx', '*4_LFCPC_*.gctx', '*4_LFCVC_*.gctx', '*MODZ.ZSPC_*.gctx', '*MODZ.ZSPC.COMBAT_*.gctx']

    data_search_patterns = [os.path.join(proj_dir, x) for x in data_pattern_list]

    gcts = [globandparse(x) for x in data_search_patterns]

    data_map = dict(zip(['mfi', 'count', 'norm', 'zspc', 'combat_zspc',
                'zsvc', 'lfcpc', 'lfcvc', 'modz', 'combat_modz'], gcts))


    print "Reading Metadata"
    inst_info = pd.read_table(glob.glob(os.path.join(proj_dir, '*inst_info*.txt'))[0], index_col='profile_id')
    sig_info = pd.read_table(glob.glob(os.path.join(proj_dir, '*sig_metrics_MODZ.ZSPC.txt'))[0], index_col='sig_id')
    cb_sig_info = pd.read_table(glob.glob(os.path.join(proj_dir, '*sig_metrics_MODZ.ZSPC.COMBAT.txt'))[0], index_col='sig_id')

    cell_info = pd.read_table(glob.glob(os.path.join(proj_dir, '*cell_info*.txt'))[0], index_col='rid')
    ssmd_info = pe.parse(glob.glob(os.path.join(proj_dir, '*ssmd*.gct'))[0]).data_df

    metadata_map = {'inst': inst_info, 'sig': sig_info, 'cell': cell_info, 'ssmd': ssmd_info, 'cb_sig': cb_sig_info}

    return data_map, metadata_map

def mk_folders(out_dir):

    if not os.path.exists(os.path.join(out_dir, 'sensitivities')):
        os.mkdir(os.path.join(out_dir, 'sensitivities'))

    if not os.path.exists(os.path.join(out_dir, 'distributions')):
        os.mkdir(os.path.join(out_dir, 'distributions'))


def mk_distributions(data_map, metadata_map,project_name, out_dir):

    plot_map = {data_map['mfi']: ['MFI', 0, 50000], data_map['count']: ['Bead Count', 0, 220],
              data_map['norm']: ['NORM', -10, 6],
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
                                      outfile=os.path.join(out_dir, 'distributions', 'histogram_{}_{}.png'.format(project_name,plot_map[df][0])),
                                      xlim=[plot_map[df][1], plot_map[df][2]])

        prism_plots.stacked_heatmap(df=df.data_df, column_metadata=meta.loc[df.data_df.columns],
                                    title='Median {} Across {}'.format(plot_map[df][0],project_name),
                                    outfile=os.path.join(out_dir, 'distributions', 'heatmap_{}_{}.png'.format(project_name,plot_map[df][0])),
                                                         lims=[plot_map[df][1], plot_map[df][2]])


def main(proj_dir, out_dir, invar=True, sense=False):
    # READ EM ALL IN

    data_map, metadata_map = read_build_data(proj_dir=proj_dir)

    mk_folders(out_dir=out_dir)

    project_name = os.path.basename(proj_dir)

    mk_distributions(data_map, metadata_map, project_name, out_dir)

    if sense == True:
        expected_sense.wtks(data_map['combat_modz'], metadata_map['sig'], os.path.join(out_dir, 'sensitivities'))

    prism_plots.sc_plot(metadata_map['sig'], os.path.join(out_dir,'sc_modz.zspc.png'))

    comp.modz_dist(data_map['combat_modz'], metadata_map['cb_sig'], [], os.path.join(out_dir, 'modz_dist.png'))

    if invar==True:
        inv.invariant_monotonicity(data_map['mfi'], metadata_map['inst'], out_dir)

    ssmd.ssmd_by_pool(metadata_map['ssmd'], metadata_map['cell'], out_dir)

    ssmd.ssmd_ecdf(GCToo.GCToo(data_df=data_map['norm'], col_metadata_df=metadata_map['inst'].loc[data_map['norm'].data_df.columns],
                               row_metadata_df=data_map['norm'].row_metadata_df),
                   GCToo.GCToo(data_df=data_map['mfi'].data_df, col_metadata_df=metadata_map['inst'].loc[data_map['mfi'].data_df.columns],
                               row_metadata_df=data_map['mfi'].row_metadata_df), 'SSMD ECDF for {}'.format(os.path.dirname(proj_dir))
                   ,os.path.join(out_dir))





    for plate in metadata_map['inst']['prism_replicate'].unique():

        print plate


        if plate in ['CBRANT008_KJ100.48H_X5']:
            continue

        plate_mfi = GCToo.GCToo(data_df=data_map['mfi'].data_df[metadata_map['inst'][metadata_map['inst']['prism_replicate'] == plate].index],
                                col_metadata_df=metadata_map['inst'][metadata_map['inst']['prism_replicate'] == plate],
                                row_metadata_df=data_map['mfi'].row_metadata_df)
        print plate_mfi.data_df.shape

        plate_count = GCToo.GCToo(data_df=data_map['count'].data_df[metadata_map['inst'][metadata_map['inst']['prism_replicate'] == plate].index],
                                col_metadata_df=metadata_map['inst'][metadata_map['inst']['prism_replicate'] == plate],
                                row_metadata_df=data_map['count'].row_metadata_df)

        plate_norm = GCToo.GCToo(data_df=data_map['norm'].data_df.loc[:,[x for x in metadata_map['inst'][metadata_map['inst']['prism_replicate'] == plate].index if x in data_map['norm'].data_df.columns]],
                                col_metadata_df=metadata_map['inst'].loc[[x for x in metadata_map['inst'][metadata_map['inst']['prism_replicate'] == plate].index if x in data_map['norm'].data_df.columns]],
                                row_metadata_df=data_map['norm'].row_metadata_df)

        plate_zscore = GCToo.GCToo(data_df=data_map['zspc'].data_df.loc[:,[x for x in metadata_map['inst'][metadata_map['inst']['prism_replicate'] == plate].index if x in data_map['norm'].data_df.columns]],
                                col_metadata_df=metadata_map['inst'].loc[[x for x in metadata_map['inst'][metadata_map['inst']['prism_replicate'] == plate].index if x in data_map['norm'].data_df.columns]],
                                row_metadata_df=data_map['zspc'].row_metadata_df)

        plate_viability = GCToo.GCToo(data_df=data_map['lfcpc'].data_df.loc[:,[x for x in metadata_map['inst'][metadata_map['inst']['prism_replicate'] == plate].index if x in data_map['norm'].data_df.columns]],
                                col_metadata_df=metadata_map['inst'].loc[[x for x in metadata_map['inst'][metadata_map['inst']['prism_replicate'] == plate].index if x in data_map['norm'].data_df.columns]],
                                row_metadata_df=data_map['zspc'].row_metadata_df)

        if not os.path.exists(os.path.join(out_dir, plate)):
            os.mkdir(os.path.join(out_dir, plate))
            for fold in ['invariants', 'distributions', 'heatmaps', 'ssmd', 'cp_strength']:
                os.mkdir(os.path.join(out_dir, plate, fold))

       # dist.distributions(plate_norm, plate_mfi, plate_count, plate_zscore, plate_viability,
        #                   plate_norm.col_metadata_df, os.path.join(out_dir, plate, 'distributions'))


        if invar == True:
            for func in [inv.invariant_monotonicity, inv.invariant_curves_plot, inv.invariant_range_distributions]:
                try:
                    func(plate_mfi, plate_mfi.col_metadata_df, os.path.join(out_dir, plate, 'invariants'))
                except:
                    print 'skipping function'


            inv.invariant_heatmap(plate_mfi, os.path.join(out_dir, plate, 'invariants', 'inv_heatmap.png'), lims=[0,25000])

        map.mk_heatmap(plate_norm.data_df, 'Heatmap of Median Norm Values',
                       os.path.join(out_dir, plate,'heatmaps', 'NORM.png'), lims=[-5,5])

        map.mk_heatmap(plate_mfi.data_df, 'Heatmap of Median MFI Values',
                       os.path.join(out_dir, plate, 'heatmaps', 'MFI.png'), lims=[0,20000], colormap='Reds')

        map.mk_heatmap(plate_count.data_df, 'Heatmap of Median COUNT Values',
                       os.path.join(out_dir, plate, 'heatmaps', 'COUNT.png'), lims=[0,80], colormap='Reds')

        map.mk_heatmap(plate_zscore.data_df, 'Heatmap of Median ZSCORE Values',
                       os.path.join(out_dir, plate, 'heatmaps', 'ZSCORE.png'),lims=[-5,5])

        map.mk_heatmap(plate_viability.data_df, 'Heatmap of Median FCPC Values',
                       os.path.join(out_dir, plate, 'heatmaps', 'FCPC.png'),lims=[-5,5])

        ssmd.norm_v_mfi_ssmd(plate_norm, plate_mfi, os.path.join(out_dir, plate, 'ssmd'))

        ssmd.ssmd_ecdf(plate_norm, plate_mfi, 'SSMD ECDF for {}'.format(plate),os.path.join(out_dir, plate, 'ssmd'))


        cp.median_ZSPC_histogram(plate_zscore, plate_zscore.col_metadata_df,
                                 os.path.join(out_dir, plate, 'cp_strength'), det=plate)

        cp.median_ZSPC_ecdf(plate_zscore, plate_zscore.col_metadata_df,
                                 os.path.join(out_dir, plate, 'cp_strength'), det=plate)



if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    logger.debug("args:  {}".format(args))

    main(proj_dir=args.build_folder, out_dir=args.qc_folder, invar=args.invar, sense=args.sensitivities)





