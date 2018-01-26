import glob
import os
import cmapPy.pandasGEXpress.parse as pe
import cmapPy.pandasGEXpress.GCToo as GCToo
import distributions as dist
import invariant_analysis as inv
import pandas as pd
import generic_heatmap as map
import ssmd_analysis as ssmd
import compound_strength as cp
import comp_strength_overview as comp
#import expected_sensitivities as sense
import count_heatmap as c_map
import matplotlib.pyplot as plt
import numpy as np
import sc_plot


def main(proj_dir, out_dir):
    # READ EM ALL IN

    print "Reading LEVEL2 Data"
    mfi_path = glob.glob(os.path.join(proj_dir, '*MFI*.gctx'))[0]
    mfi_gct = pe(mfi_path)

    count_path = glob.glob(os.path.join(proj_dir, '*COUNT*.gctx'))[0]
    count_gct = pe(count_path)

    print "Reading LEVEL3 Data"
    norm_path = glob.glob(os.path.join(proj_dir, '*NORM*.gctx'))[0]
    norm_gct = pe(norm_path)

    print "Reading LEVEL4 Data"
    zscore_path = glob.glob(os.path.join(proj_dir, '*_ZSPC_*.gctx'))[0]
    zscore_gct = pe(zscore_path)

    zsvc_path = glob.glob(os.path.join(proj_dir, '*_ZSVC_*.gctx'))[0]
    zsvc_gct = pe(zsvc_path)

    viability_path = glob.glob(os.path.join(proj_dir, '*4_LFCPC_*.gctx'))[0]
    viability_gct = pe(viability_path)

    fcvc_path = glob.glob(os.path.join(proj_dir, '*4_LFCVC_*.gctx'))[0]
    fcvc_gct = pe(fcvc_path)

    print "Reading LEVEL5 Data"
    modz_path = glob.glob(os.path.join(proj_dir, '*MODZ.ZSPC_*.gctx'))[0]
    modz_gct = pe(modz_path)

    print "Reading Metadata"
    inst_info = pd.read_table(glob.glob(os.path.join(proj_dir, '*inst_info*.txt'))[0], index_col='cid')
    sig_info = pd.read_table(glob.glob(os.path.join(proj_dir, '*MODZ.ZSPC_sig_metrics.txt'))[0], index_col='cid')

    reload(dist)
    reload(inv)
    reload(map)

    print os.path.join(out_dir, 'modz_dist.png')

    if not os.path.exists(os.path.join(out_dir, 'sensitivities')):
        os.mkdir(os.path.join(out_dir, 'sensitivities'))

    #sense.wtks(zscore_gct, inst_info, os.path.join(out_dir, 'sensitivities'))

    my_map = {mfi_gct: ['MFI', 0, 50000], norm_gct: ['NORM', -10, 6],
              zscore_gct: ['ZSPC', -10, 10], zsvc_gct: ['ZSVC', -10, 10],
              viability_gct: ['FCPC', -10, 10], fcvc_gct: ['FCVC', -10, 10],
              modz_gct: ['MODZ', -10, 10]}

    for df in my_map.keys():
        bins = np.linspace(my_map[df][1], my_map[df][2], 100)
        plt.hist(df.data_df.unstack().dropna(), bins)
        plt.title('Distribution for PCAL {}'.format(my_map[df][0]))
        plt.xlabel('{} Data'.format(my_map[df][0]))
        plt.savefig(out_dir + '/PCAL_{}.png'.format(my_map[df][0]))
        plt.clf()

    c_map.make_count_heatmap(count_gct, inst_info, 'Count Across PCALT2B', os.path.join(out_dir, 'count_map.png'))

    reload(comp)
    comp.modz_dist(modz_gct, sig_info, [], os.path.join(out_dir, 'modz_dist.png'))
    dist.distributions(norm_gct, mfi_gct, count_gct, zscore_gct, viability_gct, inst_info, os.path.join(out_dir))

    inv.invariant_monotonicity(mfi_gct, inst_info, out_dir)

    sc_plot.sc_plot(sig_info, 'SS vs CCQ75', os.path.join(out_dir,'sc_modz.zspc.png'))





    for plate in inst_info['prism_replicate'].unique():
        plate_mfi = GCToo.GCToo(data_df=mfi_gct.data_df[inst_info[inst_info['prism_replicate'] == plate].index],
                                col_metadata_df=inst_info[inst_info['prism_replicate'] == plate],
                                row_metadata_df=mfi_gct.row_metadata_df)
        print plate_mfi.data_df.shape

        plate_count = GCToo.GCToo(data_df=count_gct.data_df[inst_info[inst_info['prism_replicate'] == plate].index],
                                col_metadata_df=inst_info[inst_info['prism_replicate'] == plate],
                                row_metadata_df=count_gct.row_metadata_df)

        plate_norm = GCToo.GCToo(data_df=norm_gct.data_df.loc[:,[x for x in inst_info[inst_info['prism_replicate'] == plate].index if x in norm_gct.data_df.columns]],
                                col_metadata_df=inst_info.loc[[x for x in inst_info[inst_info['prism_replicate'] == plate].index if x in norm_gct.data_df.columns]],
                                row_metadata_df=norm_gct.row_metadata_df)

        plate_zscore = GCToo.GCToo(data_df=zscore_gct.data_df.loc[:,[x for x in inst_info[inst_info['prism_replicate'] == plate].index if x in norm_gct.data_df.columns]],
                                col_metadata_df=inst_info.loc[[x for x in inst_info[inst_info['prism_replicate'] == plate].index if x in norm_gct.data_df.columns]],
                                row_metadata_df=zscore_gct.row_metadata_df)

        plate_viability = GCToo.GCToo(data_df=viability_gct.data_df.loc[:,[x for x in inst_info[inst_info['prism_replicate'] == plate].index if x in norm_gct.data_df.columns]],
                                col_metadata_df=inst_info.loc[[x for x in inst_info[inst_info['prism_replicate'] == plate].index if x in norm_gct.data_df.columns]],
                                row_metadata_df=zscore_gct.row_metadata_df)

        if not os.path.exists(os.path.join(out_dir, plate)):
            os.mkdir(os.path.join(out_dir, plate))
            for fold in ['invariants', 'distributions', 'heatmaps', 'ssmd', 'cp_strength']:
                os.mkdir(os.path.join(out_dir, plate, fold))

        dist.distributions(plate_norm, plate_mfi, plate_count, plate_zscore, plate_viability,
                           plate_norm.col_metadata_df, os.path.join(out_dir, plate, 'distributions'))

        for func in [inv.invariant_monotonicity, inv.invariant_curves_plot, inv.invariant_range_distributions]:
            func(plate_mfi, plate_mfi.col_metadata_df, os.path.join(out_dir, plate, 'invariants'))

        inv.invariant_heatmap(plate_mfi, os.path.join(out_dir, plate, 'invariants', 'inv_heatmap.png'))

        map.mk_heatmap(plate_norm.data_df, 'Heatmap of Median Norm Values',
                       os.path.join(out_dir, plate,'heatmaps', 'NORM.png'), lims=[-3,3])

        map.mk_heatmap(plate_mfi.data_df, 'Heatmap of Median MFI Values',
                       os.path.join(out_dir, plate, 'heatmaps', 'MFI.png'))

        map.mk_heatmap(plate_count.data_df, 'Heatmap of Median COUNT Values',
                       os.path.join(out_dir, plate, 'heatmaps', 'COUNT.png'))

        map.mk_heatmap(plate_zscore.data_df, 'Heatmap of Median ZSCORE Values',
                       os.path.join(out_dir, plate, 'heatmaps', 'ZSCORE.png'))

        map.mk_heatmap(plate_viability.data_df, 'Heatmap of Median FCPC Values',
                       os.path.join(out_dir, plate, 'heatmaps', 'FCPC.png'))


        ssmd.norm_v_mfi_ssmd(plate_norm, plate_mfi, os.path.join(out_dir, plate, 'ssmd'))


        cp.median_ZSPC_histogram(plate_zscore, plate_zscore.col_metadata_df,
                                 os.path.join(out_dir, plate, 'cp_strength'), det=plate)

        cp.median_ZSPC_ecdf(plate_zscore, plate_zscore.col_metadata_df,
                                 os.path.join(out_dir, plate, 'cp_strength'), det=plate)





