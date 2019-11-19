import os
import glob
import sys
import functools
import logging

import argparse
import pandas as pd

import cmapPy.pandasGEXpress.parse as pe
import cmapPy.pandasGEXpress.concat as concat

import merino.setup_logger as setup_logger
import merino.build_summary.ssmd_analysis as ssmd_analysis
import merino.build_summary.expected_sensitivities as sense

logger = logging.getLogger(setup_logger.LOGGER_NAME)

def build_parser():

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # The following arguments are required. These are files that are necessary for assembly and which change
    # frequently between cohorts, replicates, etc.
    parser.add_argument("-proj_folder", "-pd", help="path to the pod directory you want to run card on",
                        type=str, required=True)
    parser.add_argument("-gmt_path", "-gmt", help="path to the pod directory you want to run card on",
                        type=str, required=True, default='/cmap/data/vdb/merino/sensitivity_files/CORE_D1_48H_KJ100_DN_s25_n5127.gmt')
    parser.add_argument("-outfolder", "-out", help="path to the pod directory you want to run card on",
                        type=str, required=True)
    parser.add_argument("-verbose", '-v', help="Whether to print a bunch of output", action="store_true", default=False)

    return parser

def run_sensitivities(proj_folder, gmt_path, out_folder):

    gct_list = [pe.parse(y) for y in [x for x in glob.glob(os.path.join(proj_folder,'card/*/*ZSPC.gct'))]]
    fail_gct = concat.hstack(gct_list)
    if not os.path.exists(os.path.join(args.outfolder, 'sense')):
        os.mkdir(os.path.join(args.outfolder, 'sense'))

    sense.wtks(gct=fail_gct, metadata=fail_gct.col_metadata_df,
               outfolder=os.path.join(out_folder, 'sense'), group_col='prism_replicate',
               gmt_path=gmt_path)

def mk_report(proj_folder, out_folder):
    mar_sense = pd.read_table(os.path.join(out_folder,'sense/expected_sensitivity_ranks.txt'),
                              index_col='det_plate')
    mar_sense = mar_sense / 384
    mar_sense = mar_sense * 100
    for x in pd.Series([y.split('_')[0] for y in mar_sense.columns]).unique():
        mar_sense[x] = mar_sense[[p for p in mar_sense.columns if p.startswith(x)]].median(axis=1)
    mar_sense = mar_sense[pd.Series([y.split('_')[0] for y in mar_sense.columns]).unique()]

    gct_list = [pe.parse(y) for y in
                [x for x in glob.glob(os.path.join(proj_folder,'card/*/*NORM*'))]]
    norm_gct = concat.hstack(gct_list)

    gct_list = [pe.parse(y) for y in
                [x for x in glob.glob(os.path.join(proj_folder,'assemble/*/*MEDIAN*'))]]
    mfi_gct = concat.hstack(gct_list)

    n_recovered = []
    invs = []
    beadsets = []
    plate = []
    med_rank = []
    dropouts = []
    for det_plate in mar_sense.index:
        temp = norm_gct.data_df[[x for x in norm_gct.data_df.columns if x.startswith(det_plate)]]
        dropouts.append(384 - temp.shape[1])
        sigs_recovered = mar_sense.loc[det_plate].dropna()[mar_sense.loc[det_plate].dropna() < 50].count()
        median_rank = mar_sense.loc[det_plate].median()
        temp = mfi_gct.data_df[[x for x in mfi_gct.data_df.columns if x.startswith(det_plate)]]
        median_inv = temp.loc[['c-661', 'c-662', 'c-663', 'c-664']].median(axis=1).median()
        beadset = det_plate.split('_')[-1].split(':')[0]
        n_recovered.append(sigs_recovered)
        invs.append(median_inv)
        beadsets.append(beadset)
        plate.append(det_plate)
        med_rank.append(median_rank)

    mar_df = pd.concat([pd.Series(plate).rename('det_plate'), pd.Series(n_recovered).rename('sigs_recovered_core'),
                        pd.Series(med_rank).rename('median_rank_core'), pd.Series(invs).rename('median_inv'),
                        pd.Series(dropouts).rename('n_dropouts'), pd.Series(beadsets).rename('beadset')], axis=1)
    mar_df.set_index('det_plate', inplace=True)
    mar_ssmd = ssmd_an.ssmd_matrix(
        norm_paths=glob.glob(os.path.join(proj_folder,'card/*/*NORM*')))
    mar_df = mar_df.join(mar_ssmd[mar_ssmd < 2].count().rename('ssmd_failures'))
    return mar_df

def main(args):

    #run_sensitivities(args.proj_folder, args.gmt_path, args.outfolder)

    report_df = mk_report(args.outfolder)

    report_df.to_csv(os.path.join(args.outfolder, 'qc_report.txt'), sep='\t')

if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    logger.debug("args:  {}".format(args))

    main(args)