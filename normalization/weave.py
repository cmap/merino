import os
import glob
import modz
import functools
import shear
import pandas as pd
import cmapPy.pandasGEXpress.parse as pe
import merino.setup_logger as setup_logger
import logging
import argparse
import sys
import cmapPy.pandasGEXpress.write_gct as wg
import ConfigParser


logger = logging.getLogger(setup_logger.LOGGER_NAME)

def build_parser():

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # The following arguments are required. These are files that are necessary for assembly and which change
    # frequently between cohorts, replicates, etc.
    parser.add_argument("-proj_dir", "-pd", help="path to the pod directory you want to run card on",
                        type=str, required=True)
    parser.add_argument("-search_pattern", "-sp",
                        help="Search for this string in the directory, only run plates which contain it. "
                             "Default is wildcard",
                        type=str, default='*', required=False)
    parser.add_argument("-input_folder", "-if",
                        help="Search for this string in the directory, only run plates which contain it.",
                        type=str, default='ZSPC', required=False)
    parser.add_argument("-verbose", '-v', help="Whether to print a bunch of output", action="store_true", default=False)
    parser.add_argument("-bad_wells", "-wells", help="List of wells to be excluded from processing", type=list,
                        default=[])
    parser.add_argument("-log_tf", "-log", help="True or false, if true log transform the data",
                        action="store_true", default=True)

    return parser


def weave(proj_dir, rep_set, input_folder='zscorepc'):

        print rep_set
        files = glob.glob(os.path.join(proj_dir, input_folder, rep_set + '*'))
        plate_names = [os.path.basename(x) for x in files]
        replicate_ids = [x.split("_")[2] for x in plate_names]
        short_reps = [x.split('.')[0] for x in replicate_ids]
        keep = []

        for r in set(short_reps):
            temp = [y for y in replicate_ids if y.startswith(r)]

            if len(temp) == 1:
                keep.append(temp[0])
            else:

                temp2 = [z for z in temp if "." in z]

                max_l = max([int(x[-1]) for x in temp2])
                temp3 = [b for b in temp2 if b.endswith(str(max_l))]
                keep.append(temp3[0])

        keep_perts = [x for x in plate_names if x.split('_')[2] in keep]
        keep_files = [glob.glob(path + '/*')[0] for path in files if os.path.basename(path) in keep_perts]


        pert = plate_names[0].split('_')[0] + '_' + plate_names[0].split('_')[1]



        if not os.path.exists(os.path.join(proj_dir, 'modz.{}'.format(input_folder), pert)):
            os.mkdir(os.path.join(proj_dir, 'modz.{}'.format(input_folder), pert))

        reload(modz)
        gct_list = []
        for path in keep_files:
            gct = pe(path)
            print gct.data_df.shape
            gct_list.append(gct)


        modZ_GCT, cc_q75_df, weights = modz.calculate_modz(gct_list)

        outfile = os.path.join(proj_dir, 'modz.{}'.format(input_folder), pert)

        if not os.path.exists(outfile):
            os.mkdir(outfile)

        weights[0].to_csv(os.path.join(outfile, rep_set + '_norm_weights.txt'), sep='\t')
        weights[1].to_csv(os.path.join(outfile, rep_set + '_raw_weights.txt'), sep='\t')
        wg.write(modZ_GCT, os.path.join(outfile, pert + '_MODZ.{}'.format(input_folder)))
        cc_q75_df.to_csv(os.path.join(outfile, rep_set + '_cc_q75.txt'), sep='\t')


def main(args):
    if not os.path.exists(os.path.join(args.proj_dir, 'modz.{}'.format(args.input_folder))):
        os.mkdir(os.path.join(args.proj_dir, 'modz.{}'.format(args.input_folder)))
    for x in set([os.path.basename(y).split('_')[0] for y in glob.glob(os.path.join(args.proj_dir, '{}/*'.format(args.input_folder)))]):

        weave(args.proj_dir, x, input_folder=args.input_folder)


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    logger.debug("args:  {}".format(args))

    main(args)