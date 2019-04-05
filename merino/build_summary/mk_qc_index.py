import os
import sys
import glob
import logging
import argparse
import pandas as pd

import cmapPy.pandasGEXpress.parse as pe

import merino.setup_logger as setup_logger
import make_gallery as galleries

logger = logging.getLogger(setup_logger.LOGGER_NAME)

def build_parser():

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # The following arguments are required. These are files that are necessary for assembly and which change
    # frequently between cohorts, replicates, etc.
    parser.add_argument("-qc_folder", "-qc", help="string designating the prefix to each build file eg. PCAL075-126_T2B",
                        type=str, required=True)
    parser.add_argument("-project_name", "-pn",
                        help="Code for project eg. PCAL",
                        type=str, required=False)
    parser.add_argument("-verbose", '-v', help="Whether to print a bunch of output", action="store_true", default=False)


    return parser


def qc_galleries(out_dir, proj_name):
    local_paths = pd.Series(glob.glob(os.path.join(out_dir, proj_name + '*', '*.html'))).sort_values().tolist()
    dex = [os.path.basename(os.path.dirname(x)) for x in local_paths]
    index_info = pd.concat([pd.read_table(x) for x in glob.glob(os.path.join(out_dir, '*', 'report.txt'))])
    index_info.set_index('plate', inplace=True)
    index_info = index_info.loc[dex]

    def make_url(ref, name):
        ref = os.path.relpath(ref, out_dir)
        return '<a target="_blank" href="{}">{}</a>'.format(ref, name)

    premadelinks = [make_url(x, dex[i]) for i, x in enumerate(local_paths)]

    headers = index_info.reset_index().columns
    table_info = zip(premadelinks, index_info['median SSMD'].tolist(),
                     index_info['n SSMD failures'].tolist(), index_info['pct SSMD failures'].tolist(),
                     index_info['median_invariant'].tolist(), index_info['median_sensitivity_rank'].tolist(),
                     index_info['n sensitivities recovered'].tolist(), index_info['n wells'].tolist(),
                     index_info['n dropouts'].tolist(), index_info['n active wells'].tolist(),
                     index_info['n unique perts'].tolist(),index_info['qc_status'].tolist()
                     )
    #print index_info
    made_gallery = galleries.mk_index(table_headers=headers, table_tuples=table_info, outfolder=out_dir, project_name=proj_name)
    if made_gallery:
        logger.info("successfully made QC gallery")

    import pdb
    pdb.set_trace()
    index_info.to_csv(os.path.join(out_dir, 'qc_report.txt'), sep='\t')

def main(args, out_dir, project_name):

    # Read in the data

    qc_galleries(out_dir, project_name)


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    logger.debug("args:  {}".format(args))

    main(args, out_dir=args.qc_folder, project_name=args.project_name)