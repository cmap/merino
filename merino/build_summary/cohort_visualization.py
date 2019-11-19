import argparse
import os
import sys
import logging
import glob
import string
logger = logging.getLogger()
logging.basicConfig(filename='/dev/null', level=logging.ERROR, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import cmapPy.pandasGEXpress.parse as pe
import cmapPy.pandasGEXpress.subset_gctoo as sub
from cmapPy.visualization import stratogram, scattergram, cohort_view
import metrics as new_metrics
import merino.setup_logger as setup_logger
import make_gallery as galleries
logger = logging.getLogger(setup_logger.LOGGER_NAME)



def build_parser():

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # The following arguments are required. These are files that are necessary for assembly and which change
    # frequently between cohorts, replicates, etc.
    parser.add_argument("-build_folder_path", '-build', help="path to folder containing necessary build files with appropriate naming conventions",
                        type=str, required=True)
    parser.add_argument("-out", "-o",
                        help="Our folder for results",
                        type=str, required=True)
    parser.add_argument("-dict_key", "-key",
                        help="String key corresponding to a cohort, determining how the data will be categorized. Check metrics.py meta_dict for possible choices",
                        type=str, default='PR500', required=False, choices=['PR500', 'COP', 'DP78'])
    parser.add_argument("-dont_overwrite_sig_metrics", "-dont_overwrite", help="True or False, if true overwrite given sig_metrics file with added metrics",
                        action="store_true")
    parser.add_argument("-verbose", '-v', help="Whether to print a bunch of output", action="store_true", default=False)


    return parser


data_mapper = dict(zip(['ZSPC', 'ZSVC'], ['*MODZ.ZSPC.COMBAT*.gctx', '*MODZ.ZSVC.COMBAT*.gctx']))
metadata_mapper = dict(zip(['ZSPC', 'ZSVC'], ['*MODZ.ZSPC.COMBAT*.txt', '*MODZ.ZSVC.COMBAT*.txt']))

pr500_section_list = [{'name': "Cohort Distributions of Level 5 Data", "images": ['ZSPC_modz_hist.png', 'ZSVC_modz_hist.png']},
                    {"name" : "Cohort Overview of Killing Metrics", "images" : ['ZSPC_scattergram.png','ZSVC_scattergram.png',
                                                                                          'ZSPC_stratogram.png', 'ZSVC_stratogram.png']},
                    {"name" : "Quality Control", "images" : ['vc_pc_scatter.png']}]

cop_section_list = [
        {'name': "Cohort Distributions of Level 5 Data", "images": ['ZSPC_vi_modz_hist.png', 'ZSVC_vi_modz_hist.png',
                                                                    'ZSPC_gex_modz_hist.png', 'ZSVC_gex_modz_hist.png']},
        {"name": "Cohort Overview of Killing Metrics", "images": ['ZSPC_scattergram.png', 'ZSVC_scattergram.png',
                                                                  'ZSPC_stratogram.png', 'ZSVC_stratogram.png',
                                                                  'ZSPC_cis_tas_scatter.png', 'ZSVC_cis_tas_scatter.png']},
        {"name": "Quality Control", "images": ['viability_vc_pc_scatter.png', 'gex_vc_pc_scatter.png']}]



def globandparse(search_pattern):
    """
    Search for inputs based on a search pattern and read in file
    Args:
        search_pattern: Identifier to find a given build file using glob search

    Returns: GCToo object found through glob search and read in

    """
    logger.info("globbing for : {}".format(search_pattern))
    #TODO check length of glob list
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
    data_pattern_list = [data_mapper[x] for x in data_mapper.keys()]

    # Link identifiers to path to build folder
    data_search_patterns = [os.path.join(proj_dir, x) for x in data_pattern_list]

    # Search for and reaad in each file in the build
    gcts = [globandparse(x) for x in data_search_patterns]

    # Make dictionary and link each GCToo object to a corresponding key
    data_map = dict(zip(data_mapper.keys(), gcts))

    print "Reading Metadata"
    # Read in each metadata file

    metadata_pattern_list = [metadata_mapper[x] for x in metadata_mapper.keys()]
    metadata_search_patterns = [os.path.join(proj_dir, x) for x in metadata_pattern_list]
    metadata_paths = [glob.glob(x)[0] for x in metadata_search_patterns]

    sigmetrics_list = [[pd.read_table(x, index_col='sig_id'), x] for x in metadata_paths]

    sig_map = dict(zip(metadata_mapper.keys(), sigmetrics_list))

    other_metadata_map = {}
    inst_info = pd.read_table(glob.glob(os.path.join(proj_dir, '*inst_info*.txt'))[0], index_col='profile_id')


    ssmd_info = pe.parse(glob.glob(os.path.join(proj_dir, '*ssmd*.gctx'))[0]).data_df

    other_metadata_map['inst'] = inst_info
    other_metadata_map['ssmd'] = ssmd_info

    return data_map, other_metadata_map, sig_map


def split_cop_data(data_map):
    new_data_map = {}
    for key in data_map:
        temp_gct = data_map[key]
        vi_gct = sub.subset_gctoo(temp_gct, rid=[x for x in temp_gct.data_df.index.tolist() if 'c-' in x])
        gex_gct = sub.subset_gctoo(temp_gct, rid=[x for x in temp_gct.data_df.index.tolist() if 'c-' not in x])
        new_data_map[key + '_vi'] = vi_gct
        new_data_map[key + '_gex'] = gex_gct

    return new_data_map


def mk_heatmap(df, title, outfile, colormap='coolwarm', lims=[]):
    #values = df.median(axis=0)
    values = df
    heatmap_df = pd.DataFrame(index=list(string.ascii_uppercase)[0:16])
    columns = [str.zfill(str(x), 2) for x in range(1, 25)]

    for col in columns:
        curr_column = values[[x[-2:] == col for x in values.index]]
        curr_column.index = [y[-3] for y in curr_column.index]
        try:
            heatmap_df[col] = curr_column
        except:
            print 'skipping {}'.format(title)
            return
    if len(lims) == 0:
        fig, ax = plt.subplots(figsize=(15,10))
        sns.heatmap(heatmap_df, linewidths=.1, cmap=colormap, ax=ax)
    elif len(lims) == 2:
        sns.heatmap(heatmap_df, linewidths=.1, cmap=colormap, vmin=lims[0], vmax=lims[1])
    else:
        raise Exception('Must pass exactly 2 color scale limits or none at all')
    plt.yticks(rotation=1)
    plt.title(title)
    plt.savefig(outfile)

def make_gallery(qc_dir, section_list):



    outfile = os.path.join(qc_dir, 'gallery.html')
    print 'here'

    galleries.mk_cohort_gal(section_list, outfile)



def plot_pc_vc_scatter(data_pc, data_vc, outfile):


    df = data_pc.data_df.unstack().to_frame('pc')
    df['vc'] = data_vc.data_df.unstack()

    if len(df['vc']) > 500000:
        df = df.sample(500000)

    plt.figure(figsize=(10, 10))
    with sns.axes_style('whitegrid', {'font.family': "Roboto"}):
        plt.scatter(df['pc'], df['vc'], alpha=0.1, s=2)
        plt.xlim(-10, 10)
        plt.ylim(-10, 10)
        plt.xlabel('PC', fontweight="bold", fontsize=16)
        plt.ylabel('VC', fontweight="bold", fontsize=16)
        plt.plot((-10, 10), (-2, -2), 'r')
        plt.plot((-2, -2), (-10, 10), 'r')

    plt.savefig(outfile)


def plot_cis_tas_scatter(sig, outfile):

    plt.figure(figsize=(7, 7))
    ms = 3

    plt.loglog(sig['distil_tas'], sig['cis_ltn2'], "o", color="#111111", alpha=0.1, markersize=ms, markeredgewidth=0.3)

    plt.grid(which="both", color="#aaaaaa", linestyle="--")
    plt.gca().set_aspect('equal')
    # plt.legend(loc=(1.02, 0))
    plt.xlabel('TAS', fontweight='bold')
    plt.ylabel('CIS', fontweight='bold')
    plt.xlim(1e-3, 1)
    plt.ylim(1e-3, 1)

    for counter, cas in enumerate(np.arange(0.1, 0.7, 0.1)):
        xx = np.logspace(-3, 0, 50)
        yy = cas ** 2 / xx
        plt.loglog(xx, yy, 'orange')

        plt.text(cas ** 2, 1.1, "{:.2}".format(cas), horizontalalignment="center", color="orange")
    plt.text(3e-3, 1.1, "CAS:", horizontalalignment="center", color="orange")
    plt.savefig(outfile)


def plot_modz_hist(df_sig, data_df, outfile, title=''):
    '''
    df_sig: subset of siginfo defining a class of compounds
    data_df: GCToo.data_df whose elements are to be plotted
    '''
    with sns.axes_style('ticks'):
        sns.set_palette('Set1')
        bins = np.linspace(-10, 10, 200)
        pert_type = df_sig.name
        if pert_type not in ['trt_cp', 'ctl_vehicle', 'trt_poscon', 'trt_poscon.es']:
            return
        if pert_type == "trt_cp":
            edgecolor = "black"
            lw = 2
            alpha = 0.9
        else:
            edgecolor = "white"
            lw = 1
            alpha = 0.5

        x = data_df.loc[:, df_sig.index].values.flatten()
        plt.hist(x,
                 bins=bins, histtype="stepfilled", lw=lw, alpha=alpha, normed=True,
                 label=new_metrics.dict_pert_type_labels[pert_type], log=True,
                 edgecolor=edgecolor)
        plt.ylim(1e-3, 10)
        sns.despine(trim=True)
        fontdict = dict(fontname="Lato", fontweight=600, color="#333333")
        plt.legend(loc=1, fontsize=10)
        plt.xlabel('MODZ', **fontdict)
        plt.ylabel('Density', **fontdict)
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)
        plt.title(title, **fontdict)
        plt.savefig(outfile)


def plot_stratogram(sig, outfile):

    stratogram.stratogram(
        sig,
        category_definition="category_label",
        category_label="category_label_abridged",
        category_order="category_order",
        outfile = outfile,
        metrics = ['sel_vi', 'rep_vi', 'cis_ltn2', 'mag_vi', 'cont_vi', 'spec_vi'],
        column_display_names = ['Selectivity', 'Reproducibility', 'CIS', 'Magnitude', 'Contrast', 'Specificity'],
        bins=21,figsize=(20,8),
            xtick_orientation="vertical",
            ylabel_fontsize=16,
            xlabel_fontsize=20,
            xlabel_fontcolor="#222222",
            ylabel_fontcolor="#222222",
            fontfamily="Roboto",
    )

def plot_scattergram(sig, outfile, scattergram_columns = ['sel_vi', 'rep_vi', 'mag_vi', 'cont_vi']):


    plt.figure()
    g = scattergram.scattergram(
        sig,
        columns=scattergram_columns,
        column_names=scattergram_columns,
        outfile = outfile,
        title='Scattergram',
        fontfamily="Roboto"
        )



"""Adding stratification annotations to siginfo fields that need to be set for various stratification operations to work:

"si_bucket" (see for example metrics.assign_ss_buckets_PR500_v7())
"category": Defines the stratum, and is used as the key for the dictionaries defining the other fields in this list. For trt_cp, use si_bucket, for others use their pert_type
"category_label": full string description of each bucket. For test compounds, this string should contain "Test compound".
"category_label_abridged": Shorter string category label.
"category_order" integer order in which items should appear in stratified views.
For each of these, see dictionaries defined in Biosensor.metrics"""


def add_annotations(df, dict_key='PR500', start_column = 'ss_ltn2'):
    my_dict = new_metrics.meta_dict[dict_key]
    ''' Equivalent of the version used for CoP nominations '''
    df['si_bucket'] = df[start_column].apply(lambda s:new_metrics.assign_ss_buckets(s, dict_key))
    df['category'] = df.apply(new_metrics.assign_category, bucket_column="si_bucket", axis=1)
    df['category_label_abridged'] = df.category.apply(lambda x:my_dict[x][0])
    df['category_label'] = df.category.apply(lambda x:my_dict[x][1])
    df['category_order'] = df.category.apply(lambda x:my_dict[x][2])
    return df


def main(args):
    #sig_metrics = pd.read_table(args.sig_metrics_path)

    data_map, metadata_map, sig_map = read_build_data(proj_dir=args.build_folder_path)


    for key in sig_map.keys():
        updated_sig_metrics = add_annotations(sig_map[key][0], args.dict_key)

        if args.dont_overwrite_sig_metrics == False:
            updated_sig_metrics.to_csv(sig_map[key][1], sep='\t')


    # Wrap in try except

    for key in sig_map:
        plot_stratogram(sig_map[key][0], os.path.join(args.out, '{}_stratogram.png'.format(key)))

        plot_scattergram(sig_map[key][0], os.path.join(args.out, '{}_scattergram.png'.format(key)))

        mk_heatmap(sig_map[key][0][sig_map[key][0]['ss_ltn2'] > 5]['pert_well'].value_counts(),
                   'Distribution of Active Signatures Across Plate', '{}_hit_map.png'.format(key))


    if args.dict_key == 'COP':
        data_map = split_cop_data(data_map)

        [plot_cis_tas_scatter(sig_map[x][0], os.path.join(args.out, '{}_cis_tas_scatter'.format(x))) for x in sig_map.keys()]

        plot_pc_vc_scatter(data_map['ZSPC_vi'], data_map['ZSVC_vi'], outfile=os.path.join(args.out, 'viability_vc_pc_scatter.png'))

        plot_pc_vc_scatter(data_map['ZSPC_gex'], data_map['ZSVC_gex'],
                           outfile=os.path.join(args.out, 'gex_vc_pc_scatter.png'))

    else:
        plot_pc_vc_scatter(data_map['ZSPC'], data_map['ZSVC'],
                           outfile=os.path.join(args.out, 'vc_pc_scatter.png'))

    for key in data_map:
        sig = sig_map[key.split('_')[0]][0]
        data = data_map[key]
        outfile = os.path.join(args.out, '{}_modz_hist.png'.format(key))
        plt.figure(figsize=(14, 5))
        sig.groupby('pert_type').apply(
            plot_modz_hist,
            data.data_df,
            title='Distribution of viability scores',
            outfile = outfile
        )
    if args.dict_key == 'COP':
        make_gallery(args.out, cop_section_list)
    else:
        make_gallery(args.out, pr500_section_list)


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    logger.debug("args:  {}".format(args))

    main(args)
