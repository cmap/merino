import sys
import utils.mysql_utils as mu
import utils.orm.lims_assemble_params_orm as lapo
import os
import argparse
import setup_logger
import logging

import requests


logger = logging.getLogger(setup_logger.LOGGER_NAME)

# Initialize MySQL Connection
default_config_filepath = os.path.expanduser('~/.prism_pipeline.cfg')
def build_parser():

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # The following arguments are required. These are files that are necessary for assembly and which change
    # frequently between cohorts, replicates, etc.
    parser.add_argument("-det_plate", "-dp", help="name of the prism replicate that is being processed",
                        type=str, required=True)
    parser.add_argument("-config_filepath", "-dmf", help="mapping of analytes to pools and davepools",
                        type=str, default = default_config_filepath)
    parser.add_argument("-pod_dir_path", '-pdp', help="Custom path to pod directory", type=str, default='/cmap/obelix/pod/custom')
    parser.add_argument("-custom_outfile_path", '-out', help="Custom outfile path", type=str, default=None)
    parser.add_argument("-verbose", '-v', help="Whether to print a bunch of output", action="store_true", default=False)

    return parser


def query_db(cursor, det_plate):
    # Query database for necessary information. See caldaia/utils/orm/lims_assemble_params_orm.py
    my_lapo = lapo.get_LAP(cursor, det_plate)

    return my_lapo


def construct_plate_map_path(pod_dir, my_lapo):
    '''
    Use proj_id to get pod direcotry and pert plate name to get map file name.
    '''
    #TODO Change the /cmap/obelix/pod/custom to an argument or put it in the config file. ".src" as well
    # Add option to use custom plate map filepath
    map_path = os.path.join(pod_dir, my_lapo.project_id, "map_src", my_lapo.pert_plate + ".src")

    return map_path


def construct_davepool_csv_path_pairs(pod_dir, my_lapo):
    '''
    Using project ID, prism replicate name, and list of connected davepools, construct paths to csv for each davepool.
    Concatenate davepool IDs and csv filepaths into a single string for passing to assemble as an arg.
    :param my_lapo:
    :return: dp csv pairs string
    '''
    # We're going to make a list which we will then convert to a string
    dlist = []

    for davepool in my_lapo.davepool_list:
        dlist.append(davepool)
        det_plate = my_lapo.pert_plate + '_' + davepool + '_' + my_lapo.replicate
        csv_path = os.path.join(pod_dir, my_lapo.project_id, 'lxb', det_plate, det_plate + '.csv')
        dlist.append(csv_path)

    arg_string = " ".join(dlist)

    return arg_string


def construct_outfile_path(my_lapo, custom_outfile, pod_dir_path):
    '''
    Check if assemble directory exists. If it does not, create it. The construct outfile path.
    :param my_lapo:
    :return: outfile string
    '''
    #TODO make assemble folder argument too
    if custom_outfile is not None:
        outdir = custom_outfile
    else:
         outdir = os.path.join(pod_dir_path, my_lapo.project_id, 'assemble',
                           my_lapo.prism_replicate_name)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    return outdir

# put default in arg builder

def build_args(cursor, det_plate, config_filepath, pod_dir):
    '''
    Query database, use returned values to construct all assemble args and put into a list.
    :param det_plate:
    :return: Raw arguments as dictionary. Keys in dict correspond to form names in flask app.
    '''

    my_lapo = query_db(cursor, det_plate)

    map_path = construct_plate_map_path(pod_dir, my_lapo)

    dp_csv_string = construct_davepool_csv_path_pairs(pod_dir, my_lapo)

    outfile = construct_outfile_path(my_lapo, args.custom_outfile_path, args.pod_dir_path)

    # Arguments will be passed as a map to the flask server.
    assemble_args = {"prn": my_lapo.prism_replicate_name, "pmp": map_path,
                "cfg": config_filepath, "dam": my_lapo.davepool_mapping_file, "csdf": my_lapo.cell_set_definition_file,
                "dp_csv": dp_csv_string, 'outfile': outfile, 'submit': 'SUBMIT'}

    return assemble_args

def main(args):
    
    db = mu.DB(host="localhost").db
    cursor = db.cursor()
    assemble_args = build_args(cursor, args.det_plate, args.config_filepath, args.pod_dir_path)
    print assemble_args
    rep = requests.post('http://127.0.0.1:5000/', data=assemble_args)
    print rep.text
    db.close()

if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    logger.debug("args:  {}".format(args))

    main(args)

