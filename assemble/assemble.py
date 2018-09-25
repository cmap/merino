"""

Command line script which takes the two CSVs belonging to a single PRISM replicate and combines them into a GCT file
along with all of the relevant meta data.

The meta data inputs are a plate map, a cell set definition file, a plate tracking file, and a davepool-analyte mapping.
"""
import merino.setup_logger as setup_logger
import logging
import davepool_data
import prism_metadata
import numpy
import argparse
import merino
import sys
import ConfigParser
import assemble_core
import os


logger = logging.getLogger(setup_logger.LOGGER_NAME)

_prism_cell_config_file_section = "PrismCell column headers"
_davepool_analyte_mapping_file_section = "DavepoolAnalyteMapping column headers"



def build_parser():

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # The following arguments are required. These are files that are necessary for assembly and which change
    # frequently between cohorts, replicates, etc.
    #todo: deprecate prism_replicate_name, implicit naming based on os.path.basename(args.csv_filepath)
    parser.add_argument("-prism_replicate_name", "-prn", help="name of the prism replicate that is being processed",
                        type=str, required=True)
    #todo: add choices to assay_type
    parser.add_argument("-config_filepath", "-cfg", help="path to the location of the configuration file", type=str,
                        default=merino.default_config_filepath)
    parser.add_argument("-assay_type", "-at", help="assay data comes from eg. PR500, PR300, KJ100",
                        type=str, required=True, choices=["DP78", "PR500", "PR300", "COPRO"])
    parser.add_argument("-pert_time", "-time", help="the assay time point in hours", type=str, required=True)

    # These arguments are optional. Some may be superfluous now and might be removed.
    parser.add_argument("-verbose", '-v', help="Whether to print a bunch of output", action="store_true", default=False)
    parser.add_argument("-plate_map_path", "-pmp",
                        help="path to file containing plate map describing perturbagens used", type=str, required=True)

    parser.add_argument("-analyte_mapping_file", "-dmf",
                        help="mapping of analytes to pools and davepools, overriding config file",
                        type=str, default=None, required=False)
    parser.add_argument("-cell_set_definition_file", "-csdf",
                        help="file containing cell set definition to use, overriding config file",
                        type=str, default=None, required=False)


    parser.add_argument("-ignore_assay_plate_barcodes", "-batmanify", help="list of assay plate barcodes that should be"
                        " ignored / excluded from the assemble", nargs="+", default=None)
    parser.add_argument("-outfile", "-out", help="location to write gct", type=str,
                        default='')
    parser.add_argument("-truncate_to_plate_map", "-trunc", help="True or false, if true truncate data to fit framework of platemap provided",
                        action="store_true", default=False)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-davepool_id_csv_filepath_pairs", "-dp_csv",
                        help="space-separated list of pairs of davepool_id and corresponding csv filepath for that davepool_id",
                        type=str, nargs="+", required=False)
    group.add_argument("-csv_filepath", "-csv",
                        help="space-separated list of pairs of davepool_id and corresponding csv filepath for that davepool_id",
                        type=str, nargs="+", required=False)

    return parser


def parse_location_to_well(location):
    split = location.split(",")
    right_paren_index = split[1].index(")")
    raw_well = split[1][0:right_paren_index]
    row = raw_well[0]
    raw_col = raw_well[1:]
    col = raw_col.zfill(2)
    return row + col


def read_davepool_data_objects(davepool_id_csv_list):
    '''
    create davepool objects and populate with data read from csv's
    :param davepool_id_csv_list: list of pairs of davepool ID and path to corresponding csv file for that davepool
    :return:
    '''
    r = []

    for (dp_id, csv_filepath) in davepool_id_csv_list:
        pd = davepool_data.read_data(csv_filepath)
        pd.davepool_id = dp_id
        r.append(pd)

    return r

def read_csv(csv_list, assay_type):
    r = []
    for csv_filepath in csv_list:
        pd = davepool_data.read_data(csv_filepath)
        pd.davepool_id = assay_type
        r.append(pd)

    return r


def combine_maps_with_checks(source_map, dest_map):
    source_keys = set(source_map.keys())
    dest_keys = set(dest_map.keys())

    common_keys = source_keys & dest_keys

    if len(common_keys) > 0:
        msg = "the source_map and dest_map had common_keys:  {}".format(common_keys)
        logger.error(msg)
        raise Exception("assemble combine_maps_with_checks " + msg)
    else:
        dest_map.update(source_map)



def build_davepool_id_csv_list(davepool_id_csv_filepath_pairs):
    '''
    break list of davepool_id and csv file path from input into pairs
    :param davepool_id_csv_filepath_pairs:
    :return:
    '''
    r = []

    for i in range(len(davepool_id_csv_filepath_pairs)/2):
        index = 2*i
        davepool_id = davepool_id_csv_filepath_pairs[index]
        csv_filepath = davepool_id_csv_filepath_pairs[index+1]
        r.append((davepool_id, csv_filepath))

    return r


def build_prism_cell_list(cell_set_definition_file, analyte_mapping_file):
    '''
    read PRISM cell line meta data from file specified in config file (at config_filepath), then associate with
    assay_plate based on pool ID.  Check for cell pools that are not associated with any assay plate
    :param assay_plates:
    :param cell_set_definition_file:
    :return:
    '''
    # cp = ConfigParser.RawConfigParser()
    # cp.read(config_filepath)

    # prism_cell_list_items = cp.items(_prism_cell_config_file_section)
    # analyte_mapping_items = cp.items(_davepool_analyte_mapping_file_section)

    prism_cell_list_items = [(x,x) for x in ["analyte_id", "pool_id", "davepool_id", "feature_id", "cell_iname", "minipool_id",
                             "ccle_name", "barcode_id", "cell_lineage", "cell_culture"]]
    analyte_mapping_items = [(x,x) for x in ["analyte_id", "feature_id", "cell_name", "davepool_id", "pool_id"]]

    prism_cell_list = prism_metadata.read_prism_cell_from_file(cell_set_definition_file, prism_cell_list_items)

    analyte_mapping = prism_metadata.read_prism_cell_from_file(analyte_mapping_file, analyte_mapping_items)

    cell_list_id_not_in_davepool_mapping = set()

    # Assign davepool mapping info to respective cell feature IDs

    cell_id_davepool_map = {}

    for dp in analyte_mapping:
        cell_id_davepool_map[dp.feature_id] = dp

    for pc in prism_cell_list:
        if pc.feature_id in cell_id_davepool_map.keys():
            cell_davepool = cell_id_davepool_map[pc.feature_id]
            pc.analyte_id = cell_davepool.analyte_id
            pc.davepool_id = cell_davepool.davepool_id
            if pc.pool_id != cell_davepool.pool_id:
                raise Exception ("Cell set pool id does not match davepool mapping pool id at cell id {}".format(pc.id))
        else:
            cell_list_id_not_in_davepool_mapping.add(pc.id)

    if len(cell_list_id_not_in_davepool_mapping) > 0:
        cell_list_id_not_in_davepool_mapping = list(cell_list_id_not_in_davepool_mapping)
        cell_list_id_not_in_davepool_mapping.sort()
        message2 = ("some cell ids were found in the cell set but not in the davepool mapping - IDs: {}".format(cell_list_id_not_in_davepool_mapping))
        raise Exception ("assemble build_prism_cell_list " + message2)
    return prism_cell_list


def truncate_data_objects_to_plate_map(davepool_data_objects, all_perturbagens, truncate_to_platemap):

    platemap_well_list = set([p.pert_well for p in all_perturbagens])
    for davepool in davepool_data_objects:
        if platemap_well_list == set(davepool.median_data.keys()):
            return davepool_data_objects
        elif truncate_to_platemap == True:
            for d in davepool_data_objects[0].median_data.keys():
                 if d not in platemap_well_list:
                     del davepool_data_objects[0].median_data[d]

            for c in davepool_data_objects[0].count_data.keys():
                if c not in platemap_well_list:
                    del davepool_data_objects[0].count_data[c]
        else:
            raise Exception("Assemble truncate data objects to plate map: Well lists of platemap and csv do not match")


    return davepool_data_objects


def main(args, all_perturbagens=None, assay_plates=None):
    if all_perturbagens is None:
        all_perturbagens = prism_metadata.build_perturbagens_from_file(args.plate_map_path, args.pert_time)

    args.ignore_assay_plate_barcodes = set(args.ignore_assay_plate_barcodes) if args.ignore_assay_plate_barcodes is not None else set()

    #read actual data from relevant csv files, associate it with davepool ID

    if args.davepool_id_csv_filepath_pairs is not None:
        davepool_id_csv_list = build_davepool_id_csv_list(args.davepool_id_csv_filepath_pairs)
        davepool_data_objects = read_davepool_data_objects(davepool_id_csv_list)

    elif args.csv_filepath is not None:
        davepool_id_csv_list = args.csv_filepath
        davepool_data_objects = read_csv(davepool_id_csv_list, args.assay_type)

    cp = ConfigParser.RawConfigParser()
    cp.read(args.config_filepath)

    #read PRISM cell line metadata from file specified in config file, and associate with assay_plate metadata
    cell_set_file = args.cell_set_definition_file if args.cell_set_definition_file else cp.get(args.assay_type, "cell_set_definition_file")
    analyte_mapping_file = args.analyte_mapping_file if args.analyte_mapping_file else cp.get(args.assay_type, "analyte_mapping_file")

    prism_cell_list = build_prism_cell_list(cell_set_file, analyte_mapping_file)

    logger.info("len(prism_cell_list):  {}".format(len(prism_cell_list)))

    # truncate csv to plate map size if indicated by args.truncate_to_plate_map
    truncate_data_objects_to_plate_map(davepool_data_objects, all_perturbagens, args.truncate_to_plate_map)

    # Pass python objects to the core assembly module (this is where command line and automated assembly intersect)
    assemble_core.main(args.prism_replicate_name, args.outfile, all_perturbagens, davepool_data_objects, prism_cell_list)


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    logger.debug("args:  {}".format(args))

    main(args)