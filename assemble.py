"""

Command line script which takes the two CSVs belonging to a single PRISM replicate and combines them into a GCT file
along with all of the relevant meta data.

The meta data inputs are a plate map, a cell set definition file, a plate tracking file, and a davepool-analyte mapping.
"""
import setup_logger
import logging
import validate_prism_gct
import davepool_data
import prism_metadata
import numpy
import argparse
import prism_pipeline
import sys
import os
import broadinstitute_cmap.io.GCToo.GCToo as GCToo
import pandas
import broadinstitute_cmap.io.GCToo.write_gctoo as write_gctoo
import ConfigParser
import utils.mysql_utils as mysql_utils
import utils.orm.assay_plates_orm as ap_orm




logger = logging.getLogger(setup_logger.LOGGER_NAME)

_prism_cell_config_file_section = "PrismCell column headers"
_davepool_analyte_mapping_file_section = "DavepoolAnalyteMapping column headers"


class DataByCell:
    def __init__(self, cell_data_map=None, well_list=None):
        self.cell_data_map = cell_data_map
        self.well_list = well_list
    def __str__(self):
        return "cell_data_map:  {}  well_list:  {}".format(self.cell_data_map, self.well_list)

_remove_row_annotations = ["id", "ignore"]
_remove_col_annotations = ["assay_plate_barcode"]

_null = "-666"
_NaN = "NaN"


def build_parser():

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # The following arguments are required. These are files that are necessary for assembly and which change
    # frequently between cohorts, replicates, etc.
    parser.add_argument("-prism_replicate_name", "-prn", help="name of the prism replicate that is being processed",
                        type=str, required=True)
    parser.add_argument("-davepool_mapping_file", "-dmf", help="mapping of analytes to pools and davepools",
                        type=str, required=True)
    parser.add_argument("-plate_map_path", "-pmp",
                        help="path to file containing plate map describing perturbagens used", type=str, required=True)
    parser.add_argument("-davepool_id_csv_filepath_pairs", "-dp_csv",
                        help="space-separated list of pairs of davepool_id and corresponding csv filepath for that davepool_id",
                        type=str, nargs="+", required=True)
    parser.add_argument("-cell_set_definition_file", "-csdf",
                        help="file containing cell set definition to use, overriding config file",
                        type=str, default=None, required=True)
    # These arguments are optional. Some may be superfluous now and might be removed.
    parser.add_argument("-verbose", '-v', help="Whether to print a bunch of output", action="store_true", default=False)
    parser.add_argument("-config_filepath", "-cfg", help="path to the location of the configuration file", type=str,
                        default=prism_pipeline.default_config_filepath)
    parser.add_argument("-ignore_assay_plate_barcodes", "-batmanify", help="list of assay plate barcodes that should be"
                        " ignored / excluded from the assemble", nargs="+", default=None)
    parser.add_argument("-plate_map_type", "-pmt", help="type of the plate map", choices=prism_metadata.plate_map_types,
                        default=prism_metadata.plate_map_type_CMap)
    parser.add_argument("-outfile", "-out", help="location to write gct", type=str,
                        default='')
    parser.add_argument("-machine_barcode", "-bc", help="machine barcode for querying assay plate table", type=int,
                        default=None, required=False)
    parser.add_argument("-truncate_to_plate_map", "-trunc", help="True or false, if true truncate data to fit framework of platemap provided",
                        action="store_true", default=False)

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


def build_davepool_id_to_cells_map(prism_cell_list):
    '''
    build one-to-many mapping between davepool ID and the corresponding PRISM cell lines that are within that davepool
    :param prism_cell_list:
    :return:
    '''
    r = {}
    for pc in prism_cell_list:
        if pc.davepool_id not in r:
            r[pc.davepool_id] = []
        r[pc.davepool_id].append(pc)

    return r


def build_data_by_cell(cells, davepool_data_obj):
    cell_to_median_data_map = {}
    median_wells = []
    cell_to_count_data_map = {}
    count_wells = []

    l = [(cell_to_median_data_map, median_wells, davepool_data_obj.median_headers, davepool_data_obj.median_data),
        (cell_to_count_data_map, count_wells, davepool_data_obj.count_headers, davepool_data_obj.count_data)]

    for (cell_data_map, wells, headers, data) in l:
        cell_header_map = {}
        for c in cells:
            if c.ignore == False:
                cell_header_map[c] = headers.index(c.analyte_id)
            cell_data_map[c] = []

        for d in data.keys():
            wells.append(d)
            for c in cells:
                if c in cell_header_map:
                    datum_index = cell_header_map[c]
                    datum = data[d][datum_index - 1]
                else:
                    datum = float('nan')

                cell_data_map[c].append(datum)

    return (DataByCell(cell_to_median_data_map, median_wells), DataByCell(cell_to_count_data_map, count_wells))


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


def process_data(davepool_data_objects, davepool_id_to_cells_map):
    logger.debug("davepool_id_to_cells_map.keys():  {}".format(davepool_id_to_cells_map.keys()))

    authoritative_well_list = []
    all_median_data_by_cell = DataByCell(cell_data_map={})
    all_count_data_by_cell = DataByCell(cell_data_map={})

    for dd in davepool_data_objects:
        cells = davepool_id_to_cells_map[dd.davepool_id]
        logger.debug("pools:  {}".format(cells))

        (median_data_by_cell), (count_data_by_cell) = build_data_by_cell(cells, dd)

        if len(authoritative_well_list) == 0:
            authoritative_well_list = median_data_by_cell.well_list
        else:
            assert authoritative_well_list == median_data_by_cell.well_list, (authoritative_well_list,
                                                                              median_data_by_cell.well_list)
            assert authoritative_well_list == count_data_by_cell.well_list, (authoritative_well_list,
                                                                             count_data_by_cell.well_list)

        combine_maps_with_checks(median_data_by_cell.cell_data_map, all_median_data_by_cell.cell_data_map)
        combine_maps_with_checks(count_data_by_cell.cell_data_map, all_count_data_by_cell.cell_data_map)

    all_median_data_by_cell.well_list = authoritative_well_list
    all_count_data_by_cell.well_list = authoritative_well_list

    return (all_median_data_by_cell, all_count_data_by_cell)


def build_gctoo(prism_replicate_name, perturbagen_list, data_by_cell):


    #build column metadata dataframe:
    def column_ID_builder(perturbagen):
        return prism_replicate_name + ":" + perturbagen.well_id

    col_metadata_df = prism_metadata.convert_objects_to_metadata_df(column_ID_builder, perturbagen_list,
                                                                             {"well_id":"pert_well"})
    for col_annot in _remove_col_annotations:
        if col_annot in col_metadata_df.columns:
            col_metadata_df.drop(col_annot, axis=1, inplace=True)

    col_metadata_df.sort_index(inplace=True)
    logger.info("my_gctoo.col_metadata_df.shape:  {}".format(col_metadata_df.shape))
    logger.debug("my_gctoo.col_metadata_df:  {}".format(col_metadata_df))
    ########################################

    #build row metadata dataframe:
    def row_ID_builder(prism_cell_obj):
        return prism_cell_obj.id

    row_metadata_df = prism_metadata.convert_objects_to_metadata_df(row_ID_builder,
        data_by_cell.cell_data_map.keys(), {})

    for row_annot in _remove_row_annotations:
        if row_annot in row_metadata_df.columns:
            row_metadata_df.drop(row_annot, axis=1, inplace=True)

    row_metadata_df.sort_index(inplace=True)
    logger.info("my_gctoo.row_metadata_df.shape:  {}".format(row_metadata_df.shape))
    logger.debug("my_gctoo.row_metadata_df:  {}".format(row_metadata_df))
    ########################################

    #build data dataframe - will have the cells as columns and rows as wells, and then transpose
    #start by building mapping between cell ID and corresponding data for that cell, to populate dataframe columns
    cell_id_data_map = {}
    for (c, data) in data_by_cell.cell_data_map.items():
        id = row_ID_builder(c)
        cell_id_data_map[id] = data
    logger.debug("cell_id_data_map:  {}".format(cell_id_data_map))

    #build index for the rows of the dataframe using what is actually the column ID (since it will be transposed)
    well_perturbagen_map = {}
    for p in perturbagen_list:
        well_perturbagen_map[p.well_id] = p
    data_df_column_ids = []
    for w in data_by_cell.well_list:
        p = well_perturbagen_map[w]
        data_df_column_ids.append(column_ID_builder(p))

    data_df = build_gctoo_data_df(cell_id_data_map, data_df_column_ids)
    ########################################

    my_gctoo = GCToo.GCToo(data_df=data_df, row_metadata_df=row_metadata_df, col_metadata_df=col_metadata_df)

    return my_gctoo


def build_gctoo_data_df(cell_id_data_map, data_df_column_ids):
    '''
    build the pandas dataframe that will be used for the data_df part of a gctoo object
    :param cell_id_data_map:
    :param data_df_column_ids:
    :return:
    '''
    data_df = pandas.DataFrame(cell_id_data_map, index=data_df_column_ids).T
    data_df.sort_index(axis=0, inplace=True)
    data_df.sort_index(axis=1, inplace=True)
    data_df.replace("", value=_NaN, inplace=True)
    logger.info("data_df.shape:  {}".format(data_df.shape))
    logger.debug("data_df:  {}".format(data_df))

    return data_df


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


def build_prism_cell_list(config_filepath, assay_plates, cell_set_definition_file, davepool_mapping_file):
    '''
    read PRISM cell line meta data from file specified in config file (at config_filepath), then associate with
    assay_plate based on pool ID.  Check for cell pools that are not associated with any assay plate
    :param config_filepath:
    :param assay_plates:
    :param cell_set_definition_file:
    :return:
    '''
    cp = ConfigParser.RawConfigParser()
    cp.read(config_filepath)

    prism_cell_list_items = cp.items(_prism_cell_config_file_section)
    davepool_mapping_items = cp.items(_davepool_analyte_mapping_file_section)

    prism_cell_list = prism_metadata.read_prism_cell_from_file(cell_set_definition_file, prism_cell_list_items)

    davepool_mapping = prism_metadata.read_prism_cell_from_file(davepool_mapping_file, davepool_mapping_items)

    pool_id_assay_plate_map = {}

    for ap in assay_plates:
        pool_id_assay_plate_map[ap.pool_id] = ap

    pool_id_without_assay_plate = set()
    cell_list_id_not_in_davepool_mapping = set()

    # Assign davepool mapping info to respective cell IDs
    cell_id_davepool_map = {}
    for dp in davepool_mapping:
        cell_id_davepool_map[dp.id] = dp

    for pc in prism_cell_list:
        if pc.id in cell_id_davepool_map.keys():
            cell_davepool = cell_id_davepool_map[pc.id]
            pc.analyte_id = cell_davepool.analyte_id
            pc.davepool_id = cell_davepool.davepool_id
            if pc.pool_id != cell_davepool.pool_id:
                raise Exception ("Cell set pool id does not match davepool mapping pool id at cell id {}".format(pc.id))
        else:
            cell_list_id_not_in_davepool_mapping.add(pc.id)
        if pc.pool_id in pool_id_assay_plate_map.keys():
            assay_plate = pool_id_assay_plate_map[pc.pool_id]
            pc.assay_plate_barcode = assay_plate.assay_plate_barcode
            pc.det_plate = assay_plate.det_plate
            pc.det_plate_scan_time = assay_plate.det_plate_scan_time
            pc.ignore = assay_plate.ignore

        else:
            pool_id_without_assay_plate.add(pc.pool_id)


    if len(pool_id_without_assay_plate) > 0:
        pool_id_without_assay_plate = list(pool_id_without_assay_plate)
        pool_id_without_assay_plate.sort()
        message1 = ("some pools were not found in the list of assay plates - pool_id_without_assay_plate:  {}".format(
            pool_id_without_assay_plate))
        logger.error(message1)
        #raise Exception ("assemble build_prism_cell_list " + message1)
    if len(cell_list_id_not_in_davepool_mapping) > 0:
        cell_list_id_not_in_davepool_mapping = list(cell_list_id_not_in_davepool_mapping)
        cell_list_id_not_in_davepool_mapping.sort()
        message2 = ("some cell ids were found in the cell set but not in the davepool mapping - IDs: {}".format(cell_list_id_not_in_davepool_mapping))
        raise Exception ("assemble build_prism_cell_list " + message2)
    return prism_cell_list


def build_assay_plates(davepool_data_objects, ignore_assay_plate_barcodes, machine_barcode, cursor):
    '''
    read all assay plate meta data from provided file at plates_mapping_path then remove assay plates whose det_plate
    does not match those in the davepool_data_objects.  Add additional metadata to assay_plates indicating the
    scan time and if they should be ignored (based if the assay plate barcode is in ignore_assay_plate_barcodes)
    :param plates_mapping_path: path to file that contains
    :param config_filepath:
    :param davepool_data_objects:
    :param ignore_assay_plate_barcodes:
    :return:
    '''

    assay_plates = []
    #TODO COULD HAVE automationify get machine barcode as well as davepool names
    det_plate_davepool_data_objects_map = {}
    for dpdo in davepool_data_objects:
        filename = os.path.basename(dpdo.csv_filepath)
        det_plate = ".".join(filename.split(".")[0:-1])
        det_plate_davepool_data_objects_map[det_plate] = dpdo
        #TODO get rid of this
        if machine_barcode is None:
            cursor.execute("""SELECT machine_barcode FROM plate WHERE det_plate = %s""", (det_plate,))
            machine_barcode = cursor.fetchone()
        det_plate_assay_plates = ap_orm.get_assay_plates(cursor, machine_barcode[0])
        assay_plates = assay_plates.extend(det_plate_assay_plates)

    logger.info("len(all_assay_plates):  {}".format(len(assay_plates)))

    logger.info("det_plate_davepool_data_objects_map.keys():  {}".format(det_plate_davepool_data_objects_map.keys()))

    #add scan time to assay_plate metadata, and indicate if the assay plate should be ignored
    for ap in assay_plates:
        ap.det_plate_scan_time = det_plate_davepool_data_objects_map[ap.det_plate].csv_datetime

        ap.ignore = ap.assay_plate_barcode in ignore_assay_plate_barcodes

    return assay_plates


def truncate_data_objects_to_plate_map(davepool_data_objects, all_perturbagens, truncate_to_platemap):

    platemap_well_list = set([p.well_id for p in all_perturbagens])

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


def main(args, all_perturbagens=None):

    db = mysql_utils.DB(host='localhost').db
    cursor = db.cursor()
    #TODO factor out plate_map_path
    if all_perturbagens is None:
        all_perturbagens = prism_metadata.build_perturbagens_from_file(args.plate_map_path, args.plate_map_type, args.config_filepath)

    args.ignore_assay_plate_barcodes = set(args.ignore_assay_plate_barcodes) if args.ignore_assay_plate_barcodes is not None else set()

    #read actual data from relevant csv files, associate it with davepool ID
    davepool_id_csv_list = build_davepool_id_csv_list(args.davepool_id_csv_filepath_pairs)
    davepool_data_objects = read_davepool_data_objects(davepool_id_csv_list)

    #read assay plate meta data relevant to current set of csv files / davepools
    import pdb
    pdb.set_trace()
    #TODO build this in the API
    assay_plates = build_assay_plates(davepool_data_objects, args.ignore_assay_plate_barcodes, args.machine_barcode, cursor)
    logger.info("len(assay_plates):  {}".format(len(assay_plates)))

    #read PRISM cell line metadata from file specified in config file, and associate with assay_plate metadata
    prism_cell_list = build_prism_cell_list(args.config_filepath, assay_plates, args.cell_set_definition_file, args.davepool_mapping_file)
    logger.info("len(prism_cell_list):  {}".format(len(prism_cell_list)))

    #build one-to-many mapping between davepool ID and the multiple PRISM cell lines that are within that davepool
    davepool_id_to_cells_map = build_davepool_id_to_cells_map(prism_cell_list)

    #truncate csv to plate map size if indicated by args.truncate_to_plate_map
    truncate_data_objects_to_plate_map(davepool_data_objects, all_perturbagens, args.truncate_to_plate_map)

    (all_median_data_by_cell, all_count_data_by_cell) = process_data(davepool_data_objects, davepool_id_to_cells_map)

    median_outfile = os.path.join(args.outfile, args.prism_replicate_name + "_MEDIAN.gct")
    median_gctoo = build_gctoo(args.prism_replicate_name, all_perturbagens, all_median_data_by_cell)
    write_gctoo.write(median_gctoo, median_outfile, data_null=_NaN, filler_null=_null)

    count_outfile = os.path.join(args.outfile, args.prism_replicate_name + "_COUNT.gct")
    count_gctoo = build_gctoo(args.prism_replicate_name, all_perturbagens, all_count_data_by_cell)
    write_gctoo.write(count_gctoo, count_outfile, data_null=_NaN, filler_null=_null)

    validate_prism_gct.check_headers(median_outfile)
    validate_prism_gct.check_headers(count_outfile)

    db.close()


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    logger.debug("args:  {}".format(args))

    main(args)