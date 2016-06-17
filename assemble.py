import setup_logger
import logging
import davepool_data
import prism_metadata
import numpy
import argparse
import prism_pipeline
import sys
import os
import GCToo.GCToo as GCToo
import pandas
import GCToo.write_gctoo as write_gctoo


logger = logging.getLogger(setup_logger.LOGGER_NAME)

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
    parser.add_argument("-verbose", '-v', help="Whether to print a bunch of output", action="store_true", default=False)
    parser.add_argument("-config_filepath", help="path to the location of the configuration file", type=str,
                        default=prism_pipeline.default_config_filepath)
    parser.add_argument("prism_replicate_name", help="name of the prism replicate that is being processed", type=str)
    parser.add_argument("plate_map_path", help="path to file containing plate map describing perturbagens used", type=str)
    parser.add_argument("plates_mapping_path", help="path to file containing the mapping between assasy plates and det_plates",
                        type=str)
    parser.add_argument("davepool_id_csv_filepath_pairs",
                        help="space-separated list of pairs of davepool_id and corresponding csv filepath for that davepool_id",
                        type=str, nargs="+")
    parser.add_argument("-ignore_assay_plate_barcodes", "-batmanify", help="list of assay plate barcodes that should be"
                        " ignored / excluded from the assemble", nargs="+", default=None)
    parser.add_argument("-plate_map_type", "-pmt", help="type of the plate map", choices=prism_metadata.plate_map_types,
                        default=prism_metadata.plate_map_type_CM)
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

        for d in data:
            wells.append(parse_location_to_well(d[0]))
            for c in cells:
                if c in cell_header_map:
                    datum_index = cell_header_map[c]
                    datum = d[datum_index]
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
    my_gctoo = GCToo.GCToo()

    #build column metadata dataframe:
    def column_ID_builder(perturbagen):
        return prism_replicate_name + ":" + perturbagen.well_id

    my_gctoo.col_metadata_df = prism_metadata.convert_objects_to_metadata_df(column_ID_builder, perturbagen_list,
                                                                             {"well_id":"pert_well"})
    for col_annot in _remove_col_annotations:
        if col_annot in my_gctoo.col_metadata_df.columns:
            my_gctoo.col_metadata_df.drop(col_annot, axis=1, inplace=True)

    my_gctoo.col_metadata_df.sort_index(inplace=True)
    logger.info("my_gctoo.col_metadata_df.shape:  {}".format(my_gctoo.col_metadata_df.shape))
    logger.debug("my_gctoo.col_metadata_df:  {}".format(my_gctoo.col_metadata_df))
    ########################################

    #build row metadata dataframe:
    def row_ID_builder(prism_cell_obj):
        return prism_cell_obj.id

    my_gctoo.row_metadata_df = prism_metadata.convert_objects_to_metadata_df(row_ID_builder,
        data_by_cell.cell_data_map.keys(), {})

    for row_annot in _remove_row_annotations:
        if row_annot in my_gctoo.row_metadata_df.columns:
            my_gctoo.row_metadata_df.drop(row_annot, axis=1, inplace=True)

    my_gctoo.row_metadata_df.sort_index(inplace=True)
    logger.info("my_gctoo.row_metadata_df.shape:  {}".format(my_gctoo.row_metadata_df.shape))
    logger.debug("my_gctoo.row_metadata_df:  {}".format(my_gctoo.row_metadata_df))
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

    my_gctoo.data_df = pandas.DataFrame(cell_id_data_map, index=data_df_column_ids).T
    my_gctoo.data_df.sort_index(axis=0, inplace=True)
    my_gctoo.data_df.sort_index(axis=1, inplace=True)
    logger.info("my_gctoo.data_df.shape:  {}".format(my_gctoo.data_df.shape))
    logger.debug("my_gctoo.data_df:  {}".format(my_gctoo.data_df))
    ########################################

    return my_gctoo


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


def build_perturbagen_list(plate_map_path, config_filepath, assay_plates):
    '''
    read all perturbagens in, and then keep only those whose assay_plate_barcode matches one of the already loaded
    assay plate barcodes.  Validate the remaining perturbagens.
    :param plate_map_path:
    :param config_filepath:
    :param assay_plates:
    :return:
    '''
    assay_plate_barcodes = set([x.assay_plate_barcode for x in assay_plates if x.ignore == False])

    all_perturbagens = prism_metadata.build_perturbagens_from_file(plate_map_path, prism_metadata.plate_map_type_CM,
                                                                   config_filepath)

    all_assay_plate_perts = [x for x in all_perturbagens if x.assay_plate_barcode in assay_plate_barcodes]

    validated_unique_perts = prism_metadata.validate_perturbagens(all_assay_plate_perts).values()

    return validated_unique_perts


def build_prism_cell_list(config_filepath, assay_plates):
    '''
    read PRISM cell line meta data from file specified in config file (at config_filepath), then associate with
    assay_plate based on pool ID.  Check for cell pools that are not associated with any assay plate
    :param config_filepath:
    :param assay_plates:
    :return:
    '''
    prism_cell_list = prism_metadata.read_prism_cell_from_file(config_filepath)

    pool_id_assay_plate_map = {}
    for ap in assay_plates:
        pool_id_assay_plate_map[ap.pool_id] = ap

    pool_id_without_assay_plate = set()

    for pc in prism_cell_list:
        if pc.pool_id in pool_id_assay_plate_map:
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
        logger.warning("some pools were not found in the list of assay plates - pool_id_without_assay_plate:  {}".format(
            pool_id_without_assay_plate))

    return prism_cell_list


def build_assay_plates(plates_mapping_path, config_filepath, davepool_data_objects, ignore_assay_plate_barcodes):
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
    all_assay_plates = prism_metadata.read_assay_plate_from_file(plates_mapping_path, config_filepath)

    #parse the csv filename to get the det_plate, build a map between det_plate and davepool object
    det_plate_davepool_data_objects_map = {}
    for dpdo in davepool_data_objects:
        filename = os.path.basename(dpdo.csv_filepath)
        det_plate = filename.split(".")[0]
        det_plate_davepool_data_objects_map[det_plate] = dpdo

    #only keep assay plates whose det_plate matches one of the loaded davepool
    assay_plates = [x for x in all_assay_plates if x.det_plate in det_plate_davepool_data_objects_map]

    #add scan time to assay_plate metadata, and indicate if the assay plate should be ignored
    for ap in assay_plates:
        ap.det_plate_scan_time = det_plate_davepool_data_objects_map[ap.det_plate].csv_datetime

        ap.ignore = ap.assay_plate_barcode in ignore_assay_plate_barcodes

    return assay_plates


def main(args):
    args.ignore_assay_plate_barcodes = set(args.ignore_assay_plate_barcodes) if args.ignore_assay_plate_barcodes is not None else set()

    #read actual data from relevant csv files, associate it with davepool ID
    davepool_id_csv_list = build_davepool_id_csv_list(args.davepool_id_csv_filepath_pairs)
    davepool_data_objects = read_davepool_data_objects(davepool_id_csv_list)

    #read assay plate meta data relevant to current set of csv files / davepools
    assay_plates = build_assay_plates(args.plates_mapping_path, args.config_filepath, davepool_data_objects,
                                      args.ignore_assay_plate_barcodes)
    logger.info("len(assay_plates):  {}".format(len(assay_plates)))

    #read PRISM cell line metadata from file specified in config file, and associate with assay_plate metadata
    prism_cell_list = build_prism_cell_list(args.config_filepath, assay_plates)
    logger.info("len(prism_cell_list):  {}".format(len(prism_cell_list)))

    #read in all the perturbagens but restrict to those that were on the provided assay_plates
    perturbagen_list = build_perturbagen_list(args.plate_map_path, args.config_filepath, assay_plates)
    logger.info("len(perturbagen_list):  {}".format(len(perturbagen_list)))

    #build one-to-many mapping between davepool ID and the multiple PRISM cell lines that are within that davepool
    davepool_id_to_cells_map = build_davepool_id_to_cells_map(prism_cell_list)

    (all_median_data_by_cell, all_count_data_by_cell) = process_data(davepool_data_objects, davepool_id_to_cells_map)

    median_gctoo = build_gctoo(args.prism_replicate_name, perturbagen_list, all_median_data_by_cell)
    write_gctoo.write(median_gctoo, args.prism_replicate_name + "_MEDIAN.gct", data_null=_NaN, filler_null=_null)

    count_gctoo = build_gctoo(args.prism_replicate_name, perturbagen_list, all_count_data_by_cell)
    write_gctoo.write(count_gctoo, args.prism_replicate_name + "_COUNT.gct", data_null=_NaN, filler_null=_null)


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    logger.debug("args:  {}".format(args))

    main(args)
