import setup_logger
import logging
import davepool_data
import prism_metadata
import collections
import numpy
import argparse
import prism_pipeline
import sys
import os


logger = logging.getLogger(setup_logger.LOGGER_NAME)

DataByCell = collections.namedtuple("DataByCell", "cell_data_map well_list")
MatrixAndAnnots = collections.namedtuple("MatrixAndAnnots", "sorted_unique_cells sorted_unique_wells matrix")

_gct_version = "#1.3"

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
                                                                           " ignored / excluded from the assemble",
                        nargs="+", default=None)
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
    r = []
    for (dp_id, csv_filepath) in davepool_id_csv_list:
        pd = davepool_data.read_data(csv_filepath)
        pd.davepool_id = dp_id
        r.append(pd)

    return r


def build_davepool_id_to_cells_map(prism_cell_list):
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


def combine_maps_with_checks(source_map, dest_map):#TODO use this
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

    all_median_data_by_cell = []
    all_count_data_by_cell = []

    for dd in davepool_data_objects:
        cells = davepool_id_to_cells_map[dd.davepool_id]
        logger.debug("pools:  {}".format(cells))

        (median_data_by_cell), (count_data_by_cell) = build_data_by_cell(cells, dd)

        all_median_data_by_cell.append(median_data_by_cell)
        all_count_data_by_cell.append(count_data_by_cell)

    return (all_median_data_by_cell, all_count_data_by_cell)


def generate_sorted_unique_cells_and_wells(data_by_cells):
    unique_cells = set()
    unique_wells = set()
    for dbc in data_by_cells:
        unique_cells.update(dbc.cell_data_map.keys())
        unique_wells.update(dbc.well_list)

    sorted_unique_cells = list(unique_cells)

    def try_float_id(id):
        try:
            return float(id)
        except ValueError:
            return id

    sorted_unique_cells.sort(key=lambda c: try_float_id(c.id))

    sorted_unique_wells = list(unique_wells)
    sorted_unique_wells.sort()

    return (sorted_unique_cells, sorted_unique_wells)


def build_matrix_and_annotations(data_by_cells):
    (sorted_unique_cells, sorted_unique_wells) = generate_sorted_unique_cells_and_wells(data_by_cells)

    matrix = numpy.empty((len(sorted_unique_cells), len(sorted_unique_wells)))
    matrix[:] = numpy.NaN

    for dbc in data_by_cells:
        for (cell, data) in dbc.cell_data_map.items():
            logger.debug("data:  {}".format(data))

            row_ind = sorted_unique_cells.index(cell)
            row = matrix[row_ind]

            col_inds = [sorted_unique_wells.index(x) for x in dbc.well_list]
            logger.debug("col_inds:  {}".format(col_inds))

            row[col_inds] = data

    return MatrixAndAnnots(sorted_unique_cells, sorted_unique_wells, matrix)


def write_gct_version_and_size(file_handle, example_perturbagen, matrix_and_annots):
    file_handle.write(_gct_version + "\n")

    matrix_shape = matrix_and_annots.matrix.shape
    num_row_annots = len(matrix_and_annots.sorted_unique_cells[0].__dict__) - 1
    num_col_annots = len(example_perturbagen.__dict__)
    size_line = [matrix_shape[0], matrix_shape[1], num_row_annots, num_col_annots]

    file_handle.write("\t".join([str(x) for x in size_line]) + "\n")


def generate_column_headers(prism_replicate_name, example_cell, sorted_unique_wells):
    h = ["id"]

    cell_annot_order = [str(x) for x in example_cell.__dict__ if x != "id"]
    cell_annot_order.sort()
    h.extend(cell_annot_order)

    h.extend([prism_replicate_name + ":" + w for w in sorted_unique_wells])

    return (h, cell_annot_order)


def generate_perturbagen_annotation_header_block(perturbagen_list, sorted_unique_wells, num_cell_annot):
    annot_field_set = set()
    for p in perturbagen_list:
        annot_field_set.update(p.__dict__.keys())

    annot_fields = list(annot_field_set)
    annot_fields.sort()

    well_pert_map = {}
    for p in perturbagen_list:
        well_pert_map[p.well_id] = p

    r = []
    for af in annot_fields:
        row = []
        r.append(row)

        row.append(af)
        row.extend([_null for i in range(num_cell_annot)])

        for w in sorted_unique_wells:
            p = well_pert_map[w]
            value = p.__dict__[af] if af in p.__dict__ else _null
            value = value if value is not None else _null
            row.append(value)

    return r


def generate_row_annotation_and_data_block(matrix_and_annots, cell_annot_order):
    r = []
    for (i, c) in enumerate(matrix_and_annots.sorted_unique_cells):
        row = [c.id]
        r.append(row)

        for ca in cell_annot_order:
            value = c.__dict__[ca] if c.__dict__[ca] is not None else _null
            row.append(value)

        numerical_data = [_NaN if numpy.isnan(x) else x for x in matrix_and_annots.matrix[i]]
        row.extend(numerical_data)

    return r


def write_output_gct(output_filepath, prism_replicate_name, perturbagen_list,
                    matrix_and_annots):
    f = open(output_filepath, "w")
    write_gct_version_and_size(f, perturbagen_list[0], matrix_and_annots)

    (headers, cell_annot_order) = generate_column_headers(prism_replicate_name,
        matrix_and_annots.sorted_unique_cells[0],
        matrix_and_annots.sorted_unique_wells)

    f.write("\t".join(headers) + "\n")

    perturbation_annotation_header_block = generate_perturbagen_annotation_header_block(perturbagen_list,
        matrix_and_annots.sorted_unique_wells, len(cell_annot_order))

    for header_block_row in perturbation_annotation_header_block:
        row_output_strings = []
        for x in header_block_row:
            if isinstance(x, float):
                row_output_strings.append("%.2f" % x)
            else:
                row_output_strings.append(str(x))

        f.write("\t".join(row_output_strings) + "\n")

    row_annotation_and_data_block = generate_row_annotation_and_data_block(matrix_and_annots, cell_annot_order)
    for row in row_annotation_and_data_block:
        f.write("\t".join([str(x) for x in row]) + "\n")

    f.close()


def build_davepool_id_csv_list(davepool_id_csv_filepath_pairs):
    r = []
    for i in range(len(davepool_id_csv_filepath_pairs)/2):
        index = 2*i
        davepool_id = davepool_id_csv_filepath_pairs[index]
        csv_filepath = davepool_id_csv_filepath_pairs[index+1]
        r.append((davepool_id, csv_filepath))

    return r


def build_perturbagen_list(plate_map_path, config_filepath, assay_plates):
    assay_plate_barcodes = set([x.assay_plate_barcode for x in assay_plates if x.ignore == False])

    all_perturbagens = prism_metadata.read_perturbagen_from_file(plate_map_path, config_filepath)

    perts = [x for x in all_perturbagens if x.assay_plate_barcode in assay_plate_barcodes]

    prism_metadata.validate_perturbagens(perts)

    return perts


def build_prism_cell_list(config_filepath, assay_plates):
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
    all_assay_plates = prism_metadata.read_assay_plate_from_file(plates_mapping_path, config_filepath)

    det_plate_davepool_data_objects_map = {}
    for dpdo in davepool_data_objects:
        filename = os.path.basename(dpdo.csv_filepath)
        det_plate = filename.split(".")[0]
        det_plate_davepool_data_objects_map[det_plate] = dpdo

    assay_plates = [x for x in all_assay_plates if x.det_plate in det_plate_davepool_data_objects_map]
    for ap in assay_plates:
        ap.det_plate_scan_time = det_plate_davepool_data_objects_map[ap.det_plate].csv_datetime

        ap.ignore = ap.assay_plate_barcode in ignore_assay_plate_barcodes

    return assay_plates


def main(args):
    args.ignore_assay_plate_barcodes = set(args.ignore_assay_plate_barcodes) if args.ignore_assay_plate_barcodes is not None else set()

    davepool_id_csv_list = build_davepool_id_csv_list(args.davepool_id_csv_filepath_pairs)
    davepool_data_objects = read_davepool_data_objects(davepool_id_csv_list)

    assay_plates = build_assay_plates(args.plates_mapping_path, args.config_filepath, davepool_data_objects,
                                      args.ignore_assay_plate_barcodes)
    logger.info("len(assay_plates):  {}".format(len(assay_plates)))

    prism_cell_list = build_prism_cell_list(args.config_filepath, assay_plates)
    logger.info("len(prism_cell_list):  {}".format(len(prism_cell_list)))

    perturbagen_list = build_perturbagen_list(args.plate_map_path, args.config_filepath, assay_plates)
    logger.info("len(perturbagen_list):  {}".format(len(perturbagen_list)))

    davepool_id_to_cells_map = build_davepool_id_to_cells_map(prism_cell_list)

    (all_median_data_by_cell, all_count_data_by_cell) = process_data(davepool_data_objects, davepool_id_to_cells_map)

    median_matrix_and_annots = build_matrix_and_annotations(all_median_data_by_cell)
    count_matrix_and_annots = build_matrix_and_annotations(all_count_data_by_cell)

    write_output_gct(args.prism_replicate_name + "_MEDIAN.gct", args.prism_replicate_name, perturbagen_list,
                     median_matrix_and_annots)
    write_output_gct(args.prism_replicate_name + "_COUNT.gct", args.prism_replicate_name, perturbagen_list,
                     count_matrix_and_annots)


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    logger.debug("args:  {}".format(args))

    main(args)
