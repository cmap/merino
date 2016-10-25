import logging
import prism_pipeline.setup_logger as setup_logger
import ConfigParser
import prism_pipeline.prism_metadata as prism_metadata
import prism_pipeline.parse_data as parse_data


logger = logging.getLogger(setup_logger.LOGGER_NAME)


def read_assay_plates(tsv_filepaths, config_filepath, compound_plate_col_basename, poscon_plate_col_basename):
    cp = ConfigParser.RawConfigParser()
    cp.read(config_filepath)

    has_problems = False
    r = {}

    compound_plate_cols_first = None
    for tsv_filepath in tsv_filepaths:
        (headers, data) = parse_data.read_data(tsv_filepath)

        compound_plate_cols = get_all_related_plate_cols(compound_plate_col_basename, headers)
        if compound_plate_cols_first is None:
            compound_plate_cols_first = compound_plate_cols
        else:
            assert compound_plate_cols == compound_plate_cols_first, "the compound plate columns are not the same in the provided files - compound_plate_cols_first:  {}  compound_plate_cols:  {}".format(compound_plate_cols_first, compound_plate_cols)

        poscon_plate_cols = get_all_related_plate_cols(poscon_plate_col_basename, headers)

        header_map = parse_data.generate_header_map(headers, cp.items("column headers"), True)

        for (i, row) in enumerate(data):
            apb = row[header_map["assay_plate_barcode"]]
            pool_id = row[header_map["pool_id"]]

            compound_plates = build_plates_from_row(row, header_map, compound_plate_cols)

            poscon_plates = build_plates_from_row(row, header_map, poscon_plate_cols)

            if apb not in r:
                new_ap = prism_metadata.AssayPlate(assay_plate_barcode=apb, pool_id=pool_id)
                new_ap.compound_plates = compound_plates
                new_ap.poscon_plates = poscon_plates
                r[apb] = new_ap

            ap = r[apb]

            if ap.pool_id != pool_id:
                logger.error("assay plate with 2 different pool_id - i:  {}  ap.assay_plate_barcode:  {}  ap.pool_id:  {}  "
                                "pool_id:  {}".format(i, ap.assay_plate_barcode, ap.pool_id, pool_id))
                has_problems = True

            if compound_plates != ap.compound_plates:
                logger.error("assay plate has been specified with 2 different sets of compound plates - i:  {}"
                    "ap.assay_plate_barcode:  {}  ap.compound_plates:  {}  compound_plates:  {}  "
                    "compound_plate_col_basename:  {}".format(i, ap.assay_plate_barcode, ap.compound_plates,
                    compound_plates, compound_plate_col_basename))
                has_problems = True

            if poscon_plates != ap.poscon_plates:
                logger.error("assay plate has been specified with 2 different sets of poscon plates - i:  {}"
                    "ap.assay_plate_barcode:  {}  ap.poscon_plates:  {}  poscon_plates:  {}  "
                    "poscon_plate_col_basename:  {}".format(i, ap.assay_plate_barcode, ap.poscon_plates,
                    poscon_plates, poscon_plate_col_basename))
                has_problems = True

    if has_problems:
        raise Exception("assay_plate read_assay_plates problems encountered, see error messages")

    return (r.values(), compound_plate_cols)


def get_all_related_plate_cols(col_basename, headers):
    len_basename = len(col_basename)

    r = []
    for h in headers:
        if h == col_basename:
            r.append(h)
        elif h.startswith(col_basename) and len(h) > len_basename+1:
            rest_of_col = h[len_basename:]
            if rest_of_col[0] == "_":
                try:
                    index = int(rest_of_col[1:])
                    r.append(h)
                except ValueError:
                    pass

    return r


def build_plates_from_row(row, header_map, plate_cols):
    r = []
    for pc in plate_cols:
        plate_col_ind = header_map[pc]
        plate = row[plate_col_ind]
        r.append(plate)

    return tuple(r)
