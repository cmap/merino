import logging
import prism_pipeline.setup_logger as setup_logger
import ConfigParser
import prism_pipeline.prism_metadata as prism_metadata
import prism_pipeline.parse_data as parse_data


logger = logging.getLogger(setup_logger.LOGGER_NAME)


def read_assay_plates(tsv_filepaths, config_filepath):
    cp = ConfigParser.RawConfigParser()
    cp.read(config_filepath)

    has_problems = False
    r = {}

    for tsv_filepath in tsv_filepaths:
        (headers, data) = parse_data.read_data(tsv_filepath)

        header_map = parse_data.generate_header_map(headers, cp.items("column headers"), True)

        for (i, row) in enumerate(data):
            apb = row[header_map["assay_plate_barcode"]]
            pool_id = row[header_map["pool_id"]]
            compound_plate_map_name = row[header_map["compound_plate_map_name"]]
            poscon_plate_map_name = row[header_map["poscon_plate_map_name"]]
            if apb not in r:
                new_ap = prism_metadata.AssayPlate(assay_plate_barcode=apb, pool_id=pool_id)
                new_ap.compound_plate = {compound_plate_map_name:poscon_plate_map_name}
                r[apb] = new_ap

            ap = r[apb]

            if ap.pool_id != pool_id:
                logger.error("assay plate with 2 different pool_id - i:  {}  ap.assay_plate_barcode:  {}  ap.pool_id:  {}  "
                                "pool_id:  {}".format(i, ap.assay_plate_barcode, ap.pool_id, pool_id))
                has_problems = True

            if compound_plate_map_name in ap.compound_plate:
                if ap.compound_plate[compound_plate_map_name] != poscon_plate_map_name:
                    logger.error("assay plate with 2 different compound_plate_map_name for the same compound_plate_barcode - "
                                "i:  {}  ap.assay_plate_barcode:  {}  ap.compound_plate:  {}  compound_plate_map_name:  {}  "
                                    "poscon_plate_map_name:  {}".format(i, ap.assay_plate_barcode, ap.compound_plate,
                                                                      compound_plate_map_name, poscon_plate_map_name))
                    has_problems = True
            else:
                ap[compound_plate_map_name] = poscon_plate_map_name

    if has_problems:
        raise Exception("assay_plate read_assay_plates problems encountered, see error messages")

    return r.values()

