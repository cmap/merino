import logging
import setup_logger
import argparse
import sys
import assay_plate
import ConfigParser
import prism_det_plate
import math


logger = logging.getLogger(setup_logger.LOGGER_NAME)


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("project_id", help="project_id will be used as prefix for pert_plate", type=str)
    parser.add_argument("pert_plate_start_number", help="starting number to use for the CMAP pert plate", type=int)
    parser.add_argument("cellset_id", help="CellSet ID used for these plates", type=str)
    parser.add_argument("input_files",
                        help="""input file(s) containing output of compound management Pipeline Pilot Webport
                              "CMAP_Assay_Plate_Compound_Map_Names_Cell_Line_Metadata_from_RunID""",
                        type=str, nargs="+")
    parser.add_argument("-output_file", help="file to write the output to", type=str, default="output.txt")
    parser.add_argument('-verbose', '-v', help='Whether to print a bunch of output.', action='store_true',
                        default=False)
    parser.add_argument("-config_filepath", help="path to configuration file", type=str,
                        default="register_prism_plates.cfg")
    parser.add_argument("-hostname", help="lims db host name", type=str, default="getafix-v")
    parser.add_argument("-dont_assert_plate_numbers_match", 
        help="instead of asserting that len(assay_plates) %% len(pool_id_sorted) == 0 and len(assay_plates) %% len(compound_plates_sorted) == 0 just issue warnings",
        action="store_true", default=False)
    return parser


def load_and_validate_assay_plates(input_fles, config_filepath):
    assay_plates = assay_plate.read_assay_plates(args.input_files, args.config_filepath)
    logger.info("len(assay_plates):  {}".format(len(assay_plates)))

    r = [x for x in assay_plates if len(x.compound_plate) > 1]
    logger.info("assay plates with multiple compound plates:  {}".format(r))
    assert len(r) == 0
    for ap in assay_plates:
        for x in ap.compound_plate:
            ap.compound_plate_map_name = x

    assay_plates.sort(key=lambda x: x.assay_plate_barcode)

    return assay_plates


def determine_pool_ids(assay_plates):
    r = set([x.pool_id for x in assay_plates])
    logger.info("number of pool_id found:  {}".format(len(r)))
    pool_id_sorted = list(r)
    pool_id_sorted.sort()
    logger.info("pool_id_sorted:  {}".format(pool_id_sorted))
    return pool_id_sorted


def determine_compound_plates(assay_plates):
    r = set([x.compound_plate_map_name for x in assay_plates])
    compound_plates_sorted = list(r)
    compound_plates_sorted.sort()
    logger.info("number of compound plates:  {}".format(len(compound_plates_sorted)))
    logger.info("compound_plates_sorted:  {}".format(compound_plates_sorted))
    return compound_plates_sorted


def build_pert_plate_mapping(compound_plates_sorted, pert_plate_start_number, project_id):
    pert_plate_mapping = {}
    for (i, cp) in enumerate(compound_plates_sorted):
        pert_plate = project_id + str(i + pert_plate_start_number).zfill(3)
        pert_plate_mapping[cp] = pert_plate
    logger.info("pert_plate_mapping:  {}".format(pert_plate_mapping))
    return pert_plate_mapping


def build_pool_id_to_davepool_id_mapping(config_filepath):
    config = ConfigParser.RawConfigParser()
    config.read(config_filepath)
    pool_id_to_davepool_id_map = {}
    for (pool_id, davepool_id) in config.items("pool_id to davepool_id mapping"):
        pool_id_to_davepool_id_map[pool_id.upper()] = davepool_id
    logger.debug("pool_id_to_davepool_id_map:  {}".format(pool_id_to_davepool_id_map))
    return pool_id_to_davepool_id_map


def main(args):
    assay_plates = load_and_validate_assay_plates(args.input_files, args.config_filepath)

    pool_id_sorted = determine_pool_ids(assay_plates)

    if (len(assay_plates) % len(pool_id_sorted)) != 0:
        msg = "number of assay plates is not an even multiple of number of pools - remainder len(assay_plates) % len(pool_id_sorted):  {}".format(len(assay_plates) % len(pool_id_sorted))
        if args.dont_assert_plate_numbers_match:
            logger.warning(msg)
        else:
            raise Exception("register_prism_plates main " + msg)

    compound_plates_sorted = determine_compound_plates(assay_plates)

    if (len(assay_plates) % len(compound_plates_sorted)) != 0:
        msg = "number of assay plates is not an even multiple of number of compound plates - remainder len(assay_plates) % len(compound_plates_sorted):  {}".format(len(assay_plates) % len(compound_plates_sorted))
        if args.dont_assert_plate_numbers_match:
            logger.warning(msg)
        else:
            raise Exception("register_prism_plates main " + msg)

    num_repos = None
    if args.dont_assert_plate_numbers_match:
        num_reps = int(math.ceil(float(len(assay_plates)) / float(len(compound_plates_sorted)) / float(len(pool_id_sorted))))
    else:
        num_reps = len(assay_plates) / len(compound_plates_sorted) / len(pool_id_sorted)
    logger.info("num_reps:  {}".format(num_reps))

    pert_plate_mapping = build_pert_plate_mapping(compound_plates_sorted, args.pert_plate_start_number, args.project_id)

    pool_id_to_davepool_id_map = build_pool_id_to_davepool_id_mapping(args.config_filepath)

    rep_counter = {}
    rows = []
    for ap in assay_plates:
        pert_plate = pert_plate_mapping[ap.compound_plate_map_name]

        davepool_id = pool_id_to_davepool_id_map[ap.pool_id]

        k = (pert_plate, ap.pool_id)
        if k not in rep_counter:
            rep_counter[k] = {}

        cur_rc = rep_counter[k]
        if len(cur_rc) == 0:
            cur_rc[ap.assay_plate_barcode] = 1
        else:
            if ap.assay_plate_barcode in cur_rc:
                raise Exception("assay plate barcode appeared twice in rep_counter - cur_rc:  {}  ap:  {}".format(cur_rc,
                                                                                                                  ap))
            else:
                cur_rc[ap.assay_plate_barcode] = max(cur_rc.values()) + 1

        rep_num = cur_rc[ap.assay_plate_barcode]
        rep_str = "X" + str(rep_num)
        det_plate = pert_plate + "_" + davepool_id + "_" + rep_str
        prism_replicate = pert_plate + "_" + args.cellset_id + "_" + rep_str
        true_pool_id = ap.pool_id[:4]
        rows.append((ap.assay_plate_barcode, true_pool_id, ap.compound_plate_map_name, davepool_id, det_plate,
                     prism_replicate))

    rows.sort(key=lambda x: (x[5], x[1]))

    rows.insert(0, ("assay_plate_barcode", "pool_id", "compound_plate_map_name", "davepool_id", "det_plate", "prism_replicate"))

    f = open(args.output_file, "w")
    for r in rows:
        f.write("\t".join(r) + "\n")
    f.close()


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])

    setup_logger.setup(verbose=args.verbose)

    logger.debug("args:  {}".format(args))

    main(args)
