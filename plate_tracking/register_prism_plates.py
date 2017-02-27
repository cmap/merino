import logging
from prism_pipeline import setup_logger
import argparse
import sys
import utils.mysql_utils as mysql_utils
import utils.orm.assay_plates_orm as ap_orm
import plate_tracking_metadata
import ConfigParser
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
    parser.add_argument("-compound_plate_col_basename", "-cpcb", help="base name to use when identifying the compound plate(s) " +
                        "e.g. for the default value of compound_map would use the column \"compound_map\" but would also look " +
                        " for columns compound_map_2, compound_map_3 etc.", type=str, default="compound_map")
    parser.add_argument("-poscon_plate_col_basename", "-ppcb", help="base name to use when identifying the poscon plate(s) " +
                        "e.g. for the default value of poscon_map would use the column \"poscon_map\" but would also look " +
                        " for columns poscon_map_2, poscon_map_3 etc.", type=str, default="poscon_map")
    parser.add_argument("-do_not_commit_to_db", "-dont",
                        help="Do not commit updates to the database",
                        action="store_true", default=False)
    return parser


def load_and_sort_assay_plates(input_files, config_filepath, compound_plate_col_basename, poscon_plate_col_basename):
    (assay_plates, compound_plate_cols) = plate_tracking_metadata.read_assay_plates(input_files, config_filepath, compound_plate_col_basename,
                                                             poscon_plate_col_basename)
    logger.info("len(assay_plates):  {}".format(len(assay_plates)))

    assay_plates.sort(key=lambda x: x.assay_plate_barcode)

    return (assay_plates, compound_plate_cols)


def determine_pool_ids(assay_plates):
    r = set([x.pool_id for x in assay_plates])
    logger.info("number of pool_id found:  {}".format(len(r)))
    pool_id_sorted = list(r)
    pool_id_sorted.sort()
    logger.info("pool_id_sorted:  {}".format(pool_id_sorted))
    return pool_id_sorted


def determine_compound_plates(assay_plates):
    r = set([x.compound_plates for x in assay_plates])
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

    db = mysql_utils.DB(host='localhost').db
    cursor = db.cursor()

    (assay_plates, compound_plate_cols) = load_and_sort_assay_plates(args.input_files, args.config_filepath, args.compound_plate_col_basename,
                                              args.poscon_plate_col_basename)

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

    num_reps = None
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
        pert_plate = pert_plate_mapping[ap.compound_plates]

        davepool_id = pool_id_to_davepool_id_map[ap.pool_id]

        k = (pert_plate, ap.pool_id)
        if k not in rep_counter:
            rep_counter[k] = {}

        cur_rc = rep_counter[k]
        if len(cur_rc) == 0:
            cur_rc[ap.assay_plate_barcode] = 1
        else:
            if ap.assay_plate_barcode in cur_rc:
                raise Exception("assay plate barcode appeared twice in rep_counter - cur_rc:  {}  ap:  {}".format(cur_rc, ap))
            else:
                cur_rc[ap.assay_plate_barcode] = max(cur_rc.values()) + 1

        rep_num = cur_rc[ap.assay_plate_barcode]
        rep_str = "X" + str(rep_num)
        det_plate = pert_plate + "_" + davepool_id + "_" + rep_str
        prism_replicate = pert_plate + "_" + args.cellset_id + "_" + rep_str
        true_pool_id = ap.pool_id[:4]

        cur_row = [ap.assay_plate_barcode, true_pool_id]
        cur_row.extend(ap.compound_plates)
        cur_row.extend([davepool_id, det_plate, prism_replicate])
        rows.append(cur_row)

    rows.sort(key=lambda x: (x[5], x[1]))

    headers = ["assay_plate_barcode", "pool_id"]
    headers.extend(compound_plate_cols)
    headers.extend(["davepool_id", "det_plate", "prism_replicate"])
    rows.insert(0, headers)
    dp_index = headers.index("det_plate")

    count = 0
    for r in rows[1:]:
        count += 1
        print count
        ap_orm.insert_into_plate_tracking(cursor, r[0], r[1], r[dp_index])

    if not args.do_not_commit_to_db:
        db.commit()

    f = open(args.output_file, "w")
    for r in rows:
        f.write("\t".join(r) + "\n")
    f.close()

    import pdb
    pdb.set_trace()

    db.close()

if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])

    setup_logger.setup(verbose=args.verbose)

    logger.debug("args:  {}".format(args))

    main(args)
