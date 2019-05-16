import sys
import os
import argparse
import logging
from string import digits

import merino.setup_logger as setup_logger
import merino.assemble.assemble as assemble
import merino.card.card as card
import merino.build_summary.plate_summary as plate_qc

logger = logging.getLogger(setup_logger.LOGGER_NAME)


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-plate_name', '-p', help='name of plate to process up to card and run plate qc')
    parser.add_argument("-rep_map", "-rm", help="Whether to use replicate level plate maps or not",
                        type=str, required=False, default="FALSE",choices=["TRUE", "FALSE"])

    return parser

def main(args):
    plate_entry = parse_plate_name(args.plate_name)

    (project_dir, jcsv_path, plate_map_path, assemble_out_path, qc_out_path) = build_paths(plate_entry, args.rep_map)

    assemble_args = assemble.build_parser().parse_args(['-pmp', plate_map_path, '-csv', jcsv_path, '-out', assemble_out_path])
    logger.info("Running assemble with args: {}".format(assemble_args))
    assemble.main(assemble_args)

    card_args = card.build_parser().parse_args(['-proj_dir', project_dir, '-plate_name', args.plate_name])
    logger.info("Running card with args: {}".format(card_args))
    card.main(card_args)

    gmt_path = '/cmap/data/vdb/merino/sensitivity_files/CORE_D1_48H_KJ100_DN_s25_n5127.gmt'

    plate_qc_args = plate_qc.build_parser().parse_args(['-project_folder', project_dir, '-plate_name', args.plate_name,
                                                        '-qc', qc_out_path, '-gmt', gmt_path])

    logger.info("Running PlateQC with args: {}".format(plate_qc_args))
    plate_qc.main(plate_qc_args)


def parse_plate_name(plate_name):
    pieces = plate_name.split("_")
    pert_plate = pieces[0]
    rep_id = pieces[3]
    project_code = pert_plate.translate(None, digits)

    return {"pert_plate": pert_plate, "project_code": project_code, "det_plate": plate_name, "rep_id": rep_id}

def build_paths(plate_entry, rep_map):
    project_dir = os.path.join('/cmap/obelix/pod/custom/', plate_entry['project_code'])
    jcsv_path = os.path.join(project_dir, 'lxb', plate_entry['det_plate'], plate_entry['det_plate'] + '.jcsv')

    if rep_map == "FALSE":
        plate_map_path = os.path.join(project_dir, 'map_src', plate_entry['pert_plate'] + '.src')

    elif rep_map == "TRUE":
        plate_map_path = os.path.join(project_dir, 'map_src', plate_entry['pert_plate'] + '.' + plate_entry['rep_id'] + '.src')

    assemble_out_path = project_dir
    qc_out_path = os.path.join(project_dir, 'qc')
    if not os.path.exists(qc_out_path):
        os.mkdir(qc_out_path)
    return (project_dir, jcsv_path, plate_map_path, assemble_out_path, qc_out_path)

if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])

    setup_logger.setup(verbose=True)

    main(args)