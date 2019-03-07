import sys
import os
import argparse

import caldaia.utils.orm.lims_plate_orm as lpo

import merino.assemble.assemble as assemble
import merino.card.card as card
import merino.build_summary.plate_summary as plate_qc

def build_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-plate_name', '-p', help='name of plate to process up to card and run plate qc')

    return parser

def main(args):
    plate_entry = lpo.get_by_det_plate(args.plate_name)
    (project_dir, jcsv_path, plate_map_path, assemble_out_path, qc_out_path) = build_paths(plate_entry)

    if plate_entry:
        assemble_args = assemble.build_parser().parse_args(['-pmp', plate_map_path, '-csv', jcsv_path, '-out', assemble_out_path])
        assemble.main(assemble_args)
        card_args = card.build_parser().parse_args(['-proj_dir', project_dir, '-plate_name', args.plate_name])
        card.main(card_args)
        plate_qc.build_parser().parse_args(['-project_folder', project_dir, '-plate_name', args.plate_name, '-qc', qc_out_path])

def build_paths(plate_lpo):
    project_dir = os.path.join('/cmap/obelix/pod/custom/', plate_lpo.project_code)
    jcsv_path = os.path.join(project_dir, 'lxb', plate_lpo.det_plate + '.jcsv')
    plate_map_path = os.path.join(project_dir, 'map_src', plate_lpo.pert_plate +'.src')
    assemble_out_path = os.path.join(project_dir, 'assemble', plate_lpo.det_plate)
    qc_out_path = os.path.join(project_dir, 'qc')
    return (project_dir, jcsv_path, plate_map_path, assemble_out_path, qc_out_path)

if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    main(args)