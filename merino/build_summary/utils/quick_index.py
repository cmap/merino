import sys
import os
import logging
import argparse
from string import digits


import merino.setup_logger as setup_logger
import merino.build_summary.make_gallery as galleries

base_path = '/cmap/obelix/pod/custom/'

logger = logging.getLogger(setup_logger.LOGGER_NAME)

def build_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-plate_grp_path', help='plate to grp of plates to gather into index', default=os.path.join(base_path, 'merino_plates.grp'))
    return parser


def main(args):
    list_of_plates = get_list_of_plates_from_grp(args.plate_grp_path)
    (card_dir_paths, gallery_paths, index_path) = set_up_paths(list_of_plates)

    def make_url(ref, name):
        ref = os.path.relpath(ref, index_path)
        return '<a target="_blank" href="{}">{}</a>'.format(ref, name)

    links = {list_of_plates[i]:make_url(x, list_of_plates[i]) for i, x in enumerate(gallery_paths)}
    pass_or_fail = get_pass_or_fail(card_dir_paths)

    headers = ['plate', 'pass/fail']

    table_tuples = []
    for i in list_of_plates:
        table_tuples.append((links[list_of_plates[i]], pass_or_fail[list_of_plates[i]]))

    logger.debug("info for table {}".format(table_tuples))

    made_gallery = galleries.mk_index(table_headers=headers, table_tuples=table_tuples, outfolder=index_path, project_name=extract_project_code(list_of_plates[0]))
    if made_gallery:
        logger.info("successfully made index")

def get_list_of_plates_from_grp(grp_path):
    list_of_plates = []
    with open(grp_path) as grp:
        for line in grp:
            list_of_plates.append(line.strip())
    return list_of_plates

def extract_project_code(plate_name):
    return plate_name.split("_")[0].translate(None, digits)

def set_up_paths(list_of_plates):

    card_dir_paths = { plate:os.path.join(base_path, extract_project_code(plate),'card', plate) for plate in list_of_plates}
    gallery_paths = [os.path.join(base_path, extract_project_code(plate), 'qc', plate, 'gallery.html') for plate in list_of_plates]
    index_path = os.path.join(base_path, extract_project_code(list_of_plates[0]), 'qc')
    return (card_dir_paths, gallery_paths, index_path)

def get_pass_or_fail(card_dir_paths):
    map = {}
    for plate, path in card_dir_paths.items():
        if os.path.exists(os.path.join(path, 'success.txt')):
            map[plate] = 1
        elif os.path.exists(os.path.join(path, 'failure.txt')):
            map[plate] = 2
    return map

if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])

    setup_logger.setup(verbose=True)
    main(args)
