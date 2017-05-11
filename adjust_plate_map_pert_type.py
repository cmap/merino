import logging
import setup_logger as setup_logger
import argparse
import sys
import pandas
import os.path


logger = logging.getLogger(setup_logger.LOGGER_NAME)

pert_type = "pert_type"
pert_type_poscon = "trt_poscon"


def build_parser():
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-verbose", "-v", help="Whether to print a bunch of output.", action="store_true", default=False)
	parser.add_argument("plate_map_files", nargs="+", help="plate map files to be adjusted", default=None)
	parser.add_argument("-poscon_pert_inames", nargs="+", help="pert_inames to use to identify poscons in the plate maps",
						default=["Bortezomib", "MG-132"])
	#parser.add_argument("-hostname", help="lims db host name", type=str, default="getafix-v")
	#parser.add_argument("-queue_choice", "-qc", help="which of the queues to work on - valid values are roast, brew, both", type=str,
	#	choices=["roast", "brew", "both"], default="both")
	#parser.add_argument("-add_to_queue", "-a", help="add the det_plate entries to the roast_queue", type=str, nargs="+", default=None)
	# To make --option1 and --option2 mutually exclusive, one can define mutually_exclusive_group in argparse, 
	# argparse asserts that the options added to the group are not used at the same time and throws exception if otherwise
    #mutually_exclusive_group = parser.add_mutually_exclusive_group()
    #mutually_exclusive_group.add_argument("--option1", action="store", dest="option1", help="provide argument for option1", default=None)
    #mutually_exclusive_group.add_argument("--option2", action="store", dest="option2", help="provide argument for option2", default=None)
	return parser


def read_plate_maps(plate_map_files):
	result = {}
	for pmf in plate_map_files:
		basename = os.path.basename(pmf)
		result[basename] = pandas.read_csv(pmf, sep="\t", dtype=str)

	return result

def write_plate_maps(filename_df_map):
	for (filename, df) in filename_df_map.items():
		df.to_csv(filename, sep="\t", index=False)

def replace_pert_type_for_poscons(filename_df_map, poscon_pert_inames):
	for platemap_df in filename_df_map.values():
		platemap_df.loc[platemap_df.pert_iname.isin(poscon_pert_inames), pert_type] = pert_type_poscon

def main(args):
	filename_df_map = read_plate_maps(args.plate_map_files)

	replace_pert_type_for_poscons(filename_df_map, args.poscon_pert_inames)

	write_plate_maps(filename_df_map)


if __name__ == "__main__":
	args = build_parser().parse_args(sys.argv[1:])

	setup_logger.setup(verbose=args.verbose)

	logger.debug("args:  {}".format(args))

	main(args)

