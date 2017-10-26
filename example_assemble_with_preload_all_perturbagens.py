import merino.assemble as assemble
import merino.setup_logger as setup_logger
import merino.prism_metadata as prism_metadata
import logging
import os
import glob
import shutil


logger = logging.getLogger(setup_logger.LOGGER_NAME)


def main():
	config_filepath = "2016-08-12_PCAL015-44_prism_pipeline.cfg"
	plate_mapping_filepath = "../map_src/2016-08-12_PCAL015-44_CS1.2_plates_mapping.txt"
	plate_map_filepath = "../map_src/2016-08-13_PCAL015-44_plate_maps_7159-03-A04-01-07-18.txt"

	base_assemble_dir = "../assemble/"
	base_lxb_dir = "../lxb/"

	all_perturbagens = assemble.read_all_perturbagens_from_file(plate_map_filepath, config_filepath,
																prism_metadata.plate_map_type_CM)
	
	parser = assemble.build_parser()
	args_list = ["-config_filepath", config_filepath, None, "", plate_mapping_filepath, "DP7", None, "DP8", None]

	for pert_plate_index in range(15,45):
		for rep_num in range(1,4):
			pert_plate = "PCAL0{}".format(pert_plate_index)
			replicate = "X{}".format(rep_num)
			prism_replicate = "{}_CS1.2_{}".format(pert_plate, replicate)
			dp7 = "{}_DP7_{}".format(pert_plate, replicate)
			dp8 = "{}_DP8_{}".format(pert_plate, replicate)
			logger.info((pert_plate, replicate, prism_replicate, dp7, dp8))

			assemble_dir = os.path.join(base_assemble_dir, prism_replicate)
			if not os.path.exists(assemble_dir):
				os.mkdir(assemble_dir)

			dp7_filepath = os.path.join(base_lxb_dir, dp7, dp7 + ".csv")
			assert os.path.exists(dp7_filepath)
			dp8_filepath = os.path.join(base_lxb_dir, dp8, dp8 + ".csv")
			assert os.path.exists(dp8_filepath)

			args_list[2] = prism_replicate
			args_list[6] = dp7_filepath
			args_list[8] = dp8_filepath
			logger.debug("args_list:  {}".format(args_list))

			args = parser.parse_args(args_list)
			logger.debug("args:  {}".format(args))

			assemble.main(args, all_perturbagens)

			output_files = glob.glob(prism_replicate + "_*.gct")
			for of in output_files:
				f = open(of)
				t = f.read().strip()
				len_t = len(t.split("\n"))
				print "output file of:  {}  length:  {}".format(of, len_t)
				
				dest_file = os.path.join(assemble_dir, of)
				if os.path.exists(dest_file):
					os.remove(dest_file)
				shutil.move(of, assemble_dir)


if __name__ == "__main__":
	setup_logger.setup(verbose=False)
	main()



