import argparse
import prism_pipeline
import prism_metadata
import assemble
import setup_logger
import logging
import davepool_data
import sys
import GCToo.write_gctoo as write_gctoo


logger = logging.getLogger(setup_logger.LOGGER_NAME)


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-verbose", '-v', help="Whether to print a bunch of output", action="store_true", default=False)
    parser.add_argument("-config_filepath", help="path to the location of the configuration file", type=str,
                        default=prism_pipeline.default_config_filepath)
    parser.add_argument("prism_replicate_name", help="name of the prism replicate that is being processed", type=str)
    parser.add_argument("plate_map_path", help="path to file containing plate map describing perturbagens used", type=str)
    parser.add_argument("csv_filepath", help="path to csv file containing data", type=str)
    parser.add_argument("-plate_map_type", "-pmt", help="type of the plate map", choices=prism_metadata.plate_map_types,
                        default=prism_metadata.plate_map_type_CM)
    parser.add_argument("-cell_set_definition_file", "-csdf",
                        help="file containing cell set definition to use, overriding config file", type=str, default=None)
    return parser


def main(args):
    #read actual data from relevant csv files, associate it with davepool ID
    my_davepool = davepool_data.read_data(args.csv_filepath)

    #read PRISM cell line metadata from file specified in config file, and associate with assay_plate metadata
    prism_cell_list = prism_metadata.read_prism_cell_from_file(args.config_filepath, args.cell_set_definition_file)
    logger.info("len(prism_cell_list):  {}".format(len(prism_cell_list)))

    #read in all the perturbagens but restrict to those that were on the provided assay_plates
    perturbagen_list = prism_metadata.build_perturbagens_from_file(args.plate_map_path, args.plate_map_type,
                                                                   args.config_filepath)
    logger.info("len(perturbagen_list):  {}".format(len(perturbagen_list)))

    (median_data_by_cell, count_data_by_cell) = assemble.build_data_by_cell(prism_cell_list, my_davepool)

    median_gctoo = assemble.build_gctoo(args.prism_replicate_name, perturbagen_list, median_data_by_cell)
    write_gctoo.write(median_gctoo, args.prism_replicate_name + "_MEDIAN.gct", data_null=assemble._NaN,
                      filler_null=assemble._null)

    count_gctoo = assemble.build_gctoo(args.prism_replicate_name, perturbagen_list, count_data_by_cell)
    write_gctoo.write(count_gctoo, args.prism_replicate_name + "_COUNT.gct", data_null=assemble._NaN,
                      filler_null=assemble._null)

    return (median_gctoo, count_gctoo)


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    logger.debug("args:  {}".format(args))

    main(args)
