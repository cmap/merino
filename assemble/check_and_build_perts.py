"""

Command line script for validating and separating individual plate maps from a cohort plate map.

Given a cohort plate map, this script checks that all metadata for a given pert plate is the same,
and then writes the metadata out to a separate plate map for each pert plate.

Given a plate tracking file, the script builds a map for which assay plates correspond to which pert plates.
It then uses this to map out which perturbagen objects in the plate map correspond to which pert plate,
based on their assay plate. Once there is a 'perturbagen map' for each pert plate, these perturbagen sets can
be validated (this just checks that all of the corresponding wells received the same treatment),
and written into a new plate map.

"""

import prism_metadata
import numpy
import pandas as pd
import logging
import setup_logger
import argparse
import sys
import merino

__author__ = "Evan Lemire"
__email__ = "elemire@broadinstitute.org"

logger = logging.getLogger(setup_logger.LOGGER_NAME)

def build_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Optional
    parser.add_argument("-verbose", '-v', help="Whether to print a bunch of output", action="store_true", default=False)
    # Next three required
    parser.add_argument("-config_filepath", "-cfg", help="path to the location of the configuration file", type=str,
                        default=merino.default_config_filepath)
    parser.add_argument("-plate_map_path", "-pmp", help="path to file containing plate map describing perturbagens used", type=str, required=True)
    parser.add_argument("-plate_tracking_path", "-ptp", help="path to file containing the mapping between assasy plates and det_plates",
                        type=str, required=True)
    return parser

def build_assayplate_pertplate_map(plate_tracking_file):

    # Read in plate tracking file as pandas data frame
    plate_tracking_df = pd.read_table(plate_tracking_file, header=0)

    # Initialize map
    assayplate_pertplate_map = {}

    # For each assay plate, create a key in the map and add its pert plate as the value
    for (apb, pp) in plate_tracking_df[["assay_plate_barcode", "pert_plate"]].values:
            assayplate_pertplate_map[apb] = pp

    return assayplate_pertplate_map


def build_pertplate_perturbagen_map(all_perturbagens, assayplate_pertplate_map):

    pertplate_perturbagens_map = {}

    # Read through each pert object in the cohort map, get its assay plate and use that to look up its pert plate.
    # Add the pert object to a map as a value under the key of its respective pert plate.

    missing_plates = []
    for p in all_perturbagens:
            apb = p.assay_plate_barcode
            if apb in assayplate_pertplate_map:
                    pertpl = assayplate_pertplate_map[apb]
                    if pertpl not in pertplate_perturbagens_map:
                            pertplate_perturbagens_map[pertpl] = []
                    current_perturbagens = pertplate_perturbagens_map[pertpl]
                    current_perturbagens.append(p)
            # Record the names of plate present in the plate map but not in plate tracking
            else:
                    missing_plates.append(apb)

    # If there are plates missing from plate tracking, fail.
    if len(missing_plates) > 0:
            msg = "{} assay plates in plate map were missing from plate tracking file:  {}".format(
                    len(missing_plates), set(missing_plates))
            logger.error(msg)
            raise Exception("check_and_build_perts build_pertplate_perturbagen_map " + msg)

    return pertplate_perturbagens_map


def create_plate_map_dataframes(pertplate_perturbagens_map):

        # This index builder currently serves no purpose. convert_objects_to_metadata_df requires an 'index builder'
        # as it was initially used to make GCTs and you labelled your rows and columns using pieces from multiple
        # fields in the df. In the case of the plate map we just want to make one column the index, but refactoring
        # convert_objects_to_metadata_df would complicate the code a lot, so we make a dummy index builder.
        def index_builder(pert):
                return pert.well_id

        # Use validate perturbagens to check that all perts under a give pert plate are the same and
        # write out a well map of perts for each pert plate. Turn this well map into a list of perts to pass to df
        # builder. Sort and index df, then write to tab delimited .src file.
        dataframe_map = {}
        for (pertplate, perturbagens) in pertplate_perturbagens_map.iteritems():
                validated_perturbagens = prism_metadata.validate_perturbagens(perturbagens)
                dataframe = prism_metadata.convert_objects_to_metadata_df(index_builder,
                                                                          validated_perturbagens.values())
                dataframe_indexed = dataframe.set_index("pert_well")
                dataframe_indexed.drop('assay_plate_barcode', axis=1, inplace=True)
                dataframe_indexed.sort_index(axis='index', inplace=True)
                dataframe_indexed.fillna(value='-666', inplace=True)
                dataframe_map[pertplate] = dataframe_indexed

        return dataframe_map


def write_plate_maps(dataframe_map):
        # Write plate maps to tab delimited text files
        for (pertplate, dataframe) in dataframe_map.items():
                dataframe.to_csv('{}.src'.format(pertplate), sep='\t')


def read_all_perturbagens_from_file(plate_map_path, config_filepath, plate_map_type):
    logger.info("loading all perturbagens...")
    all_perturbagens = prism_metadata.build_perturbagens_from_file(plate_map_path, plate_map_type)

    logger.info("finished loading all perturbagens")

    return all_perturbagens


def main(args):

        # Create mapping of all assay plates to their respective pert plate. There are usually multiple
        # assay plates per pert plate, e.g. for a typical experiment with CellSet 1.2 there are 23 cell pools
        # done in 3 replicates = 69 assay plates per pert plate
        assayplate_pertplate_map = build_assayplate_pertplate_map(args.plate_tracking_path)

        # Read in all perturbagens from the entire cohort from the platemap
        all_perturbagens = read_all_perturbagens_from_file(args.plate_map_path, args.config_filepath, prism_metadata.plate_map_type_CM)

        # Separate out all pert pbjects into their respective pert plates using assayplate to pertplate map
        pertplate_perturbagens_map = build_pertplate_perturbagen_map(all_perturbagens, assayplate_pertplate_map)

        # Return a pandas dataframe of perts for each pert plate
        pertplate_dataframes = create_plate_map_dataframes(pertplate_perturbagens_map)

        # Write out individual plate map for each pert plate
        write_plate_maps(pertplate_dataframes)

if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    logger.debug("args:  {}".format(args))

    main(args)
