import prism_metadata
import numpy
import pandas as pd
import logging
import setup_logger

logger = logging.getLogger(setup_logger.LOGGER_NAME)


def build_assayplate_pertplate_map(plate_tracking_file):

        plate_tacking_df = pd.read_table(plate_tracking_file, header=0, sep='\t')
        # Return sorted list of all assay plate barcodes
        assay_plates = sorted(set(plate_tacking_df['assay_plate_barcode']))
        assayplate_pertplate_map = {}

        # For each assay plate, find its entry and add its pert plate name to the map
        for (apb, pp) in plate_tacking_df[["assay_plate_barcode", "pert_plate"]].values:
                assayplate_pertplate_map[apb] = pp

        return assayplate_pertplate_map


def build_pertplate_perturbagen_map(all_perturbagens, assayplate_pertplate_map):

        pert_plate_perturbagens_map = {}

        # Read through each pert object in the cohort map, get its assay plate and use that to look up its pert plate.
        # Add the pert object to a map as a value under the key of its respective pert plate.
        for p in all_perturbagens:
                apb = p.assay_plate_barcode
                missing_plates = []
                if apb in assayplate_pertplate_map:
                        pertpl = assayplate_pertplate_map[apb]
                        if pertpl not in pert_plate_perturbagens_map:
                                pert_plate_perturbagens_map[pertpl] = []
                        current_perturbagens = pert_plate_perturbagens_map[pertpl]
                        current_perturbagens.append(p)
                else:
                        missing_plates.append(apb)

                if len(missing_plates) > 0:
                        msg = "{} assay plates in plate map were missing from plate tracking file:  {}".format(
                                len(missing_plates), missing_plates)
                        logger.error(msg)
                        raise Exception("check_and_build_perts build_pertplate_perturbagen_map " + msg)

        return pert_plate_perturbagens_map


def create_plate_map_dataframes(pert_plate_perturbagens_map):

        def index_builder(pert):
                return pert.well_id

        # Use validate perturbagens to check that all perts under a give pert plate are the same and
        # write out a well map of perts for each pert plate. Turn this well map into a list of perts to pass to df
        # builder. Sort and index df, then write to tab delimited .src file.
        dataframe_map = {}
        for (pertplate) in pert_plate_perturbagens_map:
                validated_perturbagens = prism_metadata.validate_perturbagens(pert_plate_perturbagens_map[pertplate])
                dataframe = prism_metadata.convert_objects_to_metadata_df(index_builder, validated_perturbagens.values(),
                                                           {"well_id": "pert_well"})
                dataframe_indexed = dataframe.set_index("pert_well")
                dataframe_indexed.drop('assay_plate_barcode', axis=1, inplace=True)
                dataframe_indexed.sort_index(axis='index', inplace=True)
                dataframe_map[pertplate] = dataframe_indexed

        return dataframe_map


def write_plate_maps(dataframe_map):

        for (pertplate, dataframe) in dataframe_map.items():
                dataframe.to_csv('{}.src'.format(pertplate), sep='\t')


def read_all_perturbagens_from_file(plate_map_path, config_filepath, plate_map_type):
    logger.info("loading all perturbagens...")
    all_perturbagens = prism_metadata.build_perturbagens_from_file(plate_map_path, plate_map_type, config_filepath)

    logger.info("finished loading all perturbagens")

    return all_perturbagens


def main(plate_map_path, plate_tracking_file, config_filepath):

        # Create mapping of all assay plates to their respective pert plate. Usually 69 assay plates per pert plate
        # This uses the plate tracking file to build the mapping
        assayplate_pertplate_map = build_assayplate_pertplate_map(plate_tracking_file)

        # Read in all perturbagens from the entire cohort from the platemap
        all_perturbagens = read_all_perturbagens_from_file(plate_map_path, config_filepath, prism_metadata.plate_map_type_CM)

        # Separate out all pert pbjects into their respective pert plates using assayplate to pertplate map.
        pert_plate_perturbagens_map = build_pertplate_perturbagen_map(all_perturbagens, assayplate_pertplate_map)

        # Return a pandas dataframe of perts for each pert plate
        pertplate_dataframes = create_plate_map_dataframes(pert_plate_perturbagens_map)

        # Write out individual plate map for each pert plate
        write_plate_maps(pertplate_dataframes)
