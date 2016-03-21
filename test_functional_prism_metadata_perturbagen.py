import prism_metadata
import setup_logger
import logging


logger = logging.getLogger(setup_logger.LOGGER_NAME)


if __name__ == "__main__":
    setup_logger.setup(verbose=True)

    assay_plates = prism_metadata.read_assay_plate_from_file("requirements_artifacts/REP Master Plate Key_All Pools.txt",
                                                             "prism_pipeline.cfg")
    replicates_det_plates = {"PCAL001_P7_X1", "PCAL001_P8_X1"}
    assay_plates = [x for x in assay_plates if x.det_plate in replicates_det_plates]
    logger.info("assay_plates:  {}".format(assay_plates))

    assay_plate_barcodes = set([x.assay_plate_barcode for x in assay_plates])

    perts = prism_metadata.read_perturbagen_from_file("requirements_artifacts/7159-03-A02-01-01_02-29-16_12.20.32.txt",
                                                      "prism_pipeline.cfg")

    replicate_perts = [x for x in perts if x.assay_plate_barcode in assay_plate_barcodes]

    well_pert_map = prism_metadata.validate_perturbagens(replicate_perts)

    logger.info("len(well_part_map):  {}".format(len(well_pert_map)))
    assert len(well_pert_map) == 384, len(well_pert_map)
