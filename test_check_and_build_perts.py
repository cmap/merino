import check_and_build_perts
import unittest
import setup_logger
import logging
import prism_metadata
import numpy
import os
import glob

logger = logging.getLogger(setup_logger.LOGGER_NAME)

plate_tracking = "functional_tests/test_assemble/test_build_perturbagen_list/plate_tracking_example.tsv"
config = 'functional_tests/test_assemble/test_build_perturbagen_list/prism_pipeline.cfg'

wells = {'A01': 'SCW0117638', 'B01': 'SCW0116413', 'C01': 'SCW0115438', 'D01': 'SCW0114128', 'E02': 'SCW0118663'}
all_perturbagens = []

for (well, assay_plate_barcode) in wells.items():
    pert_instance = prism_metadata.Perturbagen(well_id=well)
    pert_instance.assay_plate_barcode = assay_plate_barcode
    all_perturbagens.append(pert_instance)

class TestCheckAndBuildPerts(unittest.TestCase):
    def test_build_assayplate_pertplate_map(self):

        test_map = check_and_build_perts.build_assayplate_pertplate_map(plate_tracking)
        assayplate_expected_pert_map = {'SCW0117638': 'PCAL015', 'SCW0116413': 'PCAL015', 'SCW0115438': 'PCAL015',
                                        'SCW0114128': 'PCAL016', 'SCW0118663': 'PCAL016'}

        for (assayplate, expected_pertplate) in assayplate_expected_pert_map.items():
            assert test_map[assayplate] == expected_pertplate
        assert len(test_map) == 9
        assert type(test_map) == dict

    def test_build_pertplate_perturbagen_map(self):

        assayplate_pertplate_map = check_and_build_perts.build_assayplate_pertplate_map(plate_tracking)
        test_map = check_and_build_perts.build_pertplate_perturbagen_map(all_perturbagens, assayplate_pertplate_map)

        assert len(test_map) == 2
        assert type(test_map) == dict
        assert test_map['PCAL015'][1].well_id == 'C01'
        assert test_map['PCAL016'][1].well_id == 'D01'

    def test_write_plate_maps(self):
        for map_file in glob.glob('PCAL*.src'):
            os.remove(map_file)

        assayplate_pertplate_map = check_and_build_perts.build_assayplate_pertplate_map(plate_tracking)

        pert_plate_perturbagens_map = check_and_build_perts.build_pertplate_perturbagen_map(all_perturbagens, assayplate_pertplate_map)

        dataframes = check_and_build_perts.create_plate_map_dataframes(pert_plate_perturbagens_map)

        check_and_build_perts.write_plate_maps(dataframes)

        num_files = glob.glob("PCAL*.src")
        assert len(num_files) == 2

        for map_file in glob.glob('PCAL*.src'):
            x = os.path.getsize(map_file)
            assert x > 0
            os.remove(map_file)


if __name__ == "__main__":
    setup_logger.setup(verbose=True)


    unittest.main()
