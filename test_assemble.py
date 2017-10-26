import assemble
import unittest
import setup_logger
import logging
import prism_metadata
import davepool_data
import check_and_build_perts
import numpy
import os
import glob


logger = logging.getLogger(setup_logger.LOGGER_NAME)


class TestAssemble(unittest.TestCase):
    def test_parse_location_to_well(self):
        r = assemble.parse_location_to_well("1(1,A1)")
        assert r == "A01", r

        r = assemble.parse_location_to_well("384(1,P24)")
        assert r == "P24", r

    def test_read_davepool_data_objects(self):
        l = [(1, "requirements_artifacts/PCAL001_P1_X1.csv"), (2, "requirements_artifacts/PCAL001_P2_X1.csv")]
        r = assemble.read_davepool_data_objects(l)
        assert len(r) > 0
        logger.debug("r:  {}".format(r))

        assert r[0].davepool_id == 1, r[0].davepool_id
        assert r[1].davepool_id == 2, r[1].davepool_id


    def test_combine_maps_with_checks(self):
        a = {1:2, 3:5}
        r = {7:11, 13:17}

        assemble.combine_maps_with_checks(a, r)
        logger.debug("b:  {}".format(r))

        for (k,v) in a.items():
            assert k in r
            assert r[k] == v

        with self.assertRaises(Exception) as context:
            assemble.combine_maps_with_checks(a, r)
        assert context.exception is not None
        logger.debug("context.exception:  {}".format(context.exception))
        assert "source_map and dest_map had common_keys" in str(context.exception), str(context.exception)

    def test_build_davepool_id_csv_list(self):
        r = assemble.build_davepool_id_csv_list(["a", "1", "b", "2", "c", "3"])
        logger.debug("r:  {}".format(r))
        assert len(r) == 3, len(r)

        assert r[0][0] == "a", r[0][0]
        assert r[0][1] == "1", r[0][1]
        assert r[2][0] == "c", r[2][0]
        assert r[2][1] == "3", r[2][1]

    def test_full_functional(self):
        expected_files = ["PCAL003_CS1_X1_COUNT.gct", "PCAL003_CS1_X1_MEDIAN.gct"]
        for ef in expected_files:
            if os.path.exists(ef):
                os.remove(ef)

        config_filepath = "/Users/elemire/Workspace/prism_pipeline/functional_tests/test_assemble/full_functional_test/prism_pipeline.cfg"
        prism_replicate_name = "PCAL003_CS1_X1"
        plates_mapping_path = "functional_tests/test_assemble/full_functional_test/2016-03-22_PCAL_plate_mapping.txt"
        args1 = check_and_build_perts.build_parser().parse_args(
            ["-cfg", "functional_tests/test_assemble/full_functional_test/prism_pipeline.cfg", "-pmp",
             "functional_tests/test_assemble/full_functional_test/7159-03-A04-01-01_03-22-16_15.34.24.txt",
             "-ptp", "functional_tests/test_assemble/full_functional_test/2016-03-22_PCAL_plate_mapping.txt"])
        check_and_build_perts.main(args1)
        plate_map_path = 'PCAL003.src'
        dp7_csv_path = "functional_tests/test_assemble/full_functional_test/PCAL003_DP7_X1.csv"
        dp8_csv_path = "functional_tests/test_assemble/full_functional_test/PCAL003_DP8_X1.csv"
        csdf_path = "/Users/elemire/Workspace/prism_pipeline/requirements_artifacts/CalicoTranche1PrimaryMetaData_02252016.txt"
        dmf_path = "/Users/elemire/Workspace/prism_pipeline/requirements_artifacts/test_davepool_analyte_mapping.txt"

        args = assemble.build_parser().parse_args(["-config_filepath", config_filepath, "-prn", prism_replicate_name,
            "-pmp", plate_map_path, "-ptp", plates_mapping_path, "-dp_csv", "DP7", dp7_csv_path, "DP8", dp8_csv_path,
                                                   "-csdf", csdf_path, "-dmf", dmf_path])

        logger.debug("args:  {}".format(args))

        assemble.main(args)

        for ef in expected_files:
            assert os.path.exists(ef), ef
            os.remove(ef)

        for map_file in glob.glob('PCAL*.src'):
            x = os.path.getsize(map_file)
            assert x > 0
            os.remove(map_file)


if __name__ == "__main__":
    setup_logger.setup(verbose=True)

    unittest.main()
