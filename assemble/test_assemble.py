import assemble
import unittest
import merino.setup_logger as setup_logger
import logging
import os
import glob
import mock


logger = logging.getLogger(setup_logger.LOGGER_NAME)


class TestAssemble(unittest.TestCase):

    def test_read_davepool_data_objects(self):
        l = [(1, "../requirements_artifacts/PDOG001_P1_X1.csv"), (2, "../requirements_artifacts/PDOG001_P2_X1.csv")]
        r = assemble.read_davepool_data_objects(l)
        assert len(r) > 0
        logger.debug("r:  {}".format(r))

        assert r[0].davepool_id == 1, r[0].davepool_id
        assert r[1].davepool_id == 2, r[1].davepool_id


    def test_build_davepool_id_csv_list(self):
        r = assemble.build_davepool_id_csv_list(["a", "1", "b", "2", "c", "3"])
        logger.debug("r:  {}".format(r))
        assert len(r) == 3, len(r)

        assert r[0][0] == "a", r[0][0]
        assert r[0][1] == "1", r[0][1]
        assert r[2][0] == "c", r[2][0]
        assert r[2][1] == "3", r[2][1]

    def test_full_functional_DP78(self):
        expected_files = ["PDOG003_DP78_120H_X1/PDOG003_DP78_120H_X1_COUNT.gct", "PDOG003_DP78_120H_X1/PDOG003_DP78_120H_X1_MEDIAN.gct"]
        for ef in expected_files:
            if os.path.exists(ef):
                os.remove(ef)


        config_filepath = "../prism_pipeline.cfg"

        plate_map_path = '../functional_tests/test_data/PDOG/map_src/PDOG003.src'
        dp7_csv_path = "../functional_tests/test_data/PDOG/lxb/PDOG003_DP7_24H_X1_B1.csv"
        dp8_csv_path = "../functional_tests/test_data/PDOG/lxb/PDOG003_DP8_24H_X1_B1.csv"
        cell_set_def_file = "../functional_tests/test_data/vdb/cell_set_definitions/PRISM_DP78.CS1_definition.txt"
        analyte_mapping_file = "../functional_tests/test_data/vdb/analyte_mapping/DP78_mapping.txt"
        assay_type = "DP78"

        args = assemble.build_parser().parse_args(["-config_filepath", config_filepath,
            "-pmp", plate_map_path, "-dp_csv", "DP7", dp7_csv_path, "DP8", dp8_csv_path, "-pert_time", "120",
                                                   "-csdf", cell_set_def_file, "-amf", analyte_mapping_file,
                                                   "-at", assay_type])

        logger.debug("args:  {}".format(args))

        assemble.main(args)

        for ef in expected_files:
            assert os.path.exists(ef), ef
            os.remove(ef)
        os.rmdir("./PDOG003_DP78_120H_X1")

        for map_file in glob.glob('PDOG*.src'):
            x = os.path.getsize(map_file)
            assert x > 0
            os.remove(map_file)

    def test_full_functional_PR500(self):
        expected_files = ["PASG003_PR500_120H_X251/PASG003_PR500_120H_X251_COUNT.gct",
                          "PASG003_PR500_120H_X251/PASG003_PR500_120H_X251_MEDIAN.gct"]
        for ef in expected_files:
            if os.path.exists(ef):
                os.remove(ef)

        config_filepath = "../prism_pipeline.cfg"

        plate_map_path = '../functional_tests/test_data/PASG/map_src/PASG003.src'
        csv_path = "../functional_tests/test_data/PASG/lxb/PASG003_PR500.2_X251/PASG003_PR500.2_X251.csv"
        cell_set_def_file = "../functional_tests/test_data/vdb/cell_set_definitions/PRISM_PR500.CS5_definition.txt"
        analyte_map_file = "../functional_tests/test_data/vdb/analyte_mapping/PR500_mapping.txt"
        assay_type = "PR500"

        args = assemble.build_parser().parse_args(["-config_filepath", config_filepath, "-pmp", plate_map_path,
                                                   "-csv", csv_path,  "-pert_time", "120", "-csdf", cell_set_def_file,
                                                   "-amf", analyte_map_file, "-at", assay_type])

        logger.debug("args:  {}".format(args))

        assemble.main(args)

        for ef in expected_files:
            assert os.path.exists(ef), ef
            os.remove(ef)
        os.rmdir("./PASG003_PR500_120H_X251")

        for map_file in glob.glob('PASG*.src'):
            x = os.path.getsize(map_file)
            assert x > 0
            os.remove(map_file)


    @mock.patch("assemble.assemble_core.validate_prism_gct.column_metadata_fields")
    def test_full_functional_PR300(self, mock_col_header_validation):
        # the plate_map in test_data is CM map, and thus fails header checks within validate_prism_gct
        # mocking the column metadata fields within output validation function given
        # the real test here is analyte and cell set mapping
        mock_col_header_validation.return_value = []

        expected_files = ["PSPA001_PR300_120H_X1/PSPA001_PR300_120H_X1_COUNT.gct",
                          "PSPA001_PR300_120H_X1/PSPA001_PR300_120H_X1_MEDIAN.gct"]

        for ef in expected_files:
            if os.path.exists(ef):
                os.remove(ef)

        config_filepath = "../prism_pipeline.cfg"

        plate_map_path = '../functional_tests/test_data/PSPA/map_src/PSPA001.src'
        csv_path = "../functional_tests/test_data/PSPA/lxb/PSPA001_PR300_X1.csv"
        cell_set_def_file = "../functional_tests/test_data/vdb/cell_set_definitions/PRISM_PR300.CS1_definition.txt"
        analyte_map_file = "../functional_tests/test_data/vdb/analyte_mapping/PR300_mapping.txt"
        assay_type = "PR300"

        args = assemble.build_parser().parse_args(["-config_filepath", config_filepath, "-pmp", plate_map_path,
                                                   "-csv", csv_path,  "-pert_time", "120", "-csdf", cell_set_def_file,
                                                   "-amf", analyte_map_file, "-at", assay_type])

        logger.debug("args:  {}".format(args))

        assemble.main(args)

        for ef in expected_files:
            assert os.path.exists(ef), ef
            os.remove(ef)
        os.rmdir("./PSPA001_PR300_120H_X1")

        for map_file in glob.glob('PSPA*.src'):
            x = os.path.getsize(map_file)
            assert x > 0
            os.remove(map_file)

    def test_full_functional_COPRO(self):
        expected_files = ["PGUM001_KJ100_120H_X1/PGUM001_KJ100_120H_X1_COUNT.gct",
                          "PGUM001_KJ100_120H_X1/PGUM001_KJ100_120H_X1_MEDIAN.gct"]
        for ef in expected_files:
            if os.path.exists(ef):
                os.remove(ef)

        config_filepath = "../prism_pipeline.cfg"

        plate_map_path = '../functional_tests/test_data/PGUM/map_src/PGUM001.src'
        csv_path = "../functional_tests/test_data/PGUM/lxb/PGUM001_KJ100_X1.jcsv"
        cell_set_def_file = "../functional_tests/test_data/vdb/cell_set_definitions/PRISM_KJ100.CS5_definition.txt"
        analyte_map_file = "../functional_tests/test_data/vdb/analyte_mapping/KJ100_mapping.txt"
        assay_type = "KJ100"

        args = assemble.build_parser().parse_args(["-config_filepath", config_filepath, "-pmp", plate_map_path,
                                                   "-csv", csv_path, "-pert_time", "120", "-csdf", cell_set_def_file,
                                                   "-amf", analyte_map_file, "-at", assay_type])

        logger.debug("args:  {}".format(args))

        assemble.main(args)

        for ef in expected_files:
            assert os.path.exists(ef), ef
            os.remove(ef)
        os.rmdir("./PGUM001_KJ100_120H_X1")

        for map_file in glob.glob('PGUM*.src'):
            x = os.path.getsize(map_file)
            assert x > 0
            os.remove(map_file)

if __name__ == "__main__":
    setup_logger.setup(verbose=True)

    unittest.main()
