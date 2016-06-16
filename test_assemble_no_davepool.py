import unittest
import logging
import setup_logger
import assemble_no_davepool
import prism_metadata
import os.path


logger = logging.getLogger(setup_logger.LOGGER_NAME)


class TestAssembleNoDavepool(unittest.TestCase):
    def test_main(self):
        expected_files = ["my_prism_replicate_COUNT.gct", "my_prism_replicate_MEDIAN.gct"]
        for ef in expected_files:
            if os.path.exists(ef):
                os.remove(ef)

        config_filepath = "functional_tests/test_assemble_no_davepool/prism_pipeline.cfg"
        prism_replicate_name = "my_prism_replicate"
        plate_map_path = "functional_tests/test_assemble_no_davepool/PMEL.A001.src"
        csv_filepath = "functional_tests/test_assemble_no_davepool/PMEL.A001_CS2_X1.csv"

        args = assemble_no_davepool.build_parser().parse_args(["-config_filepath", config_filepath,
            prism_replicate_name, plate_map_path, csv_filepath, "-plate_map_type", prism_metadata.plate_map_type_CMap])

        (median_gctoo, count_gctoo) = assemble_no_davepool.main(args)

        well_expected_values_map = {"P24": [(135, "9946"), (484, "8910"), (515, "11153"), (59, "10992.5")],
                                    "A01": [(445, "12780"), (276, "3944"), (504, "10044.5"), (139, "4586")]}
        for (well_id, expected_values) in well_expected_values_map.items():
            column_name = prism_replicate_name + ":" + well_id

            for (cell_id, value) in expected_values:
                logger.debug("checking values in median_gctoo.data_df - well_id:  {}  cell_id:  {}  value:  {}".format(well_id, cell_id, value))
                assert median_gctoo.data_df[column_name].loc[cell_id] == value, \
                    (median_gctoo.data_df[column_name].loc[cell_id], value)

        for ef in expected_files:
            assert os.path.exists(ef), ef
            # os.remove(ef)


if __name__ == "__main__":
    setup_logger.setup(verbose=True)

    unittest.main()
