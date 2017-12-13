import unittest
import logging
import merino.setup_logger as setup_logger
import os.path
import glob
import plate_tracking_metadata as ptm


logger = logging.getLogger(setup_logger.LOGGER_NAME)

test_dir = "functional_tests/test_plate_tracking_metadata/"


class TestPlateTrackingMetadata(unittest.TestCase):
    def test_read_assay_plates(self):
        #happy path
        cur_test_dir = os.path.join(test_dir, "test_read_assay_plates")

        config_filepath = os.path.join(cur_test_dir, "register_prism_plates.cfg")

        search_str = os.path.join(cur_test_dir, "7159-03-A04-01-*_REPORT_dlahr.txt")
        logger.debug("search_str:  {}".format(search_str))
        tsv_filepaths = glob.glob(search_str)
        logger.debug("tsv_filepaths:  {}".format(tsv_filepaths))

        r = ptm.read_assay_plates(tsv_filepaths, config_filepath, "compound_map", "poscon_map")
        self.assertIsNotNone(r)
        r_ap = r[0]
        self.assertGreater(len(r_ap), 0)
        logger.debug("happy path result - r_ap[0]:  {}".format(r_ap[0]))

        #mismatch pool_id path
        tsv_filepaths = [os.path.join(cur_test_dir, "mismatch_pool_id.txt")]
        with self.assertRaises(Exception) as context:
            ptm.read_assay_plates(tsv_filepaths, config_filepath, "compound_map", "poscon_map")
        self.assertIsNotNone(context.exception)
        logger.debug("context.exception:  {}".format(context.exception))
        self.assertIn("problems encountered, see error messages", str(context.exception))

        #mismatch compound plate map
        tsv_filepaths = [os.path.join(cur_test_dir, "mismatch_compound_plate.txt")]
        with self.assertRaises(Exception) as context:
            ptm.read_assay_plates(tsv_filepaths, config_filepath, "compound_map", "poscon_map")
        self.assertIsNotNone(context.exception)
        logger.debug("context.exception:  {}".format(context.exception))
        self.assertIn("problems encountered, see error messages", str(context.exception))

        #mismatch poscon plate map
        tsv_filepaths = [os.path.join(cur_test_dir, "mismatch_poscon_plate.txt")]
        with self.assertRaises(Exception) as context:
            ptm.read_assay_plates(tsv_filepaths, config_filepath, "compound_map", "poscon_map")
        self.assertIsNotNone(context.exception)
        logger.debug("context.exception:  {}".format(context.exception))
        self.assertIn("problems encountered, see error messages", str(context.exception))

    def test_read_assay_plates_combination_treatments(self):
        #happy path
        cur_test_dir = os.path.join(test_dir, "test_read_assay_plates_combination_treatments")

        config_filepath = os.path.join(cur_test_dir, "register_prism_plates.cfg")

        tsv_filepaths = [os.path.join(cur_test_dir, "7159-03-A04-01-23_REPORT_dlahr.txt")]

        r = ptm.read_assay_plates(tsv_filepaths, config_filepath, "compound_plate", "poscon_map")
        self.assertIsNotNone(r)
        r_ap = r[0]
        self.assertGreater(len(r_ap), 0)
        r_ap = r_ap[0]
        logger.debug("r_ap:  {}".format(r_ap))
        self.assertEqual(2, len(r_ap.compound_plates))
        self.assertIn("AB00020382", r_ap.compound_plates)
        self.assertIn("AB00021543", r_ap.compound_plates)

    def test_read_assay_plates_mismatch_compound_plate_cols(self):
        #happy path
        cur_test_dir = os.path.join(test_dir, "test_read_assay_plates_mismatch_compound_plate_cols")

        config_filepath = os.path.join(cur_test_dir, "register_prism_plates.cfg")

        tsv_filepaths = [os.path.join(cur_test_dir, "input_file1.txt"), os.path.join(cur_test_dir, "input_file2.txt")]

        with self.assertRaises(Exception) as context:
            ptm.read_assay_plates(tsv_filepaths, config_filepath, "compound_plate", "poscon_map")
        self.assertIsNotNone(context.exception)
        logger.debug("context.exception:  {}".format(context.exception))
        self.assertIn("the compound plate columns are not the same in the provided files", str(context.exception))

    def test_get_all_related_plate_cols(self):
        headers = ["a", "a_1", "d", "a_13", "a", "a_5", "b", "c"]
        r = ptm.get_all_related_plate_cols("a", headers)
        self.assertIsNotNone(r)
        logger.debug("r:  {}".format(r))
        self.assertEqual(["a", "a_1", "a_13", "a", "a_5"], r)

    def test_build_plates_from_row(self):
        row = [2, 5, 7, 11]
        header_map = {"a":1, "a_2":0, "b":3}
        plate_cols = {"a", "a_2"}
        r = ptm.build_plates_from_row(row, header_map, plate_cols)
        self.assertIsNotNone(r)
        logger.debug("r:  {}".format(r))
        self.assertEqual((5, 2), r)
        # self.assertIn(2, r)
        # self.assertIn(5, r)


if __name__ == "__main__":
    setup_logger.setup(verbose=True)

    unittest.main()
