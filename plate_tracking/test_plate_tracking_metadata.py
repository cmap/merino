import unittest
import logging
import prism_pipeline.setup_logger as setup_logger
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

        r = ptm.read_assay_plates(tsv_filepaths, config_filepath)
        self.assertIsNotNone(r)
        self.assertGreater(len(r), 0)
        logger.debug("r[0]:  {}".format(r[0]))

        #mismatch pool_id path
        tsv_filepaths = [os.path.join(cur_test_dir, "mismatch_pool_id.txt")]
        with self.assertRaises(Exception) as context:
            ptm.read_assay_plates(tsv_filepaths, config_filepath)
        self.assertIsNotNone(context.exception)
        logger.debug("context.exception:  {}".format(context.exception))
        self.assertIn("problems encountered, see error messages", str(context.exception))

        #mismatch poscon plate map
        tsv_filepaths = [os.path.join(cur_test_dir, "mismatch_poscon_plate.txt")]
        with self.assertRaises(Exception) as context:
            ptm.read_assay_plates(tsv_filepaths, config_filepath)
        self.assertIsNotNone(context.exception)
        logger.debug("context.exception:  {}".format(context.exception))
        self.assertIn("problems encountered, see error messages", str(context.exception))


if __name__ == "__main__":
    setup_logger.setup(verbose=True)

    unittest.main()
