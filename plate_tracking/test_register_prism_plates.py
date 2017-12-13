import unittest
import logging
import merino.setup_logger as setup_logger
import register_prism_plates as rpp
import os.path
import glob
import tempfile
import merino.prism_metadata as prism_metadata


logger = logging.getLogger(setup_logger.LOGGER_NAME)

test_dir = "functional_tests/test_register_prism_plates"

class TestRegisterPrismPlates(unittest.TestCase):
    def test_main(self):
        cur_test_dir = os.path.join(test_dir, "test_main")
        config_filepath = os.path.join(cur_test_dir, "register_prism_plates.cfg")

        output_file = tempfile.NamedTemporaryFile(mode="w")
        logger.debug("output_file.name:  {}".format(output_file.name))

        arg_list = ["-config_filepath", config_filepath, "-output_file", output_file.name, "-dont", "PDOG", "45", "CS3.1"]

        search_str = os.path.join(cur_test_dir, "7159-03-A04-01-*_REPORT_dlahr.txt")
        input_files = glob.glob(search_str)
        arg_list.extend(input_files)

        args = rpp.build_parser().parse_args(arg_list)
        logger.debug("args:  {}".format(args))


        rpp.main(args)

        f = open(output_file.name)
        r = f.read().strip().split("\n")
        f.close()
        N = 5
        self.assertGreater(len(r), N)

        for i in xrange(N):
            logger.debug("r[{}]:  {}".format(i, r[i]))

    def test_main_combination_treatments(self):
        cur_test_dir = os.path.join(test_dir, "test_main_combination_treatments")
        config_filepath = os.path.join(cur_test_dir, "register_prism_plates.cfg")

        output_file = tempfile.NamedTemporaryFile(mode="w")
        logger.debug("output_file.name:  {}".format(output_file.name))

        arg_list = ["-config_filepath", config_filepath, "-output_file", output_file.name,
                    "-compound_plate_col_basename", "compound_plate", "PMEL", "35", "CS2.1", "-dont"]

        search_str = os.path.join(cur_test_dir, "7159-03-A04-01-*_REPORT_dlahr.txt")
        input_files = glob.glob(search_str)
        arg_list.extend(input_files)

        args = rpp.build_parser().parse_args(arg_list)
        logger.debug("args:  {}".format(args))

        rpp.main(args)

        f = open(output_file.name)
        r = f.read().strip().split("\n")
        f.close()
        N = 5
        self.assertGreater(len(r), N)

        for i in xrange(N):
            logger.debug("r[{}]:  {}".format(i, r[i]))

    def test_determine_pool_ids(self):
        assay_plates = [prism_metadata.AssayPlate(pool_id=i/2) for i in xrange(10)]
        r = rpp.determine_pool_ids(assay_plates)
        self.assertIsNotNone(r)
        logger.debug("r:  {}".format(r))
        self.assertEqual(r, range(5))

    def test_determine_compound_plates(self):
        #happy path just 1 compound plate per assay plate
        expected_compound_plates = []
        assay_plates = []
        for i in xrange(10):
            ap = prism_metadata.AssayPlate()
            compound_plate = (i/2,)
            expected_compound_plates.append(compound_plate)
            ap.compound_plates = compound_plate
            assay_plates.append(ap)

        expected_compound_plates = list(set(expected_compound_plates))
        expected_compound_plates.sort()

        r = rpp.determine_compound_plates(assay_plates)
        self.assertIsNotNone(r)
        logger.debug("r:  {}".format(r))
        self.assertEqual(r, expected_compound_plates)

    def test_build_pert_plate_mapping(self):
        r = rpp.build_pert_plate_mapping([(0,), (1,)], 1, "FAKEPROJ")
        self.assertIsNotNone(r)
        logger.debug("r:  {}".format(r))
        self.assertEqual(2, len(r))
        self.assertIn((0,), r)
        self.assertEqual(r[(0,)], "FAKEPROJ001")
        self.assertIn((1,), r)
        self.assertEqual(r[(1,)], "FAKEPROJ002")

if __name__ == "__main__":
    setup_logger.setup(verbose=True)

    unittest.main()
