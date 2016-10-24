import unittest
import logging
import prism_pipeline.setup_logger as setup_logger
import register_prism_plates as rpp
import os.path
import glob
import tempfile


logger = logging.getLogger(setup_logger.LOGGER_NAME)

test_dir = "functional_tests/test_register_prism_plates"

class TestRegisterPrismPlates(unittest.TestCase):
    def test_main(self):
        cur_test_dir = os.path.join(test_dir, "test_main")
        config_filepath = os.path.join(cur_test_dir, "register_prism_plates.cfg")

        output_file = tempfile.NamedTemporaryFile(mode="w")
        logger.debug("output_file.name:  {}".format(output_file.name))

        arg_list = ["-config_filepath", config_filepath, "-output_file", output_file.name, "PCAL", "45", "CS3.1"]

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


if __name__ == "__main__":
    setup_logger.setup(verbose=True)

    unittest.main()
