import unittest
import logging
import setup_logger as setup_logger
import adjust_plate_map_pert_type as apmpt
import os
import pandas


logger = logging.getLogger(setup_logger.LOGGER_NAME)


# Some notes on testing conventions (more in cuppers convention doc):
#	(1) Use "self.assert..." over "assert"
#		- self.assert* methods: https://docs.python.org/2.7/library/unittest.html#assert-methods
#       - This will ensure that if one assertion fails inside a test method,
#         exectution won't halt and the rest of the test method will be executed
#         and other assertions are also verified in the same run.  
# 	(2) For testing exceptions use:
#		with self.assertRaises(some_exception) as context:
#			[call method that should raise some_exception]
#		self.assertEqual(str(context.exception), "expected exception message")
#
#		self.assertAlmostEquals(...) for comparing floats


class TestAdjustPlateMapPertType(unittest.TestCase):
    def test_main(self):
        expected_file = "test_adjust_plate_map_pert_type.src"
        if os.path.exists(expected_file):
            os.remove(expected_file)

        args_str_list = [os.path.join("functional_tests", expected_file)]
        args = apmpt.build_parser().parse_args(args_str_list)

        apmpt.main(args)

        self.assertTrue(os.path.exists(expected_file))

        df = pandas.read_csv(expected_file, sep="\t")
        poscon_count = sum(df[apmpt.pert_type] == apmpt.pert_type_poscon)
        logger.debug("output files poscon_count:  {}".format(poscon_count))
        self.assertGreater(poscon_count, 0)

        os.remove(expected_file)

    def test_replace_pert_type_for_poscons(self):
        df = pandas.DataFrame({apmpt.pert_type:range(3), "pert_iname":["my pos 1", "my pos 2", "not a pos"]})
        logger.debug("before replacement df:  {}".format(df))

        apmpt.replace_pert_type_for_poscons({"fake filename":df}, {"my pos 1", "my pos 2"})
        logger.debug("after replacement df:  {}".format(df))

        new_pert_type = df.loc[df.pert_iname == "my pos 1", apmpt.pert_type].values[0]
        self.assertEqual(apmpt.pert_type_poscon, new_pert_type)


if __name__ == "__main__":
    setup_logger.setup(verbose=True)

    unittest.main()
