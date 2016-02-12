import assemble
import unittest
import setup_logger
import logging


logger = logging.getLogger(setup_logger.LOGGER_NAME)


class TestAssemble(unittest.TestCase):
    def test_parse_location_to_well(self):
        r = assemble.parse_location_to_well("1(1,A1)")
        assert r == "A1", r

        r = assemble.parse_location_to_well("384(1,P24)")
        assert r == "P24", r


if __name__ == "__main__":
    setup_logger.setup(verbose=True)

    unittest.main()
