from unittest import TestCase
import setup_logger
import logging
import unittest
import prism_cell as pc
import ConfigParser


logger = logging.getLogger(setup_logger.LOGGER_NAME)

test_file = "functional_tests/test_prism_cell/prism_cell_tsv.txt"

class TestPrismCell(unittest.TestCase):
    def test___init__(self):
        r = pc.PrismCell()
        assert hasattr(r, "pool_id")
        assert hasattr(r, "analyte_id")

    def test__read_data(self):
        (h, d) = pc._read_data(test_file)

        assert h is not None
        logger.debug("h:  {}".format(h))
        assert "pool_id" in h

        assert d is not None
        logger.debug("d:  {}".format(d))
        assert len(d) > 0

    def test__generate_header_map(self):
        headers = ["pool_id", "analyte", "strippedname"]

        cp = ConfigParser.RawConfigParser()
        cp.read("prism_pipeline.cfg")

        r = pc._generate_header_map(headers, cp)
        logger.debug("r:  {}".format(r))
        assert len(r) == 3, len(r)

    def test__parse_data(self):
        headers = ["pool_id", "analyte", "strippedname"]
        cp = ConfigParser.RawConfigParser()
        cp.read("prism_pipeline.cfg")
        header_map = pc._generate_header_map(headers, cp)

        data = [["1", "analyte 2", "my cell's name"], ["3", "analyte 5", "autre cell nom"]]

        r = pc._parse_data(header_map, data)
        logger.debug("r:  {}".format(r))
        assert len(r) == len(data), len(r)

    def test_read_prism_cell_from_file(self):
        r = pc.read_prism_cell_from_file(test_file, "prism_pipeline.cfg")
        assert len(r) > 0
        logger.debug("r:  {}".format(r))

if __name__ == "__main__":
    setup_logger.setup(verbose=True)

    unittest.main()
