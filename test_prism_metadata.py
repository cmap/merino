from unittest import TestCase
import setup_logger
import logging
import unittest
import prism_metadata as pm
import ConfigParser


logger = logging.getLogger(setup_logger.LOGGER_NAME)

test_file = "functional_tests/test_prism_metadata/prism_cell_tsv.txt"


class TestPrismMetadata(unittest.TestCase):
    def test___init__(self):
        r = pm.PrismCell()
        assert hasattr(r, "pool_id")
        assert hasattr(r, "analyte_id")

    def test__read_data(self):
        (h, d) = pm._read_data(test_file)

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

        r = pm._generate_header_map(headers, cp, pm._prism_cell_config_file_section)
        logger.debug("r:  {}".format(r))
        assert len(r) == 3, len(r)

    def test__parse_data(self):
        headers = ["pool_id", "analyte", "strippedname"]
        cp = ConfigParser.RawConfigParser()
        cp.read("prism_pipeline.cfg")
        header_map = pm._generate_header_map(headers, cp, pm._prism_cell_config_file_section)

        data = [["1", "analyte 2", "my cell's name"], ["3", "analyte 5", "autre cell nom"]]

        r = pm._parse_data(header_map, data, pm.PrismCell)
        logger.debug("r:  {}".format(r))
        assert len(r) == len(data), len(r)

        header_map["extra header that doesn't have data in any row"] = 10
        r = pm._parse_data(header_map, data, pm.PrismCell)

        data.append(["7", "", "blah"])
        r = pm._parse_data(header_map, data, pm.PrismCell)
        assert r[2].analyte_id is None

    def test_read_prism_cell_from_file(self):
        r = pm.read_prism_cell_from_file("prism_pipeline.cfg")
        assert len(r) > 0
        logger.debug("r:  {}".format(r))

    def test_read_perturbagen_from_file(self):
        r = pm.read_perturbagen_from_file("functional_tests/test_prism_metadata/perturbagen.txt", "prism_pipeline.cfg")
        assert len(r) > 0
        for x in r:
            logger.debug("x:  {}".format(x))

    def test__build_additional_perturbagen_info(self):
        cp = ConfigParser.RawConfigParser()
        cp.read("prism_pipeline.cfg")

        p = pm.Perturbagen()
        #worst case scenario, these are both ints
        p.compound_well_mmoles_per_liter = int(5)
        p.dilution_factor = int(2000)
        p.pert_mfc_id = "BRD-K12345678-910-11-1"
        p.pert_type = "Test"
        p.pert_mfc_desc = "my fake compound name"

        pm._build_additional_perturbagen_info(cp, [p])
        logger.debug("p:  {}".format(p))
        assert p.pert_dose == (1000.0 * float(5) / float(2000)), p.pert_dose
        assert p.pert_dose_unit == "uM", p.pert_dose_unit
        assert p.pert_idose is not None

        assert p.pert_id == "BRD-K12345678", p.pert_id
        assert p.pert_type == "trt_cp", p.pert_type
        assert p.pert_iname == "my fake compound name", p.pert_iname
        assert p.pert_time == "120", p.pert_time
        assert p.pert_time_unit == "h", p.pert_time_unit
        assert p.pert_itime == "120 h", p.pert_itime

        p.pert_type = "unrecognized"
        pm._build_additional_perturbagen_info(cp, [p])
        assert p.pert_type == "ctl_vehicle", p.pert_type

        p.pert_type = None
        pm._build_additional_perturbagen_info(cp, [p])
        assert p.pert_type == "ctl_vehicle", p.pert_type

        p.pert_mfc_id = None
        pm._build_additional_perturbagen_info(cp, [p])
        assert p.pert_id == "DMSO", p.pert_id

        p.compound_well_mmoles_per_liter = "not a valid concentration"
        with self.assertRaises(Exception) as context:
            pm._build_additional_perturbagen_info(cp, [p])
        assert context.exception is not None
        logger.debug("context.exception:  {}".format(context.exception))
        assert "the concentration or dilution factors should be numbers" in str(context.exception), str(context.exception)

        p.compound_well_mmoles_per_liter = 5.0
        del p.pert_type
        with self.assertRaises(Exception) as context:
            pm._build_additional_perturbagen_info(cp, [p])
        assert context.exception is not None
        logger.debug("context.exception:  {}".format(context.exception))
        assert "pert_type attribute is missing from perturbagen" in str(context.exception), str(context.exception)


if __name__ == "__main__":
    setup_logger.setup(verbose=True)

    unittest.main()
