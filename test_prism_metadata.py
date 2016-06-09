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

        headers = ["well_position", "compound_well_mmoles_per_liter", "dilution_factor"]
        cp = ConfigParser.RawConfigParser()
        cp.read("prism_pipeline.cfg")
        header_map = pm._generate_header_map(headers, cp, pm._perturbagen_config_file_section)

        data = [["A01", "1.010101", "2"], ["B07", "3.030303", "5"]]
        r = pm._parse_data(header_map, data, pm.Perturbagen)
        logger.debug("r:  {}".format(r))
        assert len(r) == len(data), len(r)
        assert isinstance(r[0].compound_well_mmoles_per_liter, float)
        assert r[0].compound_well_mmoles_per_liter == 1.010101, r[0].compound_well_mmoles_per_liter
        assert isinstance(r[0].dilution_factor, int)
        assert r[0].dilution_factor == 2, r[0].dilution_factor
        assert isinstance(r[1].compound_well_mmoles_per_liter, float)
        assert isinstance(r[1].dilution_factor, int)

    def test_read_prism_cell_from_file(self):
        r = pm.read_prism_cell_from_file("prism_pipeline.cfg")
        assert len(r) > 0
        logger.debug("r:  {}".format(r))

    def test_read_perturbagen_from_CM_file(self):
        r = pm.read_perturbagen_from_CM_file("functional_tests/test_prism_metadata/perturbagen.txt", "prism_pipeline.cfg")
        assert len(r) > 0
        for x in r:
            logger.debug("x:  {}".format(x))

    def test_read_assay_plate_from_file(self):
        r = pm.read_assay_plate_from_file("functional_tests/test_prism_metadata/assay_plate.txt", "prism_pipeline.cfg")
        assert len(r) > 0
        for x in r:
            logger.debug("x:  {}".format(x))
            assert x.assay_plate_barcode is not None
            assert x.det_plate is not None
            assert x.pool_id is not None

    def test__build_additional_perturbagen_info(self):
        cp = ConfigParser.RawConfigParser()
        cp.read("prism_pipeline.cfg")

        p = pm.Perturbagen()
        #worst case scenario, these are both ints
        p.compound_well_mmoles_per_liter = int(5)
        p.dilution_factor = int(2001)
        p.pert_mfc_id = "BRD-K12345678-910-11-1"
        p.pert_type = "Test"
        p.pert_mfc_desc = "my fake compound name"

        pm._build_additional_perturbagen_info(cp, [p])
        logger.debug("p:  {}".format(p))
        assert p.pert_dose == (1000.0 * float(5) / 2001.0), p.pert_dose
        assert p.pert_dose_unit == "uM", p.pert_dose_unit
        assert p.pert_idose is not None

        assert p.pert_id == "BRD-K12345678", p.pert_id
        assert p.pert_type == "trt_cp", p.pert_type
        assert p.pert_iname == "my fake compound name", p.pert_iname
        assert p.pert_time == "120", p.pert_time
        assert p.pert_time_unit == "h", p.pert_time_unit
        assert p.pert_itime == "120 h", p.pert_itime

        p.pert_mfc_desc = None
        pm._build_additional_perturbagen_info(cp, [p])
        assert p.pert_iname == p.pert_id, p.pert_iname

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

    def test__build_additional_perturbagen_info_vehicle(self):
        cp = ConfigParser.RawConfigParser()
        cp.read("prism_pipeline.cfg")

        p = pm.Perturbagen()
        p.pert_type = "Empty"
        p.pert_vehicle = "DMSO"
        p.compound_well_mmoles_per_liter = None
        p.pert_mfc_id = None

        pm._build_additional_perturbagen_info(cp, [p])
        logger.debug("p:  {}".format(p))

        assert p.pert_dose is None, p.pert_dose
        assert p.pert_dose_unit is None, p.pert_dose_unit
        assert p.pert_idose is None, p.pert_idose

    def test_validate_perturbagens(self):
        perts = []
        num_wells = 3

        for assay_plate_barcode in ["a" + str(i) for i in range(2)]:
            new_perts = [pm.Perturbagen(well_id=i) for i in range(num_wells)]
            for (i,np) in enumerate(new_perts):
                np.assay_plate_barcode = assay_plate_barcode
                np.pert_id = "BRD-K" + str(int(11.0*i))
                np.pert_idose = str(13.0*i) + " uM"
            perts.extend(new_perts)
        logger.debug("perts:  {}".format(perts))

        r = pm.validate_perturbagens(perts)
        logger.debug("r:  {}".format(r))

        assert len(r) == num_wells, (len(r), num_wells)
        for i in range(num_wells):
            assert i in r
            assert r[i].well_id == i, r[i]

        perts[num_wells].pert_id = perts[0].pert_id + " extra junk"
        with self.assertRaises(Exception) as context:
            pm.validate_perturbagens(perts)
        assert context.exception is not None
        logger.debug("context.exception:  {}".format(context.exception))
        assert "the perturbagens provided contain different compounds in the same wells of different assasy plates" in str(context.exception), str(context.exception)

    def test__parse_raw_value(self):
        r = pm._parse_raw_value("")
        assert r is None

        r = pm._parse_raw_value("6")
        assert r == 6, r

        r = pm._parse_raw_value("6.1")
        assert r == 6.1, r

        r = pm._parse_raw_value("hello world")
        assert r == "hello world", r

    def test_read_perturbagen_from_CMap_file(self):
        #1 perturbagen per well
        r = pm.read_perturbagen_from_CMap_file("functional_tests/test_prism_metadata/LJP005.src")
        assert len(r) > 0, len(r)
        logger.debug("r[0]:  {}".format(r[0]))
        logger.debug("r[6]:  {}".format(r[6]))
        assert hasattr(r[0], "pert_id")
        assert hasattr(r[0], "pert_type")
        assert r[0].well_id is not None

        #multiple perturbagens per well
        r = pm.read_perturbagen_from_CMap_file("functional_tests/test_prism_metadata/PMEL.A001.src")
        assert len(r) > 0, len(r)
        logger.debug("r[0]:  {}".format(r[0]))
        logger.debug("r[6]:  {}".format(r[6]))
        assert hasattr(r[0], "pert_id")
        assert hasattr(r[0], "pert_type")
        assert r[0].well_id is not None


if __name__ == "__main__":
    setup_logger.setup(verbose=True)

    unittest.main()
