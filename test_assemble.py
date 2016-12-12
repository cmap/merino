import assemble
import unittest
import setup_logger
import logging
import prism_metadata
import davepool_data
import check_and_build_perts
import numpy
import os
import glob


logger = logging.getLogger(setup_logger.LOGGER_NAME)


class TestAssemble(unittest.TestCase):
    def test_parse_location_to_well(self):
        r = assemble.parse_location_to_well("1(1,A1)")
        assert r == "A01", r

        r = assemble.parse_location_to_well("384(1,P24)")
        assert r == "P24", r

    def test_read_davepool_data_objects(self):
        l = [(1, "requirements_artifacts/PCAL001_P1_X1.csv"), (2, "requirements_artifacts/PCAL001_P2_X1.csv")]
        r = assemble.read_davepool_data_objects(l)
        assert len(r) > 0
        logger.debug("r:  {}".format(r))

        assert r[0].davepool_id == 1, r[0].davepool_id
        assert r[1].davepool_id == 2, r[1].davepool_id

    def test_build_davepool_to_cells_map(self):
        l = [prism_metadata.PrismCell(pool_id=str(x+20), analyte_id=str(x+10), davepool_id=str(x/2)) for x in range(6)]
        logger.debug("l:  {}".format(l))

        r = assemble.build_davepool_id_to_cells_map(l)
        assert len(r) == 3, len(r)
        logger.debug("r:  {}".format(r))

        for x in range(3):
            y = str(x)
            assert y in r, y
            assert len(r[y]) == 2, r[y]

    def test_combine_maps_with_checks(self):
        a = {1:2, 3:5}
        r = {7:11, 13:17}

        assemble.combine_maps_with_checks(a, r)
        logger.debug("b:  {}".format(r))

        for (k,v) in a.items():
            assert k in r
            assert r[k] == v

        with self.assertRaises(Exception) as context:
            assemble.combine_maps_with_checks(a, r)
        assert context.exception is not None
        logger.debug("context.exception:  {}".format(context.exception))
        assert "source_map and dest_map had common_keys" in str(context.exception), str(context.exception)

    def test_build_data_by_cell(self):
        cells = [prism_metadata.PrismCell(pool_id=str(x+20), analyte_id=str(x+10), davepool_id=str(x/2)) for x in range(2)]
        logger.debug("cells:  {}".format(cells))

        davepool_data_obj = davepool_data.DavepoolData()
        davepool_data_obj.median_headers = ["Location", "10", "11"]
        davepool_data_obj.median_headers.extend([str(x) for x in range(30,40)])

        davepool_data_obj.median_data = [["(1,A01)", 1,2], ["(1,J13)", 3,5]]
        davepool_data_obj.median_data[0].extend([x for x in range(40,50)])
        davepool_data_obj.median_data[1].extend([x for x in range(50,60)])

        davepool_data_obj.count_headers = davepool_data_obj.median_headers
        davepool_data_obj.count_data = [["(1,B03)", 7,11], ["(1,L17)", 13,17]]
        davepool_data_obj.count_data[0].extend([x for x in range(60,70)])
        davepool_data_obj.count_data[1].extend([x for x in range(70,80)])
        logger.debug("davepool_data_obj:  {}".format(davepool_data_obj))

        r = assemble.build_data_by_cell(cells, davepool_data_obj)

        r_med = r[0]
        assert r_med is not None
        logger.debug("r_med:  {}".format(r_med))
        assert r_med.well_list is not None
        assert len(r_med.well_list) == 2, len(r_med.well_list)
        assert "A01" in r_med.well_list, r_med.well_list
        assert "J13" in r_med.well_list, r_med.well_list

        assert len(r_med.cell_data_map) == len(cells), (len(r_med.cell_data_map), len(cells))
        for (i, c) in enumerate(cells):
            assert c in r_med.cell_data_map
            r_med_data = r_med.cell_data_map[c]
            for (j, d) in enumerate(davepool_data_obj.median_data):
                assert d[i+1] == r_med_data[j], (i, j, d[i+1], r_med_data[j])

        r_count = r[1]
        assert len(r_count.well_list) == 2, len(r_count.well_list)
        logger.debug("r_count.well_list:  {}".format(r_count.well_list))
        assert "B03" in r_count.well_list, r_count.well_list
        assert "L17" in r_count.well_list, r_count.well_list

        assert len(r_count.cell_data_map) == len(cells), (len(r_count.cell_data_map), len(cells))
        logger.debug("r_count.cell_data_map:  {}".format(r_count.cell_data_map))
        for (i, c) in enumerate(cells):
            assert c in r_count.cell_data_map
            r_count_data = r_count.cell_data_map[c]
            for (j, d) in enumerate(davepool_data_obj.count_data):
                assert d[i+1] == r_count_data[j], (i, j, d[i+1], r_count_data[j])

    def test_process_data(self):
        cells = [prism_metadata.PrismCell(pool_id=str(x+20), analyte_id=str(x+10), davepool_id=str(x/2)) for x in range(4)]
        logger.debug("cells:  {}".format(cells))
        cells_map = assemble.build_davepool_id_to_cells_map(cells)
        logger.debug("cells_map:  {}".format(cells_map))

        davepool_list = []
        davepool_data_obj = davepool_data.DavepoolData()
        davepool_list.append(davepool_data_obj)

        davepool_data_obj.davepool_id = "0"
        davepool_data_obj.median_headers = [str(x) for x in range(10,22)]
        davepool_data_obj.median_headers.insert(0, "Location")
        davepool_data_obj.median_data = [["(1,A01)", 1,2], ["(1,J13)", 3,5]]
        davepool_data_obj.median_data[0].extend([x for x in range(40,50)])
        davepool_data_obj.median_data[1].extend([x for x in range(50,60)])

        davepool_data_obj.count_headers = davepool_data_obj.median_headers
        davepool_data_obj.count_data = [["(1,A01)", 7,11], ["(1,J13)", 13,17]]
        davepool_data_obj.count_data[0].extend([x for x in range(60,70)])
        davepool_data_obj.count_data[1].extend([x for x in range(70,80)])
        logger.debug("davepool_data_obj:  {}".format(davepool_data_obj))

        davepool_data_obj = davepool_data.DavepoolData()
        davepool_list.append(davepool_data_obj)

        davepool_data_obj.davepool_id = "1"
        davepool_data_obj.median_headers = [str(x) for x in range(10,22)]
        davepool_data_obj.median_headers.insert(0, "Location")
        davepool_data_obj.median_data = [["(1,A01)", -1, -2, 19,23], ["(1,J13)", -3, -5, 29,31]]
        davepool_data_obj.median_data[0].extend([x for x in range(80,88)])
        davepool_data_obj.median_data[1].extend([x for x in range(90,98)])

        davepool_data_obj.count_headers = davepool_data_obj.median_headers
        davepool_data_obj.count_data = [["(1,A01)", -7, -11, 37,41], ["(1,J13)", -13, -17, 43,47]]
        davepool_data_obj.count_data[0].extend([x for x in range(100,108)])
        davepool_data_obj.count_data[1].extend([x for x in range(110,118)])
        logger.debug("davepool_data_obj:  {}".format(davepool_data_obj))

        r = assemble.process_data(davepool_list, cells_map)
        assert len(r) == 2, len(r)

        median_data_by_cell = r[0]
        assert median_data_by_cell is not None
        logger.debug("median_data_by_cell:  {}".format(median_data_by_cell))
        assert 2 == median_data_by_cell.cell_data_map[cells[1]][0]
        assert 3 == median_data_by_cell.cell_data_map[cells[0]][1]
        t_wells = median_data_by_cell.well_list
        logger.debug("t_wells:  {}".format(t_wells))
        assert "J13" in t_wells

        assert 23 == median_data_by_cell.cell_data_map[cells[3]][0]
        assert 29 == median_data_by_cell.cell_data_map[cells[2]][1]


    def test_build_gctoo(self):
        filepath = "functional_tests/test_assemble/test_write_output_gct.txt"
        if os.path.exists(filepath):
            os.remove(filepath)

        prn = "my_prism_replicate"

        pert_list = []
        p = prism_metadata.Perturbagen()
        p.well_id = "J01"
        p.pert_id = "BRD-K12345678"
        p.pert_dose = 11.0
        p.pert_dose_unit = "uM"
        pert_list.append(p)
        p = prism_metadata.Perturbagen()
        p.well_id = "M03"
        p.pert_id = "BRD-K91011121"
        p.pert_dose = 13.0
        p.pert_dose_unit = "mg/mL"
        pert_list.append(p)
        p = prism_metadata.Perturbagen()
        p.well_id = "B02"
        p.pert_id = "BRD-K31415171"
        p.pert_dose = 17.0
        p.pert_dose_unit = "nM"
        pert_list.append(p)

        cell_list = []
        c = prism_metadata.PrismCell(pool_id="fake pool", analyte_id="fake analyte",
            davepool_id="fake davepool", id=3)
        cell_list.append(c)
        c = prism_metadata.PrismCell(pool_id="fake pool 2", analyte_id="fake analyte 2",
            davepool_id="fake davepool 2", id=5)
        cell_list.append(c)
        c = prism_metadata.PrismCell(pool_id="fake pool 2", analyte_id="fake analyte 3",
            davepool_id="fake davepool 2", id=7)
        cell_list.append(c)

        data_by_cell = assemble.DataByCell({cell_list[1]:[1,2,11], cell_list[2]:[13,17,19], cell_list[0]:[23,29,31]},
                                           ["J01", "M03", "B02"])

        r = assemble.build_gctoo(prn, pert_list, data_by_cell)
        self.assertIsNotNone(r)
        logger.debug("r:  {}".format(r))
        logger.debug("r.col_metadata_df:  {}".format(r.col_metadata_df))
        logger.debug("r.row_metadata_df:  {}".format(r.row_metadata_df))
        logger.debug("r.data_df:  {}".format(r.data_df))

    def test_build_gctoo_data_df(self):
        #happy path - all numbers
        cell_id_data_map = {"cell1":[0,1,2,3], "cell2":[5,7,11,13]}
        data_df_column_ids = ["s1", "s2", "s3", "s5"]
        r = assemble.build_gctoo_data_df(cell_id_data_map, data_df_column_ids)
        self.assertIsNotNone(r)
        logger.debug("r:  {}".format(r))
        self.assertEquals(True, all(data_df_column_ids == r.columns),
                          "columns of result do not match provided columns - data_df_column_ids:  {}  r.columns:  {}".format(
                              data_df_column_ids, r.columns))
        #check that index entries are present and sorted
        self.assertEquals("cell1", r.index[0])
        self.assertEquals("cell2", r.index[1])

        #2nd happy path - blank entry present in data
        cell_id_data_map["cell1"][1] = ""
        cell_id_data_map["cell2"][2] = ""
        r = assemble.build_gctoo_data_df(cell_id_data_map, data_df_column_ids)
        self.assertIsNotNone(r)
        logger.debug("r:  {}".format(r))
        s = r["s2"].loc["cell1"]
        logger.debug("s:  {}".format(s))
        self.assertEquals(assemble._NaN, s, "expected blank empty string to be replaced with _NaN, was not")

        s = r["s3"].loc["cell2"]
        logger.debug("s:  {}".format(s))
        self.assertEquals(assemble._NaN, s, "expected blank empty string to be replaced with _NaN, was not")

    def test_build_davepool_id_csv_list(self):
        r = assemble.build_davepool_id_csv_list(["a", "1", "b", "2", "c", "3"])
        logger.debug("r:  {}".format(r))
        assert len(r) == 3, len(r)

        assert r[0][0] == "a", r[0][0]
        assert r[0][1] == "1", r[0][1]
        assert r[2][0] == "c", r[2][0]
        assert r[2][1] == "3", r[2][1]

    def test_full_functional(self):
        expected_files = ["PCAL003_CS1_X1_COUNT.gct", "PCAL003_CS1_X1_MEDIAN.gct"]
        for ef in expected_files:
            if os.path.exists(ef):
                os.remove(ef)

        config_filepath = "/Users/elemire/Workspace/prism_pipeline/functional_tests/test_assemble/full_functional_test/prism_pipeline.cfg"
        prism_replicate_name = "PCAL003_CS1_X1"
        plates_mapping_path = "functional_tests/test_assemble/full_functional_test/2016-03-22_PCAL_plate_mapping.txt"
        args1 = check_and_build_perts.build_parser().parse_args(
            ["functional_tests/test_assemble/full_functional_test/prism_pipeline.cfg",
             "functional_tests/test_assemble/full_functional_test/7159-03-A04-01-01_03-22-16_15.34.24.txt",
             "functional_tests/test_assemble/full_functional_test/2016-03-22_PCAL_plate_mapping.txt"])
        check_and_build_perts.main(args1)
        plate_map_path = 'PCAL003.src'
        dp7_csv_path = "functional_tests/test_assemble/full_functional_test/PCAL003_DP7_X1.csv"
        dp8_csv_path = "functional_tests/test_assemble/full_functional_test/PCAL003_DP8_X1.csv"
        csdf_path = "/Users/elemire/Workspace/prism_pipeline/requirements_artifacts/CalicoTranche1PrimaryMetaData_02252016.txt"
        dmf_path = "/Users/elemire/Workspace/prism_pipeline/requirements_artifacts/test_davepool_analyte_mapping.txt"

        args = assemble.build_parser().parse_args(["-config_filepath", config_filepath, "-prn", prism_replicate_name,
            "-pmp", plate_map_path, "-ptp", plates_mapping_path, "-dp_csv", "DP7", dp7_csv_path, "DP8", dp8_csv_path, '-plate_map_type', 'CMap',
                                                   "-csdf", csdf_path, "-dmf", dmf_path])

        logger.debug("args:  {}".format(args))

        assemble.main(args)

        for ef in expected_files:
            assert os.path.exists(ef), ef
            os.remove(ef)

        for map_file in glob.glob('PCAL*.src'):
            x = os.path.getsize(map_file)
            assert x > 0
            os.remove(map_file)


if __name__ == "__main__":
    setup_logger.setup(verbose=True)

    unittest.main()
