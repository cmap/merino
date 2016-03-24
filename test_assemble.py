import assemble
import unittest
import setup_logger
import logging
import prism_metadata
import davepool_data
import numpy
import os


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

    def test_build_cell_to_data_maps(self):
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
        r_med_wells = r_med[1]
        assert len(r_med_wells) == 2, len(r_med_wells)
        logger.debug("r_med_wells:  {}".format(r_med_wells))
        assert "A01" in r_med_wells, r_med_wells
        assert "J13" in r_med_wells, r_med_wells

        r_med_map = r_med[0]
        assert len(r_med_map) == len(cells), (len(r_med_map), len(cells))
        logger.debug("r_med_map:  {}".format(r_med_map))
        for (i, c) in enumerate(cells):
            assert c in r_med_map
            r_med_data = r_med_map[c]
            for (j, d) in enumerate(davepool_data_obj.median_data):
                assert d[i+1] == r_med_data[j], (i, j, d[i+1], r_med_data[j])

        r_count = r[1]
        r_count_wells = r_count[1]
        assert len(r_count_wells) == 2, len(r_count_wells)
        logger.debug("r_count_wells:  {}".format(r_count_wells))
        assert "B03" in r_count_wells, r_count_wells
        assert "L17" in r_count_wells, r_count_wells

        r_count_map = r_count[0]
        assert len(r_count_map) == len(cells), (len(r_count_map), len(cells))
        logger.debug("r_count_map:  {}".format(r_count_map))
        for (i, c) in enumerate(cells):
            assert c in r_count_map
            r_count_data = r_count_map[c]
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
        davepool_data_obj.count_data = [["(1,B03)", 7,11], ["(1,L17)", 13,17]]
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
        davepool_data_obj.count_data = [["(1,B03)", -7, -11, 37,41], ["(1,L17)", -13, -17, 43,47]]
        davepool_data_obj.count_data[0].extend([x for x in range(100,108)])
        davepool_data_obj.count_data[1].extend([x for x in range(110,118)])
        logger.debug("davepool_data_obj:  {}".format(davepool_data_obj))

        r = assemble.process_data(davepool_list, cells_map)
        assert len(r) == 2, len(r)

        med = r[0]
        assert len(med) == 2, len(med)
        t = med[0]
        t_data = t[0]
        logger.debug("t_data:  {}".format(t_data))
        assert 2 == t_data[cells[1]][0]
        assert 3 == t_data[cells[0]][1]
        t_wells = t[1]
        logger.debug("t_wells:  {}".format(t_wells))
        assert "J13" in t_wells

        t = med[1]
        t_data = t[0]
        logger.debug("t_data:  {}".format(t_data))
        assert 23 == t_data[cells[3]][0]
        assert 29 == t_data[cells[2]][1]
        t_wells = t[1]
        logger.debug("t_wells:  {}".format(t_wells))
        assert "J13" in t_wells

    def test_generate_sorted_unique_cells_and_wells(self):
        cells = [prism_metadata.PrismCell(id="11"), prism_metadata.PrismCell(id="1"), prism_metadata.PrismCell(id="2")]
        data_by_cells = [assemble.DataByCell({cells[0]:[1,2,3,5], cells[1]:[7,11,13,17]}, ["D01","B02","J03","D04"])]
        data_by_cells.append(assemble.DataByCell({cells[1]:[19,2,3,5], cells[2]:[23,11,13,17]}, ["D01","B03","K13","M21"]))

        expected_cells = [cells[1], cells[2], cells[0]]
        expected_wells = set(data_by_cells[0].well_list)
        expected_wells.update(data_by_cells[1].well_list)
        expected_wells = list(expected_wells)
        expected_wells.sort()

        r = assemble.generate_sorted_unique_cells_and_wells(data_by_cells)
        logger.debug("r:  {}".format(r))
        assert len(r) == 2, len(r)

        assert expected_cells == r[0], (expected_cells, r[0])
        assert expected_wells == r[1], (expected_wells, r[1])

    def test_build_matrix_and_annotations(self):
        cells = [prism_metadata.PrismCell(id=(x+30), pool_id=str(x+20), analyte_id=str(x+10), davepool_id=str(x/2)) for x in range(4)]
        data_by_cells = [assemble.DataByCell({cells[0]:[1,2,3,5], cells[1]:[7,11,13,17]}, ["A01","B02","C03","D04"])]
        data_by_cells.append(assemble.DataByCell({cells[1]:[19,2,3,5], cells[2]:[23,11,13,17]}, ["D01","B03","K13","M21"]))

        #expected column-well ordering:  A01 B02 B03 C03 D01 D04 K13 M21
        #expected matrix (generated by mapping cell data above onto expected column-well ordering):
        #1 2 nan 3 nan 5 nan nan
        #7 11 2 13 19 17 3 5
        #nan nan 11 nan 23 nan 13 17
        r = assemble.build_matrix_and_annotations(data_by_cells)
        r = r[2]
        logger.debug("r:  {}".format(r))

        r_shape = r.shape
        logger.debug("r_shape:  {}".format(r_shape))
        assert r_shape[0] == 3, r_shape[0]
        assert r_shape[1] == 8, r_shape[1]

        assert r[0][0] == 1, r[0][0]
        assert numpy.isnan(r[0][2]), r[0][2]
        assert numpy.isnan(r[0][4]), r[0][4]
        assert numpy.isnan(r[0][6]), r[0][6]
        assert r[1][1] == 11, r[1][1]
        assert r[1][3] == 13, r[1][3]
        assert r[1][5] == 17, r[1][5]
        assert r[1][7] == 5, r[1][7]
        assert numpy.isnan(r[2][0]), r[2][0]
        assert r[2][2] == 11, r[2][2]
        assert r[2][4] == 23, r[2][4]
        assert r[2][6] == 13, r[2][6]

    def test_write_gct_version_and_size(size):
        filepath = "functional_tests/test_assemble/test_write_gct_version_and_size.txt"
        if os.path.exists(filepath):
            os.remove(filepath)
        f = open(filepath, "w")

        p = prism_metadata.Perturbagen()
        p.extra_field_one = 1
        p.extra_field_two = 2

        c = prism_metadata.PrismCell()
        c.extra_field_one = 3

        matrix_and_annots = assemble.MatrixAndAnnots([c], None, numpy.empty((2,7)))

        assemble.write_gct_version_and_size(f, p, matrix_and_annots)
        f.close()

        f = open(filepath)
        r = f.read().strip().split("\n")
        f.close()
        logger.debug("r:  {}".format(r))

        assert len(r) == 2, len(r)
        assert r[0] == assemble._gct_version, r[0]

        r = r[1].split("\t")
        assert len(r) == 4, len(r)
        assert r[0] == "2", r[0]
        assert r[1] == "7", r[1]
        assert r[2] == "4", r[2]
        assert r[3] == "3", r[3]

    def test_generate_column_headers(self):
        prn = "my_prism_replicate"
        c = prism_metadata.PrismCell()
        wells = ["A01","B02","C03","D05"]

        expected = ["id", "analyte_id", "davepool_id", "pool_id"]
        expected.extend([prn + ":" + w for w in wells])

        r = assemble.generate_column_headers(prn, c, wells)
        logger.debug("r:  {}".format(r))
        assert len(r) == 2, len(r)
        h = r[0]
        assert len(h) == len(expected), (len(h), len(expected))
        assert h == expected, (h, expected)

        cao = r[1]
        assert len(cao) == 3, len(cao)
        assert cao == expected[1:4], (cao, expected[1:4])

    def test_generate_perturbagen_annotation_header_block(self):
        pert_list = []
        p = prism_metadata.Perturbagen()
        p.well_id = "J01"
        p.pert_id = "BRD-K12345678"
        p.pert_dose = 11.0
        p.pert_dose_unit = "uM"
        pert_list.append(p)
        p = prism_metadata.Perturbagen()
        p.well_id = "B02"
        p.pert_id = "BRD-K91011121"
        p.pert_dose = 13.0
        p.pert_dose_unit = "mg/mL"
        p.extra_annotation = "my extra annotation"
        pert_list.append(p)

        num_cell_annot = 3

        r = assemble.generate_perturbagen_annotation_header_block(pert_list, ["B02", "J01"], num_cell_annot)
        logger.debug("r:  {}".format(r))
        assert len(r) == 5, len(r)

        row = r[0]
        logger.debug("row:  {}".format(row))
        assert len(row) == 6, len(row)
        assert row[0] == "extra_annotation", row[0]
        assert row[1:4] == [assemble._null for i in range(3)], row[1:4]
        assert row[4] == "my extra annotation", row[4]
        assert row[5] == assemble._null, row[5]

        row = r[1]
        logger.debug("row:  {}".format(row))
        assert len(row) == 6, len(row)
        assert row[0] == "pert_dose", row[0]
        assert row[1:4] == [assemble._null for i in range(3)], row[1:4]
        assert row[4] == 13.0, row[4]
        assert row[5] == 11.0, row[5]

        row = r[3]
        logger.debug("row:  {}".format(row))
        assert len(row) == 6, len(row)
        assert row[0] == "pert_id", row[0]
        assert row[1:4] == [assemble._null for i in range(3)], row[1:4]
        assert row[4] == "BRD-K91011121", row[4]
        assert row[5] == "BRD-K12345678", row[5]

    def test_generate_row_annotation_and_data_block(self):
        cell_list = []
        c = prism_metadata.PrismCell(pool_id="fake pool", analyte_id="fake analyte",
            davepool_id="fake davepool", id=3)
        c.extra_annot = "my extra annotation"
        cell_list.append(c)
        c = prism_metadata.PrismCell(pool_id="fake pool 2", analyte_id="fake analyte 2",
            davepool_id="fake davepool 2", id=5)
        c.extra_annot = None
        cell_list.append(c)

        matrix = numpy.empty((2,2))
        matrix[:] = [[1,2],[7,float('nan')]]
        matrix_and_annots = assemble.MatrixAndAnnots(cell_list, ["B02", "J01"], matrix)

        r = assemble.generate_row_annotation_and_data_block(matrix_and_annots, ["analyte_id", "davepool_id", "pool_id",
                                                                                "extra_annot"])
        logger.debug("r:  {}".format(r))

        assert len(r) == 2, len(r)
        row = r[0]
        assert len(row) == 7, len(row)
        assert row[0] == 3, row[0]
        assert row[1] == "fake analyte", row[1]
        assert row[2] == "fake davepool", row[2]
        assert row[3] == "fake pool", row[3]
        assert row[4] == "my extra annotation", row[4]
        assert row[5] == 1, row[5]
        assert row[6] == 2, row[6]

        row = r[1]
        assert len(row) == 7, len(row)
        assert row[0] == 5, row[0]
        assert row[1] == "fake analyte 2", row[1]
        assert row[2] == "fake davepool 2", row[2]
        assert row[3] == "fake pool 2", row[3]
        assert row[4] == assemble._null, row[4]
        assert row[5] == 7, row[5]
        assert row[6] == assemble._NaN, row[6]

    def test_write_output_gct(self):
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
        p.well_id = "B02"
        p.pert_id = "BRD-K91011121"
        p.pert_dose = 13.0
        p.pert_dose_unit = "mg/mL"
        pert_list.append(p)

        cell_list = []
        c = prism_metadata.PrismCell(pool_id="fake pool", analyte_id="fake analyte",
            davepool_id="fake davepool", id=3)
        cell_list.append(c)
        c = prism_metadata.PrismCell(pool_id="fake pool 2", analyte_id="fake analyte 2",
            davepool_id="fake davepool 2", id=5)
        cell_list.append(c)

        matrix = numpy.empty((2,2))
        matrix[:] = [[1,2],[7,11]]
        matrix_and_annots = assemble.MatrixAndAnnots(cell_list, ["B02", "J01"], matrix)

        assemble.write_output_gct(filepath, prn, pert_list, matrix_and_annots)

        f = open(filepath)
        r = f.read().strip().split("\n")
        f.close()
        logger.debug("r:  {}".format(r))
        assert len(r) == 9, "expect 9:  version, size, header row, 4 pert annot, 2 data rows, found len(r):  {}".format(len(r))

    def test_build_davepool_id_csv_list(self):
        r = assemble.build_davepool_id_csv_list(["a","1","b","2","c","3"])
        logger.debug("r:  {}".format(r))
        assert len(r) == 3, len(r)

        assert r[0][0] == "a", r[0][0]
        assert r[0][1] == "1", r[0][1]
        assert r[2][0] == "c", r[2][0]
        assert r[2][1] == "3", r[2][1]


if __name__ == "__main__":
    setup_logger.setup(verbose=True)

    unittest.main()
