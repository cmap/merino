import unittest
import setup_logger
import logging
import assemble_automationify as aa
import caldaia.utils.mysql_utils as mu
import os

db = mu.DB(host="getafix-v-dev").db
cursor = db.cursor()

det_plate = 'PAVP.A001_DP9_X1'
my_lapo = aa.query_db(cursor, det_plate)


class TestAssembleAutomationify(unittest.TestCase):

    def test_query_db(self):
        self.assertEquals(my_lapo.prism_replicate_name, 'PAVP.A001_DP9_X1')
        self.assertEquals(my_lapo.cell_set_definition_file,
                          '/cmap/data/vdb/PRISM/cell_set_definitions/PAVP_cell_set.txt')
        self.assertEquals(my_lapo.davepool_mapping_file, '/cmap/obelix/pod/custom/PAVP/map_src/DP9_mapping.txt')
        self.assertEquals(my_lapo.plate_tracking_path, '/cmap/obelix/pod/custom/PAVP/map_src/PAVP_plate_tracking2.txt')
        self.assertEquals(my_lapo.project_id, 'PAVP')
        self.assertEquals(my_lapo.pert_plate, 'PAVP.A001')
        self.assertEquals(my_lapo.davepool_group_id, 3)
        self.assertEquals(my_lapo.davepool_list[0][0], 'DP9')


    def test_construct_plate_map_path(self):
        map_path = aa.construct_plate_map_path('/cmap/obelix/pod/custom', my_lapo)
        self.assertEquals(map_path, '/cmap/obelix/pod/custom/PAVP/map_src/PAVP.A001.src')

    def test_construct_davepool_csv_path_pairs(self):
        pairs = aa.construct_davepool_csv_path_pairs(my_lapo)
        self.assertEquals(pairs, 'DP9 /cmap/obelix/pod/custom/PAVP/lxb/PAVP.A001_DP9_X1/PAVP.A001_DP9_X1.csv')

    def test_construct_outfile_path(self):
        outfile = aa.construct_outfile_path(my_lapo)
        self.assertEquals(outfile, '/cmap/obelix/pod/custom/PAVP/assemble/PAVP.A001_DP9_X1')

    def test_build_args(self):
        default_config_filepath = os.path.expanduser('~/.prism_pipeline.cfg')
        args = aa.build_args(cursor, det_plate, default_config_filepath)
        self.assertIsInstance(args, dict)
        self.assertEquals(len(args.keys()), 9)


if __name__ == "__main__":
    setup_logger.setup(verbose=True)

    unittest.main()

