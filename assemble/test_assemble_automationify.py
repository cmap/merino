import unittest
import setup_logger
import logging
import assemble_automationify as aa
import sys
sys.path.append('/Users/elemire/Workspace/caldaia_root/caldaia')
import utils.mysql_utils as mu
import os

db = mu.DB(host="localhost").db
cursor = db.cursor()
det_plate = 'PAVP.A001_DP9_X1'
my_lapo = aa.query_db(cursor, det_plate)


class TestAssembleAutomationify(unittest.TestCase):

    def test_construct_plate_map_path(self):
        map_path = aa.construct_plate_map_path('/cmap/obelix/pod/custom', my_lapo)
        self.assertEquals(map_path, '/cmap/obelix/pod/custom/PAVP/map_src/PAVP.A001.src')

    def test_construct_davepool_csv_path_pairs(self):
        pod_dir = '/cmap/obelix/pod/custom'
        pairs = aa.construct_davepool_csv_path_pairs(pod_dir, my_lapo)
        self.assertEquals(pairs[0][0], 'DP9')
        self.assertEqual(pairs[0][1], '/cmap/obelix/pod/custom/PAVP/lxb/PAVP.A001_DP9_X1/PAVP.A001_DP9_X1.csv')

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

