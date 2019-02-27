import batch_adjust
import cmapPy.pandasGEXpress.GCToo as GCToo
import cmapPy.pandasGEXpress.parse as pe
import cmapPy.pandasGEXpress.write_gct as wg
import unittest
#import merino.setup_logger as setup_logger
import logging
import numpy
import os
import glob
import pandas as pd
import shutil


#logger = logging.getLogger(setup_logger.LOGGER_NAME)


class TestMerino(unittest.TestCase):
    def test_batch_adjust(self):
        gct_files = glob.glob(os.path.join('functional_tests/test_merino/for_batch_adjust', '*ZSCOREPC.gct'))
        gct_list = []
        for gct_path in gct_files:
            gct = pe.parse(gct_path)
            gct_list.append(gct)
                    
        adj_ds, adj_list = batch_adjust.combat_by_group(gct_list, col_group='pert_well', batch_field='pool_id', use_col_group_as_batch=True)

        wg.write(adj_ds, 'batch_adjusted_values.gct')
        for ctr, ds in enumerate(adj_list):
            print ds.src
            out_file='{}.COMBAT.gct'.format(os.path.splitext(os.path.basename(ds.src))[0])
            wg.write(ds, out_file)
            
        assert(len(adj_ds)==len(gct_list))
        
        # Check values
        #assert round(r.data_df['test_CS0_X1:A'][0], 4) == 1.415

if __name__ == "__main__":
 #   setup_logger.setup(verbose=True)

    unittest.main()
