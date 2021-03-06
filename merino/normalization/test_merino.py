import viability_normalization
import zscore
import normalize_with_prism_invariant as norm
import distil
import card
import weave
import merino
import cmapPy.pandasGEXpress.GCToo as GCToo
import cmapPy.pandasGEXpress.parse as pe
import unittest
import merino.setup_logger as setup_logger
import logging
import numpy
import os
import glob
import pandas as pd
import shutil


logger = logging.getLogger(setup_logger.LOGGER_NAME)

l = pe('functional_tests/test_merino/assemble/test_CS0_X1/test_CS0_X1_MEDIAN.gct')

class TestMerino(unittest.TestCase):
    def test_norm(self):
        r = norm.normalize(l)

        # Column 'C' should have been removed due to low median invariant
        assert len(r.data_df.columns) == 4
        # Check values
        assert round(r.data_df['test_CS0_X1:A'][0], 4) == 1.415
        assert round(r.data_df['test_CS0_X1:B'][2], 4) == 0.7776
        assert round(r.data_df['test_CS0_X1:D'][1], 4) == 3.4919

    def test_zscore(self):
        #TODO split plate control and vehicle control into two functions
        r = norm.normalize(l)
        zv = zscore.calculate_zscore(r, plate_control=False)
        zp = zscore.calculate_zscore(r, plate_control=True)

        assert round(zv.data_df['test_CS0_X1:A'][0], 4) == -10.0
        assert round(zv.data_df['test_CS0_X1:B'][2], 4) == -5.5729
        assert round(zv.data_df['test_CS0_X1:D'][1], 4) == 0.6745

        assert round(zp.data_df['test_CS0_X1:A'][0], 4) == -0.903
        assert round(zp.data_df['test_CS0_X1:B'][2], 4) == -0.7558
        assert round(zp.data_df['test_CS0_X1:D'][1], 4) == 0.8159

    def test_viability_normalization(self):
        r = norm.normalize(l)
        fv = viability_normalization.calculate_viability(r, plate_control=False)
        fp = viability_normalization.calculate_viability(r, plate_control=True)

        assert round(fv.data_df['test_CS0_X1:A'][0], 4) == 0.4837
        assert round(fv.data_df['test_CS0_X1:B'][2], 4) == 0.2704
        assert round(fv.data_df['test_CS0_X1:D'][1], 4) == 1.0556

        assert round(fp.data_df['test_CS0_X1:A'][0], 4) == 0.6208
        assert round(fp.data_df['test_CS0_X1:B'][2], 4) == 0.4294
        assert round(fp.data_df['test_CS0_X1:D'][1], 4) == 1.4367

    def test_modz(self):
        l2 = pe('functional_tests/test_merino/assemble/test_CS0_X2/test_CS0_X2_MEDIAN.gct')
        l3 = pe('functional_tests/test_merino/assemble/test_CS0_X3/test_CS0_X3_MEDIAN.gct')

        r = norm.normalize(l)
        r2 = norm.normalize(l2)
        r3 = norm.normalize(l3)

        z = zscore.calculate_zscore(r, plate_control=True)
        z2 = zscore.calculate_zscore(r2, plate_control=True)
        z3 = zscore.calculate_zscore(r3, plate_control=True)

        modz_gct, ccq74, weights = distil.calculate_modz([z, z2, z3])

        assert round(modz_gct.data_df['test_CS0X1:A'][0], 4) == -0.5816
        assert round(modz_gct.data_df['test_CS0X1:B'][2], 4) == -0.4991
        assert round(modz_gct.data_df['test_CS0X1:D'][1], 4) == 1.001


    def test_count_shear(self):
        count = GCToo.GCToo(data_df=pd.DataFrame({'test_CS0_X1:A': [40, 50, 30, 20, 10],
                                              'test_CS0_X1:B': [110, 80, 60, 40, 30],
                                              'test_CS0_X1:C': [5, 15, 4, 3, 5],
                                              'test_CS0_X1:D': [60, 90, 70, 8, 8],
                                              'test_CS0_X1:E': [75, 85, 60, 9, 10]},
                                             index=['1', '2', '3', '661', '662']),
                        row_metadata_df=pd.DataFrame(index=['1', '2', '3', '661', '662']),
                        col_metadata_df=pd.DataFrame(
                            {'pert_type': ['trt_cp', 'trt_cp', 'trt_cp', 'ctl_vehicle', 'ctl_vehicle']},
                            index=['test_CS0_X1:A', 'test_CS0_X1:B', 'test_CS0_X1:C', 'test_CS0_X1:D', 'test_CS0_X1:E']))
        shear = norm.remove_low_bead_wells(l,count)

        print shear.data_df.shape

        assert len(shear.data_df.columns) == 4

    def test_card(self):
        card.card('functional_tests/test_merino/', 'test_CS0_X2')
        assert os.path.isfile('functional_tests/test_merino/normalize/test_CS0_X2/test_CS0_X2_NORM.gct')
        assert os.path.isfile('functional_tests/test_merino/LFCPC/test_CS0_X2/test_CS0_X2_FCPC.gct')
        assert os.path.isfile('functional_tests/test_merino/LFCVC/test_CS0_X2/test_CS0_X2_FCVC.gct')
        assert os.path.isfile('functional_tests/test_merino/ZSPC/test_CS0_X2/test_CS0_X2_ZSPC.gct')
        assert os.path.isfile('functional_tests/test_merino/ZSVC/test_CS0_X2/test_CS0_X2_ZSVC.gct')

        shutil.rmtree('functional_tests/test_merino/normalize/test_CS0_X2/')
        shutil.rmtree('functional_tests/test_merino/LFCPC/test_CS0_X2/')
        shutil.rmtree('functional_tests/test_merino/LFCVC/test_CS0_X2/')
        shutil.rmtree('functional_tests/test_merino/ZSPC/test_CS0_X2/')
        shutil.rmtree('functional_tests/test_merino/ZSVC/test_CS0_X2/')


    def test_weave(self):
        if os.path.exists('functional_tests/test_merino/modz.zscorepc'):
            shutil.rmtree('functional_tests/test_merino/modz.zscorepc')
        os.mkdir('functional_tests/test_merino/modz.zscorepc')
        weave.weave('functional_tests/test_merino/', 'test')
        assert os.path.isfile('functional_tests/test_merino/modz.zscorepc/test_CS0/test_CS0_MODZ.zscorepc.gct')
        assert os.path.isfile('functional_tests/test_merino/modz.zscorepc/test_CS0/test_cc_q75.txt')
        shutil.rmtree('functional_tests/test_merino/modz.zscorepc')


    def test_merino(self):
        os.system(
            'python merino_main.py -pd functional_tests/test_merino/ -cn test_build -bf functional_tests/test_merino/build')



if __name__ == "__main__":
    setup_logger.setup(verbose=True)

    unittest.main()