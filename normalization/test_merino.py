import viability_normalization
import zscore
import normalize_with_prism_invariant as norm
import modz
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


logger = logging.getLogger(setup_logger.LOGGER_NAME)

l = pe('functional_tests/test_merino/assemble/test_CS0_X1/test_CS0_X1.gct')

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
        fv = viability_normalization.calculate_viability(r,plate_control=False)
        fp = viability_normalization.calculate_viability(r,plate_control=True)

        assert round(fv.data_df['test_CS0_X1:A'][0], 4) == 0.4837
        assert round(fv.data_df['test_CS0_X1:B'][2], 4) == 0.2704
        assert round(fv.data_df['test_CS0_X1:D'][1], 4) == 1.0556

        assert round(fp.data_df['test_CS0_X1:A'][0], 4) == 0.6208
        assert round(fp.data_df['test_CS0_X1:B'][2], 4) == 0.4294
        assert round(fp.data_df['test_CS0_X1:D'][1], 4) == 1.4367

    def test_modz(self):
        l2 = pe('functional_tests/test_merino/assemble/test_CS0_X2/test_CS0_X2.gct')
        l3 = pe('functional_tests/test_merino/assemble/test_CS0_X3/test_CS0_X3.gct')

        r = norm.normalize(l)
        r2 = norm.normalize(l2)
        r3 = norm.normalize(l3)

        modz_gct, ccq74, weights = modz.calculate_modz([r,r2,r3])

        print modz_gct.data_df


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
        card.card('functional_tests/test_merino/assemble/test_CS0_X2')



if __name__ == "__main__":
    setup_logger.setup(verbose=True)

    unittest.main()