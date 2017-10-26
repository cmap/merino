import glob
import os
import cmapPy.pandasGEXpress.parse as pe
import distributions as dist
import pandas as pd


def main(proj_dir, out_dir):
    reload(dist)
    # READ EM ALL IN
    norm_path = glob.glob(os.path.join(proj_dir, '*NORM*.gctx'))[0]
    norm_gct = pe(norm_path)

    mfi_path = glob.glob(os.path.join(proj_dir, '*MFI*.gctx'))[0]
    mfi_gct = pe(mfi_path)

    count_path = glob.glob(os.path.join(proj_dir, '*COUNT*.gctx'))[0]
    count_gct = pe(count_path)

    zscore_path = glob.glob(os.path.join(proj_dir, '*ZSPC*.gctx'))[0]
    zscore_gct = pe(zscore_path)

    viability_path = glob.glob(os.path.join(proj_dir, '*FCPC*.gctx'))[0]
    viability_gct = pe(viability_path)

    modz_path = glob.glob(os.path.join(proj_dir, '*MODZ*.gctx'))[0]
    modz_gct = pe(modz_path)

    inst_info = pd.read_table(os.path.join(proj_dir, 'inst_info.txt'), index_col='cid')

    dist.distributions(norm_gct, mfi_gct, count_gct, zscore_gct, viability_gct, modz_gct, inst_info, os.path.join(out_dir))


