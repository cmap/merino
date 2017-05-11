import normalize_with_prism_invariant as norm
import viability_normalization as viability
import zscore
import os
import glob

def normalize_all(proj_dir):
    #Invariant Normalize everything in project directory

    for folder in glob.glob(os.path.join(proj_dir, 'assemble/*')):
        name = os.path.basename(folder)

        if not os.path.exists(os.path.join(proj_dir, 'normalize', name)):
            os.mkdir(os.path.join(proj_dir, 'normalize', name))

            norm.normalize(os.path.join(folder, name + '_MEDIAN.gct'), os.path.join(proj_dir, 'normalize', name, name + '_NORM.gct'))


def viability_all(proj_dir):
    # Calculate viability for everything in project directory
    for folder in glob.glob(os.path.join(proj_dir, 'normalize/*')):
        name = os.path.basename(folder)

        if not os.path.exists(os.path.join(proj_dir, 'viability', name)):
            os.mkdir(os.path.join(proj_dir, 'viability', name))

        viability.calculate_viability(os.path.join(folder, name + '_NORM.gct'), write=True, outfile=os.path.join(proj_dir, 'viability', name, name + '_VIABILITY.gct'))


def zscore_all(proj_dir):
    # Calculate zscore with vehicle controlfor everything in project directory
    for folder in glob.glob(os.path.join(proj_dir, 'normalize/*')):
        name = os.path.basename(folder)

        if not os.path.exists(os.path.join(proj_dir, 'zscore', name)):
            os.mkdir(os.path.join(proj_dir, 'zscore', name))

        zscore.calculate_zscorevc(os.path.join(folder, name + '_NORM.gct'), write=True, outfile=os.path.join(proj_dir, 'zscore', name, name + '_ZSCORE.gct'))

def zscorepc_all(proj_dir):
    # Calculate zscore with plate control for everything in project directory
    for folder in glob.glob(os.path.join(proj_dir, 'normalize/*')):
        name = os.path.basename(folder)

        if not os.path.exists(os.path.join(proj_dir, 'zscorepc', name)):
            os.mkdir(os.path.join(proj_dir, 'zscorepc', name))

        zscore.calculate_zscorepc(os.path.join(folder, name + '_NORM.gct'), write=True, outfile=os.path.join(proj_dir, 'zscorepc', name, name + '_ZSCOREPC.gct'))