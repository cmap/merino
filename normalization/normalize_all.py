import normalize_with_prism_invariant as norm
import viability_normalization as viability
import zscore
import os
import glob
import modz

def log_normalize_all(proj_dir, search_pattern = '*', assemble_folder='assemble', out_folder='normalize'):
    #Invariant Normalize everything in project directory

    for folder in glob.glob(os.path.join(proj_dir, assemble_folder, search_pattern)):
        name = os.path.basename(folder)

        if not os.path.exists(os.path.join(proj_dir, out_folder, name)):
            os.mkdir(os.path.join(proj_dir, out_folder, name))

        norm.normalize(os.path.join(folder, name + '_MEDIAN.gct'), os.path.join(proj_dir, out_folder, name, name + '_NORM.gct'))

def old_norm_all(proj_dir, assemble_folder='assemble', out_folder='normalize'):
    #Invariant Normalize everything in project directory

    for folder in glob.glob(os.path.join(proj_dir, assemble_folder, '*')):
        name = os.path.basename(folder)

        if not os.path.exists(os.path.join(proj_dir, out_folder, name)):
            os.mkdir(os.path.join(proj_dir, out_folder, name))

            norm.old_normalize(os.path.join(folder, name + '_MEDIAN.gct'), os.path.join(proj_dir, out_folder, name, name + '_NORM.gct'))


def viability_all(proj_dir, input_folder='normalize', out_folder='viabilityvc', search_pattern1='*', search_pattern2='_NORM.gct'):
    # Calculate viability for everything in project directory
    for folder in glob.glob(os.path.join(proj_dir, input_folder, search_pattern1)):
        name = os.path.basename(folder)

        if not os.path.exists(os.path.join(proj_dir, out_folder, name)):
            os.mkdir(os.path.join(proj_dir, out_folder, name))

        viability.calculate_viability(os.path.join(folder, name + search_pattern2), write=True, outfile=os.path.join(proj_dir, out_folder, name, name + '_FCVC.gct'))

def viability_pc_all(proj_dir, input_folder='normalize', out_folder='viabilitypc', search_pattern1='*', search_pattern2='_NORM.gct'):
    # Calculate viability for everything in project directory
    for folder in glob.glob(os.path.join(proj_dir, input_folder, search_pattern1)):
        name = os.path.basename(folder)

        if not os.path.exists(os.path.join(proj_dir, out_folder, name)):
            os.mkdir(os.path.join(proj_dir, out_folder, name))

        viability.calculate_viability(os.path.join(folder, name + search_pattern2), write=True, outfile=os.path.join(proj_dir, out_folder, name, name + '_FCPC.gct'), plate_control=True)


def viability_from_mfi(proj_dir):
    # Calculate viability for everything in project directory
    for folder in glob.glob(os.path.join(proj_dir, 'assemble/*')):
        name = os.path.basename(folder)

        if not os.path.exists(os.path.join(proj_dir, 'viability_mfi', name)):
            os.mkdir(os.path.join(proj_dir, 'viability_mfi', name))

        viability.calculate_viability(os.path.join(folder, name + '_MEDIAN.gct'), write=True, outfile=os.path.join(proj_dir, 'viability_mfi', name, name + '_MFI.FCVC.gct'))

def viabilitypc_from_mfi(proj_dir):
    # Calculate viability for everything in project directory
    for folder in glob.glob(os.path.join(proj_dir, 'assemble/*')):
        name = os.path.basename(folder)

        if not os.path.exists(os.path.join(proj_dir, 'viabilitypc_mfi', name)):
            os.mkdir(os.path.join(proj_dir, 'viabilitypc_mfi', name))

        viability.calculate_viability(os.path.join(folder, name + '_MEDIAN.gct'), write=True, outfile=os.path.join(proj_dir, 'viabilitypc_mfi', name, name + '_MFI.FCPC.gct'), plate_control=True)


def zscore_all(proj_dir, input_folder='normalize', search_pattern1='*', search_pattern2='_NORM.gct', out_folder='zscorevc'):
    # Calculate zscore with vehicle controlfor everything in project directory
    for folder in glob.glob(os.path.join(proj_dir, input_folder, search_pattern1)):
        name = os.path.basename(folder)

        if not os.path.exists(os.path.join(proj_dir, out_folder, name)):
            os.mkdir(os.path.join(proj_dir, out_folder, name))

        zscore.calculate_zscorevc(os.path.join(folder, name + search_pattern2), write=True, outfile=os.path.join(proj_dir, out_folder, name, name + '_ZSCORE.gct'))

def zscorepc_all(proj_dir, input_folder='normalize', search_pattern1='*', search_pattern2 = '_NORM.gct', out_folder='zscorepc'):
    # Calculate zscore with plate control for everything in project directory
    for folder in glob.glob(os.path.join(proj_dir, input_folder, search_pattern1)):
        name = os.path.basename(folder)

        if not os.path.exists(os.path.join(proj_dir, out_folder, name)):
            os.mkdir(os.path.join(proj_dir, out_folder, name))

        zscore.calculate_zscorepc(os.path.join(folder, name + search_pattern2), write=True, outfile=os.path.join(proj_dir, out_folder, name, name + '_ZSCOREPC.gct'))


def modz_all(proj_dir):
    for x in set([os.path.basename(y).split('_')[0] for y in glob.glob(os.path.join(proj_dir, 'zscorepc/*'))]):
        print x
        files = glob.glob(os.path.join(proj_dir, 'zscorepc', x + '*'))
        plate_names = [os.path.basename(x) for x in files]
        replicate_ids = [x.split("_")[2] for x in plate_names]
        short_reps = [x[:2] for x in replicate_ids]
        keep = []

        for r in set(short_reps):
            temp = [y for y in replicate_ids if y.startswith(r)]

            if len(temp) == 1:
                keep.append(temp[0])
            else:

                temp2 = [z for z in temp if "." in z]

                max_l = max([int(x[-1]) for x in temp2])
                temp3 = [b for b in temp2 if b.endswith(str(max_l))]
                keep.append(temp3[0])

        keep_perts = [x for x in plate_names if x.split('_')[2] in keep]
        keep_files = [glob.glob(path + '/*')[0] for path in files if os.path.basename(path) in keep_perts]

        pert = plate_names[0].split('_')[0] + '_' + plate_names[0].split('_')[1]

        if not os.path.exists(os.path.join(proj_dir, 'modZ', pert)):
            os.mkdir(os.path.join(proj_dir, 'modZ', pert))

        modz.calculate_modz(keep_files, proj_dir)