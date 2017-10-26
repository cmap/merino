import modz
import glob
import os


def weave(proj_dir, rep_set, input_folder='zscorepc'):

        print rep_set
        files = glob.glob(os.path.join(proj_dir, input_folder, rep_set + '*'))
        plate_names = [os.path.basename(x) for x in files]
        replicate_ids = [x.split("_")[2] for x in plate_names]
        short_reps = [x.split('.')[0] for x in replicate_ids]
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

        if not os.path.exists(os.path.join(proj_dir, 'modz', pert)):
            os.mkdir(os.path.join(proj_dir, 'modz', pert))

        reload(modz)
        modz.calculate_modz(keep_files, proj_dir)


def weave_all(proj_dir):

    if not os.path.exists(os.path.join(proj_dir, 'modz')):
        os.mkdir(os.path.join(proj_dir, 'modz'))
    for x in set([os.path.basename(y).split('_')[0] for y in glob.glob(os.path.join(proj_dir, 'zscorepc/*'))]):
        weave(proj_dir, x)