import os


def cut_l1(original_list):
    new_list = []

    for pert in set(['_'.join(os.path.basename(y).split('_')[:-2]) for y in original_list]):
        print pert
        files = [x for x in original_list if os.path.basename(x).startswith(pert + '_')]
        plate_names = [os.path.basename(x) for x in files]
        replicate_ids = [x.split("_")[-2] for x in plate_names]
        #short_reps = [x[:2] for x in replicate_ids]
        keep = []

        for r in set(replicate_ids):
            temp = [y for y in replicate_ids if y.startswith(r)]

            if len(temp) == 1:
                keep.append(temp[0])
            else:

                temp2 = [z for z in temp if ".L" in z]

                if len(temp2) == 0:
                    import pdb
                    pdb.set_trace()
                max_l = max([int(x[-1]) for x in temp2])
                temp3 = [b for b in temp2 if b.endswith(str(max_l))]
                keep.append(temp3[0])

        keep_perts = [x for x in plate_names if x.split('_')[-2] in keep]
        keep_files = [path for path in files if os.path.basename(path) in keep_perts]
        [new_list.append(x) for x in keep_files]

    return new_list