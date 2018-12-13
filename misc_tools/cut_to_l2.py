import os


def cut_l1(list_of_plate_paths):
    curated_plate_path_list = []
    all_replicate_set_names = set(['_'.join(os.path.basename(replicate_path).split('_')[:-2]) for replicate_path in list_of_plate_paths])

    for replicate_set_name in all_replicate_set_names:
        print replicate_set_name
        replicate_set_paths = [path for path in list_of_plate_paths if os.path.basename(path).startswith(replicate_set_name)]
        replicate_set_plate_names = [os.path.basename(path) for path in replicate_set_paths]
        replicate_nums = [replicate.split("_")[-2] for replicate in replicate_set_plate_names]
        base_replicate_nums = [replicate[:2] for replicate in replicate_nums]

        keep = []
        for base_replicate_num in set(base_replicate_nums):
            all_possible_replicates = [replicate for replicate in replicate_nums if replicate.startswith(base_replicate_num)]

            if len(all_possible_replicates) == 1:
                keep.append(all_possible_replicates[0])

            else:
                replicateLXs = [replicate for replicate in all_possible_replicates if ".L" in replicate]

                if len(replicateLXs) == 0:
                    print base_replicate_num
                    raise Exception(UnableToDifferentiateReplicates)

                max_l = max([int(replicate[-1]) for replicate in replicateLXs])
                desired_replicate = [replicate for replicate in replicateLXs if replicate.endswith(str(max_l))]
                keep.append(desired_replicate[0])

        keep_replicates = [replicate for replicate in replicate_set_plate_names if replicate.split('_')[-2] in keep]
        keep_files = [path for path in replicate_set_paths if os.path.basename(path) in keep_replicates]
        [curated_plate_path_list.append(replicate) for replicate in keep_files]

    return curated_plate_path_list

class UnableToDifferentiateReplicates(Exception):
    pass