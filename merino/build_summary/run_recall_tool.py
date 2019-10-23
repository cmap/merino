import pandas as pd
import glob
import os

def mk_ds_list(proj_dir, input_type):
    filepaths = []
    group_ids = []
    search = proj_dir + '*/*{}.gct'.format(input_type)
    print proj_dir
    print search
    for path in glob.glob(search):
        filepaths.append(path)
        group_ids.append(os.path.basename(os.path.dirname(path)).rsplit('_', 2)[0])
    ds_list = pd.DataFrame({'group_id': group_ids, 'file_path': filepaths})
    ds_list['file_path'] = ['/' + x for x in ds_list['file_path']]
    return ds_list