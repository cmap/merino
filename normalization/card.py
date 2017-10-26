import normalize_with_prism_invariant as norm
import viability_normalization as viability
import zscore
import os
import glob
import modz
import prism_pipeline.cut_to_l2 as cut
import cmapPy.pandasGEXpress.write_gct as wgx
import cmapPy.pandasGEXpress.parse as pe
import functools
import prism_pipeline.flip_cells as flip
import shear
import pandas as pd


def reader_writer(input_file, output_file, function, check_size=False):
    plate_failure = False
    # Read in input file
    gctoo = pe(input_file)
    # Call normalizing function on gctoo
    new_gctoo = function(gctoo)

    # If told to, check size of new_gctoo and flag if too small
    if new_gctoo.data_df.shape[1] <= 349 and check_size==True:
        print '{} Plate Failure With {} Failed Wells'.format(os.path.basename(os.path.dirname(input_file)),
                                                             384 - new_gctoo.data_df.shape[1])
        plate_failure = True

    # write out new gctoo
    wgx.write_gct(new_gctoo, out_fname=output_file)
    print output_file

    return plate_failure


def card(proj_dir, plate_name, log_tf, bad_wells=[]):

    if not os.path.exists(os.path.join(proj_dir, 'normalize')):
        os.mkdir(os.path.join(proj_dir, 'normalize'))
    # Get path to raw mfi
    assemble_path = os.path.join(proj_dir, 'assemble', plate_name, plate_name + '_MEDIAN.gct')
    # Get path to beadcount values
    count_path = os.path.join(proj_dir, 'assemble', plate_name, plate_name + '_COUNT.gct')
    # Get path to LEVEL3 norm values
    norm_path = os.path.join(proj_dir, 'normalize', plate_name, plate_name + '_NORM.gct')
    # Set plate_failure variable to false
    plate_failure=False

    if not os.path.exists(os.path.join(proj_dir, 'normalize', plate_name)):
        # Create norm folder if it doesn't exist already
        os.mkdir(os.path.join(proj_dir, 'normalize', plate_name))
        # Create norm file
        reader_writer(assemble_path, norm_path, functools.partial(norm.normalize, log=log_tf))
        # Read in count file
        count_gctoo = pe(count_path)
        # Remove low bead count wells and check GCT size, if too many wells have been stripped it will qualify as a failure
        plate_failure = reader_writer(norm_path, norm_path, functools.partial(norm.remove_low_bead_wells, count_gct=count_gctoo), check_size=True)
        # Shear predetermined bad wells (if any exist)
        reader_writer(norm_path, norm_path, functools.partial(shear.shear, bad_wells=bad_wells))


    # Map denoting each type of LEVEL4 data, its folder name, the function to create it, and the file ending.
    lvl4_card_map = {'zscorevc': [zscore.calculate_zscore, '_ZSVC.gct'],
                     'zscorepc': [functools.partial(zscore.calculate_zscore, plate_control=True), '_ZSPC.gct'],
                     'viabilitypc': [functools.partial(viability.calculate_viability, plate_control=True), '_FCPC.gct'],
                     'viabilityvc': [viability.calculate_viability, '_FCVC.gct']}

    # Loop through this map to output all level 4 data
    for x in lvl4_card_map.keys():
        if not os.path.exists(os.path.join(proj_dir, x)):
            os.mkdir(os.path.join(proj_dir, x))
        if not os.path.exists(os.path.join(proj_dir, x, plate_name)):
            os.mkdir(os.path.join(proj_dir, x, plate_name))

            output_path = os.path.join(proj_dir, x, plate_name, plate_name + lvl4_card_map[x][1])

            reader_writer(input_file=norm_path, output_file=output_path, function=lvl4_card_map[x][0])

    # Return status of plate failure

    return plate_failure


def all(proj_dir, search_pattern='*', log_tf=True, bad_wells=[]):
    failure_list = []
    for folder in glob.glob(os.path.join(proj_dir, 'assemble', search_pattern)):
        name = os.path.basename(folder)
        plate_failure = card(proj_dir, name, log_tf = log_tf, bad_wells=bad_wells)
        if plate_failure == True:
            failure_list.append(name)

    print failure_list
    pd.Series(failure_list).to_csv(os.path.join(proj_dir, 'failed_plates.txt'), sep='\t')

