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


def reader_writer(input_file, output_file, function):

    gctoo = pe(input_file)

    new_gctoo = function(gctoo)

    if new_gctoo.data_df.shape[1] <= 349:
        print '{} Plate Failure With {} Failed Wells'.format(os.path.basename(os.path.dirname(input_file)),
                                                             384 - new_gctoo.data_df.shape[1])

    wgx.write_gct(new_gctoo, out_fname=output_file)
    print output_file


def card(proj_dir, plate_name):

    if not os.path.exists(os.path.join(proj_dir, 'normalize')):
        os.mkdir(os.path.join(proj_dir, 'normalize'))

    assemble_path = os.path.join(proj_dir, 'assemble', plate_name, plate_name + '_MEDIAN.gct')

    count_path = os.path.join(proj_dir, 'assemble', plate_name, plate_name + '_COUNT.gct')
    count_gctoo = pe(count_path)

    norm_path = os.path.join(proj_dir, 'normalize', plate_name, plate_name + '_NORM.gct')

    if not os.path.exists(os.path.join(proj_dir, 'normalize', plate_name)):
        os.mkdir(os.path.join(proj_dir, 'normalize', plate_name))

    #norm_path=assemble_path

    reader_writer(assemble_path, norm_path, functools.partial(norm.normalize, log=False))

    reader_writer(norm_path, norm_path,
                    functools.partial(norm.remove_low_bead_wells, count_gct=count_gctoo))

    #reader_writer(norm_path, norm_path, functools.partial(flip.flip, cell_info='/Users/elemire/Workspace/projects/PREP/cell_info.txt'))





    lvl4_card_map = {'zscorevc': [zscore.calculate_zscore, '_ZSVC.gct'],
                     'zscorepc': [functools.partial(zscore.calculate_zscore, plate_control=True), '_ZSPC.gct'],
                     'viabilitypc': [functools.partial(viability.calculate_viability, plate_control=True), '_FCPC.gct'],
                     'viabilityvc': [viability.calculate_viability, '_FCVC.gct']}

    for x in lvl4_card_map.keys():
        if not os.path.exists(os.path.join(proj_dir, x)):
            os.mkdir(os.path.join(proj_dir, x))
        if not os.path.exists(os.path.join(proj_dir, x, plate_name)):
            os.mkdir(os.path.join(proj_dir, x, plate_name))

        output_path = os.path.join(proj_dir, x, plate_name, plate_name + lvl4_card_map[x][1])

        reader_writer(input_file=norm_path, output_file=output_path, function=lvl4_card_map[x][0])


def all(proj_dir, search_pattern='*'):
    for folder in glob.glob(os.path.join(proj_dir, 'assemble', search_pattern)):
        name = os.path.basename(folder)
        card(proj_dir, name)