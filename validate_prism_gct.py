'''
Script which checks that a GCT file has all the fields that we would expect.
May add more methods in the future.
'''
import cmapPy.pandasGEXpress.parse as parse

__author__ = "Evan Lemire"
__email__ = "elemire@broadinstitute.org"

def check_headers(filepath):

    my_gctoo = parse(filepath)

    row_metadata_fields = ['ccle_name', 'davepool_id', 'assay_plate_barcode', 'analyte_id', 'det_plate', 'det_plate_scan_time',
                           'lua', 'minipool_id', 'name', 'pool_id']
    column_metadata_fields = ['pert_id', 'pert_dose', 'pert_type', 'pert_well', 'pert_dose_unit', 'pert_time',
                              'pert_time_unit', 'pert_vehicle', 'pert_iname']
    for header in row_metadata_fields:
        if header not in my_gctoo.row_metadata_df:
            raise Exception ("{} missing from row metadata".format(header))

    for header in column_metadata_fields:
        if header not in my_gctoo.col_metadata_df:
            raise Exception ("{} field missing from column metadata".format(header))
