'''
Script which checks that a GCT file has all the fields that we would expect.
May add more methods in the future.
'''
import cmapPy.pandasGEXpress.parse as pe

__author__ = "Evan Lemire"
__email__ = "elemire@broadinstitute.org"

row_metadata_fields = ['ccle_name', 'davepool_id', 'analyte_id', 'minipool_id', 'barcode_id', 'cell_iname', 'pool_id']

column_metadata_fields = ['pert_id', 'pert_dose', 'pert_type', 'pert_well', 'pert_dose_unit', 'pert_vehicle',
                          'pert_iname']

def check_headers(filepath):

    my_gctoo = pe.parse(filepath)


    for header in row_metadata_fields:
        if header not in my_gctoo.row_metadata_df:
            raise Exception ("{} missing from row metadata".format(header))

    for header in column_metadata_fields:
        if header not in my_gctoo.col_metadata_df:
            raise Exception ("{} field missing from column metadata".format(header))
