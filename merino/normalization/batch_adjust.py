"""
batch_adjust.py

Batch adjustment methods for Prism data.

"""

import os
import sys
import glob
import argparse
import logging
import multiprocessing as mp
import numpy as np

import merino.normalization.combat as combat
import merino.setup_logger as setup_logger
import cmapPy.pandasGEXpress.parse as pe
import cmapPy.pandasGEXpress.write_gct as wg
import cmapPy.pandasGEXpress.concat as cg
import cmapPy.set_io.grp as grp
# handle cmapPy version differences for subsetting gcts
try:
    from cmapPy.pandasGEXpress.subset_gctoo import subset_gctoo as gct_slice
except ImportError:
    import cmapPy.pandasGEXpress.slice_gctoo as gct_slice


LOGGER = logging.getLogger(setup_logger.LOGGER_NAME)

def data_splitter(all_ds, col_group, batch_field, use_col_group_as_batch):
    """ Split dataset into replicate chunks and reshape for ComBat
    """
    col_groups = all_ds.col_metadata_df.groupby(col_group).groups
    LOGGER.info(len(col_groups))
    batches = all_ds.row_metadata_df[batch_field]
    nbatch = len(batches.unique())
    LOGGER.info('Splitting dataset by %s into %d groups' % (col_group, len(col_groups)))
    LOGGER.info('Batch field %s has %d levels' % (batch_field, nbatch))
    
    for _, key in enumerate(sorted(col_groups)):
        LOGGER.info('key_{}'.format(key))
        this_gp = all_ds.data_df[col_groups[key]].copy()
        this_gp[batch_field] = batches
        # save the row-ids
        this_gp['row_id'] = this_gp.index
        # unpivot the matrix to a column vector
        df_long = this_gp.melt(id_vars=['row_id', batch_field])
        # append a dummy feature to allow ComBat to execute
        df_long['dummy'] = np.random.randn(df_long.shape[0], 1)
        if use_col_group_as_batch is True:
            df_long['batch'] = df_long[batch_field] + ':' + df_long['cid']
        else:
            df_long['batch'] = df_long[batch_field]
        yield df_long.transpose()

def combat_worker(df_long):
    """ Apply ComBat to a chunk of data and return a reshaped
    dataframe adjusted values """

    # handle missing values

    to_use = df_long.loc['value'].notnull().values
    vals = df_long.loc[['value', 'dummy'], to_use]
    batch = df_long.loc[['batch'], to_use].squeeze()

    adj = combat.combat(vals, batch)
    df_long.loc[['value'], to_use] = adj.loc[['value']]
    # reshape results
    return df_long.T.pivot(index='row_id', columns='cid', values='value')

def combat_by_group(gct_list, col_group='pert_well', batch_field='pool_id',
                    use_col_group_as_batch=True):
    """ Applies ComBat batch adjustment algorithm to a list of input
    GCT objects grouped by the specified column grouping.  The method
    first concatenates the input list of GCT objects, splits columns
    by specified col_group and applies ComBat on the unwrapped matrix
    values using the either batch_field + replicate id if
    use_col_group_as batch is True or batch_field alone otherwise.

    Args:
        gct_list: list of GCT objects
        col_group: Column metadata field(s) to group columns by
        batch_field: Row metadata field used to specify the batches
        use_col_group_as_batch: if True the the column identity is
        appended to the batch vector such that the number of batches
        for a group = Number of unique batch_field entries * Number of
        columns in the group

    Returns:
        all_ds: Concatenated Combat adjusted values
        adj_list: list of gct objects of the adjusted values of the
        subsets of all_ds that match gct_list
    """
    # concatenate replicate datasets by column
    LOGGER.info("now running ComBat batch adjustment")
    fields_to_remove = [x for x in gct_list[0].row_metadata_df.columns if
                        x in ['det_plate', 'det_plate_scan_time', 'assay_plate_barcode']]
    all_ds = cg.hstack(gct_list, remove_all_metadata_fields=False,error_report_file=None, fields_to_remove=fields_to_remove)


    # column groups
    #col_groups = all_ds.col_metadata_df.groupby(col_group).groups
    #row_groups = all_ds.row_metadata_df[row_group]

    pool = mp.Pool(processes=mp.cpu_count())

    chunks = data_splitter(all_ds, col_group, batch_field, use_col_group_as_batch)
    LOGGER.info('Here 1')

    adjusted_data = pool.map(combat_worker, chunks)
    LOGGER.info('Here 2')

    for res in adjusted_data:
        all_ds.data_df[res.columns] = res
    LOGGER.info('Here 3')

    combat_adjusted_gcts = []
    for _, input_ds in enumerate(gct_list):
        this_ds = gct_slice(all_ds, rid=input_ds.data_df.index.tolist(), cid=input_ds.data_df.columns.tolist())
        this_ds.src = input_ds.src
        this_ds.data_df = this_ds.data_df.astype(float)
        combat_adjusted_gcts.append(this_ds)
    LOGGER.info('Here 4')

    return all_ds, combat_adjusted_gcts


def build_parser():
    """Build argument parser."""

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required args
    parser.add_argument("--input_glob", '-g',
                        help="Glob pattern to input GCT files")
    parser.add_argument("--input_list", '-i',
                        help="Text file listing filpaths to input GCT files")
    parser.add_argument("--col_group", nargs="+",
                        help="Column metadata field(s) to group columns by")
    parser.add_argument("--batch_field",
                        help="Row metadata field used to specify the batches")
    parser.add_argument("--use_col_group_as_batch", action='store_true',
                        default=False,
                        help="Use column identity as a batch if True")
    parser.add_argument("--verbose", "-v", action="store_true", default=False,
                        help="whether to increase the # of messages reported")
    return parser

def load_data(gct_files):
    """ Read a list of GCT files and returns a list
    """
    gct_list = []
    for gct_path in gct_files:
        LOGGER.info('Reading {}'.format(gct_path))
        gct = pe.parse(gct_path)
        gct_list.append(gct)
    return gct_list

def save_data(adj_ds, adj_list):
    """ Write batch-adjusted data to files
    """
    wg.write(adj_ds, 'batch_adjusted_values.gct')
    for ctr, this_ds in enumerate(adj_list):
        if this_ds.src is not None:
            out_file = '{}.COMBAT.gct'.format(os.path.splitext(os.path.basename(this_ds.src))[0])
        else:
            out_file = 'batch_adjusted_values_X{}.gct'.format(ctr)
        wg.write(this_ds, out_file)

def main():
    """Main function
    """

    # Get args
    parser = build_parser()
    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    if args.input_glob:
        gct_files = glob.glob(args.input_glob)
    elif args.input_list:
        gct_files = grp.read(args.input_list)
    else:
        raise Exception('input_glob or input_list must be specified')

    gct_list = load_data(gct_files)

    setup_logger.setup(verbose=args.verbose)
    adj_ds, adj_list = combat_by_group(gct_list,
                                       col_group=args.col_group,
                                       batch_field=args.batch_field,
                                       use_col_group_as_batch=args.use_col_group_as_batch)
    save_data(adj_ds, adj_list)

if __name__ == "__main__":
    main()
