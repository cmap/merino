import pandas as pd
import cmapPy.pandasGEXpress.parse as pe
import cmapPy.pandasGEXpress.GCToo as GCToo
import merino.build_summary.ssmd_analysis as ssmd_an
import pestle.cmap.io.gmt as gmt


def longform_table(siginfo_path, row_meta_path, pcl_connect_path, pert_connect_path, pcl_gmt_path='/Volumes/cmap_obelix/pod/custom/ASG/sas/ASG_C_V3/pcl_connectivity/ASG_PCL_xicon_n106.gmt'):

    siginfo = pd.read_table(siginfo_path)
    row_meta = pd.read_table(row_meta_path,index_col='rid')
    row_meta.index = [str(x) for x in row_meta.index]
    row_meta['minipool'] = ['P' + x.split(' ')[2] for x in row_meta['minipool_id']]
    siginfo.set_index('sig_id', inplace=True)
    siginfo['pert_well'] = siginfo['det_well']


    diff_tas_dict = {}
    for x in siginfo['cell_id'].unique():
        print ' '
        print x
        temp_col = siginfo[siginfo['cell_id'] == x]

        data = GCToo.GCToo(data_df=pd.DataFrame(temp_col['distil_tas']).T, col_metadata_df=temp_col,
                           row_metadata_df=pd.DataFrame(index=['distil_tas']))
        diff_tas_dict[x] = {}
        for y in siginfo[siginfo['pert_type'] == 'trt_cp']['pert_iname'].unique():
            ssmd = ssmd_an.get_ssmd(df=data, pos_field='pert_iname', pos_val=y)
            diff_tas_dict[x][y] = ssmd.values[0]

    diff_tas = pd.DataFrame(diff_tas_dict)

    dicto = {}
    for comp in siginfo.index:
        if siginfo.loc[comp, 'cell_id'].split('.')[1] in row_meta['pool_id'].values:
            vals = siginfo.loc[comp, ['pert_id', 'pert_idose']].values
            if vals[1] == '0.1 um':
                vals[1] = '100 nm'
            if vals[0] == 'DMSO':
                continue
            dicto[comp] = {}
            dicto[comp]['tas'] = siginfo.loc[comp, 'distil_tas']
            dicto[comp]['pool_id'] = siginfo.loc[comp, 'cell_id'].split('.')[1]
            dicto[comp]['dose'] = vals[1]
            dicto[comp]['ssmd'] = diff_tas.loc[siginfo.loc[comp, 'pert_iname'], siginfo.loc[comp, 'cell_id']]


    master_table = pd.DataFrame(dicto).T

    pcl_map = gmt.read(pcl_gmt_path)
    pcl_map_df = pd.DataFrame(pcl_map)
    pcl_map_df.set_index('id', inplace=True)
    pcl = pe(pcl_connect_path)
    pert_ps = pe(pert_connect_path)
    pcl_rankz = pcl.data_df.rank(method='min', ascending=False)
    pert_rankz = pert_ps.data_df.rank(method='min', ascending=False)

    dicto = {}
    for x in pcl.data_df.columns:
        pert_id = siginfo.loc[x, 'pert_id']
        pert_iname = siginfo.loc[x, 'pert_iname']
        dose = siginfo.loc[x, 'pert_idose']
        if pert_id not in pcl_map_df.index.values:
            continue
        lookup = pcl_map_df.loc[pert_id, 'sig']
        for thing in lookup:
            ident = x + ' ' + thing
            dicto[ident] = {}
            dicto[ident]['sig_id'] = x
            dicto[ident]['pert_id'] = pert_id
            dicto[ident]['pert_iname'] = pert_iname
            dicto[ident]['dose'] = dose
            dicto[ident]['pcl_score'] = pcl.data_df.loc[thing, x]
            dicto[ident]['pcl_rank'] = pcl_rankz.loc[thing, x]
            dicto[ident]['pool_id'] = siginfo.loc[x, 'cell_id'].split('.')[1]
    pcl_scrank_df = pd.DataFrame(dicto).T
    pcl_scrank_df['PCLs'] = [x.split(' ')[1] for x in pcl_scrank_df.index]
    pcl_scrank_df.index = [x.split(' ')[0] for x in pcl_scrank_df.index]

    dicto = {}
    for x in pert_ps.data_df.columns:
        pert_id = siginfo.loc[x, 'pert_id']
        pert_iname = siginfo.loc[x, 'pert_iname']
        dose = siginfo.loc[x, 'pert_idose']
        if pert_id not in pert_ps.data_df.index.values:
            continue
        dicto[x] = {}
        dicto[x]['pert_id'] = pert_id
        dicto[x]['pert_iname'] = pert_iname
        dicto[x]['dose'] = dose
        dicto[x]['cp_score'] = pert_ps.data_df.loc[pert_id, x]
        dicto[x]['cp_rank'] = pert_rankz.loc[pert_id, x]
        dicto[x]['pool_id'] = siginfo.loc[x, 'cell_id'].split('.')[1]
    pert_scrank_df = pd.DataFrame(dicto).T

    master_table = master_table.join(pert_scrank_df[['cp_rank', 'cp_score']])
    master_table = master_table.join(pcl_scrank_df[['pcl_rank', 'pcl_score', 'PCLs']])

    return master_table