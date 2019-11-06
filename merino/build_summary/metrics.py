
'''
This module contains functions for computing various signature-level metrics
for a givem CoP build. The objective is to compute a siginfo-life table with
all our standard as well as other signature metrics
'''


from scipy.stats import binom, binom_test
from sklearn.metrics import jaccard_similarity_score
import pandas as pd
import numpy as np
from cmapPy.pandasGEXpress import parse, write_gctx, GCToo

import logging
logger = logging.getLogger()

dict_metric_labels = {
    'ss_ltn2' : 'Signature strength',
    'cc_q75' : 'Replicate correlation',
    'overlap_pvalue_minus_log10' : 'Replicate overlap significance',
    'killed_jaccard' : 'Replicate overlap (Jaccard)'
}


dict_pert_type_labels = {
    'ctl_untrt' : 'Untreated Well',
    'ctl_vehicle' : 'DMSO',
    'trt_cp' : 'Test compound',
    'trt_poscon' : 'Bortezomib',
    'trt_poscon.es' : 'Selective PosCon'
}

# logger.setLevel(logging.INFO)

# Used in CoP nominations (PRIME, GEX active, subtle) for JSC July 2019
dict_category_cop = {
    'ctl_vehicle' : ["DMSO", "Vehicle control - DMSO", 0],
    "trt_poscon" : ["Bortezomib", "Positive control - Bortezomib", 1],
    "trt_poscon.es" : ["Asgard", "Reference controls - Asgard", 2],
    None : [None, None, 0],
    (0,0) : ["Inert (0%)", "Test subset - inert (0%)", 3],
    (1,1) : ["Low (1%)", "Test subset - Low (1%)",4],
    (2,5) : ["Low (2-5%)", "Test subset - Low (2-5%)", 5],
    (6,10) : ["Mid (6-10%)", "Test subset - Mid (6-10%)", 6],
    (11,20) : ["Mid (11-20%)", "Test subset - Mid (11-20%)", 7],
    (21,35) : ["High (21-35%)", "Test subset - High (21-35%)", 8],
    (36,99) : ["High (>35%)", "Test subset - High (>35%)", 9],
    (100,100) : ["All Killed (100%)", "Test subset - All (100%)", 10]}



# To match buckets used in CoP nominations (PRIME, GEX active, subtle) for JSC July 2019
dict_category_PR500= {
    'ctl_vehicle' : ["DMSO", "Vehicle control - DMSO", 0],
    "trt_poscon" : ["Bortezomib", "Positive control - Bortezomib", 1],
    "trt_poscon.es" : ["Asgard","Reference controls - Asgard", 1.5],
    'trt_ref' : ["Reference compounds", "Reference compounds", 2],
    None : [None, None, 0],
    (0,0) : ["Inert: 0 (0%)", "Test subset - Inert (0%)", 3],
    (1,4) : ["Low: 1-4 (1%)", "Test subset - Low (1%)", 4],
    (5,24) : ["Low: 5-24 (2-5%)", "Test subset - Low (2-5%)", 5],
    (25,48) : ["Mid: 25-48 (6-10%)", "Test subset - Mid (6-10%)", 6],
    (49,97) : ["Mid: 49-97 (11-20%)", "Test subset - Mid (11-20%)", 7],
    (98,171) : ["High: 98-171 (21-35%)", "Test subset - High (21-35%)", 8],
    (172,484) : ["High: 172-484 (>35%)", "Test subset - High (>35%)", 9],
    (485,489) : ["All Killed: 485-489 (100%)", "Test subset - All (100%)", 10]
}


# Here total is 578 cell lines
dict_category_DP7_8 = {
    'ctl_vehicle' : ["DMSO", "Vehicle control - DMSO", 0],
    "trt_poscon" : ["Bortezomib", "Positive control - Bortezomib", 1],
    "trt_poscon.es" : ["Asgard","Reference controls - Asgard", 2],
    None : None,
    (0,0) : ["Inert: 0 (0%)", "Test subset - Inert (0%)", 3],
    (1,29) : ["Low: 1-29 (1-5%)", "Test subset - Low (1-5%)", 4],
    (30,115) : ["Mid: 30-115 (6-20%)", "Test subset - Mid (6-20%)", 5],
    (116,289) : ["High: 116-289 (21-50%)", "Test subset - High (21-50%)", 6],
    (290,577) : ["Non-Selective (>50%)", "Test subset - Non-Selective (50-99%)", 7],
    (578,578) : ["All Killed: 485-489 (100%)", "Test subset - All (100%)", 9]
}

meta_dict = {'PR500': dict_category_PR500, 'COP': dict_category_cop, 'DP78': dict_category_DP7_8}

def assign_ss_buckets(x, dict_key='PR500'):
    '''
    Assign a bucket to a ss (number of cell lines killed) value . According to
    the new scheme Aravind asked for for PBRANT cycle 1
    '''
    d = meta_dict[dict_key]
    for pair, _ in d.iteritems():
        if type(pair) != tuple:
            continue
        i, j = pair
        if (x >=i) and (x <=j):
            return (i,j)



dict_category_order = {
    'ctl_vehicle' : 0,
    "trt_poscon" : 1,
    "trt_poscon.es" : 2,
    None : 0,
    (0,0) : 3,
    (1,1) : 4,
    (2,5) : 5,
    (6,10) : 6,
    (11,20) : 7,
    (21,50) : 8,
    (51,99) : 9,
    (100,100) : 10
}





dict_category_order_PR500_v2 = {
    'ctl_vehicle' : 0,
    "trt_poscon" : 1,
    "trt_poscon.es" : 2,
    None : 0,
    (0,0) : 3,
    (1,5) : 4,
    (6,25) : 5,
    (26,50) : 6,
    (51,100) : 7,
    (101,250) : 8,
    (251,488) : 9,
    (489,489) : 10
}








# for backward compatibility
dict_category_order_PR500 = dict_category_order_PR500_v2

dict_min_bucket_sample_size = {
    '(1, 5)' : 150,
    '(6, 10)' : 100,
    '(11, 15)' : 50,
    '(16, 20)' : 50,
}


def assign_category(row, bucket_column="ss_bucket_v2"):
    '''
    Given a row of the sigstats table, compute category labels.
    If a test compound, simply return the ss_bucket_v2. Otherwise, 
    return pert_type
    '''
    if row['pert_type'] in ['ctl_vehicle', 'trt_poscon', 'trt_poscon.es', 'trt_ref']:
        return row['pert_type']
    elif row['pert_type'] == "trt_cp":
        return row[bucket_column]
    else:
        return None




def assign_ss_buckets_CoP_v2(x):
    '''
    Assign a bucket to a ss (number of cell lines killed) value . According to
    the new scheme Aravind asked fo on May 16.
    '''
    if (x == 0):
        return (0,0)
    if (x == 1):
        return (1,1)
    if (x == 100):
        return (100,100)
    
    if (x>=2 and x<=5):
        return (2,5)
    if (x>=6 and x<=10):
        return (6,10)
    if (x>=11 and x<=20):
        return (11, 20)
    if (x>=21 and x<=50):
        return (21, 50)
    if (x>=51 and x<=99):
        return (51, 99)
    

def assign_ss_buckets_CoP_v3(x):
    '''
    Assign a bucket to a ss (number of cell lines killed) value . According to
    the new scheme Aravind asked fo on June 11 for PRIME definitions.
    '''
    if (x == 0):
        return (0,0)
    if (x == 1):
        return (1,1)
    if (x == 100):
        return (100,100)
    
    if (x>=2 and x<=5):
        return (2,5)
    if (x>=6 and x<=10):
        return (6,10)
    if (x>=11 and x<=20):
        return (11, 20)
    if (x>=21 and x<=35):
        return (21, 35)
    
    if (x>=36 and x<=99):
        return (36, 99)  

    
def assign_ss_buckets_CoP_v5(x):
    '''
    Assign a bucket to a ss (number of cell lines killed) value . According to
    the new scheme Aravind asked for for PBRANT cycle 1
    '''
    if (x == 0):
        return (0,0)
    if (x == 1):
        return (1,1)
   
    if (x>=2 and x<=20):
        return (2,20)
    if (x>=21 and x<=50):
        return (21,50)
    if (x>=51):
        return (51, 100)
    
    
    
def assign_ss_buckets_PR500_v2(x):
    '''
    Assign a bucket to a ss (number of cell lines killed) value . According to
    the new scheme Aravind asked fo on May 16.
    '''
    if (x == 0):
        return (0,0)

    if (x == 489):
        return (489, 489)
    
    if (x>=1 and x<=5):
        return (1,5)
    if (x>=6 and x<=25):
        return (6,25)
    if (x>=26 and x<=50):
        return (26, 50)
    if (x>=51 and x<=100):
        return (51, 100)
    if (x>=101 and x<=250):
        return (101, 250)
    if (x>=251 and x<=488):
        return (251, 488)
    
        
def assign_ss_buckets_PR500_v3(x):
    '''
    Assign a bucket to a ss (number of cell lines killed) value . According to
    the new scheme Rajiv asked for June 7.
    '''
    
    if (x == 0):
        return (0,0)

    if (x == 489):
        return (489, 489)
    
    if (x>=1 and x<=25):
        return (1,25)
    if (x>=26 and x<=100):
        return (26,100)
    if (x>=101 and x<=250):
        return (101, 250)
    if (x>=251 and x<=488):
        return (251, 488)

def assign_ss_buckets_PR500_v5(x):
    '''
    Assign a bucket to a ss (number of cell lines killed) value . According to
    the new scheme Aravind asked for for PBRANT cycle 1
    '''
    if (x == 0):
        return (0,0)
    if (x>=1 and x<=5):
        return (1,5)
    if (x>=6 and x<=100):
        return (6,100)
    if (x>=101 and x<=250):
        return (101, 250)
    if (x>=251):
        return (251, 489)


    
    
def assign_ss_buckets_PR500_PRIMARY_v1(x):
    '''
    Assign a bucket to a ss (number of cell lines killed) value . 
    PRISM PRIMARY has a total of 578 cell lines.
    '''
    if (x == 0):
        return (0,0)

    if (x == 578):
        return (578, 578)
    
    if (x>=1 and x<=29):
        return (1,29)
    if (x>=30 and x<=115):
        return (30,115)
    if (x>=116 and x<=289):
        return (116, 289)
    if (x>=290 and x<=577):
        return (290, 577)

    
    

    
    
def assign_ss_buckets_CoP(x, x_max=100):
    '''
    Assign a bucket to a ss (number of cell lines killed) value 
    '''
    if (x == 0):
        return (0,0)
    if (x>=1 and x<=4):
        return (1,4)
    if (x>=5 and x<=20):
        return (5,20)
    if (x>=21 and x<=50):
        return (21, 50)
    if (x>=51 and x<=100):
        return (51, 100)
    
dict_ss_bucket_labels_CoP = {
    0:'inert (0)',
    4:'borderline (1-4)',
    20:'low (5-20)',
    50:'mid (21-50)',
    100:'high (51-100)'
}

def assign_tas_buckets_CoP(x, x_max=1):
    '''
    Assign a bucket to a ss (number of cell lines killed) value 
    '''
    if (x == 0):
        return (0,0)
    if (x>0 and x<=0.1):
        return (0,0.1)
    if (x>0.1 and x<=0.2):
        return (0.1,0.2)
    if (x>0.2 and x<=0.3):
        return (0.2, 0.3)
    if (x>0.3 and x<=1):
        return (0.3, 1)
    
dict_tas_bucket_labels_CoP = {
    0:'inert (0)',
    0.1:'borderline (0-0.1)',
    0.2:'low (0.1-0.2)',
    0.3:'mid (0.2-0.3)',
    1:'high (0.3-1)'
}

def assign_buckets_int(x, x_max=100, width=5):
    if x == x_max:
        x = x_max - 1
    if x == 0:
        return (0, 0)
    n = int(np.ceil(float(x)/width) * width)
    return (n - width + 1, n)
    return n


def count_common_killed_cells(row, df_killed):
    c1 = row[0]
    c2 = row[1]
    
    try:
        x, y = df_killed.loc[:,c1], df_killed.loc[:,c2]
    except:
        return np.nan
    return (x * y).sum()


def compute_jaccard(row, df_killed=None):
    c1 = row[0]
    c2 = row[1]
    
    try:
        x, y = df_killed.loc[:,c1], df_killed.loc[:,c2]
    except:
        return np.nan
    n_union = (x + y)
    n_union = len(n_union[n_union >=1])
    if n_union == 0:
        return 0
    n_common = (x * y).sum()
    return float(n_common) / n_union


def compute_pvalue(row, df_killed=None, N=100):
    c1 = row[0]
    c2 = row[1]
    
    assert df_killed is not None
    
    try:
        x, y = df_killed.loc[:,c1], df_killed.loc[:,c2]
    except:
        return np.nan
    assert len(x) == len(y) == N
    n1, n2 = x.sum(), y.sum()
    m = (x * y).sum()
    p12 = 1. *  (n1 * n2) / N/N
    pvalue = binom_test(m, N, p12, alternative="greater")
    return pvalue


def flag_high_cc_zero_replicate_concordance(row):
    cc = row['cc_q75']
    rc = row['num_common_frac']
    return 1 if (cc > 0.2 )  and (rc == 0) else 0


def flag_low_cc_high_replicate_concordance(row):
    cc = row['cc_q75']
    rc = row['num_common_frac']
    return 1 if (cc <= 0.2 )  and (rc >=0.75) else 0


def flag_low_ss_high_replicate_concordance(row):
    ss = row['ss_ltn2']
    rc = row['num_common_frac']
    return 1 if (ss <= 20 )  and (rc >=0.75) else 0


def flag_low_ss_low_replicate_concordance(row):
    ss = row['ss_ltn2']
    rc = row['num_common_frac']
    return 1 if (ss <= 20 )  and (rc ==0) else 0



'''
def compute_category_labels(df, ss_column="ss_ltn2", which="CoP", 
                            fcn_assign_bucket=assign_ss_buckets_CoP_v2,
                            dict_category_labels=dict_category_labels,
                            dict_category_labels_abridged=dict_category_labels_abridged,
                            dict_category_order=dict_category_order
                           ):
    
    Assign labels for the 12 signature categories. Two columns will be added:
    'category_label' and 'category_label_abrdiged'
    
    bucket_column = "ss_bucket"
    
    assert which in ['CoP', "PR500"]
    
    logger.info('Computing comprehensive category labels, spanning pert_types and v2 SS buckets.')
    
    logger.info('Computing ss_bucket for all signatures (not just test compounds)')
    try:
        df[bucket_column] = df['ss_ltn2'].apply(fcn_assign_bucket)
    except Exception as e:
        logger.error("Unable to assign buckets : {}".format(str(e)))
        
    
    logger.info('Computing the "category" field')
    try:
        df['category'] = df.apply(assign_category, axis=1, bucket_column=bucket_column)
    except Exception as e:
        logger.error("Unable to assign categories : {}".format(str(e)))
        
    logger.info('Computing full category labels')
    try:
        df['category_label'] = df['category'].apply(lambda x : dict_category_labels[x])
    except Exception as e:
        logger.error("Unable to assign labels : {}".format(str(e)))
        
    logger.info('Computing abridged category labels')
    try:
        df['category_label_abridged'] = df['category'].apply(lambda x : dict_category_labels_abridged[x])
    except Exception as e:
        logger.error("Unable to assign abridged labels : {}".format(str(e)))
        
    logger.info('Assiging display order to each category')
    try:
        df['category_order'] = df['category'].apply(lambda x: dict_category_order[x])
    except Exception as e:
        logger.error("Unable to assign category order : {}".format(str(e)))
        
    
    return df
    '''

def compute_sig_metric_flags(df):
    '''
    Given all the coputed sigmetrcs, compute all flags, bucket labels and category labels
    '''
    required_fields = ['ss_ltn2', 'cc_q75', 'killed_jaccard', 'overlap_pvalue_minus_log10']
    for c in required_fields:
        assert c in df.columns, "{} not in df.columns".format(c)
        
    # Compute the ss buckets used in nomination
    logger.info('Computing the SS buckets used in nomination')
    df['ss_bucket'] = df['ss_ltn2'].apply(lambda x: assign_buckets_int(x, x_max=100, width=5)[1])
    df['ss_bucket_label'] = df['ss_ltn2'].apply(lambda x:str(assign_buckets_int(x, x_max=100, width=5)) )

    logger.info('Computing the coarse SS buckets used for presentation purposes')
    df['ss_bucket_coarse'] = df['ss_ltn2'].apply(lambda x: assign_ss_buckets_CoP(x, x_max=100)[1])
    df['ss_bucket_coarse_label'] = df['ss_bucket_coarse'].apply(lambda x: dict_ss_bucket_labels_CoP[x] )
    
    # Compute labels for the comprehensive category list
    df = compute_category_labels(df)

    df['is_high_cc'] = (df['cc_q75'] > 0.2) + 0
    df['is_high_replicate_jaccard'] = (df['killed_jaccard'] > 0.3) + 0
    df['is_high_replicate_overlap_significance'] = (df['overlap_pvalue_minus_log10'] > 2) + 0
    return df






def compute_build_sig_stats(gct_level4, siginfo, assay="COP", killing_threshold=-2, num_reps_min=2):
    assert assay in ['COP', 'PR500']
    dict_replicate_field = {
        'COP' : ('distil_id', '|'),
        'PR500' : ('profile_ids', ",")
    }
    
    replicate_field_name, replicate_field_delimiter = dict_replicate_field[assay]
    # Old version
    # siginfo_2rep = siginfo[(siginfo['nprofile'] == num_reps) & (siginfo['distil_nsample'] == num_reps)]
    
    # New version, May 16. Only care about nprofile.
    siginfo_2rep = siginfo[(siginfo['nprofile'] >= num_reps_min)]
    if 'distil_nsample' in siginfo.columns:
        siginfo_2rep = siginfo_2rep[(siginfo_2rep['distil_nsample'] >= 2)]
#     siginfo_2rep = siginfo[(siginfo['nprofile'] >= num_reps_min)]
    logger.info('Number of all signatures in siginfo: {}'.format(siginfo.shape[0]))
    logger.info('Number of all signatures in siginfo with at least 2 GEX replicates and at least {} KJ100 replicates: {}'.format(num_reps_min, siginfo_2rep.shape[0]))

    
    
    
    df_reps = siginfo_2rep[replicate_field_name].apply(lambda s:pd.Series(s.split(replicate_field_delimiter)))
    df_reps.index = siginfo_2rep.index

    df_reps.loc[:, 'pert_type'] = siginfo_2rep['pert_type']
    df_reps.loc[:, 'sig_id'] = siginfo_2rep['sig_id']
    df_reps.loc[:, 'pert_id'] = siginfo_2rep['pert_id']
    

    fields = ['distil_cc_q75', 'cc_q75', 'ss_ltn2', 'cis_ltn2', 'pert_dose']
    for field in fields:
        try:
            df_reps.loc[:, field] = siginfo_2rep[field]
        except:
            print "field '{}' not in columns".format(field)

    
    logger.info('Computing the number of cell lines killed for all signatures...')
    df_killed = (gct_level4.data_df< killing_threshold) + 0
    logger.info('Done.')

    num_killed = df_killed.sum(0)

    
    logger.info('Computing all metrics...')
    logger.warning('Only processing the first 2 replictes.')
    df_reps.set_index(0, inplace=True)
        
    df_reps['n_killed_0'] = num_killed.loc[df_reps.index]
    df_reps = df_reps.reset_index().set_index(1)

    df_reps['n_killed_1'] = num_killed.loc[df_reps.index]
    df_reps.reset_index(inplace=True)

    logger.info('Computing replicate consistency stats')
    df_reps['num_common'] = df_reps.apply(count_common_killed_cells, axis=1, df_killed=df_killed)
#     df_reps['num_common_frac'] = 1.0 * df_reps['num_common'] / df_reps['ss_ltn2']
    df_reps['killed_jaccard'] = df_reps.apply(compute_jaccard, axis=1, df_killed=df_killed)
    df_reps['overlap_pvalue'] = df_reps.apply(compute_pvalue, axis=1, df_killed=df_killed)
    df_reps['overlap_pvalue_minus_log10'] = - np.log10(df_reps['overlap_pvalue'])

    # Range flags
#     df_reps['high_cc_low_rc'] = df_reps.apply(flag_high_cc_zero_replicate_concordance, axis=1)
#     df_reps['low_cc_high_rc'] = df_reps.apply(flag_low_cc_high_replicate_concordance, axis=1)
#     df_reps['low_ss_high_rc'] = df_reps.apply(flag_low_ss_high_replicate_concordance, axis=1)
#     df_reps['low_ss_low_rc'] = df_reps.apply(flag_low_ss_low_replicate_concordance, axis=1)
    logger.info('Done.')
    
    logger.info("Computing sig metric flags")
    df_reps = compute_sig_metric_flags(df_reps)
    return df_reps




def compute_and_export_stats_CoP_CORE():
    filename = "/cmap/projects/biosensor/datasets/CORE.D1/sample_info.txt"
    siginfo = pd.read_csv(filename, sep="\t")
    siginfo['nprofile'] = 2
    siginfo['distil_nsample'] = 2
    logger.warning("siginfo missing nprofile and distil_nsample, setting both to 2 for all signatures")

    filename = "/cmap/projects/biosensor/datasets/CORE.D1/CORE.D1_48h_Level4.ZSPC_n10752x872.gctx"
    gct_level4 = parse.parse(filename)

    df_reps = compute_build_sig_stats(gct_level4, siginfo, killing_threshold=-2)
    df_reps.to_csv('./datasets/CORE_viability_replicate_agreement_stats.csv', sep="\t")
    

def compute_and_export_stats_CBRANT_ALL():
    filename = "/cmap/obelix/pod/custom/CBRANT/build/CBRANT_ALL/CBRANT_ALL_VIABILITY_LEVEL4_ZSPC.COMBAT_n271674x100.gctx"
    gct_level4 = parse.parse(filename)
    filename = "/cmap/obelix/pod/custom/CBRANT/build/CBRANT_ALL/CBRANT_ALL_sig_metrics.txt"
    siginfo = pd.read_csv(filename, sep="\t")
    df_reps = compute_build_sig_stats(gct_level4, siginfo, killing_threshold=-2)
    outfile = './datasets/CBRANT_C_viability_replicate_agreement_stats.csv'
    logger.info('Exported CBRANT_C sigstats to {}'.format(outfile))
    df_reps.to_csv(outfile, sep="\t")
    


def compute_and_export_stats_CBRANT_C():
    filename = "/cmap/obelix/pod/custom/CBRANT/build/CBRANT_C_VIABILITY/build/CBRANT_C_VIABILITY_LEVEL4_ZSPC.COMBAT_n231680x100.gctx"
    gct_level4 = parse.parse(filename)
    filename = "/cmap/obelix/pod/custom/CBRANT/build/CBRANT_C_MERGED/build/CBRANT_C_sig_metrics.txt"
    siginfo = pd.read_csv(filename, sep="\t")
    df_reps = compute_build_sig_stats(gct_level4, siginfo, killing_threshold=-2)
    outfile = './datasets/CBRANT_C_viability_replicate_agreement_stats.csv'
    logger.info('Exported CBRANT_C sigstats to {}'.format(outfile))
    df_reps.to_csv(outfile, sep="\t")
    

def compute_and_export_stats_CBRANT_B():
    '''
    WARNING: This doesn't really make sens since CBRANT_B was 
    profiled in triplicate and none of the new metrics make sense
    with 3 replicates. Until I figure out what to do with this, we
    can't really apply the same metrics to CBRANT_B
    '''
    filename = "/cmap/obelix/pod/custom/CBRANT/build/CBRANT-B_VIABILITY/CBRANT_B_VIABILITY_LEVEL4_ZSPC.COMBAT_n39994x100.gctx"
    gct_level4 = parse.parse(filename)
    filename = "/cmap/obelix/pod/custom/CBRANT/build/CBRANT_B/merged_build/CBRANT_B_sig_metrics.txt"
    siginfo = pd.read_csv(filename, sep="\t")
    logger.warning(
        'For CBRANT_B, the distil_ids in the merged siginfo do not match the cids in the viability build. Modifying them in place, hoping they match up'
    )
    siginfo['distil_id'] = siginfo['distil_id'].apply(lambda s:s.replace("_C3", ""))
    df_reps = compute_build_sig_stats(gct_level4, siginfo, killing_threshold=-2, num_reps_min=2)
    outfile = './datasets/CBRANT_B_viability_replicate_agreement_stats.csv'
    logger.info('Exported CBRANT_B sigstats to {}'.format(outfile))
    df_reps.to_csv(outfile, sep="\t")

def combine_and_export_stats_CBRANT_ALL():
    '''
    Load complete sigstats for both
    CBRANT_B and CBRANT_C, concatenate them and export 
    to a file
    '''
    
    filename = './datasets/CBRANT_B_viability_replicate_agreement_stats.csv'
    sigstats_B = pd.read_csv(filename, sep="\t", index_col=0)
    
    filename = './datasets/CBRANT_C_viability_replicate_agreement_stats.csv'
    sigstats_C = pd.read_csv(filename, sep="\t", index_col=0)
    
    logger.info('Concatenating sigstat files for CBRANT_B and CBRANT_C')
    sigstats_all = pd.concat(
        [sigstats_B, sigstats_C], axis=0, sort=True
    )
    outfile = './datasets/CBRANT_all_viability_replicate_agreement_stats.csv'
    sigstats_all.to_csv(outfile, sep="\t")
    logger.info('Exported merged CBRANT sigstats to {}'.format(outfile))



def plot_metrics_by_pert_type( siginfo_filename, outfile):
    pal = sns.color_palette('Set1')
    sigstats = pd.read_csv(siginfo_filename, sep="\t",index_col=0)
    list_metrics = ['ss_ltn2', 'cc_q75', 'overlap_pvalue_minus_log10', 'killed_jaccard']
    list_metrics = [x for x in list_metrics if x in sigstats.columns]
    print list_metrics

    def break_lines(s):
        x = s.split(" ")
        if len(x) == 1:
            return s
        return "{}\n{}".format(x[0], " ".join(x[1:]))

    list_metric_labels = [break_lines(dict_metric_labels[x]) for x in list_metrics]
    list_pert_types = ['ctl_untrt', 'ctl_vehicle', 'trt_cp', 'trt_poscon', 'trt_poscon.es']
    list_pert_type_labels = [break_lines(dict_pert_type_labels[x]) for x in list_pert_types]

    list_xlims = [(0,100), (0,1), (0,10), (0,1)]
    n1, n2 = len(list_pert_types), len(list_metrics)
    counter = 1

    plt.figure(figsize=(2 * len(list_metrics),5))

    def remove_ticks():
        plt.tick_params(
            reset=True,
            axis='both',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            right='off',         # ticks along the top edge are off
            top='off',      # ticks along the bottom edge are off
            left="off",
            labelleft='off',
            labelbottom="off",
            color="#555555"
            
        #     labelbottom=False
        ) 
    font="Latowe"
    fontcolor="#333333"
    fontweight=None
    # plt.rcParams['font.sans-serif'] = "Lato"
    # pal = sns.color_palette('Set2')
    sns.set(style="ticks" )
    for row, pert_type in enumerate(list_pert_types):
        for col, (metric, xlim) in enumerate(zip(list_metrics, list_xlims)):
            metric_label = list_metric_labels[col]
            pert_type_label = list_pert_type_labels[row]
            plt.subplot(n1,n2,counter)
            if col == 0:
                plt.ylabel(pert_type_label + "  ", rotation='horizontal', fontweight=fontweight, horizontalalignment="right", fontname=font, color=fontcolor)
            counter += 1
            bins = np.linspace(xlim[0], xlim[1], 31)
            plt.hist(sigstats[sigstats['pert_type'] == pert_type][metric].dropna(), bins, histtype="stepfilled", normed=True, color=pal[col], lw=1, edgecolor="#555555")
            sns.despine(top=True, right=True, trim=True, offset=0)

            plt.xlim(xlim[0], xlim[1])
            if row == 0:
                plt.title("{}\n".format(metric_label), fontweight=fontweight, fontname=font, color=fontcolor)
            if row != 4:
                plt.xticks([])
            remove_ticks()
            plt.grid(False)
        plt.subplots_adjust(wspace=0.1, hspace=0)

        plt.savefig(outfile, dep=600, bbox_inches = 'tight')
        
        
def plot_selected_points_among_all_OLD(df_BC, df_selected, additional_sig_ids=None, hue_dict=None, label=""):
#     sns.reset_orig()
    plt.rcParams['font.family'] = "Lato"
    sig_ids = df_selected['sig_id'].values
    df = df_BC.copy()
    df = df[df['pert_type'] == "trt_cp"]
    df = df.set_index('sig_id')[['cc_q75', 'ss_ltn2', 'killed_jaccard', 'overlap_pvalue_minus_log10']]
    
    columns = ['Replicate\ncorrelation', 'Selectivity', 'Overlap (Jaccard)', 'Overlap (significance)']
    
    # Adjust the metric columns
    df['ss_ltn2_normed'] = df['ss_ltn2'] / 100.
    df['selectivity'] = 1. - df['ss_ltn2_normed'] 
    df['overlap_significance'] = df['overlap_pvalue_minus_log10'] / 5.
    df = df[df['ss_ltn2'] > 0]
    
    df = df[['cc_q75', 'selectivity', 'killed_jaccard', 'overlap_significance']]
    # rename the columns
    df.columns = columns 

    # Add selection flag
    df['Nominated'] = "No"
    df.loc[sig_ids, 'Nominated'] = "Yes"
    
    if additional_sig_ids is not None:
        df.loc[additional_sig_ids, 'is_selected'] = "by recall"
        
    df.sort_values('is_selected', inplace=True)
    df = df.dropna()
    print df.shape
    
    colors = df['is_selected'].apply(lambda x:(0,0,0,0.05) if x == 0 else (1,0,0,0.7))
    sizes = df['is_selected'].apply(lambda x: 1 if x == 0 else 40)
    
    # Plot
#     plt.figure(figsize=(30,30))
#     g = sns.PairGrid(df[['cc_q75', 'ss_ltn2', 'killed_jaccard', 'overlap_pvalue_minus_log10']])
#     g.map(plt.scatter, s=sizes, color=colors, lw=0)
    with sns.axes_style('ticks') as c1, sns.plotting_context(font_scale=5) as c2:
        sns.plotting_context(font_scale=0.5)
        g = sns.pairplot(df, hue="is_selected", vars=columns,
#                     palette = "YlOrRd", 
                    palette = hue_dict if hue_dict else "Greys", 
                    plot_kws=dict(s=10, lw=0), 
#                     diag_kws=dict(gridsize=80),
                     diag_kws=dict(normed=True, alpha=0.5, bins=np.arange(0,1,0.05), histtype="stepfilled", edgecolor="#333333" ),
                     diag_kind="hist"
                    )
        for i, j in zip(*np.triu_indices_from(g.axes, 1)):
            g.axes[i, j].set_visible(False)
#         g._legend.remove()

        

def _select_from_each_bucket(d, frac=0.5, dict_min_sizes=None):
    df = d.sort_values(['cc_q75', 'is_high_replicate_overlap_significance', 'killed_jaccard'], ascending=False)
    columns = [c for c in df.columns if c not in ['ss_bucket', 'is_high_replicate_overlap_significance']]
    N = len(df)
    n = int(frac * N)
    if dict_min_sizes:
        # ss_bucket_label
        try:
            print d.name
            min_size = dict_min_sizes[d.name]
            n = max(n, min_size)
            print "adjusting bucket size"
        except:
            print "buket label not found" , d.name
            pass
    return df[columns][:n]


def nominate_compounds(df, sample_frac=0.35, dict_min_sizes=None):
    logger.info('Nominating compounds from a total of {:,}'.format(len(df)))
    df_high_reproducibility = df[
            (df['ss_bucket'] > 0) & (df['pert_type'] == "trt_cp") & (df['is_high_cc'] == 1) & (df['is_high_replicate_jaccard'] == 1) 
    ]
    logger.info("Number of highly reproducible compounds : {}".format(len(df_high_reproducibility)))
    df_selected = (
       df_high_reproducibility
        .groupby(['ss_bucket_label'])
        .apply(_select_from_each_bucket, frac=sample_frac, dict_min_sizes=dict_min_sizes)
    )
    print "Number of selected compounds: {}".format(len(df_selected))
    return df_selected






    
def nominate():
    filename = "./datasets/CBRANT_all_viability_replicate_agreement_stats.csv"
    logger.info('Loading full sigstats for all CBRANT from file')
    df_BC = pd.read_csv(filename, sep="\t", index_col=0)
    
    logger.info('Running nominate_compounds()')
    df_selected = nominate_compounds(df_BC)
    otufile = "./nominated_compounds_v2.csv"
    
    logger.info('Exporting nominated CBRANT LITMUS compounds to file {}'.format(outfile))
    df_selected.to_csv(outfile, sep="\t")
    
    imagefile = "./nominated_compounds_paird_plot.png"
    logger.info('Saveing sigstats grid plot for CBRANT nominations to file {}'.format(imagefile))
    plot_selected_points_among_all(df_BC, df_selected, outfile=iamgefile)
    
def main():
    
#     compute_and_export_stats_CBRANT_B()
#     compute_and_export_stats_CBRANT_C()
    #combine_and_export_stats_CBRANT_ALL()

    compute_and_export_stats_CBRANT_ALL()
    return


#     nominate()
#     return
    #compute_and_export_stats_CBRANT_B()
    #compute_and_export_stats_CBRANT_C()
    compute_and_export_stats_CoP_CORE()
    return

    filenameC = './datasets/CBRANT_C_viability_replicate_agreement_stats.csv'
    filenameB = './datasets/CBRANT_B_viability_replicate_agreement_stats.csv'
    filenameAll = './datasets/CBRANT_all_viability_replicate_agreement_stats.csv'
    #plot_metrics_by_pert_type(siginfo_filename=filenameC, outfile="CBRANT_C_metrics_distributions_by_pert_type.png")
    #plot_metrics_by_pert_type(siginfo_filename=filenameB, outfile="CBRANT_B_metrics_distributions_by_pert_type.png")
    plot_metrics_by_pert_type(siginfo_filename=filenameAll, outfile="CBRANT_all_metrics_distributions_by_pert_type.png")
     
                          
if __name__ == "__main__":
    import matplotlib
    matplotlib.use('Agg')
    import seaborn as sns
    import matplotlib.pyplot as plt  
    logging.basicConfig(level=logging.INFO)
    main()
else:
    import seaborn as sns
    import matplotlib.pyplot as plt  





