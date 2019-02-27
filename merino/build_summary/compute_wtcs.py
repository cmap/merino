import pandas as pd
import numpy as np
import operator


def compute_wtcs(s, qi):
    '''

    :param s: This is the vector of data within which you are querying for enrichment
    :param qi: These are the indices of your query within the vector, s
    :return:
    '''
    #TODO make variables more self documenting (cost vector, signature, query_index etc.)
    #TODO add a test
    # Length of sig
    N = len(s)
    # Length of query
    G = len(qi)
    # initialize running sum with cost of miss
    rs = pd.Series(np.ones(N) * (-1.0 / (N - G)))
    #rank signature in descending order
    rall = s.rank(ascending=False, method='first').astype(int) - 1
    # update with running sum with hit costs
    r = rall.iloc[qi]
    x = abs(s.iloc[qi])
    hits = x / sum(x)
    rs.loc[r] = hits.values
    # calculate running sum
    csum = np.cumsum(rs)
    # get connectivity score (farthest distance from 0

    max_index = np.argmax(abs(csum))

    wtcs = csum[max_index]


    return wtcs, csum
