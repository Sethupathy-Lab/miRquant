#!/usr/bin/python2

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind as ttest


def open_condition(fi):
    '''
    Open conditions file containing file name and which group it falls into.
    Create dict with key = condition and value = list of samples of condition.
    '''
    cond_di = {}
    with open(fi, 'r') as f:
        f.next()
        for l in f:
            samp, cond = l.rstrip().split()
            try:
                cond_di[cond].append(samp)
            except KeyError:
                cond_di[cond] = [samp]
    return cond_di      


def open_RPMMM(fi):
    '''
    Read RPMMM file in as dataframe
    '''
    return pd.read_csv(fi, 
                       sep=',', 
                       header=0,
                       index_col=0)


def reorder_by_condition(df, cond_di):
    '''
    Calculates the average value for treatment and adds to dataframe.
    Orders the dataframe.
    '''
    order = []
    for k, v in cond_di.iteritems():
        print k
        df['AVG_{}'.format(k)] = df[v].mean(axis=1)
        order += sorted(v)
        order.append('AVG_{}'.format(k))
    return df[order]


def make_comparisons(df, fi, cond_di):
    '''
    Open comparison file, where each line is the comparison to do, for example,
    to compare treatmentA to control, the comparison file line would be:

    treatmentA    control

    Comparison divides avgerage treatmentA value by avgerage control value
    '''
    with open(fi, 'r') as f:
        li = [l for l in f.read().split('\n') if l]
    tdf = pd.DataFrame()
    for comp in li:
        cdf = pd.DataFrame()
        n, d = comp.split()
        name = 'vs'.join([n, d])
        # Calculate fold change
        cdf[name] = df['AVG_{}'.format(n)] / df['AVG_{}'.format(d)]
        # Run T-Test to get p-value
        cdf['{}_pVal'.format(name)] = ttest(df[cond_di[n]], 
                                            df[cond_di[d]], 
                                            axis = 1, 
                                            equal_var=True)[1]
        cdf = cdf.replace([np.inf, np.nan], 'NA')
        tdf = pd.concat([tdf, cdf], axis = 1)
        cdf['AVG_{}'.format(n)] = df['AVG_{}'.format(n)]
        cdf['AVG_{}'.format(d)] = df['AVG_{}'.format(d)]
        cdf = cdf.sort_index(by=['{}_pVal'.format(name)], ascending=[True])
        cdf.to_csv(path_or_buf='{}.csv'.format(name),
                  sep =',')

    df = pd.concat([tdf, df], axis = 1)
    return df


def write_csv(df):
    '''
    Write output as a csv
    '''
    df.to_csv(path_or_buf='statistics.csv',
              sep =',')


def main():
    cond_di = open_condition(sys.argv[1])
    df = open_RPMMM(sys.argv[2])
    df = reorder_by_condition(df, cond_di)
    df = make_comparisons(df, sys.argv[3], cond_di)
    write_csv(df)


if __name__ == '__main__':
    main()
