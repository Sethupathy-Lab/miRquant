#!/usr/bin/python2

import sys
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt
from scipy.stats import ttest_ind as ttest


usage = '''
Calculates the statistics for the final report.
These include: Average RPMMM across replicates
               Fold changes between conditions
               p-values between conditions

The main output is statistics.csv, which has the fold changes and pvalues for
the comparisons, the RPMMMs for the individual samples, and the average RPMMM
across conditions.  Additionally, there will be an output for each comparison
that contains the miRs, fold change, p-value, and average RPMMMs for each
condition.  These files are sorted by p-value.
'''


def open_RPMMM(fi):
    '''
    Read RPMMM file in as dataframe
    '''
    return pd.read_csv(fi, 
                       sep=',', 
                       header=0,
                       index_col=0)


def open_condition(fi):
    '''
    Open conditions file containing file name and which group it falls into.
    Create dict with key = condition and value = list of samples of condition.
    '''
    cond_di = {}
    with open(fi, 'r') as f:
        f.next()
        for l in f:
            cond_li = l.rstrip().split(',')
            samp, cond = cond_li[0], cond_li[1]
            try:
                cond_di[cond].append(samp)
            except KeyError:
                cond_di[cond] = [samp]
    return cond_di      


def check_input(samples, conditions, comparisons):
    '''
    holder
    '''
#    print samples
#    print conditions
    with open(comparisons) as f:
        comp = list(set([c for l in f.read().split(',') for c in l.split() if l]))
#    print comp
    samp_li = [s for k in conditions for s in conditions[k]]
    samp_issue = [s for s in samp_li if s not in samples]
    comp_issue = [c for c in comp if c not in conditions]
    if len(samp_issue) > 0:
        print '\nSample name in conditions file, but not in RPMMM file:'
        print '{}\n'.format('\n'.join(samp_issue))
        print 'Please check discrepency and execute again.\n'
        sys.exit()
    if len(comp_issue) > 0:
        print '\nCondition name in comparisons file, but not in conditions file:'
        print '{}\n'.format('\n'.join(comp_issue))
        print 'Please check discrepency and execute again.\n'
        sys.exit()


def reorder_by_condition(df, cond_di):
    '''
    Calculates the average value for treatment and adds to dataframe.
    Orders the dataframe.
    '''
    order = []
    for k, v in cond_di.iteritems():
        df['AVG_{}'.format(k)] = df[v].mean(axis=1)
        order += sorted(v)
        order.append('AVG_{}'.format(k))
    return df[order]


def log_df(df, cond_di, out_path):
    '''
    Make a log2 RPMMM_miRs_over_50 file and dataframe.
    '''
    df += 1
    df = df.applymap(np.log2)
    df = reorder_by_condition(df, cond_di)
    df = df[df.columns.drop(list(df.filter(regex='AVG_')))]
    df.to_csv(path_or_buf='{}/{}'.format(out_path, 'RPMMM_mirs_over_50_log2.csv', sep =','))
    return df


def make_comparisons(df, df_log, fi, cond_di, out_path):
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
        n, d = comp.split(',')
        name = 'vs'.join([n, d])
        # Calculate fold change
        cdf['foldChange'] = df['AVG_{}'.format(n)] / df['AVG_{}'.format(d)]
        # Run T-Test to get p-value
        #  - will spit errors if all zeros for expression, but can ignore
        cdf['pValue'.format(name)] = ttest(df[cond_di[n]], 
                                            df[cond_di[d]], 
                                            axis = 1, 
                                            equal_var=True)[1]
        cdf['pValue_log2'.format(name)] = ttest(df_log[cond_di[n]], 
                                                df_log[cond_di[d]], 
                                                axis = 1, 
                                                equal_var=True)[1]
        cdf = cdf.replace([np.inf, np.nan], 'NA')
        tdf = pd.concat([tdf, cdf], axis = 1)
        cdf['AVG_{}'.format(n)] = df['AVG_{}'.format(n)]
        cdf['AVG_{}'.format(d)] = df['AVG_{}'.format(d)]
        cdf = cdf.sort_values(by=['pValue'.format(name)], ascending=[True])
        cdf.to_csv(path_or_buf='{}/{}.csv'.format(out_path, name),
                   sep =',')

    df = pd.concat([tdf, df], axis = 1)
    return df


def write_csv(df, out_path):
    '''
    Write output as a csv
    '''
    df.to_csv(path_or_buf='{}/statistics.csv'.format(out_path),
              sep =',')


def main(RPMMM, cond, comp, out_path = './'):
    cond_di = open_condition(cond)
    df = open_RPMMM(RPMMM)
    check_input(list(df), cond_di, comp)
    df_log = log_df(df, cond_di, out_path)
    df = reorder_by_condition(df, cond_di)
    df = make_comparisons(df, df_log, comp, cond_di, out_path)
    write_csv(df, out_path)


if __name__ == '__main__':
    RPMMM = sys.argv[1]
    condition = sys.argv[2]
    comparison = sys.argv[3]
    if sys.argv[4]:
        out_path = sys.argv[4]
    main(RPMMM, condition, comparison, out_path)
