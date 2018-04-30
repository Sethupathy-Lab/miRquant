#!/usr/bin/python2

import os

def runDESeq(outPath):
    out_dir = os.path.dirname(__file__)
    counts = '{}/raw_miR_counts.csv'.format(outPath)
    cond = '{}/conditions.csv'.format(outPath)
    comp = '{}/comparisons.csv'.format(outPath)
    os.system('Rscript --vanilla {}/miRDESeq.R {} {} {}'.format(out_dir, 
                                                                counts,
                                                                cond,
                                                                comp))
