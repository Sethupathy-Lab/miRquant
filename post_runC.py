#!/usr/bin/python2

import sys
import os

summary_li = ['3p_summary', '5p_summary', 'shift_summary', 'lenDist_summary', 'ed_summary']
for sample in sys.argv[1:]:
    g1_dir = '{}/IntermediateFiles/g1Results/'.format(sample)
    try:
        os.makedirs('{}trash'.format(g1_dir))
    except OSError:
        pass
    for summ in summary_li:
        summ_loc = '{}/{}'.format(g1_dir, summ)
        os.system('bsub -J combine "cat {0}/* >> {0}.txt; mv {0} {1}/trash"'.format(summ_loc, g1_dir))
