#!/usr/bin/python2

import sys
import os
from bin.utils import return_sample_results_directories, \
                      load_sys_config_file, \
                      build_job

summary_li = ['3p_summary', '5p_summary', 'shift_summary', 'lenDist_summary', 'ed_summary']

sample_res = return_sample_results_directories(sys.argv[1])
scfg = load_sys_config_file('./bin/configuration')
job = build_job(scfg['job'])
for sample in sample_res:
    print sample
    g1_dir = '{}/IntermediateFiles/g1Results/'.format(sample)
    try:
        os.makedirs('{}trash'.format(g1_dir))
    except OSError:
        pass
    for summ in summary_li:
        summ_loc = '{}/{}'.format(g1_dir, summ)
        os.system('{0} "cat {1}/* >> {1}.txt; mv {1} {2}/trash"'.format(job, summ_loc, g1_dir))
