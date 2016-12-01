#!/usr/bin/python2

import sys
import glob
import os
from bin.utils import load_mirquant_config_file, \
                      load_sys_config_file, \
                      build_job, \
                      sample_output_paths, \
                      return_sample_results_directories


sample_res = return_sample_results_directories(sys.argv[1])
for sample in sample_res:
    cfg = load_mirquant_config_file('./bin/configuration/')
    scfg = load_sys_config_file('./bin/configuration/')
    job = build_job(scfg['job_quick'])
    out_di = sample_output_paths(cfg['paths']['output'], os.path.basename(sample)[:-1])
    g1Res_path = '{}/IntermediateFiles/g1Results/'.format(sample)
    files = glob.glob('{}CHR*.results'.format(g1Res_path))
    log_loc = '{}collect_results_logs/'.format(out_di['log'])
    try:
        os.makedirs(log_loc)
    except:
        pass
    for f in files:
        chr_dir = f.replace('.results', '')
        os.system('{} python ./bin/shrimp_collectRes.py {}'.format(job, chr_dir))
