#!/usr/bin/python2

import sys
import glob
import os
from bin.utils import load_mirquant_config_file, \
                      sample_output_paths

for sample in sys.argv[1:]:
    cfg = load_mirquant_config_file('./bin/configuration/')
    out_di = sample_output_paths(cfg['paths']['output'], os.path.basename(sample)[:-1])
    g1Res_path = '{}/IntermediateFiles/g1Results/'.format(sample)
    files = glob.glob('{}CHR*.results'.format(g1Res_path))
    log_loc = '{}collect_results_logs/'.format(out_di['log'])
    os.makedirs(log_loc)
    for f in files:
        chr_dir = f.replace('.results', '')
        log_name = '{}/lsf_CollectLog_{}'.format(log_loc, os.path.basename(chr_dir))
        os.system('bsub -o {} -q idle python ./bin/shrimp_collectRes.py {}'.format(log_name, chr_dir))
