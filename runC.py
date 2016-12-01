#!/usr/bin/python2

import sys
import glob
import os
import time
from bin.utils import load_mirquant_config_file, \
                      load_sys_config_file, \
                      build_job, \
                      sample_output_paths, \
                      return_sample_results_directories

def combine_result_files(sample, cfg, job, temp_fi):
    '''
    Submits the shrimp_collectRes script, which collects the results from
    both bowtie and SHRiMP and annotates the various edits.
    '''
    out_di = sample_output_paths(cfg['paths']['output'], os.path.basename(sample)[:-1])
    g1Res_path = '{}/IntermediateFiles/g1Results/'.format(sample)
    files = glob.glob('{}CHR*.results'.format(g1Res_path))
    log_loc = '{}collect_results_logs/'.format(out_di['log'])
    os.makedirs(log_loc)
    for f in files:
        chr_dir = f.replace('.results', '')
        temp_name = '{}_collectRes.temp'.format(os.path.basename(chr_dir))
        temp_path = '{}{}'.format(out_di['temp'], temp_name)
        with open(temp_path, 'w') as f:
            f.write('If not removed, error in shrimp_collectRes.py\n')
        temp_fi.append(temp_path)
        os.system('{} python ./bin/shrimp_collectRes.py {} {}'.format(job, chr_dir, temp_name))
    return temp_fi


def combine_chromosome_result_files(sample_res, job):
    '''
    Combine the multiple, chromosomal result files for each of the summary
    files (seen in summary_li below) into a single result file.
    '''
    summary_li = ['3p_summary', '5p_summary', 'shift_summary', 'lenDist_summary', 'ed_summary']

    for sample in sample_res:
        g1_dir = '{}/IntermediateFiles/g1Results/'.format(sample)
        os.makedirs('{}trash'.format(g1_dir))
        for summ in summary_li:
            summ_loc = '{}/{}'.format(g1_dir, summ)
            os.system('{0} "cat {1}/* >> {1}.txt; mv {1} {2}/trash"'.format(job, summ_loc, g1_dir))


def wait_for_collect_res(temp_fi, sample_res, job):
    '''
    Waits until all result collection jobs have finished before concatenating
    the chromosomal results
    '''
    c = 0
    while True:
        if len(temp_fi) == 0:
            combine_chromosome_results_files(sample_res, job)
            break
        elif c >= (60 * 24):
            print 'Run failed, too long running'
        else:
            temp_fi = [f for f in temp_fi if os.path.exists(f)]
            time.sleep(60) 
        c += 1


def main():
    sample_res = return_sample_results_directories(sys.argv[1])
    cfg = load_mirquant_config_file('./bin/configuration/')
    scfg = load_sys_config_file('./bin/configuration/')
    job = build_job(scfg['job'])
    temp_fi = []
    for sample in sample_res:
        temp_fi = combine_result_files(sample, cfg, job, temp_fi)
    wait_for_collect_res(temp_fi, sample_res, job)
