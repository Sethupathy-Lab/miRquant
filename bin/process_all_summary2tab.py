#!/usr/bin/python2

usage = '''

 Usage: python process_all_summary2tab.py path/to/sample

  Where: path/to/sample = Full path to any unique file in g1Results directory
          -eg: /proj/seth_lab/users/USER/smrnapipeline/SAMPLE.

'''

import sys
import os
import argparse
import subprocess as sp
from utils import load_mirquant_config_file, \
                  sample_output_paths, \
                  return_sample_results_directories


def summary_to_tab(sum2tab, file, config_path):
    '''
    Takes the semi-final output and submits it to summary2tab.py
    for further formatting.
    '''
    print "Running summary2tab.py on {}".format(file)
    os.system('python {} {} {}'.format(sum2tab, file, config_path))


def summary_3p_of_subtype(term, fi_term):
    '''
    Creates a 3p addition summary file for only the miRs or tRNAs (fi_term). 
    Sum all counts for total term count for stats file.
    This will be the input for the RPMMM calculations script.
    '''
    tab_fi = 'TAB_3p_summary.txt'
    trna_fi = 'TAB_3p_summary_{}.txt'.format(fi_term)
    count = 0
    with open(tab_fi, 'r') as f, open(trna_fi, 'w') as fo:
        fo.write(f.readline())
        li = []
        for l in f:
            if term in l:
                l = l.split('\t')
                li.append(l)
                count += float(l[5])
        li = sorted(li, key=lambda x: -float(x[6]))
        for l in li:
            fo.write('\t'.join(l))
    return count


def write_total_and_miR_mapped_to_stats(samp_path, samp_name, miRc, tRNAc, yRNAc):
    '''
    Opens the miR TAB_3p_summary file and sums all of the counts to get
    total number of miRs mapped
    '''
    with open('TAB_lenDist_summary.txt', 'r') as f:
        null = f.readline()
        count = f.readline().split('\t')[5]

    with open('{}/{}.stats'.format(samp_path, samp_name), 'a+') as fo:
        fo.write('Mapped: {}\n'.format(count))
        fo.write('miRMapped: {}\n'.format(miRc))
        fo.write('tRNAMapped: {}\n'.format(int(tRNAc)))
        fo.write('yRNAMapped: {}\n'.format(int(yRNAc)))
        

def run_summary2Tab_clust(cfg, spec, sample, conf):
    '''
    Takes the sample output directory as the input further processes
    the compiled results output by collectRes.py.  Moves the output
    to the top level of the sample output directory
    '''
    sum2tab = '{}bin/summary2Tab_clust.py'.format(cfg['mirquant'])
    sum_file_path = '/IntermediateFiles/g1Results/'
        
    samp_name = os.path.basename(sample)[:-1]
    sum_dir = '{}{}'.format(sample, sum_file_path)
    os.chdir(sum_dir)

    summary_to_tab(sum2tab, 'lenDist_summary.txt', conf)
    summary_to_tab(sum2tab, '3p_summary.txt', conf)
    summary_to_tab(sum2tab, 'ed_summary.txt', conf)

    miRc = summary_3p_of_subtype(spec, 'miR')
    tRNAc = summary_3p_of_subtype('tRNA', 'tRNA')
    yRNAc = summary_3p_of_subtype('yRNA', 'yRNA')

    write_total_and_miR_mapped_to_stats(sample, samp_name, miRc, tRNAc, yRNAc)

        
def move_files_to_out_dir(out_di, sample, samp_name):
    '''
    Move files to sample output directory.
    '''
    os.system('mv TAB*.txt Shrimp_results.bed {}'.format(out_di['output']))
    os.system('cp {}/{}.stats {}'.format(sample, samp_name, out_di['output']))


def write_summary_table(sample):
    '''
    Prints out the stats for the run
    '''
    samp_name = os.path.basename(sample)
    with open('{}/{}stats'.format(sample, samp_name), 'r') as f:
        for l in f:
            print l.rstrip()


def main(conf):
    os.chdir('./bin')
    cfg = load_mirquant_config_file(conf) 
    samples = return_sample_results_directories(cfg['paths']['project'])
    for sample in samples:
        print 'Processing sample {}...'.format(sample)
        samp_name = os.path.basename(sample[:-1])
        out_di = sample_output_paths(cfg['paths']['output'], samp_name)
        run_summary2Tab_clust(cfg['paths'], cfg['parameters']['species'], sample, conf)
        move_files_to_out_dir(out_di, sample, samp_name)
        write_summary_table(sample)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
             description=usage)
    parser.add_argument(
            'conf',
            action='store',
            help='Path to configuration directory')
    arg = parser.parse_args()
    main(arg.conf)
