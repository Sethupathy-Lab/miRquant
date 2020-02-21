#!/usr/bin/python2

usage = '''

 This script is called by the chainSubmission wrapper.  The purpose of this 
 script is to parse the trimmed reads into different size fastqs.

 Usage: python size_separate_reads.py minRNAlen maxRNAlen outfile_basename

'''

import argparse
import os
import sys
import logging
from itertools import islice
from utils import check_input
                  


def write_log_info(minRNA, maxRNA, lib_file):
    '''
    Write logging information about run
    '''
    logging.debug('script = {}'.format(sys.argv[0]))
    logging.info('Minimum read length = {}'.format(minRNA))
    logging.info('Maximum read length = {}'.format(maxRNA))
    logging.debug('Base name for file = {}'.format(lib_file))


def make_outfile_dict(minRNA, maxRNA, lib_file):
    '''
    Open the size-specific output files
    '''
    file_di, len_tot_di = {}, {}
    for length in range(int(minRNA), int(maxRNA) + 1):
        file_di[length] = open('{}_{}.fq'.format(lib_file, str(length)), 'w')
        len_tot_di[length] = 0
    file_di['NP'] = open('{}_notProc.fq'.format(lib_file), 'w')
    len_tot_di['NP'] = 0
    return file_di, len_tot_di


def separate_reads_by_length(lib_file, output_di, min, max, len_tot_di):
    '''
    Separate the cutadapt trimmed fastq reads by size to proper outfile
    '''
    tot_reads = 0
    with open('{}.fq'.format(lib_file), 'r') as f:
        while True:
            read = list(islice(f, 4))
            if not read:
                break
            tot_reads += 1
            read_len = len(read[1].rstrip())
            if int(min) <= read_len <= int(max):
                output_di[read_len].write(''.join(read))
                len_tot_di[read_len] += 1
            else:
                output_di['NP'].write(''.join(read))
                len_tot_di['NP'] += 1
#                print 'Out of size range: {}'.format(read[1].rstrip())
    return len_tot_di, tot_reads


def close_output_files(output_di):
    '''
    Close the output files
    '''
    for file in output_di:
        output_di[file].close()


def convert_histogram_to_values(len_tot_di, tot_reads):
    '''
    Takes total reads at each length and divids by total overall length
    '''
    if len_tot_di['NP'] != 0:
        print 'WARNING// Not all reads processed, {} reads out of size range'.format(len_tot_di['NP'])
    del len_tot_di['NP']
    for read_len, count in len_tot_di.iteritems():
        percent = float(count) / tot_reads * 100
        len_tot_di[read_len] = '{0:.2f}'.format(percent)
    return len_tot_di


def write_length_histogram(len_tot_di, basename, output_loc):
    '''
    Make sure all read lengths were in range and write length histogram
    '''
    directory, name = os.path.split(basename)
    with open('{}read_length_histo_{}.txt'.format(output_loc, name), 'w') as fo:
        fo.write('Sample\t{}\n'.format(name))
        for read_length, count in sorted(len_tot_di.iteritems()):
            fo.write('{}\t{}\n'.format(read_length, count))


def main(minRNA, maxRNA, basename, output_loc):
    check_input()
    write_log_info(minRNA, maxRNA, basename)
    output_di, len_tot_di = make_outfile_dict(minRNA, maxRNA, basename)
    len_tot_di, tot_reads = separate_reads_by_length(basename, output_di, minRNA, maxRNA, len_tot_di)
    len_tot_di = convert_histogram_to_values(len_tot_di, tot_reads)
    close_output_files(output_di)
    write_length_histogram(len_tot_di, basename, output_loc)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='''miRquant - analysis of small RNA sequencing data''',
            formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        'minRNA',
        type=int,
        action='store',
        help='Minimum read length')
    parser.add_argument(
        'maxRNA',
        type=int,
        action='store',
        help='Maximum read length')
    parser.add_argument(
        'basename',
        action='store',
        help='Basename for <input>.fq and outputs')
    parser.add_argument(
        'output_loc',
        action='store',
        help='Output location for histogram')
    arg = parser.parse_args()
    main(arg.minRNA, arg.maxRNA, arg.basename, arg.output_loc)
