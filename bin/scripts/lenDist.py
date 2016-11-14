#!/usr/bin/python2

usage='''
 Usage: python lenDist.py [--image] /path/to/sample_output_directory

 Output saved as lenDist.csv, and possibly lenDistHistograme.png

 Description:
   Calculates the read length distribution post-trimming across all samples.
   If the --image flag is provided, will out put a read length histogram for
   each sample.  For small RNA-seq, we expect to see a peak around 21-22, as
   these are the miRNAs.

'''
 

import sys
import argparse
import os
from itertools import islice
from os.path import splitext, basename, isfile
import f_utils


def set_path_to_fastqs(samples):
    '''
    Add the location of the <Sample>.fq to basePath
    '''
    fastq_path = []
    for s in samples:
        fi = '{}/IntermediateFiles/{}fq'.format(s, os.path.basename(s))
        fastq_path.append(f_utils.check_file(fi))
    return [fi for fi in fastq_path if fi]


def get_read_lengths(samples):    
    '''
    Make a dictionary of file names, read lengths, and total reads
    '''
    file_di = {}
    len_di = {}
    tot = {}
    for file in samples:
        base = os.path.basename(file)
        tot[base] = 0
        file_di[base] = {}
        with open(file, 'r') as f:
            for l in islice(f, 1, None, 4):
                read_len = len(l.rstrip())
                len_di[read_len] = 1
                try:
                    file_di[base][read_len] += 1
                except KeyError:
                    file_di[base][read_len] = 1
                tot[base] += 1
    return file_di, len_di, tot


def write_lenDist_output(file_di, len_di, tot, outPath):
    '''
    Calculate read length ratio and write to output fore each sample
    '''
    output_name = '{}{}'.format(outPath, 'length_distribution.csv')
    with open(output_name, "w") as f:
        f.write(',{}\n'.format(','.join(sorted(file_di))))
        for some_len in sorted(len_di):
            f.write(str(some_len))
            for sample in sorted(file_di):
                re=0
                if some_len in file_di[sample]:
                    re=file_di[sample][some_len]/float(tot[sample]) * 100
                try:
                    f.write(',{0:.1f}'.format(re))
                except:
                    print "WARNING: Writing failed for {}".format(sample)
            f.write('\n')
    return output_name


def create_length_dist_image(out_dir, out_name):
    '''
    Using R, creates a length distribution image for each sample
    '''
    os.system('Rscript --vanilla {}/lenDistGraph.R {}'.format(out_dir, out_name))
        

def main(basePath, outPath):
    samples = f_utils.get_sample_basePath(basePath)
    samples = set_path_to_fastqs(samples)
    file_di, len_di, tot = get_read_lengths(samples)
    out_name = write_lenDist_output(file_di, len_di, tot, outPath)
    create_length_dist_image(os.path.dirname(__file__), out_name)


if __name__ == '__main__':
    f_utils.check_for_input(sys.argv, usage)
    parser = argparse.ArgumentParser(
             description='Analyzed the length distribution of trimmed reads')
    parser.add_argument(
             'basePath', 
             action='store', 
             help='Path to where the sample output folders are located')
    parser.add_argument(
             'outPath', 
             action='store', 
             help='Path to where the output file will be located')
    arg = parser.parse_args()
    main(arg.basePath, arg.outPath)
