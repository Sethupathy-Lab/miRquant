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
import f_utils


def read_lengths_dict(samples):
    '''
    Load the length distribution data into a dictionary
    '''
    out_di = {}
    lengths = {}
    for file in samples:
        with open(file, 'r') as fi:
            name = fi.readline().rstrip().split('\t')[1]
            out_di[name] = {}
            for l in fi:
                a, b = l.rstrip().split('\t')
                out_di[name][a] = b
                lengths[a] = 1
    return out_di, sorted(lengths)


def write_length_distribution(outPath, len_di, lengths):
    '''
    Calculate read length ratio and write to output fore each sample
    '''
    output_name = '{}{}'.format(outPath, 'length_distribution.csv')
    with open(output_name, "w") as f:
        f.write('Length,{}\n'.format(','.join(sorted(len_di))))
        for len in lengths:
            f.write(len)
            for sample in sorted(len_di):
                try:
                    f.write(',{}'.format(len_di[sample][len]))
                except KeyError:
                    f.write(',0')
            f.write('\n')
    return output_name


def create_length_dist_image(out_dir, out_name):
    '''
    Using R, creates a length distribution image for each sample
    '''
    os.system('Rscript --vanilla {}/lenDistGraph.R {}'.format(out_dir, out_name))
        

def main(outPath, samples):
    samples = f_utils.set_path_to_files_glob(samples, 'ead_length_histo')
    len_di, lengths = read_lengths_dict(samples)
    out_name = write_length_distribution(outPath, len_di, lengths)
    create_length_dist_image(os.path.dirname(__file__), out_name)


if __name__ == '__main__':
    f_utils.check_for_input(sys.argv, usage)
    parser = argparse.ArgumentParser(
             description='Analyzed the length distribution of trimmed reads')
    parser.add_argument(
             'outPath', 
             action='store', 
             help='Path to where the output file will be located')
    parser.add_argument(
             'samples', 
             action='store', 
             nargs='+',
             help='Path to where the sample output folders are located')
    arg = parser.parse_args()
    main(arg.outPath, arg.samples)
