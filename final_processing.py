#!/usr/bin/python2

usage='''
 Usage:
   python final_processing.py /path/to/samples_directory

    where: samples_directory = location of sample fastqs

 Description:

 This is a wrapper script for the various scripts generating the outputs for the
 final report.  The scripts that are wrapped here are:
 
 generate_mapping_info.py - generates the mapping statistics for all samples
 lenDist.py - generates the read length histogram for all samples
 lenDistGraph.R - generates read length histogram image for all samples
 smRNAseq_correlation.R - generates correlation table and image for sample dists
 genNormalRPMMM.py - generates the reads per million miRs mapped for all samples
 genNormalRPMM.py - generates the reads per million mapped for all samples
 
 For more information on what the various scripts do, please consult the script
 file.
 
'''

import sys
import os
import glob
import argparse
import datetime
from bin.utils import load_mirquant_config_file
from bin.final_analysis import f_utils, \
                        generate_mapping_info, \
                        lenDist, \
                        generate_normalized_counts, \
                        statistics, \
                        assemble_xls


def create_output_folder(outPath):
    '''
    Create a uniquely named output folder for final results
    '''
    now = datetime.datetime.now()
    c = 1
    while True:
        d = '{}/{}_{}_{}_miRquant_{}/'.format(outPath, now.year, now.month, now.day, c)
        try:
            os.makedirs(d)
            break
        except OSError:
            c += 1
    return d


def sample_input_location(basePath, outPath):
    '''
    Get the sample input location for the final output scripts
    '''
    samps = [s[:-1] for s in glob.glob('{}/*.'.format(basePath))]
    return ['{}{}/output/'.format(outPath, os.path.basename(s)) for s in samps]


def mapping_stats(outPath, samples):
    '''
    Run mapping_stats.py on all the samples
    '''
    print "Calculating mapping statistics... "
    generate_mapping_info.main(outPath, samples)
    print "DONE!\n"


def length_distribution(outPath, samples):
    '''
    Run lenDist.py on all the samples, with the --image flag
    '''
    print "Running length distribution script..."
    lenDist.main(outPath, samples)
    print "DONE!\n"


def RPMMandRPMMM(species, outPath, samples):
    '''
    Run RPMM and RPMMM on all the samples
    '''
    print "Generating normalized counts tables..."
    generate_normalized_counts.main(species, outPath, samples)
    print "DONE!\n"


def calculate_statistics(basePath, outPath):
    '''
    If conditions.txt exists, calculate statistics.
    If comparisons.txt exists, calculate the pair-wise comparisons.
    '''
    if os.path.exists('{}/conditions.txt'.format(basePath)):
        os.system('cp {}/conditions.txt {}/'.format(basePath, outPath))
        os.system('cp {}/comparisons.txt {}/'.format(basePath, outPath))
        RPMMM = '{}/RPMMM_mirs_over_50.csv'.format(outPath)
        cond = '{}/conditions.txt'.format(outPath)
        comp = '{}/comparisons.txt'.format(outPath)
        statistics.main(RPMMM, cond, comp, outPath)


def main(conf):
    cfg = load_mirquant_config_file(conf)
    outPath = create_output_folder(cfg['paths']['output'])
    samples = sample_input_location(cfg['paths']['project'], cfg['paths']['output'])
    mapping_stats(outPath, samples)
    length_distribution(outPath, samples)
    RPMMandRPMMM(cfg['parameters']['species'], outPath, samples)
    calculate_statistics(cfg['paths']['project'], outPath)


if __name__ == '__main__':
    f_utils.check_for_input(sys.argv, usage)
    parser = argparse.ArgumentParser(
             description='Runs all the final summary scripts')
    parser.add_argument(
            'conf',
            action='store',
            help='Path to configuration directory')
    arg = parser.parse_args()
    main(arg.conf)
