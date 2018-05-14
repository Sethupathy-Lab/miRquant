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
                        generate_normalized_RPMYM_counts, \
                        statistics, \
                        assemble_xls, \
                        autoDESeq


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


def mapping_stats(basePath, outPath, samples):
    '''
    Run mapping_stats.py on all the samples
    '''
    print "Calculating mapping statistics... "
    generate_mapping_info.main(basePath, outPath, samples)
    print "DONE!\n"


def length_distribution(outPath, samples):
    '''
    Run lenDist.py on all the samples, with the --image flag
    '''
    print "Running length distribution script..."
    lenDist.main(outPath, samples)
    print "DONE!\n"


def RPMMandRPMMM(species, base_path, outPath, samples):
    '''
    Run RPMM and RPMMM on all the samples
    '''
    print "Generating normalized counts tables..."
    generate_normalized_counts.main(species, outPath, base_path, samples)
    if species == 'hsa':
        generate_normalized_RPMYM_counts.main(species, outPath, base_path, 'RPMYM', 100, samples)
    print "DONE!\n"


def calculate_statistics(basePath, outPath):
    '''
    If conditions.txt exists, calculate statistics.
    If comparisons.txt exists, calculate the pair-wise comparisons.
    '''
    if os.path.exists('{}/conditions.csv'.format(basePath)):
        print "Calculating statistics..."
        os.system('cp {}/conditions.csv {}/'.format(basePath, outPath))
        os.system('cp {}/comparisons.csv {}/'.format(basePath, outPath))
        RPMMM = '{}/RPMMM_mirs_over_50.csv'.format(outPath)
        cond = '{}/conditions.csv'.format(outPath)
        comp = '{}/comparisons.csv'.format(outPath)
        statistics.main(RPMMM, cond, comp, outPath)
        print "DONE!\n"


def DESeq(basePath, outPath, D):
    '''
    If DESeq flag given and conditions file exists, run DESeq2 on 
    the raw counts file.
    '''
    if os.path.exists('{}/conditions.csv'.format(basePath)) and D:
        print "Running DESeq..."
        autoDESeq.runDESeq(outPath)
        print "DONE!\n"


def main(arg):
    cfg = load_mirquant_config_file(arg.conf)
    outPath = create_output_folder(cfg['paths']['output'])
    samples = sample_input_location(cfg['paths']['project'], cfg['paths']['output'])
    mapping_stats(cfg['paths']['project'], outPath, samples)
    length_distribution(outPath, samples)
    RPMMandRPMMM(cfg['parameters']['species'], cfg['paths']['project'], outPath, samples)
    calculate_statistics(cfg['paths']['project'], outPath)
    DESeq(cfg['paths']['project'], outPath, arg.DESeq)
    assemble_xls.main(outPath)


if __name__ == '__main__':
    f_utils.check_for_input(sys.argv, usage)
    parser = argparse.ArgumentParser(
             description='Runs all the final summary scripts')
    parser.add_argument(
            'conf',
            action='store',
            help='Path to configuration directory')
    parser.add_argument(
            '-d','--DESeq',
            action='store_true',
            help='Run DESeq2 on raw counts file (requires conditions file)')
    arg = parser.parse_args()
    main(arg)
