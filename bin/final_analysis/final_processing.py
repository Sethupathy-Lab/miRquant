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
import f_utils


def mapping_stats(basePath):
    '''
    Run mapping_stats.py on all the samples
    '''
    print "Calculating mapping statistics... "
    cmd = 'python generate_mapping_info.py {}'.format(basePath)
    os.system(cmd)
    print "DONE!\n"


def length_distribution(basePath):
    '''
    Run lenDist.py on all the samples, with the --image flag
    '''
    print "Running length distribution script..."
    cmd = 'python lenDist.py {} --image'.format(basePath)
    os.system(cmd)
    print "DONE!\n"


def RPMMandRPMMM(basePath, species):
    '''
    Run RPMM and RPMMM on all the samples
    '''
    print "Generating normalized counts tables..."
    cmd = 'python genNormalRPMM.py -sp {} {}'.format(species, basePath)
    os.system(cmd)
    cmd = 'python genNormalRPMMM.py -sp {} {}'.format(species, basePath)
    os.system(cmd)
    print "DONE!\n"
       

def sample_correlation(file):
    '''
    Run the sample correlation R script on the RPMM or RPMMM file
    '''
    print 'Calculating sample correlations...'
    cmd = 'Rscript sample_correlation.R {}'.format(file)
    os.system(cmd)
    print 'DONE!\n'


def main():
    f_utils.check_for_input(sys.argv, usage)
    mapping_stats(sys.argv[1])
    length_distribution(sys.argv[1])
    RPMMandRPMMM(sys.argv[1], sys.argv[2])
    sample_correlation('RPMM_mirs_over_100.csv')


if __name__ == '__main__':
    main()

        

