#!/usr/bin/python2

import sys
import glob
import os


def check_for_input(arguments, usage):
    '''
    Check to make sure script is provided input
    '''
    if len(arguments) < 2:
        print usage
        sys.exit()
        

def get_sample_basePath(basePath):
    '''
    Returns the paths to the sample output folder.
    '''
    return glob.glob('{}/*.'.format(basePath))


def check_file(fi):
    '''
    Check if file exists, and if so, return.
    If file does not exist, raise warning and omit from analysis.
    '''
    if os.path.isfile(fi):
        return(fi)
    else:
        print 'WARNING: {} does not exist! Omitting sample.'.format(os.path.basename(fi))


def set_path_to_files(basePaths, path_to_file):
    '''
    Set the path to the files necessary for the script.  The basePath is the 
    location of the sample output directory, the path_to_file is the relative
    path to the necessary files from the basePath.
    '''
    fastq_path = []
    for fi in basePaths:
        path = '{}/{}'.format(fi, path_to_file)
        fastq_path.append(check_file(path))
    if len([f for f in fastq_path if f]) == 0:
        print 'No files found, exiting script: {}'.format(sys.argv[0])
    return [path for path in fastq_path if path]


def set_path_to_files_basename(basePaths, path_to_file, ext):
    '''
    Set the path to the files necessary for the script.  The basePath is the 
    location of the sample output directory, the path_to_file is the relative
    path to the necessary files from the basePath.
    '''
    fastq_path = []
    for fi in basePaths:
        FOI = '{}/{}{}'.format(path_to_file, os.path.basename(fi), ext)
        path = '{}/{}'.format(fi, FOI)
        fastq_path.append(check_file(path))
    if len([f for f in fastq_path if f]) == 0:
        print 'No files found, exiting script: {}'.format(sys.argv[0])
    return [path for path in fastq_path if path]


def check_for_one_glob_file(s, term):
    '''
    Check to make sure there is only one file returned for each sample
    '''
    files = glob.glob('{}/*{}*'.format(s, term))
    if len(files) == 1:
        return files[0]
    elif len(files) < 1:
        print 'No files for {}'.format(os.path.basename(s))
        sys.exit()
    elif len(files) > 1:
        print 'Too many files for {}'.format(os.path.basename(s))
        print 'Should return 1, returned these: {}'.format(', '.join(files))
        sys.exit()


def set_path_to_files_glob(samples, term):
    '''
    Search directory for file and add to sample list.
    '''
    return [check_for_one_glob_file(s, term) for s in samples]
