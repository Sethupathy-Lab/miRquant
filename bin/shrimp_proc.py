#!/usr/bin/python2

usage = '''
 Wrapper for the SHRiMP alignment section of miRquant pipeline.  The SHRiMP
 alignment step is necessary for the identification of the various edits to
 the reads.  The unmapped reads from bowtie (those which contained mismatches)
 are aligned to the previously assembled genomic windows.

 Usage:
     python shrimp_proc.py size windows base_name log_path

 Where:
     size = read length
     windows = windows.fasta file (*LIB.fa)
     base_name = base name for naming of files
     log_path = location where the logs will be saved

'''

import sys
import os
import math
import glob
import logging
import argparse
import shrimp_postProcGS
from utils import initiate_logging, \
                  load_mirquant_config_file


def make_out_dir(lib, size):
    '''
    Check if output directory exists.  If so, delete directory
    and create new, empty directory. Change into directory.
    '''
    out_dir = '{}/readSize_{}'.format(os.path.dirname(lib), size)
    try:
        os.makedirs(out_dir)
    except OSError:
        os.system('rm -r {}'.format(out_dir))
        os.makedirs(out_dir)
    os.chdir(out_dir)


def make_base_seed(size):
    '''
    Make seeds.  Seeds look something like this, 110111111100, and indicate
    where mismatches are allowed (0 allows mismatch, 1 doesn't).  Since we are
    looking for internal edits and 3' non-translated additions, the end will
    be zeros and there will be one internal zero
    '''
    logging.info('### Running SHRiMP for read size {} ###'.format(size))
    size = int(size)
    if size > 23:
        nzeros = 3
    elif size > 19:
        nzeros = 2
    elif size > 15:
        nzeros = 1
    else:
        nzeros = 0
    
    baseSeed = '1' * (size - nzeros)
    baseSeed = '0' + baseSeed[1:]

    seedList = []
    for i in range(0, len(baseSeed)):
        seed = baseSeed[-i:] + baseSeed[:-i] + '0' * nzeros
        seedList.append(seed)
    logging.info('Number of seeds: {}'.format(len(seedList)))
    return seedList


def split_db(pypath, shrimpFolder, prefix, seedGroup, lib, shl_dir):
    '''
    Split-db command from SHRiMP.
    '''
    script = '{}utils/split-db.py'.format(shrimpFolder)
    log = '{}01_splitdb_{}.log'.format(shl_dir, prefix)
    cmd = '{} {} --ram-size 46 --prefix {} --h-flag --seed {} {} &> {}'.format(
                pypath, script, prefix, seedGroup, lib, log)
    os.system(cmd)
    

def project_db(pypath, shrimpFolder, seedGroup, groupName, shl_dir, prefix):
    '''
    Project-db command from SHRiMP
    '''
    script = '{}utils/project-db.py'.format(shrimpFolder)
    log = '{}02_proj_{}.log'.format(shl_dir, prefix)
    cmd = '{} {} --shrimp-mode ls --h-flag --seed {} {} &> {}'.format(
        pypath, script, seedGroup, groupName, log)
    os.system(cmd)


def gmapper_ls(shrimpFolder, seedName, reads, qual, prefix, shl_dir):
    '''
    Gmapper-ls command from SHRiMP
    '''
    script = '{}bin/gmapper-ls'.format(shrimpFolder)
    gmpr_log = '{}03_gmapper_{}.err'.format(shl_dir, prefix)
    cmd = '{} -L {} {} -Q -N 16 -F -q -100 -g -100 -e -10 -f -10 -n 1 --qv-offset {} --shrimp-format > {}.out 2> {}'.format(
                script, seedName, reads, qual, prefix, gmpr_log)
    os.system(cmd)


def shrimp_submission(pypath, shrimp_cfg, seedList, lib, shl_dir, reads):
    '''
    These are the various SHRiMP submission scripts, including the spliting
    of the genome, and the alignment
    '''
    os.environ['SHRIMP_FOLDER'] = shrimp_cfg['path']
    shrimpFolder = shrimp_cfg['path']
    qual = shrimp_cfg['quality']
    pypath = ''
    nSeeds = len(seedList)
    nGroups = int(math.ceil(float(nSeeds) / 16))
    logging.info('Number of groups: {}\n'.format(nGroups))
    group = []
    for i in range(0, nGroups):
        prefix = 'group{}'.format(i)
        rev_seedList = seedList[::-1]
        group = rev_seedList[i * 16 : i * 16 + 16]
        logging.info('Seeds in group{}:\n{}\n'.format(i, ',\n'.join(group)))
        seedGroup = ','.join(group[::-1])
         
        split_db(pypath, shrimpFolder, prefix, seedGroup, lib, shl_dir)

        nSplits = len(glob.glob('{}*.fa'.format(prefix)))
        temp = '_'.join(['12']*len(group))

        groupName = '{}-46gb-{}seeds-{}of{}.fa'.format(prefix, temp, nSplits, nSplits)

        project_db(pypath, shrimpFolder, seedGroup, groupName, shl_dir, prefix)

        for j in range(1, nSplits + 1):
            seedName = '{}-46gb-{}seeds-{}of{}-ls'.format(prefix, temp, j, nSplits)
            gmapper_ls(shrimpFolder, seedName, reads, qual, prefix, shl_dir)
            cmd = 'rm {}*seed*'.format(seedName)
            os.system(cmd)
        cmd = 'rm {} *.genome'.format(groupName)
        os.system(cmd)
    

def main(size, lib, base, log_dir):
    cfg = load_mirquant_config_file()
    initiate_logging(log_dir, 'SHRiMP_{}.log'.format(size))

    reads = '{}{}.noHit'.format(base, size)
    working_dir = os.getcwd()

    make_out_dir(lib, size)
    seedList = make_base_seed(size)
    group_Name = shrimp_submission('pypath', cfg['shrimp'], seedList, lib, log_dir, reads)
    shrimp_postProcGS.main(working_dir, os.getcwd())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='''SHRiMP alignment of reads unable to be perfectly
            aligned with bowtie''',
            formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('size',
                        action='store',
                        help='Read length')
    parser.add_argument('lib',
                        action='store',
                        help='Windows fasta file')
    parser.add_argument('base',
                        action='store',
                        help='Basename for outputs')
    parser.add_argument('log_dir',
                        action='store',
                        help='Location of log outputs')
    arg =parser.parse_args()
    main(arg.size, arg.lib, arg.base, arg.log_dir)
