#!/usr/bin/python2

import os
import sys
import argparse
import glob
from bin.utils import load_mirquant_config_file, \
                      load_sys_config_file, \
                      build_job


usage='''
miRquant is a pipeline for quantification and annotation of small RNA
sequencing data, specifically miRNAs.
'''


def check_config_path(conf_path):
    '''
    Make sure that configuration directory path exists
    '''
    if not os.path.isdir(conf_path):
        print "Configuration directory does not exist at specified location"
        sys.exit()


def get_fastqs(proj_dir):
    '''
    Get all files ending in .fastq or .fq and passes fastq to miRquant
    '''
    fqs = []
    for type in ['*.fq', '*.fastq']:
        fqs.extend(glob.glob('{}/{}'.format(proj_dir, type)))
    if len(fqs) == 0:
        print "No fastqs in project directory"
        print "If fastqs exist, must end with extension .fq or .fastq"
        sys.exit()
    else:
        return fqs


def main(args):
    check_config_path(args.conf)
    cfg = load_mirquant_config_file(args.conf)
    scfg = load_sys_config_file(args.conf)
    job = build_job(scfg['job'])
    fastqs = get_fastqs(cfg['paths']['project'])
    for sample in fastqs:
        print 'Running for sample: {}'.format(sample)
        os.system('{} python ./bin/chainSubmission.py {} {}'.format(job, args.conf, sample))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description=usage,
            formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('conf',
                        action='store',
                        help='Path to configuration directory')
    main(parser.parse_args())
