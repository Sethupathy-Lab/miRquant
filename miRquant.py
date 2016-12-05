#!/usr/bin/python2

import os
import sys
import argparse
from bin.utils import load_sys_config_file, \
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


def main(args):
    check_config_path(args.conf)
    scfg = load_sys_config_file(args.conf)
    job = build_job(scfg['job'])
    for sample in args.samples:
        print 'Running for sample: {}'.format(sample)
        os.system('{} python ./bin/chainSubmission.py {} {}'.format(job, args.conf, sample))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description=usage,
            formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('conf',
                        action='store',
                        help='Path to configuration directory')
    parser.add_argument('samples',
                        nargs='*',
                        action='store',
                        help='Path to samples')
    main(parser.parse_args())
