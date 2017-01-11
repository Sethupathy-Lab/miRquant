#!/usr/bin/python2

usage='''
This determines the adaptor sequence for a sample based on the sample's file
name, either automatically or by user input.
'''

import argparse
import sys
import os
import re
import logging
from os.path import isfile, \
                    join
from itertools import islice
from utils import check_input, \
                  load_mirquant_config_file


def check_for_adapter_file(file):
    '''
    Checks if adapter file exists for file, if it doesn't,
    returns file name and location
    '''
    logging.info('\n### Generating adapter file ###')
    dir = '{}/'.format(os.path.dirname(file))
    name = os.path.basename(file).replace('.fastq', '')
    if os.path.isfile('{}{}.adaptor'.format(dir, name)):
        logging.info('Existing adapter file to be used for {}'.format(name))
        return 'null', 'null', False
    return dir, name, True


def check_for_barcode_file(dirc, name, cfg):
    '''
    Checks for optional file containing barcodes
    '''
    if isfile(join(dirc, 'barcodes.txt')):
        barcode_di = {}
        with open(join(dirc, 'barcodes.txt'), 'r') as f:
            for l in f:
                file, barcode = l.rstrip().split()
                barcode_di[file] = barcode
        try:
            barcode = barcode_di['{}.fastq'.format(name)]
        except KeyError:
            logging.error('ERROR: No file/barcode in barcodes.txt for {}'.format(file))
            sys.exit()
        adapter = create_adapter(barcode, cfg['cutadapt']['adapter'])
        write_adapter_file(dirc, name, adapter)
        return True
    else:
        return False

        
def write_barcode_distribution(barcode_di, keys, output):
    '''
    Writes the statistics on the distribution of identified barcodes.
    The top barcode is the true barcode and will have much more than
    others.
    '''
    out_name = join(output, 'barcode_distribution.log')
    with open(out_name, 'w') as f:
        f.write('barcode\tcount\n')
        for k in keys:
            f.write('{}\t{}\n'.format(k, barcode_di[k]))
            

def scan_fastq_for_barcode(file, name, log_dir):
    '''
    Scans first 50,000 reads and takes the barcode from the 
    readname line.  Determines the barcode with the highest
    frequency to be the true barcode.
    '''
    c = 0
    barcode_di = {}
    with open(file, 'r') as f:
        while c < 50000:
            li = list(islice(f, 4))
            barcode = li[0].rstrip().split(':')[-1]
            try:
                barcode_di[barcode] += 1
            except KeyError:
                barcode_di[barcode] = 1
            c += 1
    keys = sorted(barcode_di, key = barcode_di.get, reverse = True)
    write_barcode_distribution(barcode_di, keys, log_dir)
    return keys[0]


def validate_adapter(adapter):
    '''
    Checks to make sure adapter is only made of nucleotides
    '''
    nucs = ['A', 'T', 'G', 'C', 'N']
    for n in adapter:
        if n not in nucs:
            print 'ERROR: Adapter contains non-nucleotides'
            print 'Adapter = {}'.format(adapter)
            sys.exit()


def create_adapter(barcode, adapter):
    '''
    Inserts the barcode into the adapter in correct location
    '''
    if 'X' in adapter:
        adapter = adapter.split('X')
        adapter = '{}{}{}'.format(adapter[0], barcode, adapter[-1])
    validate_adapter(adapter)
    return adapter
    

def write_adapter_file(dirc, name, adapter):
    '''
    Write the adapter sequence to the adapter file.
    '''
    with open('{}{}.adaptor'.format(dirc, name), 'w') as f:
        logging.info('Adapter generation completed for:')
        logging.info('Sample: {}'.format(name))
        logging.info('Output saved as: {}'.format(f.name))
        logging.info('Adapter sequence: {}\n'.format(adapter))
        f.write('{}\n'.format(adapter))
                

def main(file, log_dir,  conf):
    check_input()
    cfg = load_mirquant_config_file(conf)
    dirc, name, need_adapt = check_for_adapter_file(file)
    if need_adapt == True:
        if not check_for_barcode_file(dirc, name, cfg):
            barcode = scan_fastq_for_barcode(file, name, log_dir)
            adapter = create_adapter(barcode, cfg['cutadapt']['adapter'])
            write_adapter_file(dirc, name, adapter)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
             description=usage
             )
    parser.add_argument(
            'file',
            action='store',
            help='Path to fastq file')
    parser.add_argument(
            '-l',
            action='store',
            dest='log_dir',
            default='.',
            help='Path to log directory')
    parser.add_argument(
            'conf',
            action='store',
            help='Path to configuration file')
    arg = parser.parse_args()
    main(arg.file, arg.log_dir, arg.conf)
