#!/usr/bin/python2

usage = '''
  Takes the bowtie results and removes duplicates and assigns positional
  weights to those reads that map to multiple locations in the genome.

  Usage:
   python bt_postProcEMParallel.py MERGE BOWTIE

      MERGE = windows bed file
      BOWTIE = all the bowtie results for all read lengths

'''

import argparse
import sys
import os
import subprocess as sp
import logging


def load_bowtie_hits(hitName):
    '''
    Loads bowtie hits in first loop and counts total alignments and unique
    alignments.  In the second loop, assigns a proportional count based on
    the number of locations a read
    '''
    loc = 0
    hits = {}
    counts = {}
    tags = {}
    bowtie_aln_c = 0
    with open(hitName, 'r') as f:
        for l in f:
            bowtie_aln_c += 1
            chrB, locB, endB, readTag1, gs, strB = l.rstrip().split('\t')
            readTag, seq, qual, null = readTag1.split(' ') 
    
    # tags (bowtie alignment name) dict has the tag name as the key and the number of times it appears in the results as the value
            try:
                tags[readTag] += 1
            except KeyError:
                tags[readTag] = 1
        unique_bowtie_aln = len(tags)
        logging.info("Unique bowtie alignments: {}".format(unique_bowtie_aln))
        logging.info("Total number of bowtie alignments: {}\n".format(bowtie_aln_c))
    
    # This does a positional assignment of the reads mapping to multiple locations in genome
        f.seek(0)
        c = 0
        for l in f:
            chrB1, locB, endB, readTag1, gs, strB = l.rstrip().split('\t')
            readTag, seq, qual, null = readTag1.split(' ') 
    
            denom = tags[readTag]
            posInfo = '-'.join([locB, endB])
            chrB = ':'.join([chrB1, strB])
            c += float(1) / tags[readTag]
            if chrB not in hits:
                hits[chrB] = {}
                counts[chrB] = {}
            try:
                counts[chrB][posInfo] += float(1) / tags[readTag]
            except KeyError:
                hits[chrB][posInfo] = ':'.join([chrB1, locB, endB, strB])
                counts[chrB][posInfo] = float(1) / tags[readTag]
        
        logging.debug('Total proportional assignment count = {}\n'.format(int(c)))
    return hits, counts


def bowtie_aln_overlapping_windows(loc, tname, LibBed):
    '''
    Check to see if bowtie read falls within a window.  Creates bed file of
    bowtie alignment on chromosome and strand basis.
    windowBed is a bedtools program which compares a bedLine region from one 
    file (-a) overlaps with any bedline region in second file (-b) ON THE SAME 
    STRAND (-sm).  -w is a buffer on either side to include, which is set to
    zero, requiring an exact overlap.
    '''
    with open(tname, 'w') as f:
        for location in loc:
            chN, locB, endB, strB = loc[location].split(':')
            f.write('{}\n'.format('\t'.join([chN, locB, endB, location, '1', strB])))
    cmd = 'windowBed -a {} -b {} -sm -w 0'.format(tname, LibBed)
    readWin = sp.Popen(cmd, stdout=sp.PIPE, shell = True).communicate()[0]
    return readWin.split('\n')[:-1]


def write_bowtie_results_file(dirName, chr_str, output_li):
    '''
    Write the results to a chromosome#.results file in the g1Results directory
    '''
    chromosome = chr_str.split(':')[0].upper()
    fname = '{}{}.results'.format(dirName, chromosome)
    with open(fname, 'a') as fo:
        fo.write('{}\n'.format('\n'.join(output_li)))


def bowtie_chromosome_results(bowtie_res, hits, counts, windows, temp_d):    
    '''
    Checks for an overlap between the bowtie result and a window region.
    Adds all counts that fall within that window region.
    Write output to chromosome-specific results file.
    '''
    c = 0
    baseDirName = os.path.split(os.path.abspath(bowtie_res))[0]
    dirName = '{}/g1Results/'.format(baseDirName)
    os.system('mkdir -p {}'.format(dirName))
    tname = '{}bowtie_postprocessing_temp.txt'.format(temp_d)
    for chr, loc in sorted(hits.iteritems()):
        readWin = bowtie_aln_overlapping_windows(loc, tname, windows)
        outArray = []
        for l in readWin:
            parts = l.rstrip().split('\t')
            location = parts[3]
            Hchr = parts[6].upper()
            HSt = int(parts[7])
            HEd = int(parts[8])
            HStr = parts[11]
            winName = '{}:{}-{}({})'.format(Hchr, HSt, HEd, HStr)
            chN, locB, endB, strB = hits[chr][location].split(':')
            chN = chN.upper()
            locB = int(locB)
            readSize = int(endB) - int(locB) + 1
            if strB == '-':
                mystart = HEd - locB - readSize + 1
                myEnd = mystart + readSize - 1
            else:
                mystart = locB - HSt
                myEnd = mystart + readSize - 1
            outLine = '\t'.join(map(str, [chN, mystart, myEnd, winName, counts[chr][location], strB]))
            outArray.append(outLine)
            c += counts[chr][location]
        
        logging.info("Total counts = {0:.2f}, after adding counts from {1}".format(c, chr))
        write_bowtie_results_file(dirName, chr, outArray)


def main(windows, bowtie_res, temp_d):
    logging.info('\n### Running bowtie post-processing ###')
    hits, counts = load_bowtie_hits(bowtie_res)
    bowtie_chromosome_results(bowtie_res, hits, counts, windows, temp_d)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='''miRquant - analysis of small RNA sequencing data''',
            formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        'windows',
        action='store',
        help='Windows file (*merge.bed)')
    parser.add_argument(
        'bowtie_results',
        action='store',
        help='Combined bowtie results (*allGS.bed)')
    parser.add_argument(
        'temp_dir',
        action='store',
        help='Define folder where temporary file will be saved')
    arg = parser.parse_args()
    main(arg.windows, arg.bowtie_results, arg.temp_dir)
