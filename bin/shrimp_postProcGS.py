#!/usr/bin/python

usage = '''
 Post-processing of SHRiMP output. Separates the output into chromosome folders
 and windows files.

 Usage:
  python shrimp_postProcGS.pl mirquant_dir shrimp_dir

 Where:
  mirquant   =  location to miRquant bin directory
  shrimp_dir =  location of the shrimp results for a given read length

'''

import sys
import os
import glob
import re
import logging
from utils import check_input, \
                  load_mirquant_config_file, \
                  resource_paths, \
                  sample_output_paths, \
                  ftoi


def load_mir_info(mir_fi):
    '''
    Load the information from the miR resource file into dictionaries
    '''
    logging.info('Loading miRNA resource file...')
    mirList = {}
    mirStrand = {}
    with open(mir_fi, 'r') as f:
        for l in f:
            mName, mchr, mSt, mEd, mStr, mSeq, mHp = l.rstrip().split('\t')
            try:
                ind = mHp.index(mSeq)
            except ValueError:
                logging.info('{}: miR seq not in hairpin sequence'.format(mName))
                logging.info('miR seq: {}'.format(mSeq))
                logging.info('hairpin seq: {}\n'.format(mHp))
                ind = 1
            if mStr == '+':
                loc = int(mSt) + ind
            else:
                loc = int(mEd) - ind - len(mSeq) + 1
            mchr = mchr.upper()
            try:
                mirList[mchr][loc] = mName
                mirStrand[mchr][loc] = mStr
            except KeyError:
                mirList[mchr] = {loc : mName}
                mirStrand[mchr] = {loc : mStr}
    return mirList, mirStrand
    

def load_SHRiMP_res(shrimp_dir):
    '''
    Load the output of the SHRiMP alignment into dictionaries.  Count how
    many unique alignments there are.  If string edit is only a digit, count
    read as perfect match.
    '''
    tagCount = 0
    
    files = glob.glob('{}/*.out'.format(shrimp_dir))

    hits = {}
    maps = {}
    perMatch = {}
    for fi in files:
        with open(fi, 'r') as f:
            header = f.readline()
            for l in f:
                line = l.strip("\n")
                readTag, window, strand, cstart, cend, rstart, rend, rlen, score, estr = l.rstrip().split('\t') 
                perfect_Match = estr.isdigit()
                if readTag not in maps: 
                    maps[readTag] = {window : int(score)}
                    tagCount += 1
                if window not in hits:
                    hits[window] = {readTag : line}
    
                if window not in maps[readTag]:
                    maps[readTag][window] = int(score)
                    tagCount += 1
                if readTag not in hits[window]:
                    hits[window][readTag] = line
    
                if int(score) > maps[readTag][window]:
                    maps[readTag][window] = int(score)
                    hits[window][readTag] = line
            numTags = len(maps)
    logging.info('Number of unique reads after loading SHRiMP results: {}'.format(numTags))
    return tagCount, hits, maps


def alter_window_name(window):
    '''
    Replaces (+), (-), (R+), and (R-) in name, then splits into constitutive
    parts, specifically chromosome, start, stop, strand.
    '''
    window = window.replace('(+)', ':P')
    window = window.replace('(-)', ':M')
    window = window.replace('(R+)', ':RP')
    window = window.replace('(R-)', ':RM')
    return re.split(r'[:-]', window)


def determine_max_score(windows):
    '''
    For each window, dertermine what the best alignment score is.
    If there are equal max scores, prefer the tRNA read.
    '''
    maxScore = 0
    maxStr = ''
    for win, score in windows.iteritems():
        chr, loc, str1, null = re.split(r'[:()]', win)
        if float(score) > maxScore:
            maxScore = float(score)
            maxStr = str1
        elif float(score) == maxScore:
            if 'R' in str1:
                maxStr = str1
    return maxScore, maxStr
        

def proportional_count_assignment(windows, tags, tag, maps, pCount):
    '''
    Proportional assignment of reads that map equally well to multiple
    locations within the genome.
    '''
    denom = 0
    for l in windows:
        chr, start, end, window, gs, strand = l.rstrip().split('\t')
        denom += 1
        try:
            tags[tag][window] = 1
        except KeyError:
            tags[tag] = {window : 1}
    for win in maps[tag]:
        if win not in tags[tag]:
            logging.ERROR("ERROR: No tag recorded for keys: {} {}".format(tag, win))
            sys.exit()
        tags[tag][win] /= float(denom) 
        pCount += tags[tag][win]
    
    numTags = len(maps)
    logging.info('NumTags 3: {}'.format(numTags))
    return tags, pCount


def get_best_alignments(hits, maps):
    '''
    Checks that reads are the best reads.
    maps = [readTag][window][score]
    hits = [window][readTag][res_line]
    '''
    pCount = 0
    tags = {}
    for tag, wind in maps.iteritems():
        WindArr = []
        maxScore, maxStr = determine_max_score(wind)
        for win in maps[tag].keys():
            chr, loc, str1, null = re.split(r'[:()]', win)
            if maps[tag][win] < maxScore:
                del maps[tag][win]
                del hits[win][tag]
            elif (maps[tag][win] == maxScore and 'R' in maxStr) and 'R' not in str1:
                del maps[tag][win]
                del hits[win][tag]
            else:
                chr, posA, posB, strand = alter_window_name(win)
                readTag, window, str2, cSt2, cEnd2, rSt2, rEnd2, rLen2, score2, estr2, null = hits[win][tag].split('\t')
                if 'P' in strand:
                    strand = '+'
                    t1 = int(posA)
                    posA = t1 + int(cSt2)
                    posB = t1 + int(cEnd2)
                else:
                    strand = '-'
                    t1 = int(posB)
                    posA = t1 - int(cEnd2)
                    posB = t1 - int(cSt2)
                chr = chr.replace('CHR', 'chr')
                bedline = '\t'.join(map(str, [chr, posA, posB, win, 1, strand]))
                WindArr.append(bedline)
         
        tags, pCount = proportional_count_assignment(WindArr, tags, tag, maps, pCount)
    return hits, tags, pCount


def check_if_mir(chr, pos, winStr, mirList, mirStrand):
    '''
    Check if read overlaps with miRNA, and thus, is a miRNA
    '''
    mir = 'NA'
    try:
        for p in mirList[chr]:
            if mirStrand[chr][p] == winStr:
                d = abs(int(p) - pos)
                if d < 9:
                    return mirList[chr][p]
        return mir
    except KeyError:
        logging.info('Chromosome {} has no recorded miRs'.format(chr))
        return mir


def write_processed_shrimp_output(hits, tags, mirLi, mirStr, shrimp_dir):
    '''
    Write the processed shrimp output to chromosomal directory in g1Results.
    '''
    dname = '{}/g1Results'.format(os.path.dirname(shrimp_dir))
    read_size = os.path.basename(shrimp_dir)
    for win, tag in hits.iteritems():
        chr, st, sp, str = alter_window_name(win)
        dirName = '{}/{}'.format(dname, chr)
        if not os.path.exists(dirName):
            os.makedirs(dirName)

        fname = '{}/{}.results_{}'.format(dirName, win, read_size)
        if len(tag) > 0:
            with open(fname, 'w') as fo:
                for h in hits[win]:
                    readTag, window, strand, cSt, cEn, rSt, rEn, rLen, score, estr, null = hits[win][h].split('\t')
                    if 'M' in str:
                        pos = int(sp) - int(cEn)
                        winStr = '-'
                    else:
                        pos = int(st) + int(cSt)
                        winStr = '+'
                    mir = check_if_mir(chr, pos, winStr, mirLi, mirStr)
                    fo.write('{}{}\t{}\n'.format(hits[win][h], mir, ftoi(tags[h][win])))


def remove_temp_file(readSize, tempDir):
    '''
    Remove temporary file
    '''
    print '{}{}_SHRiMPwait.txt'.format(tempDir, readSize)
    os.system('rm {}{}_SHRiMPwait.txt'.format(tempDir, readSize))


def main(conf, shrimp_dir):
    print 'Shrimp post-processing'
    check_input()
    cfg = load_mirquant_config_file(conf)
    res_li = resource_paths(cfg['parameters']['species'], cfg['paths'], cfg['parameters'])
    sample = os.path.basename(shrimp_dir.split('./IntermediateFiles/')[0])
    out_di = sample_output_paths(cfg['paths']['output'], sample) 
    logging.info('\n\n### Processing SHRiMP results ###\n')
    print 'a'
    mir_fi = res_li[1]
    print 'b'
    mirList, mirStrand = load_mir_info(mir_fi)
    print 'c'
    tagCount, hits, maps = load_SHRiMP_res(shrimp_dir)    
    print 'd'
    hits, tags, pCount = get_best_alignments(hits, maps)
    print 'e'
    write_processed_shrimp_output(hits, tags, mirList, mirStrand, shrimp_dir)
    print 'f'
    remove_temp_file(os.path.basename(shrimp_dir).split('_')[1], out_di['temp'])
    print 'g'
    logging.info('Total SHRiMP alignments = {}; proportional count = {}'.format(tagCount, pCount))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description=usage,
            formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('conf',
                        action='store',
                        help='Path to configuration directory')
    parser.add_argument('shrimp_dir',
                        action='store',
                        help='Path to readSize folder with SHRiMP results')
    arg =parser.parse_args()
    main(arg.conf, arg.shrimp_dir)
