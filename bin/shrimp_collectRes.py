#!/usr/bin/python2

import sys
import os
import glob
import re
import tempfile
import logging
import subprocess as sp
from utils import load_mirquant_config_file,\
                  resource_paths,\
                  ftoi,\
                  sample_output_paths,\
                  initiate_logging


def print_run_info(dirc, bt, species, outDir, windows_fa):
    logging.info('Run Info:')
    logging.info('Chromosome # and path: {}'.format(dirc))
    logging.info('Bowtie results file: {}'.format(bt))
    logging.info('Species genome: {}'.format(species))
    logging.info('Output directory: {}'.format(outDir))
    logging.info('Window fasta: {}\n'.format(windows_fa))


def result_file_dict(bt, dirc):
    '''
    Open bowtie results for given chromosome.  Load bowtie counts for for each
    window (how many times a read maps exactly to that window / number of
    places of perfect alignment in the genome) into a dictionary.  Put window
    in the bowtie windows and shrimp windows into a windows dictionary.
    '''
    files = glob.glob('{}/*.results'.format(dirc))
    files = [os.path.basename(f) for f in files]
    btWins = {}     # bowtie windows
    EM = {}         # exact matchs
    with open(bt, 'r') as f:
        for l in f:
            try:
                chrB, offB, endB, windB, countB, strB = l.rstrip().split('\t')
            except:
                print bt
                print l
                sys.exit()
            oName = ':'.join([offB, endB])
            if windB not in EM:
                EM[windB] = {}
            try:
                EM[windB][oName] += float(countB)
            except KeyError:
                EM[windB][oName] = float(countB)
            winFile = '{}.results'.format(windB)
            btWins[winFile] = 1
        
    # Add the SHRiMP results files to the bowtie window files
    for f in files:
        btWins[f] = 1
    return btWins, EM


def get_miR_info(mir_file):
    '''
    Load miR information from miRNA resource file.
    mirList is chromosome, location, name
    mirStrand is chromosome, location, strand
    mirLoc is name and chr:loc:strand
    '''
    mirs = {'li' : {}, 'str' : {}, 'loc' : {}}
    
    with open(mir_file, 'r') as f:
        for mir in f:
            mName, mchr, mSt, mEd, mStr, mSeq, mHp = mir.rstrip().split('\t')
            try:
                ind = mHp.index(mSeq)
            except:
                logging.warning('Mature not in hairpin')
                logging.warning(mName)
                logging.warning('Hairpin: {}'.format(mHp))
                logging.warning('Mature seq {}\n'.format(mSeq))
                ind = -1
            if mStr == "+":
                loc = int(mSt) + ind - 1
            else:
                loc = int(mEd) - ind - 1
            mchr = mchr.upper()
    
            try:
                mirs['li'][mchr][loc] = mName
                mirs['str'][mchr][loc] = mStr
            except KeyError:
                mirs['li'][mchr] = {loc : mName}
                mirs['str'][mchr] = {loc : mStr}
            mirs['loc'][mName] = "{}:{}:{}".format(mchr, loc, mStr)
    return mirs


def windowBed_temp_file(dirc, line, tRNA_resource, mirquant_output):
    '''
    Check if window overlaps with tRNA by calling bedtools windowBed function.
    '''
    tmp = '{}{}_ShrimpCollectResTempFile'.format(mirquant_output['temp'], os.path.basename(dirc))
    with open(tmp, 'w') as f:
        f.write('{}\n'.format(line))
    cmd = 'windowBed -a {} -b {} -sm -w 40'.format(tmp, tRNA_resource)
    return sp.Popen(cmd, stdout=sp.PIPE, shell = True).communicate()[0]


def make_output_file_directories(outDir):
    '''
    Make the output directories for the 5p, 3p, edit, shift, and lenDist dictionaries
    '''
    directories = ['5p_summary', '3p_summary', 'ed_summary', 'shift_summary', 'lenDist_summary']
    for d in directories:
        if not os.path.isdir('{}/{}'.format(outDir, d)):
            cmd = 'mkdir {}/{}'.format(outDir, d)
            os.system(cmd)

    
def write_output_file(outDir, res, name, e_di, sum, exMat):
    '''
    Write the output of edits and bowtie files.
    '''
    fname = '{}/{}/{}_{}.txt'.format(outDir, name, res, name)
    with open(fname, 'a+') as f:
        for k, v in e_di.iteritems():
            f.write('{}\tCount:{}\tEM:{}'.format(k, ftoi(sum[k]), ftoi(exMat[k])))
            for edit in sorted(v):
                item = '{}:{}'.format(edit, ftoi(v[edit]))
                f.write('\t{}'.format(item))
            f.write('\n')
    with open(fname, 'r') as f:
        logging.info(fname)
        logging.info(f.read())


def open_shrimp_file(filename, eCounts):
    '''
    Opens the processed SHRiMP alignment files and load into a dictionary
    where the line is the key and the count is the value
    '''
    with open(filename, 'r') as f:
        for l in f:
            readTag, window, strand, cst, cen, rst, ren, rlen, score, estr, mir, pcount = l.rstrip().split('\t')
            readkey = '\t'.join([window, strand, cst, cen, rst, ren, rlen, score, estr, mir])
            try:
                eCounts[readkey] += float(pcount)
            except KeyError:
                eCounts[readkey] = float(pcount)
    return eCounts


def set_tRNA_lookup_key(tRNAlist, lkey, tRNAlu, sSTR, STR):
    '''
    Check output from windowBed to see if there was overlap with tRNA.
    If there is a tRNA, format the value for the tRNAlu dictionary.
    '''
    if len(tRNAlist) > 0:
        parts = tRNAlist.split('\t')
        tlen = int(parts[8]) - int(parts[7])
        if parts[5] == '-':
            tst = int(parts[8]) - int(parts[2])
            ted = int(parts[8]) - int(parts[1]) - 1
        else:
            tst = int(parts[1]) - int(parts[7])
            ted = int(parts[2]) - int(parts[7]) - 1
        if int(parts[15]) > 1: 
            blocksz = parts[16].split(',')
            blockst = parts[17].split(',')
            if 'R' in STR:
                if tst > int(blocksz[0]):
                    tst += int(blockst[1]) - int(blocksz[0])
                if ted > int(blocksz[0]):
                    ted += int(blockst[1]) - int(blocksz[0])
        tRNAlu[lkey] = '_'.join(map(str, [parts[9], tst, ted, tlen, sSTR]))
    else:
        tRNAlu[lkey] = 'NA'
    return tRNAlu


def counts_over_thresh(mName, out_di, edits, lens, pCountArray, REFString, key, tRNAi, kkey):
    '''
    Add counts to results dictionaries if expression of 10 is met
    '''
    if mName not in out_di['exact']:
        out_di['exact'][mName] = 0
        out_di['lenDist'][mName] = {}
        out_di['p5'][mName] = {}
        out_di['p3'][mName] = {}
        out_di['E'][mName] = {}
        out_di['summary'][mName] = 0
        out_di['shift'][mName] = {}
    for i, e in enumerate(edits):
        countE = float(pCountArray[i])

                # recording length distribution for each mName
        l = lens[i]
        try:
            out_di['lenDist'][mName][l] += countE
        except KeyError:
            out_di['lenDist'][mName][l] = countE
               
        if e == 'EM' or e.isdigit():
            out_di['exact'][mName] += countE
            if 'E' not in out_di['p5'][mName]:
                out_di['p5'][mName]['E'] = 0
            if 'E' not in out_di['p3'][mName]:
                out_di['p3'][mName]['E'] = 0
            if 'E' not in out_di['E'][mName]:
                out_di['E'][mName]['E'] = 0
        else:
            misMatches = re.split('\d+', e)
            nmisMatches = len(misMatches)

            if misMatches[0].isalpha():
                try:
                    out_di['p5'][mName][misMatches[0]] += countE
                except KeyError:
                    out_di['p5'][mName][misMatches[0]] = countE
            else:
                try:
                    out_di['p5'][mName]['E'] += countE
                except KeyError:
                    out_di['p5'][mName]['E'] = countE

            if misMatches[-1].isalpha():
                try:
                    out_di['p3'][mName][misMatches[-1]] += countE
                except KeyError:
                    out_di['p3'][mName][misMatches[-1]] = countE
            else:
                try:
                    out_di['p3'][mName]['E'] += countE
                except KeyError:
                    out_di['p3'][mName]['E'] = countE
                        
            editsX = [c for c in re.split('(\d+)', e) if c]
            c = 0
            for x in editsX[:-1]:
                if x.isdigit():
                    c += int(x)
                else:
                    for p in x:
                        ed_loc = c + 1
                        refN = REFString[int(key) + c - 1]
                        edN = p
                        eName = '>'.join(map(str, [ed_loc, refN, edN]))
                               
                        if tRNAi[kkey] != 'NA':
                            logging.info('{}: {}: {},{}, {}: {}\t{}'.format(e, kkey, key, p, eName, tRNAi[kkey], REFString))
                        try:
                            out_di['E'][mName][eName] += countE # Add Pcount val
                        except KeyError:
                            out_di['E'][mName][eName] = countE
                        c += 1
            if len(misMatches) <= 2: 
                try:
                    out_di['E'][mName]['E'] += countE
                except KeyError:
                    out_di['E'][mName]['E'] = countE
         # Counting exact matches
        out_di['summary'][mName] += countE  # Add Pcount val
        l = '0'
        if l not in out_di['shift'][mName]:
            out_di['shift'][mName][l] = 0
        out_di['shift'][mName][l] += countE
    return out_di


def mainChunk(res, counters, bedFile, dirc, EM, TRNAfile, miRi, outDir, mirquant_output):
    # Create all these dicts for each window
    maps = {}
    ExMat = {} 
    mirs = {}             
    tRNAi = {}
    counts = {}
    mirLens = {}
    countList = {}
    eCounts = {}
    
    # Replace +- with PM for splitting, extract information about window from file name
    tmp = res.replace('(+)', ':P:')
    tmp = tmp.replace('(-)', ':M:')
    tmp = tmp.replace('(R+)', ':RP:')
    tmp = tmp.replace('(R-)', ':RM:')
    CHR, START, STOP, STR, EXT = re.split('[:-]', tmp)

    if 'P' in STR:
        sSTR = STR.replace('P', '+')
    elif 'M' in STR:
        sSTR = STR.replace('M', '-')
    else:
        logging.warning('String direction now known for: {}'.format(tmp))

    filename = '{}/{}'.format(dirc, res)
    if os.path.isfile(filename): # Check is SHRiMP result
        tRNAlu = {}
        open_shrimp_file(filename, eCounts) 

        for readkey in eCounts:
            window, strand, cstart, cend, rstart, rend, rlen, score, estr, mir = readkey.split('\t')
            pcount = eCounts[readkey]

    # Combine results file name info and line info into bed line for bedfile output
            if 'M' in STR:
                genSTART = int(STOP) - int(cend)
                genEND = int(STOP) - int(cstart) + 1
                strSYM = '-'
            else:
                genSTART = int(START) + int(cstart) - 1
                genEND = int(START) + int(cend)
                strSYM = '+'
            chrLC = CHR
            chrLC = re.sub('CHR', 'chr', chrLC)
            bedFile.append([chrLC, genSTART, genEND, estr, pcount, strSYM])
            line = '\t'.join(map(str, [chrLC, genSTART, genEND, mir, pcount, strSYM]))
    # lkey = location key (start end strand) within window
            lkey = '{}-{}-{}'.format(cstart, cend, STR)

    # Check to see if result is a tRNA; put in look up table if so, else set as NA
            if lkey not in tRNAlu:
                tRNAlist = windowBed_temp_file(dirc, line, TRNAfile, mirquant_output)
                tRNAlu = set_tRNA_lookup_key(tRNAlist, lkey, tRNAlu, sSTR, STR)


            loc = int(cstart)
            if 'M' in STR:
                pos = int(STOP) - int(cstart) + 1
            else:
                pos = int(START) + int(cstart) - 1
            leng = int(cend) - int(cstart) + 1
            tName = tRNAlu[lkey]
            if tName != 'NA':
                loc = '{}:{}:{}'.format(cstart, cend, sSTR)
            if mir == 'NA':
                mir = ':'.join(map(str, [CHR, pos, STR]))
            try:
                maps[loc] = ','.join([maps[loc], estr])
                mirLens[loc] = ','.join([mirLens[loc], str(leng)])
                countList[loc] = ','.join([countList[loc], str(pcount)])
                counts[loc] += pcount
            except KeyError:
                maps[loc] = estr
                mirLens[loc] = str(leng)
                countList[loc] = str(pcount)
                counts[loc] = float(pcount)
            if loc not in mirs:
                mirs[loc] = mir
            if loc not in tRNAi:
                tRNAi[loc] = tName
    

    # Add bowtie hits to lists; EM dict contain exact matches from bowtie
    window = res.split('.')[0]
    try:
        for item in EM[window]:
            a, b = item.split(':')
            a = int(a) + 1
            leng = int(b) - int(a) + 1

            if STR == 'M':
                pos = int(STOP) - int(b)
                winStr = '-'
                endPos = pos + leng - 1
            else:
                pos = int(START) + a - 1
                winStr = '+'
                endPos = pos + leng - 1
            mir = ':'.join(map(str, [CHR, pos, STR]))

            chrLC = CHR.replace('CHR', 'chr')
            line = '{}\n'.format('\t'.join(map(str, [chrLC, pos, endPos, mir, 1, winStr])))
            bedFile.append([chrLC, pos, endPos, 'Exact', EM[window][item], winStr])

    # Check to see if window corresponds to a tRNA
            tRNAlist = windowBed_temp_file(dirc, line, TRNAfile, mirquant_output)
            tName = 'NA'

    # This checks to see if there was an output from windowBed above, will be it tRNA exists in window
            if len(tRNAlist) > 0: # Check to see if the window overlaps with tRNA window
                parts = tRNAlist.split('\t')
                tlen = int(parts[8]) - int(parts[7])
                if parts[5] == '-':
                    tst = int(parts[8]) - int(parts[2])
                    ted = int(parts[8]) - int(parts[1]) - 1
                else:
                    tst = int(parts[1]) - int(parts[7])
                    ted = int(parts[2]) - int(parts[7]) - 1
                tName = '_'.join(map(str, [parts[9], tst, ted, tlen, winStr]))
                a = '{}:{}:{}'.format(a, b, winStr)

    # Check to see if read is a miRNA
            if CHR not in miRi['li']:
                miRi['li'] = {CHR : []}
            for p in miRi['li'][CHR]:
                if miRi['str'][CHR][p] == winStr:
                    d = abs(int(p) - int(pos))
                    if d < 9:
                        mir = miRi['li'][CHR][p]
                        break

            if a not in maps:
                maps[a] = 'EM'
                mirLens[a] = str(leng)
                mirs[a] = mir
                tRNAi[a] = tName
                counts[a] = EM[window][item]     #  Adding EM from Bowtie to Shrimp counts
                countList[a] = EM[window][item]
            else:
                maps[a] = ','.join([maps[a], 'EM'])
                mirLens[a] = ','.join(map(str, [mirLens[a], str(leng)]))
                counts[a] += EM[window][item]
                countList[a] = ','.join(map(str, [countList[a], EM[window][item]]))
    except KeyError:
        print window
        pass
    # End part adding bowtie hits to list


    # Sort hit locations by count values... Process each potential hit from max
    # Cluster all hits around max and mark as visited
    keys = sorted(counts)
    # If there is a bowtie hit, there will be window, if there is window, get sequence of window
    if len(keys) > 0:
        filesLib = glob.glob('{}/../*LIB.fa'.format(os.path.dirname(dirc)))[0]
        searchName = res.split('.')[0]
        cmd = 'grep -A1 "{}" {}'.format(searchName, filesLib)
        grepResults = sp.Popen(cmd, stdout=sp.PIPE, shell = True).communicate()[0]
        try:
            REFString = grepResults.split('\n')[1]
        except:
            print searchName
            print filesLib
            print grepResults.split('\n')
            raise

    out_di = {'exact' : ExMat, 'summary' : {}, 'p5' : {}, 'p3' : {}, 'E' : {}, 'shift' : {}, 'lenDist' : {}} 
    visited = {}
    counters['dis'] = 0

    for kkey in keys:
        if tRNAi[kkey] == 'NA':
            key = kkey
        else:
            key, endKey, kstr = kkey.split(':')

        counters['tot'] += counts[kkey]

    # kkey is the offset in the window
        if kkey not in visited:
            visited[kkey] = 1
            offset = 'NA'
            if 'M' in STR:
                key_loc = int(STOP) - int(key)
            else:
                key_loc = int(key) + int(START) - 1
            minDistance = 'NA'
            minDistanceS = 'NA'
            # If no miR, this will stay chromosome, key location, strand
            mNameMir = '{}:{}:{}'.format(CHR, key_loc, STR)
            if tRNAi[kkey] != 'NA':
                mNameMir = tRNAi[kkey]

            if CHR not in miRi['li']:
                miRi['li'] = {CHR : []}
            for p in sorted(miRi['li'][CHR]):
                if STR == 'P':
                    winStr = '+'
                    dS = key_loc - int(p)
                else:
                    winStr = '-'
                    dS = int(p) - key_loc
    # Assign miRNA name
                if miRi['str'][CHR][p] == winStr:
                    d = abs(int(p) - key_loc)
                    if minDistance == 'NA':
                        minDistance = d
                        minDistanceS = dS
                    else:
                        if d < minDistance:
                            minDistance = d
                            minDistanceS = dS
                    if d < 9:
                        mNameMir = miRi['li'][CHR][p]
                        break

            offset = minDistanceS
            if len(miRi['li'][CHR]) == 0:
                offset = 999  # if there is no miRNAs on the CHR, we don't need to modify the miRNA name to be relative to the offset
            if mNameMir in miRi['loc']:
                mirC, mirL, mirS = miRi['loc'][mNameMir].split(':')
                if STR == 'M':
                    offset = int(mirL) - key_loc
                else:
                    offset = key_loc - int(mirL)
            mNameMir2 = mNameMir
            if offset != 0:
                try:
                    if abs(int(offset)) < 40:
                        if offset < 0:
                            mNameMir2 = '{}_-_{}'.format(mNameMir, abs(int(offset)))
                        else:
                            mNameMir2 = '{}_+_{}'.format(mNameMir, abs(int(offset)))
                except:
                    print mNameMir2


            mName = '{},{},{}'.format(mNameMir2, tRNAi[kkey], offset)
            lens = mirLens[kkey].split(',')
            edits = maps[kkey].split(',')
            pCountArray = str(countList[kkey]).split(',')
            if counts[kkey] > 10:
                counters['tot_over_thres'] += counts[kkey]
                out_di = counts_over_thresh(mName, out_di, edits, lens, pCountArray, REFString, key, tRNAi, kkey)
            else:
                counters['dis'] += counts[kkey]

    counters['tot_dis'] += counters['dis']


    ExMat = out_di['exact']
    SUMMARY = out_di['summary']
    p5 = out_di['p5']
    p3 = out_di['p3']
    E = out_di['E']
    Shift = out_di['shift']
    lenDist = out_di['lenDist']

    make_output_file_directories(outDir)
    write_output_file(outDir, res, '5p_summary', p5, SUMMARY, ExMat)
    write_output_file(outDir, res, '3p_summary', p3, SUMMARY, ExMat)
    write_output_file(outDir, res, 'ed_summary', E, SUMMARY, ExMat)
    write_output_file(outDir, res, 'shift_summary', Shift, SUMMARY, ExMat)
    write_output_file(outDir, res, 'lenDist_summary', lenDist, SUMMARY, ExMat)
    return bedFile, counters


def write_shrimp_results_bed(out_name, bedFile):
    '''
    Write the combined results file
    '''
    fname = '{}_Shrimp_results.bed'.format(out_name)
    with open(fname, 'w') as fo:
        for l in bedFile:
            fo.write('{}\n'.format('\t'.join(map(str, l))))
        

def initialize_counters():
    '''
    Initialize various counters in a counter dictionary.
    '''
    return {'tot' : 0, 'tot_dis' : 0, 'tot_over_thres' : 0, 'dis' : 0 }


def write_summary_to_log(counters):
    '''
    Write how mand counts were included or discarded
    '''
    logging.info('Total counts: {}'.format(counters['tot']))
    logging.info('Total discard count: {}'.format(counters['tot_dis']))
    logging.info('Counts over 10 threshold: {}'.format(counters['tot_over_thres']))
    

def main():
    os.chdir('./bin') 
    dirc = sys.argv[1]
    cfg = load_mirquant_config_file(sys.argv[3])
    name = os.path.basename(dirc.split('./IntermediateFiles/')[0])
    mirquant_output = sample_output_paths(cfg['paths']['output'], name)
    logName  = '{}_collectRes.log'.format(os.path.basename(dirc))
    initiate_logging('{}/collect_results_logs/'.format(mirquant_output['log']), logName)
    bt = '{}.results'.format(sys.argv[1])
    res_li = resource_paths(cfg['parameters']['species'], cfg['paths'], cfg['parameters'])
    mir_file, TRNAfile = res_li[1], res_li[5]
    outDir = os.path.dirname(bt)
    filesLib = glob.glob('{}/../*LIB.fa'.format(os.path.dirname(dirc)))[0] 
    print_run_info(dirc, bt, cfg['parameters']['species'], outDir, filesLib)
    btWins, EM = result_file_dict(bt, dirc)
    mirs =  get_miR_info(mir_file)
    counters = initialize_counters()

    bedFile = [] 
    for res in sorted(btWins):
        bedFile, counters = mainChunk(res, counters, bedFile, dirc, EM, TRNAfile, mirs, outDir, mirquant_output)

    write_summary_to_log(counters)
    write_shrimp_results_bed(sys.argv[1], bedFile)
    os.system('rm {}'.format(sys.argv[2])) 


if __name__ == '__main__':
    main()
