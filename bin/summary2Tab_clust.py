#!/usr/bin/python2

import os
import yaml
import re
import subprocess as sp
import sys
from utils import load_mirquant_config_file, \
                  resource_paths


def mirANDtrna_to_bed(miR_file, tRNA_file):
    '''
    Opens the resource files SPEC_tableL.bed and SPEC_tRNA12.bed and reformats
    there lines into a different type of bedline containing chromosome, start,
    end, and strand.
    '''
    bedInfo = {}
    with open(miR_file, 'r') as f:
        for l in f:
            chr, st, ed, name, score, strand = l.rstrip().split('\t')
            mir, seq = name.split(':')
            bedInfo[mir] = '\t'.join([chr, st, ed, strand])
    
    with open(tRNA_file, 'r') as f:
        for l in f:
            chr, st, ed, name, score, strand = l.rstrip().split('\t')
            bedInfo[name] = '\t'.join([chr, st, ed, strand])
    return bedInfo


def load_features_into_di(res_line):
    '''
    Will load the different features in the result line into a features
    dictionary.  Add the all the counts together in a total counts 
    dictionary.  Add the counts to the counts dictionary as well.
    '''
    for item in res_line[1:]:
        L, V = item.split(':')
        if N not in features:
            features[N] = {}
        features[N][L] = V
        try:
            tot[L] += float(features[N][L])
        except KeyError:
            tot[L] = float(features[N][L])
    return features, tot



res = sys.argv[1]
cfg = load_mirquant_config_file(sys.argv[2])
species = cfg['parameters']['species']
res_li = resource_paths(species, cfg['paths'], cfg['parameters'])
genome, mmuFile, tRNAFile, refAnn = res_li[0], res_li[2], res_li[4], res_li[6]

bedInfo = mirANDtrna_to_bed(mmuFile, tRNAFile)

expression = {}
baseExp = {}
baseKeys = {}
bedLine = {}
features = {}
Novelarr = {}
tot = {}
tmp = res
tfile_name = 'mktemp'
with open(res, 'r') as f, open(tfile_name, 'w') as fo:
    for l in f:
        info = l.rstrip().split('\t')
        N = info[0]
        en, et, off = N.split(',')
        Nbase = re.sub('_[+-]_\d', '', en)
        C, nC = info[1].split(':')
        expression[N] = float(nC)  # holds total counts for that line

# if microRNA
        if species in N:
            try:
                baseExp[Nbase] += float(nC)
            except KeyError:
                baseExp[Nbase] = float(nC)
            # these are list of miR variants (shifts)
            try:
                baseKeys[Nbase].append(N)
            except KeyError:
                baseKeys[Nbase] = [N]
            c, s, e, r = bedInfo[Nbase].split('\t')
            s, e, off = map(int, [s, e, off])
            if r == '+':
                bedLine[N] = '\t'.join(map(str, [c, s + off, s + off + 7, N, 1, r]))
            else:
                bedLine[N] = '\t'.join(map(str, [c, e - off - 8, e - off - 1, N, 1, r]))

# if tRNA
        elif 'tRNA' in N or 'yRNA' in N:
            parts = N.split('_')
            tNbase = parts[0]
            o1 = parts[1]
            try:
                baseExp[tNbase] += float(nC)
            except KeyError:
                baseExp[tNbase] = float(nC)
            try:
                baseKeys[tNbase].append(N)
            except KeyError:
                baseKeys[tNbase] = [N]
            c, s, e, r = bedInfo[tNbase].split('\t')
            s, e, o1 = map(int, [s, e, o1])
            if r == '+':
                bedLine[N] = '\t'.join(map(str, [c, s + o1 + 1, s + o1 + 8, N, 1, r]))
            else:
                bedLine[N] = '\t'.join(map(str, [c, e - o1 - 8, e - o1 - 1, N, 1, r]))
# if chromo region
        else:
            c, s, r = Nbase.split(':')
            s = int(s)
            c = re.sub('CHR', 'chr', c)
            c = re.sub('UN', 'Un', c)
            c = re.sub('RANDOM', 'random', c)
            c = re.sub('.chr', 'chr', c)
            c = re.sub('HET', 'Het', c)
            c = re.sub('EXTRA', 'extra', c)
            R = r[0]
            if R == 'P':
                bedLine[N] = '\t'.join(map(str, [c, s + 1, s + 8, N, 1, '+']))
            else:
                bedLine[N] = '\t'.join(map(str, [c, s - 8, s - 1, N, 1, '-']))

            key = ':'.join([c, R])
            if key not in Novelarr:
                Novelarr[key] = {}
            Novelarr[key][N] = 1


        fo.write('{}\n'.format(bedLine[N]))
        features, tot = load_features_into_di(info)

# Not sure what the point of this code is
for k, v in Novelarr.iteritems():
    hits = sorted(v)
    NovelMax = {}
    for h in hits:
        en, et, off = h.split(',')
        Nbase = en
        Nbase = re.sub('_[+-]_\d', '', Nbase)
        c, s, r = Nbase.split(':')
        R = r[0]
        NbaseNov = ':'.join([c, s, R])
        nkeys = len(NovelMax)
        if nkeys < 1:
            keylist = NovelMax.keys()
        else:
            found = 0
            for item in NovelMax:
                distance = abs(int(s) - int(item))
                if distance < 9:
                    NbaseNov = ':'.join([c, item, R])
                    found = 1
                    break
            if found == 0:
                NovelMax[s] = 1

        try:
            baseExp[NbaseNov] += expression[h]
        except KeyError:
            baseExp[NbaseNov] = expression[h]
        try:
            baseKeys[NbaseNov].append(h)
        except KeyError:
            baseKeys[NbaseNov] = [h]

seed = {}
gAnn = {}
gAS = {}

tfile_name2 = 'mktemp2'
cmd = 'fastaFromBed -fi {} -s -name -tab -bed {} -fo {}'.format(genome, tfile_name, tfile_name2)
os.system(cmd)

with open(tfile_name2, 'r') as f:
    for l in f:
        N, seq = l.rstrip().split('\t')
        seed[N] = seq.upper()
os.system('rm mktemp2')

cmd = 'windowBed -w 0 -sm -a {} -b {}'.format(tfile_name, refAnn)
ann = sp.Popen(cmd, stdout=sp.PIPE, shell = True).communicate()[0]
ann = [l for l in ann.split('\n') if l]
for A in ann:
    parts = A.split('\t')
    gAnn[parts[3]] = parts[9]
    gAS[parts[3]] = parts[9]

    
keys = sorted(tot, key=tot.get, reverse = True)

expKeys = sorted(expression, key=expression.get, reverse = True)
expBaseKeys = sorted(baseExp, key=baseExp.get, reverse = True)

fOUT = 'TAB_{}'.format(res)
with open(fOUT, 'w') as fo:
    # Write header line
    fo.write('Name\tAnnotation\tmiRbaseOffset\tSeed\tPercent')
    for k in keys:
        fo.write('\t{}'.format(k))
    fo.write('\n')
    # Write total line
    fo.write('Total\t\t\t\t')
    for k in keys:
        fo.write('\t{}'.format(tot[k]))
    fo.write('\n')
    # Write feature line
    for e1 in expBaseKeys:
        expKeys1 = sorted(baseKeys[e1], key=baseExp.get, reverse = True)
        sumE = 0
        for e in expKeys1:
            sumE += float(features[e]['Count'])
        for e in expKeys1:
            en, et, off = e.split(',')
            frac = float(features[e]['Count']) / sumE
            if et == 'NA':
                if e in gAnn:
                    et = gAnn[e]
                else:
                    if e in gAS:
                        et = 'AS:{}'.format(gAS[e])
            if e not in seed:
                seed[e] = ''
            fo.write('\t'.join(map(str, [en, et, off, seed[e], frac])))
            for k in keys:
                if k not in features[e]:
                    features[e][k] = 0
                fo.write('\t{}'.format(features[e][k]))
            fo.write('\n')

#print '{}: {}\t{}: {}\n'.format(keys[0], tot[keys[0]], keys[1], tot[keys[1]])
