#!/usr/bin/env python2.7

from Bio import SeqIO
from collections import defaultdict
import re
import sys
import time


def fasta_di(fasta, species):
    records = (r for r in SeqIO.parse(fasta, 'fasta') if species in r.id)
    di = {}
    p_di = {}
    for r in records:
        seq_id = r.description.split()[1]
        seq = str(r.seq).replace('U', 'T')
        di[seq_id] = seq
        if r.id.endswith('p'):
            p_di[seq_id] = r.id.split('-')[-1]
    return di, p_di
        

def gff_parse(gff, m_di, h_di, p_di):
    miRNA_di = defaultdict(list)
    gff_di = defaultdict(dict)
    with open(gff) as f:
        for l in f:
            if l.startswith('#'):
                continue
            c,d,t,start,stop,d2,st,d3,info = l.rstrip().split('\t')
            if t == 'miRNA':
                id, alias, name, derive = info.split(';')
                id = id.split('_')[0].replace('ID=', '')
                derive = derive.split('=')[1]
                miRNA_di[derive].append(id)
            else:
                id, alias, name = info.split(';')
                key = id.split('_')[0].replace('ID=', '')
                name = name.replace('Name=', '')

                gff_di[key]['N'] = name
                gff_di[key]['C'] = c
                gff_di[key]['S'] = start
                gff_di[key]['E'] = stop
                gff_di[key]['D'] = st
                gff_di[key]['H'] = h_di[key]

    ids_to_remove = []
    for id in gff_di.copy():
        for miR_id in miRNA_di[id]:
            if miR_id not in p_di.keys():
                gff_di[id]['M'] = m_di[miR_id]
            elif miR_id in p_di.keys():
                n_id = '{}_{}'.format(id, p_di[miR_id])
                gff_di[n_id] = gff_di[id].copy()
                gff_di[n_id]['N'] = '{}-{}'.format(gff_di[id]['N'], p_di[miR_id])
                gff_di[n_id]['M'] = m_di[miR_id]
                ids_to_remove.append(id)

    for id in set(ids_to_remove):
        gff_di.pop(id, None)

    return gff_di


def write_table_output(gff_di, species, basename):
    with open('{}_table.txt'.format(basename), 'w') as fo:
        for k, v in gff_di.iteritems():
            try:
                fo.write('{}\t\n'.format('\t'.join([v['N'],
                                              v['C'],
                                              v['S'],
                                              v['E'],
                                              v['D'],
                                              v['M'],
                                              v['H']])))
            except:
                print k, v


def process_coordinates(gff_di, species, basename):
    with open('{}_tableL.bed'.format(basename), 'w') as fo:
        mirseqs = defaultdict(list)
        mirHP = {}
        for k, v in gff_di.iteritems():
            i = v['H'].index(v['M'])

            if v['D'] == '+':
                locA = int(v['S']) + i
                locB = int(v['S']) + i + len(v['M']) - 1
            else:
                locA = int(v['E']) - i - len(v['M']) + 1
                locB = int(v['E']) - i

            name = '{}:{}'.format(v['N'], v['M'])

            fo.write('{}\n'.format('\t'.join(map(str, [v['C'], locA, locB, name, 1, v['D']]))))

            baseName = re.sub(r"-.p", "", v['N'])
            mirseqs[baseName].append(v['M'])
            mirHP[baseName] = v['H']


        for k, v in mirHP.iteritems():
            if len(mirseqs[k]) == 2:
                seq1, seq2 = mirseqs[k]
                ind1 = v.index(seq1)
                ind2 = v.index(seq2)

                if ind1 < ind2:
                    sep = ind2 - (ind1 + len(seq1) - 1)
                    l = len(seq1) - 1
                    fo.write('{}\n'.format('\t'.join(map(str, [k, sep, ind2, ind1, l]))))
                else:
                    sep = ind1 - (ind2 + len(seq2) - 1)
                    l = len(seq2) - 1
                    fo.write('{}\n'.format('\t'.join(map(str, [k, sep, ind1, ind2, l]))))


def create_empty_tRNA_files(sp):
    open('{}_mature_tRNA_LIB.fa'.format(sp), 'a').close()
    open('{}_tRNA12.bed'.format(sp), 'a').close()
    open('{}_tRNA.bed'.format(sp), 'a').close()
    open('{}_tRNAlu.bed'.format(sp), 'a').close()



def main():
    sp = sys.argv[1]
    basename = sys.argv[3]
    mature_di, p_di = fasta_di('mature.fa', sp)
    hairpin_di, null = fasta_di('hairpin.fa', sp)
    gff_di = gff_parse(sys.argv[2], mature_di, hairpin_di, p_di)
    write_table_output(gff_di, sp, basename)

    process_coordinates(gff_di, sp, basename)
    
    create_empty_tRNA_files(basename)


if __name__ == '__main__':
    main()
