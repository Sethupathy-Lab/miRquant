#!/usr/bin/python2

import sys
import glob

ls = glob.glob('CHR*')

fi_di = {}
for l in ls:
    try:
        fi_di[l.split('_')[0]].append(l)
    except KeyError:
        fi_di[l.split('_')[0]] = [l]


for k, v in fi_di.iteritems():
    with open(k, 'w') as fo:
        for file in v:
            with open(file, 'r') as fi:
                for line in fi:
                    fo.write(line)
