#!/usr/bin/python2

import os
import sys

for sample in sys.argv[1:]:
    print 'Running for sample: {}'.format(sample)
    os.system('bsub python ./bin/chainSubmission.py {}'.format(sample))
