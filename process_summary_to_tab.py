#!/usr/bin/python2

import sys
from bin import process_all_summary2tab
from bin.utils import return_sample_results_directories

conf = sys.argv[1]
samples = return_sample_results_directories(sys.argv[2])

process_all_summary2tab.main(conf, samples)
