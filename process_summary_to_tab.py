#!/usr/bin/python2

import sys
from bin import process_all_summary2tab
from bin.utils import return_sample_results_directories

samples = return_sample_results_directories(sys.argv[1])

process_all_summary2tab.main(samples)

