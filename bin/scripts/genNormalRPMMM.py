#!/usr/bin/python2

usage='''

 Usage: python script.py species path_to_files

   species: species abbreviation (mmu, hsa, rno, cast)
   path_to_files: full path to TAB_3p_summary_mir.txt

 Outputs:
   RPMMM.csv: reads per million miR mapped 
   RPMMM_over_50.csv: miRs with one sample that has a RPMMM over 50

 Description:
   This script calculates normalizes the data across samples by calculating the 
   reads assigned per million miRs mapped reads.
'''

import sys
import os
import argparse
import f_utils


def check_species_input(species):
    '''
    This checks that a the species input makes sense.  If species input
    known to be correct, but not in list below, add species abbreviation
    '''
    if species not in ['mmu', 'hsa', 'rno', 'cast']:
        print '\nSpecies not recognized, known species: mmu, hsa, rno, cast'
        print 'If species abbreviation known to be correct, edit script.\n'
        print 'Usage: python genNormalRPMM.py species path_to_files > output_name\n'
        sys.exit()


def get_data_from_file(samples):
    '''
    Brings in data and counts from TAB_lenDist file
    '''
    datout, window, total_counts = {}, {}, {}
    for file in samples:
        with open(file, 'r') as f:
            path=os.path.dirname(file)
            file = os.path.basename(path)
            header = f.readline().split('\t')
            f.readline()
            count_i = header.index('Count')
            total_counts[file] = 0
            datout[file] = {}
            for l in f:
                l = l.split('\t')
                total_counts[file] += float(l[count_i])
                datout[file][l[0]] = l[count_i]
                window[l[0]] = 1
    return datout, window, total_counts
    

def calculate_RPMMM_and_seperate_miRs(datout, window, total_counts, species):
    '''
    Calculates the RPMM for each line.  Separates the miRs from non-mir windows,
    and makes a list of miRs that are over a threshold (default 100) for at
    least one of the samples
    '''
    sorted_files = sorted(datout)
    sample_wind = {f : {} for f in sorted_files}
    miRs_over_thresh = {}
    for wind in window:
        for file in sorted_files:
            if wind in datout[file]:
                re = 1000000*(float(datout[file][wind])/total_counts[file])
                if re >= 50 and species in wind:
                    miRs_over_thresh[wind] = 1
            else:
                re = 0
            sample_wind[file][wind] = re
    return sample_wind, miRs_over_thresh


def mirs_over_50_only(sample_miRs, miRs_over_thresh):
    '''
    Remove mirs for which there are no samples with an RPMM over 100
    '''
    mirs_over_50 = {}
    for sample in sample_miRs:
        mirs_over_50[sample] = {}
        for mir in miRs_over_thresh:
            mirs_over_50[sample][mir] = sample_miRs[sample][mir]
    return mirs_over_50


def write_output(sample_dict, output_name, outPath):
    '''
    Writes miRs where RPMM is greater than threshold for at least one sample
    to an output file (called RPMM_miRs_over_(threshold).tsv
    '''
    output_name = '{}/{}'.format(outPath, output_name)
    print output_name
    sample_list = sorted(sample_dict)
    window_list = sorted(sample_dict[sample_list[0]])
    with open(output_name, 'w') as f:
        f.write(',{}\n'.format(','.join(sample_list)))
        for window in window_list:
            f.write(window)
            for sample in sample_list:
                f.write(',{0:.2f}'.format(sample_dict[sample][window]))
            f.write('\n')
    

def main(basePath, outPath, species):
    samples = f_utils.get_sample_basePath(basePath)
    samples = f_utils.set_path_to_files(samples, 'TAB_3p_summary_miR.txt')
    datout, window, total_counts = get_data_from_file(samples)
    sample_miRs, miRs_over_thresh = calculate_RPMMM_and_seperate_miRs(datout, window, total_counts, species)
    mirs_over_50 = mirs_over_50_only(sample_miRs, miRs_over_thresh)
    write_output(sample_miRs, 'RPMMM_all.csv', outPath)
    write_output(mirs_over_50, 'RPMMM_mirs_over_50.csv', outPath)
    

if __name__ == '__main__':
    f_utils.check_for_input(sys.argv, usage)
    parser = argparse.ArgumentParser(
             description='Calculates the RPMMM across samples')
    parser.add_argument(
             'basePath', 
             action='store', 
             help='Path to where the sample output folders are located')
    parser.add_argument(
             'outPath', 
             action='store', 
             help='Path to where the output file will be located')
    parser.add_argument(
             '-sp', 
             default='sp',
             action='store', 
             help='Species used in this study')
    arg = parser.parse_args()
    main(arg.basePath, arg.outPath, arg.sp)
