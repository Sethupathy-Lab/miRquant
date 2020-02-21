#!/usr/bin/python2

usage='''
 Usage: python script.py species path_to_files

   species: species abbreviation (mmu, hsa, rno, cast)
   path_to_files: full path to TAB_lenDist_summary.txt

 Outputs:
   RPMM_all.tsv: RPMM for all contigs
   RPMM_miRs_only.tsv: RPMM for only the miRs
   RPMM_miRs_over_100.tsv: RPMM for miRs in which RPMM was over 100 for at
                           least one sample
 
 Description:
   This script normalizes the data across samples by calculating the 
   reads assigned per million mapped reads.


'''

import sys
import os
import argparse
import f_utils


def get_data_from_file(samples, spec):
    '''
    Brings in data and counts from TAB_lenDist file
    '''
    datout, window, total_c = {}, {}, {}
    mirs_dat, mirs, mir_c = {}, {}, {}
    for file in samples:
        with open(file, 'r') as f:
            file = file.split('/')[-3]
            count_i = f.readline().split('\t').index('Count')
            f.next()

            datout[file], mirs_dat[file] = {}, {}
            total_c[file], mir_c[file] = 0, 0

            for l in f:
                l = l.split('\t')
                total_c[file] += float(l[count_i])
                datout[file][l[0]] = float(l[count_i])
                window[l[0]] = 1
                if spec in l[0]:
                    mir_c[file] += float(l[count_i])
                    mirs_dat[file][l[0]] = float(l[count_i])
                    mirs[l[0]] = 1
    return datout, window, total_c, mirs_dat, mirs, mir_c


def write_raw_counts_table(datout, window, out_path):
    '''
    Write table of non-normalized count data
    '''
    sorted_files = sorted(datout)
    with open('{}raw_miR_counts.csv'.format(out_path), 'w') as f:
        f.write('miR,{}\n'.format(','.join(sorted_files)))
        for wind in window:
            if 'mir' not in wind and 'let' not in wind and 'miR' not in wind:
                continue
            f.write('{}'.format(wind))
            for fi in sorted_files:
                f.write(',{}'.format(datout[fi].get(wind, 0)))
            f.write('\n')
    

def calc_normalized_count(di, file, wind, c_di):
    '''
    Calculate the normalized count for each window.
    '''
    try:
        return 1000000 * ( di[file][wind] / c_di[file] )
    except KeyError:
        return 0


def windows_to_norm_counts(datout, window, tot_c):
    '''
    Calculates the RPMM for each line.  Separates the miRs from non-mir windows,
    and makes a list of miRs that are over a threshold (default 100) for at
    least one of the samples
    '''
    sorted_files = sorted(datout)
    norm_wind = {f : {} for f in sorted_files}
    for wind in window:
        for fi in sorted_files:
            norm_wind[fi][wind] = calc_normalized_count(datout, fi, wind, tot_c)
    return norm_wind


def mirs_over_thresh(data, thresh, keys, spec):
    '''
    Remove mirs for which there are no samples with an RPMM over 100
    '''
    mirs_over_thresh = {}
    for k in keys:
        li = [f for f in data if data[f][k] > thresh]
        if li and spec in k:
            for fi in data:
                try:
                    mirs_over_thresh[fi][k] = data[fi][k]
                except KeyError:
                    mirs_over_thresh[fi] = {k : data[fi][k]}
    return mirs_over_thresh


def write_output(sample_dict, output_name, outPath):
    '''
    Writes miRs where RPMM is greater than threshold for at least one sample
    to an output file (called RPMM_miRs_over_(threshold).tsv
    '''
    output_name = '{}{}'.format(outPath, output_name)
    sample_list = sorted(sample_dict)
    window_list = sorted(sample_dict[sample_list[0]])
    with open(output_name, 'w') as f:
        f.write(',{}\n'.format(','.join(sample_list)))
        for window in window_list:
            f.write(window)
            for sample in sample_list:
                f.write(',{0:.2f}'.format(sample_dict[sample][window]))
            f.write('\n')
    return output_name


def create_boxcox_trans(out_name):
    '''
    Creates a boxcox transformation of the RPMMM_mir_50.
    '''
    cmd = 'Rscript {}/boxcox.R {}'.format(os.path.dirname(__file__), out_name)
    os.system(cmd)
    

def main(species, outPath, base_path, samples):
    samples = f_utils.set_path_to_files_glob(samples, 'TAB_lenDist_summary.txt')
    datout, window, tot_c, mirs_dat, mirs, mirs_c = get_data_from_file(samples, species)

    write_raw_counts_table(datout, window, outPath)
    
    all_wind = windows_to_norm_counts(datout, window, tot_c)
    RPMM_mir_100 = mirs_over_thresh(all_wind, 100, window, species)
    out_name = write_output(all_wind, 'RPMM_all.csv', outPath)
    out_name = write_output(RPMM_mir_100, 'RPMM_mirs_over_100.csv', outPath)

    mir_wind = windows_to_norm_counts(mirs_dat, mirs, mirs_c)
    RPMMM_mir_50 = mirs_over_thresh(mir_wind, 50, mirs, species)
    out_name = write_output(mir_wind, 'RPMMM_all.csv', outPath)
    out_name = write_output(RPMMM_mir_50, 'RPMMM_mirs_over_50.csv', outPath)
    print 'DONE!\n'
    create_boxcox_trans(out_name)
    if os.path.exists('{}/conditions.csv'.format(base_path)):
        print 'Generating sample correlation heatmap...'
        cmd = 'Rscript {}/sample_correlation.R {} {}'.format(os.path.dirname(__file__), out_name, '{}/conditions.csv'.format(base_path))
        os.system(cmd)
        print 'DONE!\n'
        print 'Generating PCA...'
        cmd = 'Rscript {}/pca.R {} {}'.format(os.path.dirname(__file__), out_name, '{}/conditions.csv'.format(base_path))
        os.system(cmd)
    else:
        print 'Generating sample correlation heatmap...'
        cmd = 'Rscript {}/sample_correlation.R {}'.format(os.path.dirname(__file__), out_name)
        os.system(cmd)
        print 'DONE!\n'
        print 'Generating PCA...'
        cmd = 'Rscript {}/pca.R {}'.format(os.path.dirname(__file__), out_name)
        os.system(cmd)



if __name__ == '__main__':
    f_utils.check_for_input(sys.argv, usage)
    parser = argparse.ArgumentParser(
             description='Calculates the RPMM across samples')
    parser.add_argument(
             'sp', 
             action='store', 
             help='Species used in this study')
    parser.add_argument(
             'outPath', 
             action='store', 
             help='Path to where the output file will be located')
    parser.add_argument(
             'base_path', 
             action='store', 
             help='Path to where the output file will be located')
    parser.add_argument(
             'samples', 
             action='store', 
             nargs='+',
             help='Path to where the sample output folders are located')
    arg = parser.parse_args()
    main(arg.sp, arg.outPath, arg.base_path, arg.samples)
