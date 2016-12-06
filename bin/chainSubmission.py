#!/usr/bin/python

usage = '''
 Wrapper to run the small RNA-seq pipeline.

 Usage:
     python chainSubmission.sh -o OVERLAP -e ERR -s SPEC -n NoGS -i PATH_TO_FASTQS

 Options:
     OVERLAP - cutadapt parameter, how much overlap allowed (Default = 10)
     ERR - cutadapt parameter, how many errors allowed (Default = 1)
     SPEC - which specie is this being run for (hsa, mmu, rno)
     QUAL - phred quality score requirement (Default = 33)
     NoGS - For gro-seq data, use NoGS (Default = NoGS)

'''

import sys
import argparse
import logging
import datetime
import subprocess
import os
from os.path import dirname, basename, splitext
import glob
import generate_adapter_files
import size_separate_reads
import bt_postProcEM
import reduce_shrimp
import time
from utils import load_mirquant_config_file, \
                  load_sys_config_file, \
                  build_job, \
                  initiate_logging, \
                  resource_paths, \
                  sample_output_paths, \
                  remove_if_exists


def pipe_to_logger(cmd):
    '''
    write external programs stdout to log
    '''
    proc = subprocess.Popen(cmd, 
            shell = True, 
            stdout=subprocess.PIPE , 
            stderr=subprocess.STDOUT)
    for line in proc.stdout:
        logging.info(line.strip())


def define_input_varibles(para):
    date = datetime.datetime.now().isoformat()
    MINrna = para['Minimum Read Length']
    GENOME = para['species']

    logging.info('Date: {}\n\nRun inputs:'.format(date))
    logging.info('Minimum read length = {}'.format(MINrna))
    logging.info('Genome = {}'.format(GENOME))
    return MINrna, GENOME


def set_up_output_folder(fi, out_path):
    '''
    Get file name components and create IntermediateFiles directory and
    output directories
    '''
    dr, fi_name = os.path.dirname(fi), os.path.basename(fi)
    fi_base, fi_ext = splitext(fi_name)
    dr_i = '{}/{}./IntermediateFiles'.format(dr, fi_base)
    if not os.path.isdir(dr_i):
        os.makedirs(dr_i)
    out_dir = '{}{}'.format(out_path, fi_base)
    remove_if_exists(out_dir)
    for out_loc in ['output', 'log', 'temp']:
        os.makedirs('{}/{}'.format(out_dir, out_loc))
    return dr, dr_i, fi_base


def get_maxRNA_length(fi, cutadapt):
    '''
    Get the readlength from the fastq file
    '''
    with open(fi, 'r') as f:
        next(f)
        return len(next(f).rstrip()) - cutadapt['overlap'] + cutadapt['error']


def set_lib(dir_i, fi_base):
    '''
    Set basename for most outputs in IntermediateFiles folder
    '''
    return '{}/{}'.format(dir_i, fi_base)


def cutadapt_cmd(fi, lib, cutadapt):
    '''
    Cutadapt command to submit reads for trimming
    '''
    overlap = cutadapt['overlap']
    error = cutadapt['error']
    error_rate = float(error) / overlap
    min_read_length = cutadapt['Minimum_Read_Length']
    output = '{}.fq'.format(lib)
    untrimmed = '{}_UT.fq'.format(lib)
    too_short = '{}_ST.fq'.format(lib)
    with open('{}adaptor'.format('/'.join(lib.split('/')[:-2])), 'r') as f:
        adapter = f.read().rstrip()

    cut_adapt_cmd = 'cutadapt -a {} -e {} -O {} -m {} \
            --untrimmed-output={} --too-short-output={} -o {} {}'.format(
                    adapter, error_rate, overlap, min_read_length,
                    untrimmed, too_short, output, fi)
    logging.info('\n### Trimming adapters ###')
    logging.debug('adapter = {}'.format(adapter))
    logging.debug('error rate = {}'.format(error_rate))
    logging.debug('overlap = {}'.format(overlap))
    logging.debug('minRNA = {}'.format(min_read_length))
    logging.debug('untrimmed = {}'.format(untrimmed))
    logging.debug('too short = {}'.format(too_short))
    logging.debug('output = {}'.format(output))
    logging.debug('input = {}'.format(fi))
    pipe_to_logger(cut_adapt_cmd)


def separate_by_read_length(MINrna, MAXrna, lib, output_loc):
    '''
    Separates read lengths by size into individual fastqs
    '''
    logging.info('\n### Separating reads by length ###')
    size_separate_reads.main(MINrna, MAXrna, lib, output_loc)


def bowtie(fi, length, BI, bowtie):
    '''
    Use bowtie to get reads that perfectly align to genome
    '''
    logging.info('\n### Aligning perfectly to genome ###')
    base = splitext(fi)[0]
    unaligned = '{}.noHit'.format(base)
    aligned = '{}.hits'.format(base)
    logging.info('\nRead length = {}'.format(length))
    logging.debug('Input file = {}'.format(fi))
    logging.debug('Bowtie indexes location = {}'.format(BI))
    logging.debug('Aligned reads file = {}'.format(aligned))
    logging.debug('Unaligned reads file = {}'.format(unaligned))
    if bowtie['quality'] == 33:
        bt_cmd = 'bowtie -q -nomaqround -a -m 20 --phred33-quals \
                -n 0 -e 70 -l {} -p 8 --seed=197 --un {} {} {} {}'.format(
                        length, unaligned, BI, fi, aligned)
    else:
        bt_cmd = 'bowtie -q -nomaqround -m 20 --solexa1.3-quals \
            -n 0 -e 70 -l {} -p 8 --seed=197 --un {} {} {} {}'.format(
                length, unaligned, BI, fi, aligned)
    logging.info('Bowtie command = {}'.format(bt_cmd))
    pipe_to_logger(bt_cmd)


def combine_all_bowtie_results(minRNA, maxRNA, lib):
    '''
    Combine all bowtie results and re-format result line
    '''
    logging.info('Combining all bowtie hits...')
    remove_if_exists('{}_allGS.bed'.format(lib))
    for length in range(minRNA, maxRNA + 1):
        OUTgs = '{}_allGS.bed'.format(lib)
        A1 = '{}_{}.hits'.format(lib, length)
        cmd = 'awk -v N={}-1 -F"\\t" \'{{print $3"\\t"$4"\\t"$4+N"\\t"$1" "$5" "$6"\\t1\\t"$2}}\' {} >> {}'.format(length, A1, OUTgs)
        os.system(cmd)
    return OUTgs


def window_creation(MINrna, MAXrna, lib, BI, tRNA, tmRNA):
    '''
    create genomic window for alignment of unaligned reads using SHRiMP
    '''
    logging.info('\n### Creating genomic windows ###')
    OUTgs = combine_all_bowtie_results(MINrna, MAXrna, lib)

    OUTm = '{}_merge.bed'.format(lib)
    OUTt = '{}_tmpMK.bed'.format(lib)
    OUTt2 = '{}_tmpMK2.bed'.format(lib)
    OUTgse = '{}_allGSE.bed'.format(lib)
    OUTgseu = '{}_allGSEu.bed'.format(lib)
    OUTl = '{}_LIB.fa'.format(lib)

    logging.info('slopBed...')
    os.system('slopBed -b 0 -i {} -g {}.chromSizes >> {}'.format(
        OUTgs, BI, OUTgse))

    os.system('awk -F"\\t" \'{{print $1"\\t"$2"\\t"$3"\\tNAME\\t1\\t"$6}}\' {} | uniq > {}'.format(
        OUTgse, OUTt))

    os.system('sort -k1,1 -k2,2n {} | uniq > {}'.format(OUTt, OUTgseu))

    logging.info('Add tRNA transcripts to windows with slopBed...')
    os.system('slopBed -b 40 -i {} -g {}.chromSizes >> {}'.format(tRNA, BI, OUTgseu))
    os.system('sort -k1,1 -k2,2n {} > {}'.format(OUTgseu, OUTt2))

    logging.info('Merging into windows bed...')
    os.system('mergeBed -d 65 -s -c 1 -o count -i {} > {}'.format(OUTt2, OUTt))
    os.system('''awk '{{print $1"\\t"$2"\\t"$3"\\tNAME\\t"$5"\\t"$4}}' {} > {}'''.format(
        OUTt, OUTm))
    os.system('slopBed -b 5  -i {} -g {}.chromSizes > {}'.format(OUTm, BI, OUTt))                             # add 5nt to each side of merged bed file (this is the windows to align to)
    os.system('mv {} {}'.format(OUTt, OUTm))

    logging.info('Make fasta from bed file...')
    os.system('fastaFromBed -s -fi {}.fa -bed {} -fo {}'.format(BI, OUTm, OUTt))
    os.system("tr '[:lower:]' '[:upper:]' < {} > {}".format(OUTt, OUTl))
    os.system('cat {} >> {}'.format(tmRNA, OUTl))
    os.system('rm {} {}'.format(OUTt, OUTt2))


def file_line_count(fi):
    '''
    calculates read counts by counting lines and dividing by 4
    '''
    with open(fi) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def mapping_statistics(ARG, lib, dr, fi_base):
    '''
    calculate various statistics for run
    '''
    lines = file_line_count(ARG) / 4
    tr_lines = file_line_count('{}.fq'.format(lib)) / 4
    untr_lines = file_line_count('{}_UT.fq'.format(lib)) /4
    cmd = "awk '{{print $4}}' {}_allGS.bed | sort | uniq | wc -l".format(lib)
    align_reads = subprocess.check_output(cmd, shell=True).rstrip()
    cmd = "awk '{{if (NR%4==1) print $0}}' {}*noHit | sort | uniq | wc -l".format(lib)
    unalign_reads = subprocess.check_output(cmd, shell=True).rstrip()
    
    with open('{0}/{1}./{1}.stats'.format(dr, fi_base), 'w') as f:
        f.write('file:{}\n'.format(ARG))
        f.write('TotReads:{}\n'.format(lines))
        f.write('TrimmReads:{}\n'.format(tr_lines))
        f.write('ShortReads:{}\n'.format(lines - tr_lines - untr_lines))
        f.write('EMhits:{}\n'.format(align_reads))
        f.write('EMmiss:{}\n'.format(unalign_reads))


def create_SHRiMP_tempFiles(MINrna, MAXrna, tmpDir):
    '''
    Create temporary file that needs to be removed before the SHRiMP
    result files are concatenated
    '''
    for length in range(MINrna, MAXrna + 1):
        with open('{}{}_SHRiMPwait.txt'.format(tmpDir, length), 'w') as fo:
            fo.write('Waiting for SHRiMP alignment to finish...')


def run_shrimp_alignment(MINrna, MAXrna, lib, log_dir, tmpDir, job, conf):
    '''
    Run shrimp and bowtie post-processing
    '''
    logging.info('\n### Align mismatched reads (SHRiMP alignment) ###')
    logging.info('SHRiMP job submissions:')
    create_SHRiMP_tempFiles(MINrna, MAXrna, tmpDir)
    for length in range(MINrna, MAXrna + 1):
        os.system('rm {}_{}.fq'.format(lib, length))
        shrimp_log_dir = '{}SHRiMP/{}/'.format(log_dir, length)
        if not os.path.isdir(shrimp_log_dir):
            os.makedirs(shrimp_log_dir)
        cmd = '{} python shrimp_proc.py {} {}_LIB.fa {}_ {} {}'.format(
                    job, length, lib, lib, conf, shrimp_log_dir)
        logging.info('SHRiMP submission: {}'.format(cmd))
        os.system(cmd)


def reduce_shrimp_res(temp_dir, dr_i, job):
    '''
    Check to see if jobs have finished running
    '''
    c = 0
    logging.info('Waiting for SHRiMP to finish...')
    while True:
        files = glob.glob('{}*SHRiMPwait.txt'.format(temp_dir))
        if len(files) == 0:
            logging.info('\nReducing SHRiMP size files to single file')
            os.system('{} python reduce_shrimp.py {}'.format(job, dr_i))
            break
        elif c >= (60 * 24):
            logging.warning('SHRiMP alignment did not finish in alotted time')
            logging.warning('Check logs and try running miRquant again')
            sys.exit()
        elif c % 5 == 0:
            logging.info('{} minutes elapsed'.format(c))
            time.sleep(60)
        else:
            time.sleep(60)
        c += 1


def main(arg):
    os.chdir('./bin')
    cfg = load_mirquant_config_file(arg.conf)
    scfg = load_sys_config_file(arg.conf)
    job = build_job(scfg['job'])
    dr, dr_i, fi_base = set_up_output_folder(arg.sample, cfg['paths']['output'])
    out_di = sample_output_paths(cfg['paths']['output'], fi_base)
    initiate_logging(out_di['log'], 'chainSubmission.log')
    res_li = resource_paths(cfg['parameters']['species'], cfg['paths'], cfg['parameters'])
    tRNA, tmRNA, BI = res_li[4], res_li[3], res_li[7]
    MINrna, GENOME = define_input_varibles(cfg['parameters'])
    generate_adapter_files.main(arg.sample, out_di['log'], arg.conf)
    MAXrna = get_maxRNA_length(arg.sample, cfg['cutadapt'])
    lib = set_lib(dr_i, fi_base)
    cutadapt_cmd(arg.sample, lib, cfg['cutadapt'])
    separate_by_read_length(MINrna, MAXrna, lib, out_di['output'])
    for length in range(MINrna, MAXrna + 1):
        fi = '{}_{}.fq'.format(lib, length)
        bowtie(fi, length, BI, cfg['bowtie'])
    window_creation(MINrna, MAXrna, lib, BI, tRNA, tmRNA)
    mapping_statistics(arg.sample, lib, dr, fi_base)
    run_shrimp_alignment(MINrna, MAXrna, lib, out_di['log'], out_di['temp'], job, arg.conf)
    bt_postProcEM.main('{}_merge.bed'.format(lib), '{}_allGS.bed'.format(lib), out_di['temp'])
    reduce_shrimp_res(out_di['temp'], dr_i, job)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='''Wrapper for cutadapt, bowtie, window generation, and SHRiMP''',
            formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('conf',
                        action='store',
                        help='Path to configuration directory')
    parser.add_argument('sample',
                        action='store',
                        help='Path to sample fastq')
    main(parser.parse_args())
