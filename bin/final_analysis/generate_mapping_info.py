#!/usr/bin/python2

usage='''
 Usage: python generate_mapping_info_table.py /path/to/sample.stats
    eg: python generate_mapping_info_table.py .../smallRNA/PROJECT/*/*.stats

 Output saved as MappingInfoTable.csv

 Description:
   Calculates the mapping statistics for the miRquant run.

'''

import sys
import os
from os.path import basename, dirname, abspath
import glob
import f_utils


def mapping_stats_dict(samples):
    '''
    For each sample input, opens the SAMPLE./SAMPLE.stats file and imports
    the contents into an output dictionary.
    '''
    out_di = {}
    for file in samples:
        name = os.path.splitext(os.path.basename(file))[0]
        out_dir = '/'.join(file.split('/')[:-2])
        out_di[name] = {"File" : out_dir}
        with open(file, 'r') as fi:
            for l in fi:
                a, b = l.rstrip().split(":")
                if "." in b:
                    b = b.split(".")[0]
                out_di[name][a] = b
    return out_di


def calculate_additional_stats(out_di):
    '''
    For each sample input, calculates trimming and mapping efficencies.
    Loads the calculated efficency into the output dictionary.
    '''
    for name in out_di.keys():
        tot_r = float(out_di[name]["TotReads"])
        trim_r = float(out_di[name]["TrimmReads"])
        short_r = float(out_di[name]["ShortReads"])
        EM_r = float(out_di[name]["EMhits"])
        Emiss_r = float(out_di[name]["EMmiss"])
        map_r = float(out_di[name]["Mapped"])
        miRmap_r = float(out_di[name]["miRMapped"])
        tRNAmap_r = float(out_di[name]["tRNAMapped"])
        yRNAmap_r = float(out_di[name]["yRNAMapped"])
        out_di[name]["TrimReadPer"] = "{0:.2f}".format(trim_r / tot_r * 100)
        out_di[name]["ShortReadPer"] = "{0:.2f}".format(short_r / tot_r * 100)
        out_di[name]["EMReadPer"] = "{0:.2f}".format(EM_r / trim_r * 100)
        out_di[name]["EMissReadPer"] = "{0:.2f}".format(Emiss_r / trim_r * 100)
        out_di[name]["MapReadPer"] = "{0:.2f}".format(map_r / trim_r * 100)
        out_di[name]["miRMapPer"] = "{0:.2f}".format(miRmap_r / map_r * 100)
        out_di[name]["tRNAMapPer"] = "{0:.2f}".format(tRNAmap_r / map_r * 100)
        out_di[name]["yRNAMapPer"] = "{0:.2f}".format(yRNAmap_r / map_r * 100)
    return out_di


def output_line_headers():
    '''
    Make list of list of out_di keys and line header values
    '''
    return [["File","File"], 
            ["TotReads","Total Reads"], 
            ["TrimmReads","Trimmed Reads"],
            ["TrimReadPer","Percent Trimmed Reads"],
            ["ShortReads","Too Short Reads"],
            ["ShortReadPer","Percent Too Short"],
            ["EMhits","Exact Match Reads"],
            ["EMReadPer","Percent Exact Matches"],
            ["EMmiss","Mismatch Reads"],
            ["EMissReadPer","Percent Mismatched"],
            ["Mapped","Mapped Reads"],
            ["MapReadPer","Percent Mapped"],
            ["miRMapped","miR Mapped Reads"],
            ["miRMapPer","Percent miR Mapped"],
            ["tRNAMapped","tRNA Mapped Reads"],
            ["tRNAMapPer","Percent tRNA Mapped"],
            ["yRNAMapped","yRNA Mapped Reads"],
            ["yRNAMapPer","Percent yRNA Mapped"]]
        
        
def write_mapping_file(out_di, out_dir, line_head_li):
    '''
    Writes the output dictionary to a file, called MappingInfoTable.csv
    '''
    with open('{}/Mapping_Statistics.csv'.format(out_dir), 'w') as f:
        header = [h[1] for h in line_head_li]
        key_li = [h[0] for h in line_head_li]
        print header
        f.write('Sample_name,{}\n'.format(','.join(header)))
        for sample in sorted(out_di):
            f.write('{}'.format(sample))
            for item in key_li:
                f.write(',{}'.format(out_di[sample][item]))
            f.write('\n')


def main(outPath, samples):
    samples = f_utils.set_path_to_files_glob(samples, 'stats')
    out_di = mapping_stats_dict(samples)
    out_di = calculate_additional_stats(out_di)
    line_head_li = output_line_headers()
    write_mapping_file(out_di, outPath, line_head_li)


if __name__ == '__main__':
    f_utils.check_for_input(sys.argv, usage)
    parser = argparse.ArgumentParser(
             description=usage)
    parser.add_argument(
             'outPath', 
             action='store', 
             help='Path to where the output file will be located')
    parser.add_argument(
             'samples', 
             nargs = '+',
             action='store', 
             help='Path to where the sample output folders are located')
    arg = parser.parse_args()
    main(arg.outPath, arg.samples)
