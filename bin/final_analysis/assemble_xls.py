#!/usr/bin/python2

import sys
import glob
import xlsxwriter

def verify_outputs_exist(dir, files):
    '''
    Verify the necessary output files exist in directory.
    If not, state which files are missing.
    '''
    req_files = ['Mapping_Statistics.csv', 
                 'length_distribution.csv',
                 'sample_correlation_values.csv',
                 'RPMMM_all.csv',
                 'RPMMM_mirs_over_50.csv',
                 'statistics.csv']

    file_exist = []
    for f in req_files:
        if f in files:
            file_exist.append(f)
        else:
            print 'Warning: Missing output file -> {}'.format(f)
    vs_files = glob.glob('*vs*')
    if len(vs_files) == 0:
        print '\nNo statistic comparison files found'
    return file_exist, vs_files


def insert_image(file, wb):
    '''
    Insert image in file at indicated location.
    '''
    name = file.split('.')[0].replace('_', ' ')
    image = file.replace('.csv', '.png')
    try:
        ws = wb.get_worksheet_by_name(name)
        with open(file) as f:
            for i, l in enumerate(f):
                pass
        ws.insert_image('A{}'.format(i + 3), image)
        return wb
    except AttributeError:
        return wb
    except IOError:
        return wb


def add_comparison_sheets(wb, vs_files):
    '''
    For comparison files, insert into workbook as own sheet.
    Highlight cells with p-value <= .05 as yellow.
    '''
    if len(vs_files) > 0:
        format1 = wb.add_format({'bg_color': '#FFFF00'})
        for fi in vs_files:
            name = fi.split('.')[0]
            ws = wb.add_worksheet(name)
            with open(fi) as f:
                li = [l.split(',') for l in f.read().split('\n') if l]
            for i, row in enumerate(li):
                ws.write_row(i, 0, row)
            ws.conditional_format('C2:C{}'.format(i), {'type': 'cell', 
                                                       'criteria': '<=', 
                                                       'value': .05, 
                                                       'format': format1})
        return wb
    else:
        return wb


def create_workbook(files, vs_files):
    '''
    Create excel workbook with sheets corresponding to output files.
    Insert images from outputs in proper location.
    '''
    wb = xlsxwriter.Workbook('miRquant2.0_results.xlsx', {'strings_to_numbers' : True})
    for f in files:
        name = f.split('.')[0].replace('_', ' ')
        ws = wb.add_worksheet(name)
        with open(f) as fi:
            li = [l.split(',') for l in fi.read().split('\n') if l]
        for i, row in enumerate(li):
            ws.write_row(i, 0, row)
    wb = insert_image('length_distribution.csv', wb)
    wb = insert_image('sample_correlation_values.csv', wb)
    wb = add_comparison_sheets(wb, vs_files)
    wb.close()


def main():
    files, vs_files = verify_outputs_exist('.', glob.glob('*'))
    create_workbook(files, vs_files)


if __name__ == '__main__':
    main()

