#!/usr/bin/python2

import sys
import glob
import shutil
import os

def main(chr):
   print chr
   fi_di = {}
   files = glob.glob('{}/*readSize*'.format(chr))
   for f in files:
       try:
           fi_di[f.split('_readSize_')[0]].append(f)
       except KeyError:
           fi_di[f.split('_readSize_')[0]] = [f]
 
 
   for outfi, in_files in fi_di.iteritems():
       with open(outfi, 'w') as fo:
           for fi in in_files:
               with open(fi, 'r') as rf:
                   shutil.copyfileobj(rf, fo)
   
   
   for file in files:
       os.system('rm "{}"'.format(file))


if __name__ == '__main__':
    main(sys.argv[1])
