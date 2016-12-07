#miRquant               
Last update to README: 12/6/16

## miRquant setup
####Check out a copy of the smRNA pipeline code:
If you don’t already have a directory:

```
$ mkdir /path/to/miRquant
$ cd /path/to/miRquant
$ module load git
$ git clone https://github.com/Sethupathy-Lab/miRquant.git
```

Now you will have a miRquant directory containing all the miRquant scripts.
All of the code is run from this directory!

####Install required programs

cutadapt v1.0  
bedtools v2.25.0  
bowtie v1.1.0  
python v2.7.6  
SHRiMP v2.2.2  
R v3.2.2  

####Download relevant genome fasta files and generate Bowtie indexes

miRquant is currently set up to work with human, mouse and rat, with fruitfly support coming.

The specific genome releases used in miRquant are:

human - [hg19](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/)  
mouse - [mm9](http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/)  
rat - [rn4](http://hgdownload.cse.ucsc.edu/goldenPath/rn4/bigZips/)  

Download the appropriate genomes, and the chromosome sizes (<release>.chrom.sizes) 

Change the genome fasta name to [prefix].fa and the chromosome sizes file to [prefix].chromSizes.  The prefixes for each species is as follows:

human - hg19  
mouse - mm9  
rat - rn4  

Generate genome indexes.


##Running miRQuant

####Load proper modules and environmental variables:
```
$ cd /proj/seth_lab/users/ONYEN/miRquant
$ source uncENV.sh
```

####Enter run parameters into the miRquant configuration file

Copy the configuration directory to the directory containing the small RNA-seq fastqs.
```
cp -r /path/to/miRquant/configuration /path/to/fastq_containing_directory
```
The configuration directory contains two configuration files.

1. conf_miRquant.yml
  * Configuration that will be edited on a project by project basis
2. conf_system.yml
  * Configuration file for the cluster you are working on, currently filled out for lsf job scheduler.

The miRquant configuration file (conf_miRquant.yml) is as follows:
```
# Directory locations
paths:
    genome:
        /path/to/genome/
    mirquant:
        /path/to/miRquant/
    output:
        /path/to/output/
    resources:
        /path/to/miRquant/bin/resources/
    project:
        /path/to/fastqs/

# Location of necessary files
parameters:
    genome_release:
        mm9           <- prefix for genome release (eg, mouse release 9)
    species:
        mmu           <- species, currently hsa, mmu, or rno
    Minimum Read Length:
        14
# Load in options for cutAdapt
cutadapt:
    adapter:
        'TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACXXXXXXATCTCGTATGCCGTCTTCTGCTTG'
    overlap:
        10
    error:
        1
    Minimum_Read_Length:
        14
# Load in options for bowtie
bowtie:
    quality:
        33
# Load in options for SHRiMP
shrimp:
    path:
        /proj/.test/roach/miRNA/SHRiMP_2_2_2/
    dependencies:
        python:
            /proj/.test/roach/miRNA/lib/python/
    quality:
        33
```

####Run the chain submission script:
From the miRquant directory:
```
$ cd /path/to/miRquant

$ python miRquant.py path/to/configuration/

```

####Once all jobs have finished:
Once the chain submission has finished, check for any errors in the log file.  Each sample will have a multiple output directories setup in the output directory specified in the miRquant configuration file, in the following structure.
```
SAMPLE_NAME
  -logs
  -output
  -temp
```

The logs from each step of miRquant will be in the logs directory.  

In your project directory, there will be a directory for each \<SAMPLE\>.fastq (called \<SAMPLE\>.)

In that directory, there will be an IntermediateFiles subdirectory and a \<SAMPLE\>.stats file.

The \<SAMPLE\>.stats will have statistics to inform on the degree of trimming and aligning, type:
```
$ cat /path/to/fastqs/*/*.stats

file:/path/to/fastqs/<SAMPLE>.fastq
TotReads:100000
TrimmReads:90000
ShortReads:8500
EMhits:40000
EMmiss:45000

TotReads = Total number of reads for this file
TrimmReads = # of reads successfully trimmed of 3’ adapter
ShortReads = # of reads too short after trimming (< 13 )
EMhits =  # of reads with an exact alignment to the genome
EMmiss = # of reads that fail to exactly align to genome
```

####Run the next stage to collect results:
From your pipeline directory (/path/to/miRquant):
```
$ python runC.py path/to/configuration
```

####Run the next stage to generate TAB separated files:
```
$ python process_summary_to_tab.py path/to/configuration
```
After run finishes, you should see:
```
$ cd path/to/fastqs
$ cat */*.stats

file:path/to/fastqs/<SAMPLE>.fastq
TotReads:6149484.00000000000000000000
TrimmReads:3730081.00000000000000000000
ShortReads:1938220.00000000000000000000
EMhits:2509867
EMmiss:1220214
Mapped: 2881057.09082651
miRMapped: 1104952.30513136
```

Mapped and miRMapped indicate the number of reads mapped to the genome and to miRNAs respectively.

The log file above will also contain a table that you can use to put together the mapping stats for your project.

Outputs:
For each Sample:
  TAB_3p_summary.txt       -   3'-end differences
  TAB_3p_summary_miR.txt   -   3'-end differences (miRNA loci only)
  TAB_ed_summary.txt       -   central differences
  TAB_lenDist_summary.txt  -   length differences
  Shrimp_results.bed       -   bed file containing all results

##Final processing
To produce the final summary files, run:
```
$ python final_processing.py path/to/configuration
```
This will produce the mapping statistics, read length distribution, expression correlation heatmap, reads per million mapped (RPMM), and reads per million miRs mapped (RPMMM).

These final outputs will be in the output folder specified in the configuration file, in a directory named year_month_day_miRquant_num, where the year, month, and day refer to the date and the num will correspond to how many times miRquant had been run on that day.
