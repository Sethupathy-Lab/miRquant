#miRquant Tutorial
Last update to README: 12/6/16

#miRquant 2.0  

1. Information
2. Installation
  1. Requirements
  2. Setup
3. Steps
  1. Adapter trimming
  2. Alignment
    1. Bowtie alignment
    2. Genomic window construction
    3. SHRiMP alignment
  3. Annotation
  4. Final analysis
4. Output

##Tutorial info
This tutorial is for two mouse samples.  The mouse samples are were prepared by TruSeq.  We will
Brief overview of miRquant here.

##miRquant Installation  
###Requirements
#####Software
miRquant 2.0 can be downloaded as a zip file or cloned from the [miRquant GitHub page](https://github.com/Sethupathy-Lab/miRquant).  

In addition to these scripts, miRquant 2.0 requires the following software for various steps of the pipeline.

* python v2.7.6
* pip 
* bedtools v2.25.0  
* bowtie v1.1.0  
* SHRiMP v2.2.2  
* R v3.2.2 

Install these programs with in your system path.

#####Resources

miRquant is currently set up to work with human, mouse and rat, with fruitfly support coming.

The specific genome releases used in miRquant are:

human - [hg19](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/)  
mouse - [mm9](http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/)  
rat - [rn4](http://hgdownload.cse.ucsc.edu/goldenPath/rn4/bigZips/)  

Download the appropriate genomes and the chromosome sizes for that genome release (\<release\>.chrom.sizes) 

###Setup

Once python/2.7.6 and pip are installed, change to the miRquant directory and type:

```
pip install -r requirements.txt
```

Change the genome fasta name to \<prefix\>.fa and the chromosome sizes file to \<prefix\>.chromSizes.  The prefixes for each species is as follows:

human - hg19  
mouse - mm9  
rat - rn4  

Store the genomes and the chromosome size files *in the same location*.

Build Bowtie genome indexes for each genome.  Information on this can be found in the [Bowtie tutorial](http://bowtie-bio.sourceforge.net/tutorial.shtml).


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

# miRquant parameters
parameters:
    genome_release:
        mm9           <- prefix for genome release (eg, mouse release 9)
    species:
        mmu           <- species, currently hsa, mmu, or rno
# Load in options for cutAdapt
cutadapt:
    adapter:
        'TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACXXXXXXATCTCGTATGCCGTCTTCTGCTTG'           <- TruSeq adapter (Xs denote barcode)
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

*How to fill out the miRquant configuration file:*
* Directory locations
  - genome - path to the genome and chromosome sizes directory
  - mirquant - path to the miRquant program directory
  - output - path to where outputs will be saved.  Directory will be created if it does not exist.
  - resources - path to resources directory.  This should remain in the miRquant bin directory.
  - project - location of the small RNA sequencing results.  All fastqs in location will be processed by miRquant.
* miRquant parameters
  - genome release - appropriate prefix for the species (see prefixes in setup)
  - species - hsa, mmu, or rno for human, mouse, or rat, respectively.
  - Minimum_Read_Length - minimum read length to be included in the analysis, should match cutadapt below
* cutadapt options
  - adapter - 3' adapter sequence; if barcode present, replace with Xs; if degenerate bases present at 5' end, add as Ns
  - overlap - number of overlapping nucleotides for trimming to occur
  - error - error tolerance (1 = 0.1 or 10%)
  - minimum read length - minumum length of trimmed read to be included
* Bowtie options
  - quality - quality cutoff for Bowtie alignment
* SHRiMP options
  - path - path to SHRiMP executables
  - quality - quality cutoff for SHRiMP

For the tutorial, only the paths section of conf_mirquant.yml will have to be altered.

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
The outputs will be saved in the output directory.
Temporary files generated during a miRquant run will be stored in the temp directory.

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
