#miRquant 2.0  

1. Introduction
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

##miRquant introduction  
Brief overview of miRquant here.

##miRquant Installation  
###Requirements
#####Software
miRquant 2.0 can be downloaded as a zip file or cloned from the [miRquant GitHub page](https://github.com/Sethupathy-Lab/miRquant).  

In addition to these scripts, miRquant 2.0 requires the following software for various steps of the pipeline.

* cutadapt v1.0  
* bedtools v2.25.0  
* bowtie v1.1.0  
* python v2.7.6  
* SHRiMP v2.2.2  
* R v3.2.2 

Install these programs and add their locations to the system path.

#####Resources

miRquant is currently set up to work with human, mouse and rat, with fruitfly support coming.

The specific genome releases used in miRquant are:

human - [hg19](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/)  
mouse - [mm9](http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/)  
rat - [rn4](http://hgdownload.cse.ucsc.edu/goldenPath/rn4/bigZips/)  

Download the appropriate genomes and the chromosome sizes for that genome release (\<release\>.chrom.sizes) 

###Setup

Change the genome fasta name to \<prefix\>.fa and the chromosome sizes file to \<prefix\>.chromSizes.  The prefixes for each species is as follows:

human - hg19  
mouse - mm9  
rat - rn4  

Store the genomes and the chromosome size files *in the same location*.

Build Bowtie genome indexes for each genome.  Information on this can be found in the [Bowtie tutorial](http://bowtie-bio.sourceforge.net/tutorial.shtml).

To test that everything was installed correctly, follow the tutorial.

##miRquant steps
###Adapter trimming
###Alignment
#####Bowtie alignment
#####Genomic window construction
#####SHRiMP alignment
###Annotation
###Final analysis

##miRquant output
