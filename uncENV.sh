# set path
#
export PATH=/proj/.test/roach/miRNA/bin:./:$PATH
export PYTHONPATH=/proj/.test/roach/miRNA/lib/python/
export PATH=$PATH:/proj/.test/roach/miRNA/SHRiMP_2_2_2/bin
#export SHRIMP_FOLDER=/proj/.test/roach/miRNA/SHRiMP_2_2_2

#export JB_GENOME_DIR=/proj/seth_lab/projects/genome/


module unload bedtools
module unload perl 
module unload git
module unload python

module load r/3.2.2
module load bowtie/1.1.0
module load bedtools/2.25.0
module load bamtools/1.0.2
module load perl/5.12.0
module load python/2.7.6

