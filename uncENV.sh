# set path

export PATH=~/.local/bin:$PATH
export PYTHONPATH=~/.local/lib:$PYTHONPATH

module unload python

module load r/3.2.2
module load bowtie/1.1.0
module load bedtools/2.26.0
module load python/2.7.6
