## miRquant improvements

#### Fixed file locking issue
-Previously the file locking would cause jobs to start and not finish
-Removed file locking and introduced a map / reduce workflow

#### Implementing modular code
-Code was initially in perl
-Code is now in python, with long scripts being broken into smaller functions
--Advantages: 
1. Code is easily readable, with description of the function allowing for easy comprehension
2. Code is re-usable, with each function being able to be called multiple times in various scripts

#### Logging
-Logging for all samples were combined into a single file, and not divided into sections.  Logs were buried in output folders, making it hard to access.
Logging improvements:
1. Logs are output to sample-specific log directories.
2. Logs are cleaner.  Logs clearly presented, with sections of code explicitly stated.
3. Logs are more verbose, and verbosisity is adjustable.

#### Code is easier to run
A configuration file is used to hold run parameters.
The use of a configuration file is useful due to the following reasons.
1. Can easily go back and check what the run parameters were for sample group
2. Parameters don't need to be passed to code at various steps of the pipeline

#### Final output improvements
Final scripts improvements:
1. Created wrapper for final scripts, simplifying assembly of final output.
2. Length distribution is now reported in percentage
3. Length distribution produced as a histogram to easily visual distribution.
4. Sample correlations is visualized in a hierarchical clustered heatmap.
5. Expression is reported as reads per million miRs mapped (RPMMM)

#### Run multiple samples in parallel
Samples were run in sesquential order previously,
now samples can be run at the same time, decreasing overall run time as samples number increases.
