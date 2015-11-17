ChIP-seq Pipeline Information
==========================

Requirements:
------------
  * Python (2.7+)
  * [SICER][SICER]
  * [MACS2][MACS2] (script coming soon)
  * [Trimmomatic][trim]
  * [STAR][STAR]
  * [phantompeakqualtools][phantompeakqualtools]

Usage:
------
In order the pipeline should be run:
  1. trim_adapters
  2. make_indices
  3. map_reads
  4. run_sicer or run_sicer-rb

Each of the python files has its own usage on the command line with the -h option.

  * Make sure to change any default paths in each of the scripts (or include the optional commandline arguments).
      * For example: Paths to trimmomatic, java command, and fasta adapters in trim_adapters.py
  * make_indices only needs to be run once

Wrapper scripts:
------------------------
  * Python wrapper script: run_chip-seq.py
  * Python to shell wrapper script: run_chip-seq_bash.py
  * Both scripts require a config file -- sample provided in the git.
  * Check both scripts configuration parameters with the -h (python run_chip-seq.py -h)
      * Python wrapper script is a bit more polished.
      * Shell can be output to either bash (piped) or to a file to see which commands are run in what order.
      * The scripts use a method of glob discovery to find files.
          * All of the extension requirements can be changed within the wrapper script in the glob.glob settings.
          * trim_adapters: Files must be in the same directory and have .fastq as their extension
          * map_reads: Files must be in the same directory and have .fastq.filtered (output of trim_adapters) as their extension
          * bed: Files must be in separate subdirectories under one master directory and have .sortedByCoord.out.bam (output of map_reads) as their extension
          * bed_chr: Files must be in the same directory and have .bed (output of bed) as their extension
          * sicer or sicer-rb: Files must be in the same directory and have .sicer.bed (output of bed_chr) as their extension.
              * If using controls or samples make sure sample files have "sample" and "RepXX" (i.e. Rep03) and control files have "control" and "RepXX" (i.e. Rep03 for the same control for the same sample.) Replicate numbers must match for a sample and control to be run together otherwise nothing is run.
  

Individual Steps:
-----------------

#### QC Steps ####
  * After mapping reads to genome, check the cross-correlation plot for each bam file with phantompeakqualtools.
  * `Rscript run_spp.R -c=$bamfile -savp -out=$bamfile.spp_out`
  * $bamfile is the path to your bamfile.
  * Run multiple at a time with parallel or a variation of the bash for loop given below in SICER steps.
  * There will be two output files with run_spp.R: 1) text file with statistics 2) pdf file with cross-correlation plot

#### 1. trim_adapters ####
  * Make sure to provide path to Trimmomatic, java command, and fasta adapters in command or edit script.
  * The default java heapsizes are (8g, 15g) change if necessary.

#### 2. make_indices ####
  * Only necessary if STAR indices are not already made.
  * Ability to use GFF or GTF (-usegtf and pass path to gtf with -gff)

#### 3. map_reads  ####
  * Provide STAR indices directory (from previous step or if already made)
  * Relative paths are allowed.
  * Creates a directory within input directory with out_file_dir as name and basepre as prefix for STAR output files.
      * (directory of input_file)/out_file_dir/basepre/basepre

#### 4. run_sicer/run_sicer-rb ####
  * Provide SICER.sh or SICER-rb.sh filepath in script or in command.
  * SICER requires Bed files to run and chromosomes to be named chr1-X.
      * Use `bam2bed < bam_file > bed_output` to convert.
      * Convert several bam files to bed:
          * `for x in $bam/*.bam;do y=$(basename ${x%.*});bam2bed < $x > $output/$y.bed;done`
              * One at a time.
          * `for x in $bam/*.bam;do y=$(basename ${x%.*});bam2bed < $x > $output/$y.bed & done;wait;echo "Converting to Bed done"`
              * All at once.
          * $bam is the bam file directory for all bams, $output is output directory, $y is the bam trading .bam for .bed extension
      * To convert to several bed files to chr[1-X]:
          *  `sed -- 's/$original/chr/g' $infile > $outfile`
          *  $original is your non-chr chromosome name prefix i.e. Bd for Bd1-5 or scaffold for scaffold12345
      * Options are self-explanatory.
      * run_sicer is to be used with control files & input files.
      * run_sicer-rb is to be used only with input files.
          * Check the SICER documentation for more information on sicer-rb/sicer.

[SICER]: http://home.gwu.edu/~wpeng/Software.htm "SICER"

[trim]: http://www.usadellab.org/cms/?page=trimmomatic "Trimmomatic"

[STAR]: https://github.com/alexdobin/STAR "STAR"

[MACS2]: https://github.com/taoliu/MACS/ "MACS2"

[phantompeakqualtools]: https://code.google.com/p/phantompeakqualtools/ "Phantom Peak Qualtools"


