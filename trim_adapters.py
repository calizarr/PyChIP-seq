#!/usr/bin/python
from __future__ import print_function
import argparse
import subprocess
import os
import sys


def options():
    parser = argparse.ArgumentParser(description="Trim ChiP/ATAC-seq reads",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-nThreads", help="Number of threads to use for Trimmomatic", default="8")
    parser.add_argument("-cutoff_qual", help="Quality for nucleotide cutoff from head and tail", default="30")
    parser.add_argument("-minlen", help="Minimun length of read before drop off", default="30")
    parser.add_argument("-phred", help="Phred scores of the reads", default="33")
    parser.add_argument("-input_file", help="Read file for adapter trimming and filtering")
    parser.add_argument("-fa", "--fasta_adapters",
                        default = "/nfs4shares/bioinfosw/installs_current/Trimmomatic-0.32/adapters/TruSeq3-SE.fa",
                        help="Adapters file for trimmming adapters")
    parser.add_argument("-trim", "--trimmomatic",
                        default = "/nfs4shares/bioinfosw/installs_current/Trimmomatic-0.32/trimmomatic-0.32.jar",
                        help="Path to Trimmomatic jar file.")
    parser.add_argument("-jc", "--java_command",
                        default = "java",
                        help="Name of shell java command or path")
    parser.add_argument("-mih", "--min_heap",
                        default = "8g",
                        help="Minimum size of java command RAM usage")
    parser.add_argument("-mah", "--max_heap",
                        default = "15g",
                        help="Maximum size of java command RAM usage")
    parser.add_argument("-config",
                        default = None,
                        help = "Path to config file if used.")
    args = parser.parse_args()
    return args


def main():
    args = options()
    # Either using a Config file or straight command line parameters.
    if args.config:
        if sys.version_info[0] < 3:
            import ConfigParser
            Config = ConfigParser.ConfigParser()
        else:
            import configparser
            Config = configparser.ConfigParser()
        Config.read(args.config)
        nThreads = Config.get("ALL", "nThreads")
        phred = Config.get("TRIM", "phred")
        input_file = args.input_file
        fasta_adapters = Config.get("TRIM", "fasta_adapters")
        cutoff_qual = Config.get("TRIM", "cutoff")
        minlen = Config.get("TRIM", "minlen")
        java_command = Config.get("JAVA", "java_command")
        min_heap = Config.get("JAVA", "min_heap")
        max_heap = Config.get("JAVA", "max_heap")
        trimmomatic = Config.get("TRIM", "trimmomatic")
    elif not args.input_file:
        sys.exit("Provide input_file please -input_file [filepath] \nCheck script defaults as well( -h).")
    else:
        nThreads = args.nThreads
        phred = args.phred
        input_file = args.input_file
        fasta_adapters = args.fasta_adapters
        cutoff_qual = args.cutoff_qual
        minlen = args.minlen
        java_command = args.java_command
        min_heap = args.min_heap
        max_heap = args.max_heap
        trimmomatic = args.trimmomatic
    # Getting the absolute path to the file if a relative path was given (../file.ext)
    input_file = os.path.abspath(input_file)
    # Getting just the directory name
    dirname = os.path.dirname(input_file)
    print("Dirname is: {0}".format(dirname))
    # Getting just the file name
    basename = os.path.basename(input_file)
    print("Basename is: {0}".format(basename))
    output_base = os.path.join(dirname, "Filtered")
    print("Attempting to create directory base: {0}".format(output_base))
    try:
        os.makedirs(output_base)
    except OSError:
        if not os.path.isdir(output_base):
            raise
    # Creating the output filename.
    output_file = os.path.join(output_base, basename + ".filtered")
    call = "{0} -Xms{1} -Xmx{2} -XX:+UseG1GC -XX:+UseStringDeduplication -jar {3}".format(java_command, min_heap, max_heap, trimmomatic)
    cmd = "{0} SE -threads {1} -phred{2} {3} {4} ILLUMINACLIP:{5}:2:30:10 LEADING:{6} TRAILING:{7} MINLEN:{8}" \
          .format(call, nThreads, phred, input_file, output_file, fasta_adapters,
                  cutoff_qual, cutoff_qual, minlen)
    print("Running Trimmommatic with this command:")
    print(cmd, "\n")
    # subprocess.call(cmd, shell=True)

if __name__ == "__main__":
    main()


