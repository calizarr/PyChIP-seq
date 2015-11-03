#!/usr/bin/python
from __future__ import print_function
import argparse
import subprocess
import os


def options():
    parser = argparse.ArgumentParser(description="Trim ChiP/ATAC-seq reads",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("nThreads", help="Number of threads to use for Trimmomatic")
    parser.add_argument("cutoff_qual", help="Quality for nucleotide cutoff from head and tail")
    parser.add_argument("minlen", help="Minimun length of read before drop off")
    parser.add_argument("phred", help="Phred scores of the reads")
    parser.add_argument("input_file", help="Read file for adapter trimming and filtering")
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
    args = parser.parse_args()
    return args


def main():
    args = options()
    # Getting the absolute path to the file if a relative path was given (../file.ext)
    args.input_file = os.path.abspath(args.input_file)
    # Getting just the directory name
    dirname = os.path.dirname(args.input_file)
    print("Dirname is: {0}".format(dirname))
    # Getting just the file name
    basename = os.path.basename(args.input_file)
    print("Basename is: {0}".format(basename))
    # Creating the output filename.
    output_file = os.path.join(dirname, "Filtered", basename + ".filtered")
    call = "{0} -Xms8g -Xmx15g -XX:+UseG1GC -XX:+UseStringDeduplication -jar {1}".format(args.java_command, args.trimmomatic)
    cmd = "{0} SE -threads {1} -phred{2} {3} {4} ILLUMINACLIP:{5}:2:30:10 LEADING:{6} TRAILING:{7} MINLEN:{8}" \
          .format(call, args.nThreads, args.phred, args.input_file, output_file, args.fasta_adapters,
                  args.cutoff_qual, args.cutoff_qual, args.minlen)
    print("Running Trimmommatic with this command:")
    print(cmd, "\n")
    subprocess.call(cmd, shell=True)

if __name__ == "__main__":
    main()


