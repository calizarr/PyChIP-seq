#!/usr/bin/python
from __future__ import print_function
import argparse
import subprocess
import os


def options():
    parser = argparse.ArgumentParser(description="Map ChIP/ATAC-seq reads with STAR",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("nThreads", help="Number of threads to use for STAR")
    parser.add_argument("input_file", help="Read file to be mapped with STAR")
    parser.add_argument("out_file_dir", help="Directory for Output Directories")
    parser.add_argument("basepre", help="Base prefix for output file names in Output Directory")
    parser.add_argument("-st",
                        metavar="Star Indices",
                        default = "/shares/tmockler_share/clizarraga/Brachypodium_distachyon/Phytozome/v3.1/assembly/STAR_indices_masked",
                        help="Directory for STAR indices")
    args = parser.parse_args()
    return args


def main():
    args = options()
    # Getting the absolute path to file if relative path given.
    args.input_file = os.path.abspath(args.input_file)
    # Getting the directory name
    dirname = os.path.dirname(args.input_file)
    # Getting just the filename
    basename = os.path.basename(args.input_file)
    if not args.basepre:
        basepre = '.'.join(basename.split('.')[:2])
    else:
        basepre = args.basepre
    output_base = os.path.join(dirname, args.out_file_dir, basepre)
    print("Attempting to create directory base: {0}".format(output_base))
    try:
        os.makedirs(output_base)
    except OSError:
        if not os.path.isdir(output_base):
            raise
    # Creating the output file path.
    output_file = os.path.join(output_base, basepre+".")
    # Generate Indices
    cmd = "STAR --runThreadN {0} --genomeDir {1} --readFilesIn {2} --outFileNamePrefix {3} --alignIntronMax 1 --alignEndsType EndToEnd --quantMode TranscriptomeSAM --outSAMtype BAM SortedByCoordinate" \
          .format(args.nThreads, args.st, args.input_file, output_file)
    print("Running cmd: ")
    print(cmd)
    subprocess.call(cmd, shell=True)

if __name__ == "__main__":
    main()
