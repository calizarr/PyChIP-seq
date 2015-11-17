#!/usr/bin/python
from __future__ import print_function
import argparse
import subprocess
import os
import sys


def options():
    parser = argparse.ArgumentParser(description="Map ChIP/ATAC-seq reads with STAR",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-nThreads", help="Number of threads to use for STAR", default="8")
    parser.add_argument("-input_file", help="Read file to be mapped with STAR")
    parser.add_argument("-out_file_dir", help="Directory for Output Directories")
    parser.add_argument("-basepre", help="Last n characters to cut from filename for prefix", default=15)
    parser.add_argument("-st",
                        metavar="Star Indices",
                        default = "/shares/tmockler_share/clizarraga/Brachypodium_distachyon/Phytozome/v3.1/assembly/STAR_indices_masked",
                        help="Directory for STAR indices")
    parser.add_argument("-config",
                        default = None,
                        help = "Path to config file if used.")
    args = parser.parse_args()
    return args


def main():
    args = options()
    if args.config:
        if sys.version_info[0] < 3:
            import ConfigParser
            Config = ConfigParser.ConfigParser()
        else:
            import configparser
            Config = configparser.ConfigParser()
        Config.read(args.config)
        nThreads = Config.get("ALL", "nThreads")
        input_file = args.input_file
        out_file_dir = Config.get("ALL", "out_file_dir")
        basecut = Config.getint("STAR", "cut")
        basepre = os.path.basename(input_file[:-basecut])
        st = Config.get("STAR", "star_indices")
    elif not args.input_file:
        sys.exit("Provide input_file please -input_file [filepath] \nCheck script defaults as well( -h).")
    elif not args.out_file_dir:
        sys.exit("Provide output_file_dir please -out_file_dir [dirpath] \nCheck script defaults as well.")
    else:
        nThreads = args.nThreads
        input_file = args.input_file
        out_file_dir = args.out_file_dir
        basecut = args.basepre
        basepre = os.path.basename(input_file[:-basecut])
        st = args.st
    # Getting the absolute path to file if relative path given.
    input_file = os.path.abspath(input_file)
    # Getting the directory name
    # dirname = os.path.dirname(input_file)
    # Getting just the filename
    basename = os.path.basename(input_file)
    if not basepre:
        basepre = basename[:-basecut]
    else:
        basepre = basepre
    output_base = os.path.join(out_file_dir, basepre)
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
          .format(nThreads, st, input_file, output_file)
    print("Running cmd: ")
    print(cmd)
    subprocess.call(cmd, shell=True)

if __name__ == "__main__":
    main()
