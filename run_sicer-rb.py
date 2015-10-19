#!/usr/bin/python
from __future__ import print_function
import argparse
import subprocess
import os


def options():
    parser = argparse.ArgumentParser(description="Run SICER on Treatment vs Random Background Noise")
    parser.add_argument("InputDir", help="Directory where input files are located")
    parser.add_argument("bed_file", help="Treatment bed file for sample size")
    parser.add_argument("OutputDir", help="Directory where output files go")
    parser.add_argument("--species", dest="species", help="Species in SICER", default="bdist")
    parser.add_argument("--rdthresh", dest="rdthresh", help="SICER redundancy threshold", default="1")
    parser.add_argument("--winsize", dest="winsize", help="SICER window size", default="200")
    parser.add_argument("--fragsize", dest="fragsize", help="SICER fragment size", default="150")
    parser.add_argument("--egf", dest="egf", help="SICER Effective Genome Fraction", default="0.95")
    parser.add_argument("--gap_size", dest="gap_size", help="SICER gap size", default="400")
    parser.add_argument("--eval", dest="eval", help="SICER RB E-value", default="100")
    args = parser.parse_args()
    return args


def main():
    args = options()
    print(args)
    # Getting the absolute path to file if relative path given.
    # args.input_file = os.path.abspath(args.bed_file)
    # args.control_file = os.path.abspath(args.control_file)
    # Getting just the filename
    args.OutputDir = os.path.abspath(args.OutputDir)
    args.InputDir = os.path.abspath(args.InputDir)
    if "index" in args.bed_file:
        pre = args.bed_file.split('_')
        basepre = '_'.join(pre[:5])
    else:
        pre = args.bed_file.split('-')
        basepre = '-'.join(pre[:6])
    pre = pre[-1].split('.')[0]
    basepre = '-'.join([basepre, pre])+"-rb"
    OutputDir = os.path.join(args.OutputDir, basepre)
    print(OutputDir)
    try:
        os.makedirs(OutputDir)
    except OSError:
        if not os.path.isdir(OutputDir):
            raise
    # Creating the output file path.
    SICER = "/shares/tmockler_share/clizarraga/usr/local/SICER_V1.1/SICER/SICER-rb.sh"
    # Generate Indices
    # sh DIR/SICER.sh ["InputDir"] ["bed file"] ["control file"] ["OutputDir"]
    # ["Species"] ["redundancy threshold"] ["window size (bp)"] ["fragment size"]
    # ["effective genome fraction"] ["gap size (bp)"] ["FDR"]
    cmd = "sh {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10}" \
          .format(SICER, args.InputDir, args.bed_file, OutputDir, args.species,
                  args.rdthresh, args.winsize, args.fragsize, args.egf, args.gap_size, args.eval)
    print("Running cmd: ")
    print(cmd)
    subprocess.call(cmd, shell=True)

if __name__ == "__main__":
    main()
