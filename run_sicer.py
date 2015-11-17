#!/usr/bin/python
from __future__ import print_function
import argparse
import subprocess
import os
import sys


def options():
    parser = argparse.ArgumentParser(description="Run SICER on Treatment vs Control",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-InputDir", help="Directory where input files are located")
    parser.add_argument("-sample", help="Treatment bed file for sample size")
    parser.add_argument("-control", help="Control bed file for sample size")
    parser.add_argument("-OutputDir", help="Directory where output files go")
    parser.add_argument("--species", dest="species", help="Species in SICER", default="bdist")
    parser.add_argument("--rdthresh", dest="rdthresh", help="SICER redundancy threshold", default="1")
    parser.add_argument("--winsize", dest="winsize", help="SICER window size", default="200")
    parser.add_argument("--fragsize", dest="fragsize", help="SICER fragment size", default="150")
    # Maybe efg for Brachy should be: 0.84 or 0.98 depending on how you look at it.
    parser.add_argument("--egf", dest="egf", help="SICER Effective Genome Fraction", default="0.95")
    parser.add_argument("--gap_size", dest="gap_size", help="SICER gap size", default="600")
    parser.add_argument("--FDR", dest="FDR", help="SICER False Discovery Rate", default="1E-2")
    parser.add_argument("--sicer",
                        default = "/shares/tmockler_share/clizarraga/usr/local/SICER_V1.1/SICER/SICER.sh",
                        help = "SICER.sh absolute path")
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
        InputDir = Config.get("ALL", "bed_chr_dir")
        OutputDir = Config.get("ALL", "sicer_dir")
        sample = args.sample
        control = args.control
        species = Config.get("SICER", "species")
        rdthresh = Config.get("SICER", "rdthresh")
        winsize = Config.get("SICER", "winsize")
        fragsize = Config.get("SICER", "fragsize")
        egf = Config.get("SICER", "egf")
        gap_size = Config.get("SICER", "gap_size")
        FDR = Config.get("SICER", "FDR")
        sicer = Config.get("SICER", "sicer")
    elif not args.sample:
        sys.exit("Provide sample please -sample [filepath] \nCheck script defaults as well( -h).")
    elif not args.control:
        sys.exit("Provide control please -control [filepath] \nCheck script defaults as well( -h).")
    elif not args.InputDir:
        sys.exit("Provide directory path for bed files -InputDir [dirpath] \nCheck script defaults as well( -h).")
    elif not args.OutputDir:
        sys.exit("Provide directory path for output files -OutputDir [dirpath] \nCheck script defaults as well( -h).")
    else:
        InputDir = args.InputDir
        OutputDir = args.OutputDir
        sample = args.sample
        species = args.species
        rdthresh = args.rdthresh
        winsize = args.winsize
        fragsize = args.fragsize
        egf = args.egf
        gap_size = args.gap_size
        FDR = args.FDR
        sicer = args.sicer
    # Getting the absolute path to file if relative path given.
    # args.input_file = os.path.abspath(args.sample)
    # args.control = os.path.abspath(args.control)
    # Getting just the filename
    OutputDir = os.path.abspath(OutputDir)
    InputDir = os.path.abspath(InputDir)
    # sample and control should be basename no path.
    sample = os.path.basename(sample)
    control = os.path.basename(control)
    # Creating the output file path.
    print("Attempting to make base output directory: {0}".format(OutputDir))
    try:
        os.makedirs(OutputDir)
    except OSError:
        if not os.path.isdir(OutputDir):
            raise
    # Making base prefix for specific output directory.
    basepre = '.'.join(sample.split('.')[:-5])
    spec_OutputDir = os.path.join(OutputDir, basepre)
    print("Attempting to make specific file output directory: {0}".format(spec_OutputDir))
    try:
        os.makedirs(spec_OutputDir)
    except OSError:
        if not os.path.isdir(spec_OutputDir):
            raise
    # Generate Indices
    # sh DIR/SICER.sh ["InputDir"] ["bed file"] ["control file"] ["OutputDir"]
    # ["Species"] ["redundancy threshold"] ["window size (bp)"] ["fragment size"]
    # ["effective genome fraction"] ["gap size (bp)"] ["FDR"]
    cmd = "sh {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11}" \
          .format(sicer, InputDir, sample, control, spec_OutputDir, species,
                  rdthresh, winsize, fragsize, egf, gap_size, FDR)
    print("Running cmd: ")
    print(cmd)
    subprocess.call(cmd, shell=True)

if __name__ == "__main__":
    main()
