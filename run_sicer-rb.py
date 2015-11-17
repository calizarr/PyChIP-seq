#!/usr/bin/python
from __future__ import print_function
import argparse
import subprocess
import os
import sys


def options():
    parser = argparse.ArgumentParser(description="Run SICER on Treatment vs Random Background Noise",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-InputDir", help="Directory where input files are located")
    parser.add_argument("-bed_file", help="Treatment bed file for sample size")
    parser.add_argument("-OutputDir", help="Directory where output files go")
    parser.add_argument("--species", dest="species", help="Species in SICER", default="bdist")
    parser.add_argument("--rdthresh", dest="rdthresh", help="SICER redundancy threshold", default="1")
    parser.add_argument("--winsize", dest="winsize", help="SICER window size", default="200")
    parser.add_argument("--fragsize", dest="fragsize", help="SICER fragment size", default="150")
    # Maybe efg for Brachy should be: 0.84 or 0.98 depending on how you look at it.
    parser.add_argument("--egf", dest="egf", help="SICER Effective Genome Fraction", default="0.95")
    parser.add_argument("--gap_size", dest="gap_size", help="SICER gap size", default="600")
    parser.add_argument("--evalue", dest="evalue", help="SICER RB E-value", default="100")
    parser.add_argument("--sicer",
                        default = "/shares/tmockler_share/clizarraga/usr/local/SICER_V1.1/SICER/SICER-rb.sh",
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
        OutputDir = Config.get("ALL", "sicer_rb_dir")
        bed_file = args.bed_file
        species = Config.get("SICER", "species")
        rdthresh = Config.get("SICER", "rdthresh")
        winsize = Config.get("SICER", "winsize")
        fragsize = Config.get("SICER", "fragsize")
        egf = Config.get("SICER", "egf")
        gap_size = Config.get("SICER", "gap_size")
        evalue = Config.get("SICER", "evalue")
        sicer = Config.get("SICER", "sicer_rb")
    elif not args.bed_file:
        sys.exit("Provide bed_file please -bed_file [filepath] \nCheck script defaults as well( -h).")
    elif not args.InputDir:
        sys.exit("Provide directory path for bed files -InputDir [dirpath] \nCheck script defaults as well( -h).")
    elif not args.OutputDir:
        sys.exit("Provide directory path for output files -OutputDir [dirpath] \nCheck script defaults as well( -h).")
    else:
        InputDir = args.InputDir
        OutputDir = args.OutputDir
        bed_file = args.bed_file
        species = args.species
        rdthresh = args.rdthresh
        winsize = args.winsize
        fragsize = args.fragsize
        egf = args.egf
        gap_size = args.gap_size
        evalue = args.evalue
        sicer = args.sicer
    # Getting the absolute path to file if relative path given.
    # args.input_file = os.path.abspath(args.bed_file)
    # args.control_file = os.path.abspath(args.control_file)
    # Getting just the filename
    OutputDir = os.path.abspath(OutputDir)
    InputDir = os.path.abspath(InputDir)
    # bed_file should only be file basename no path.
    bed_file = os.path.basename(bed_file)
    # if "index" in bed_file:
    #     pre = bed_file.split('_')
    #     basepre = '_'.join(pre[:5])
    # else:
    #     pre = bed_file.split('-')
    #     basepre = '-'.join(pre[:6])
    # pre = pre[-1].split('.')[0]
    # basepre = '-'.join([basepre, pre])+"-rb"
    print("Attempting to make base output directory: {0}".format(OutputDir))
    try:
        os.makedirs(OutputDir)
    except OSError:
        if not os.path.isdir(OutputDir):
            raise
    # Making base prefix for specific output directory.
    basepre = bed_file.split('.')[:-5]
    basepre[-1] = basepre[-1]+"-rb"
    basepre = '.'.join(basepre)
    spec_OutputDir = os.path.join(OutputDir, basepre)
    print("Attempting to make specific file output directory: {0}".format(spec_OutputDir))
    try:
        os.makedirs(spec_OutputDir)
    except OSError:
        if not os.path.isdir(spec_OutputDir):
            raise
    # Creating the output file path.
    # sh DIR/SICER.sh ["InputDir"] ["bed file"] ["control file"] ["OutputDir"]
    # ["Species"] ["redundancy threshold"] ["window size (bp)"] ["fragment size"]
    # ["effective genome fraction"] ["gap size (bp)"] ["FDR"]
    cmd = "sh {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10}" \
          .format(sicer, InputDir, bed_file, spec_OutputDir, species,
                  rdthresh, winsize, fragsize, egf, gap_size, evalue)
    print("Running cmd: ")
    print(cmd)
    subprocess.call(cmd, shell=True)

if __name__ == "__main__":
    main()
