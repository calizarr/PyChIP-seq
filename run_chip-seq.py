#!env/bin/python
from __future__ import print_function
import argparse
import subprocess
import os
import sys
import datetime
import glob
import pdb

def options():
    parser = argparse.ArgumentParser(description="Create bash script to run ChIP-seq pipeline",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-trim", action="store_true", help="Run Trim Adapters step")
    parser.add_argument("-indices", action="store_true", help="Make STAR indices step")
    parser.add_argument("-map_reads", action="store_true", help="Map reads to reference")
    parser.add_argument("-sicer_rb", action="store_true", help="Run SICER without controls")
    parser.add_argument("-sicer", action="store_true", help="Run SICER with controls")
    parser.add_argument("-config", default = None, help="Path to config file")
    parser.add_argument("-outdir", help="Path to output directory (log files")
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
    cmdlist = []
    now = datetime.datetime.now().strftime("%Y-%m-%d")
    log_out = "{0}_ChIP-seq.out.log".format(now)
    log_err = "{0}_ChIP-seq.err.log".format(now)
    # logcommand = "(time {cmd} | tee -a $outdir/$LOG_OUT) 3>&1 1>&2 2>&3 | tee -a $outdir/$LOG_ERR"
    logcommand = "(time {insert} | tee -a {0}/{1}) 3>&1 1>&2 2>&3 | tee -a {0}/{2}".format(args.outdir, log_out, log_err, insert="{insert}")
    if args.trim:
        # trim_adapters.py nThreads cutoff_qual minlen phred input_file -fa -trim -jc -mih -mah
        # More info available by running python trim_adapters.py -h
        input_dir = Config.get("ALL", "fastq_dir")
        files = glob.glob(os.path.join(input_dir,"*.fastq"))
        # files = Config.get("TRIM", "input_file").split(',')
        # files = [f.strip() for f in files]
        for f in files:
            cmd = "python trim_adapters.py -input_file {0} -config {1}".format(f, args.config)
            # cmd = logcommand.format(insert=cmd)
            cmdlist.append(cmd)
    if args.indices:
        # make_indices.py nThreads read_length -st -mg -gff -usegtf
        # More info available by running python make_indices.py -h
        cmd = "python make_indices.py -config {0}".format(args.config)
        # cmd = logcommand.format(insert=cmd)
        cmdlist.append(cmd)
    if args.map_reads:
        # map_reads.py nThreads input_file out_file_dir basepre -st
        # More info available by running python map_reads.py -h
        map_dir = input_dir+"Filtered/"
        files = glob.glob(os.path.join(map_dir,"*.filtered"))
        for f in files:
            cmd  = "python map_reads.py -input_file {0} -config {1}".format(f, args.config)
            # cmd = logcommand.format(insert=cmd)
            cmdlist.append(cmd)
    if args.sicer_rb or args.sicer:
        bam_dir = os.path.join(Config.get("ALL", "map_dir"), "STAR")
        files = glob.glob(os.path.join(bam_dir,"*/*Aligned.out.bam"))
        
    if args.sicer_rb:
        # run_sicer-rb.py InputDir bed_file OutputDir --species --rdthresh --winsize --fragsize --egf --gap_size --eval --sicer
        # More info available by running python run_sicer-rb.py -h
        cmd = "python run_sicer-rb.py -config {0}".format(args.config)
        # cmd = logcommand.format(insert=cmd)
        cmdlist.append(cmd)

    print('\n'.join(cmdlist))

if __name__ == "__main__":
    main()
