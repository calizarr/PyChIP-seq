#!env/bin/python
from __future__ import print_function
import argparse
import subprocess
import os
import sys
import datetime
import glob
import pdb
import re

def options():
    parser = argparse.ArgumentParser(description="Create bash script to run ChIP-seq pipeline",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-trim", action="store_true", help="Run Trim Adapters step")
    parser.add_argument("-indices", action="store_true", help="Make STAR indices step")
    parser.add_argument("-map_reads", action="store_true", help="Map reads to reference")
    parser.add_argument("-bed", action="store_true", help="Convert bams to beds for SICER")
    parser.add_argument("-bed_chr", action="store_true", help="Convert beds to bed_chr for SICER")
    parser.add_argument("-sicer_rb", action="store_true", help="Run SICER without controls")
    parser.add_argument("-sicer", action="store_true", help="Run SICER with controls")
    parser.add_argument("-config", default = None, help="Path to config file")
    parser.add_argument("-outdir", help="Path to output directory (log files")
    args = parser.parse_args()
    return args

def main():
    # Getting the command line options and reading config.
    args = options()
    if args.config:
        if sys.version_info[0] < 3:
            import ConfigParser
            Config = ConfigParser.ConfigParser()
        else:
            import configparser
            Config = configparser.ConfigParser()
        Config.read(args.config)
    else:
        sys.exit("Wrapper script needs a Config file to continue. Please provide -config.\nA sample file has been provided in the git.")
    # List of commands to be printed at end.
    cmdlist = []
    # Getting the current date and setting log files.
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
            cmdlist.append(cmd)
            cmd = logcommand.format(insert=cmd)
            print(cmd)
            subprocess.call(cmd, shell=True)
    if args.indices:
        # make_indices.py nThreads read_length -st -mg -gff -usegtf
        # More info available by running python make_indices.py -h
        cmd = "python make_indices.py -config {0}".format(args.config)
        cmdlist.append(cmd)
        cmd = logcommand.format(insert=cmd)
        print(cmd)
        subprocess.call(cmd, shell=True)
    if args.map_reads:
        # map_reads.py nThreads input_file out_file_dir basepre -st
        # More info available by running python map_reads.py -h
        map_dir = Config.get("ALL", "map_dir")
        files = glob.glob(os.path.join(map_dir,"*fastq.filtered"))
        for f in files:
            cmd  = "python map_reads.py -input_file {0} -config {1}".format(f, args.config)
            cmdlist.append(cmd)
            cmd = logcommand.format(insert=cmd)
            print(cmd)
            subprocess.call(cmd, shell=True)
    if args.bed:
        # Converting BAM files to Bed files to Chromosome numbered Bed files for SICER
        print("echo Converting BAM files to Bed files to Chromosome numbered Bed files for SICER")
        bam_dir = Config.get("ALL", "out_file_dir")
        files = glob.glob(os.path.join(bam_dir,"*/*sortedByCoord.out.bam"))
        output_base = os.path.join(bam_dir, "Bed")
        print("echo Attempting to create directory base: {0}".format(output_base))
        try:
            os.makedirs(output_base)
        except OSError:
            if not os.path.isdir(output_base):
                raise
        bam2bed = {'bam_dir': bam_dir, 'outdir': output_base}
        cmd = "for x in {bam_dir}/*/*.sortedByCoord.out.bam;do y=$(basename ${{x%.*}});bam2bed < $x > {outdir}/$y.bed;done;echo Converting to Bed done".format(**bam2bed)
        cmdlist.append(cmd)
        cmd = logcommand.format(insert=cmd)
        print(cmd)
        subprocess.call(cmd, shell=True)
    if args.bed_chr:
        output_bed = Config.get("ALL", "bed_chr_dir")
        output_base = Config.get("ALL", "bed_dir")
        print("echo Attempting to create directory base: {0}".format(output_bed))
        try:
            os.makedirs(output_bed)
        except OSError:
            if not os.path.isdir(output_bed):
                raise
        bed2bed = {'bed_dir': output_base, 'original': "Bd", 'outdir':output_bed}
        cmd = "for x in {bed_dir}/*.bed;do y=$(basename ${{x%.*}});sed -- 's/{original}/chr/g' $x > {outdir}/$y.sicer.bed;done;echo Converting to Bed chromosome done".format(**bed2bed)
        cmdlist.append(cmd)
        cmd = logcommand.format(insert=cmd)
        print(cmd)
        subprocess.call(cmd, shell=True)
    if args.sicer_rb:
        # run_sicer-rb.py InputDir bed_file OutputDir --species --rdthresh --winsize --fragsize --egf --gap_size --eval --sicer
        # More info available by running python run_sicer-rb.py -h
        bed_chr_dir = Config.get("ALL", "bed_chr_dir")
        files = (glob.glob(os.path.join(bed_chr_dir,"*.sicer.bed")))
        for f in files:
            dic = {"input":f, "conf":args.config}
            cmd = "python run_sicer-rb.py -bed_file {input} -config {conf}".format(**dic)
            cmdlist.append(cmd)
            cmd = logcommand.format(insert=cmd)
            print(cmd)
            subprocess.call(cmd, shell=True)
    if args.sicer:
        # run_sicer.py InputDir bed_file control_file OutputDir --species --rdthresh --winsize --fragsize --egf --gap_size --FDR --sicer
        # More info available by running python run_sicer.py -h
        bed_dir = Config.get("ALL", "bed_chr_dir")
        samples = glob.glob(os.path.join(bed_dir, "*sample*.sicer.bed"))
        controls = glob.glob(os.path.join(bed_dir, "*control*.sicer.bed"))
        # if len(controls) == len(samples):
        #     # Assuming equal number of samples and controls properly labelled.
        #     joined = zip(samples, controls)
        #     for j in joined:
        #         dic = {"sample":j[0], "control":j[1], "conf":args.config}
        #         cmd = "python run_sicer.py -sample {sample} -control {control} -config {conf}".format(**dic)
        #         cmdlist.append(cmd)
        #         cmd = logcommand.format(insert=cmd)
        #         print(cmd)
        #         subprocess.call(cmd, shell=True)
        # else:
        while len(samples) > 0:
            sample = samples.pop()
            control = None
            pattern = re.compile(r'Rep[0-9][0-9]')
            result = pattern.search(sample)
            for index in range(len(controls)):
                c_result = pattern.search(controls[index])
                if result.group() == c_result.group():
                    control = controls[index]
                    del controls[index]
                    break
            if control:
                dic = {"sample": sample, "control": control, "conf": args.config}
                cmd = "python run_sicer.py -sample {sample} -control {control} -config {conf}".format(**dic)
                cmdlist.append(cmd)
                cmd = logcommand.format(insert=cmd)
                print(cmd)
                subprocess.call(cmd, shell=True)
                
if __name__ == "__main__":
    main()
