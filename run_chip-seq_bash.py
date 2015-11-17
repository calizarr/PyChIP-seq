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
        dic = {'input_dir':input_dir, 'config':args.config}
        cmd = "for f in {input_dir}/*.fastq\ndo\n python trim_adapters.py -input_file $f -config {config}\ndone\necho Done trimming reads...".format(**dic)
        cmd = logcommand.format(insert=cmd)
        print(cmd)
    if args.indices:
        # make_indices.py nThreads read_length -st -mg -gff -usegtf
        # More info available by running python make_indices.py -h
        cmd = "python make_indices.py -config {0}".format(args.config)
        cmd = logcommand.format(insert=cmd)
        # cmdlist.append(cmd)
        print(cmd)
    if args.map_reads:
        # map_reads.py nThreads input_file out_file_dir basepre -st
        # More info available by running python map_reads.py -h
        map_dir = Config.get("ALL", "map_dir")
        dic = {'map_dir': map_dir, 'config': args.config}
        cmd = "for f in {map_dir}/*.fastq.filtered\ndo\n python map_reads.py -input_file $f -config {config}\ndone\necho Done mapping reads...".format(**dic)
        cmd = logcommand.format(insert=cmd)
        print(cmd)
    if args.bed:
        # Converting BAM files to Bed files to Chromosome numbered Bed files for SICER
        print("echo Converting BAM files to Bed files to Chromosome numbered Bed files for SICER")
        bam_dir = Config.get("ALL", "out_file_dir")
        output_base = os.path.join(bam_dir, "Bed")
        print("echo Attempting to create directory base: {0}".format(output_base))
        cmd = "[ -d {0} ] && echo 'Directory exists' || mkdir -p {0}".format(output_base)
        print(cmd)
        bam2bed = {'bam_dir': bam_dir, 'outdir': output_base}
        cmd = "for x in {bam_dir}/*/*.sortedByCoord.out.bam\ndo\n y=$(basename ${{x%.*}})\n bam2bed < $x > {outdir}/$y.bed\ndone\necho Converting to Bed done".format(**bam2bed)
        cmd = logcommand.format(insert=cmd)
        # cmdlist.append(cmd)
        print(cmd)
        # subprocess.call(cmd, shell=True)
    if args.bed_chr:
        output_bed = Config.get("ALL", "bed_chr_dir")
        output_base = Config.get("ALL", "bed_dir")
        print("echo Attempting to create directory base: {0}".format(output_bed))
        cmd = "[ -d {0} ] && echo 'Directory exists' || mkdir -p {0}".format(output_bed)
        print(cmd)
        bed2bed = {'bed_dir': output_base, 'original': "Bd", 'outdir':output_bed}
        cmd = "for x in {bed_dir}/*.bed\ndo\n y=$(basename ${{x%.*}})\n sed -- 's/{original}/chr/g' $x > {outdir}/$y.sicer.bed\ndone\necho Converting to Bed chromosome done".format(**bed2bed)
        cmd = logcommand.format(insert=cmd)
        # cmdlist.append(cmd)
        print(cmd)
        # subprocess.call(cmd, shell=True)
    if args.sicer_rb:
        # run_sicer-rb.py InputDir bed_file OutputDir --species --rdthresh --winsize --fragsize --egf --gap_size --eval --sicer
        # More info available by running python run_sicer-rb.py -h
        bed_chr_dir = Config.get("ALL", "bed_chr_dir")
        dic = {'bed_chr_dir': bed_chr_dir, 'config': args.config}
        cmd = "for f in {bed_chr_dir}/*.sicer.bed\ndo\n python run_sicer-rb.py -bed_file $f -config {config}\ndone\necho SICER RB Done...".format(**dic)
        cmd = logcommand.format(insert=cmd)
        print(cmd)
    if args.sicer:
        # run_sicer.py InputDir bed_file control_file OutputDir --species --rdthresh --winsize --fragsize --egf --gap_size --FDR --sicer
        # More info available by running python run_sicer.py -h
        bed_chr_dir = Config.get("ALL", "bed_chr_dir")
        cmd = "samples=( $(find {bed_chr_dir} -name '*sample*.sicer.bed') )".format(bed_chr_dir=bed_chr_dir)
        print(cmd)
        cmd = "controls=( $(find {bed_chr_dir} -name '*control*.sicer.bed') )".format(bed_chr_dir=bed_chr_dir)
        print(cmd)
        dic = {'config': args.config}
        testcmd = 'if [ "${#controls[@]}" -eq "${#samples[@]}" ]; then'
        print(testcmd)
        cmd = ("for ((i=0;i<${{#controls[@]}};++i)); do\n     python run_sicer.py -sample ${{samples[i]}} -control ${{controls[i]}} -config {config}\n   done".format(**dic))
        cmd = logcommand.format(insert=cmd)
        print("   "+cmd)
        elsecmd = 'else\n regex="(Rep[0-9][0-9])"'
        print(elsecmd)
        forcmd = """ for ((i=0;i<${{#samples[@]}};++i)); do\n     [[ ${{samples[i]}} =~ $regex ]]\n     result="${{BASH_REMATCH[1]}}"\n     for ((j=0;j<${{#controls[@]}};++j));do \n         [[ ${{controls[i]}} =~ $regex ]]\n         c_result="${{BASH_REMATCH[1]}}"\n         if [[ $result = $c_result ]]; then\n          python run_sicer.py -sample ${{samples[i]}} -control ${{controls[j]}} -config {config}\n         fi\n     done\n done\nfi """.format(**dic)
        print(forcmd)
        print("unset samples")
        print("unset controls")
        print("unset sample")
        print("unset control")
        print("unset result")
        print("unset c_result")

if __name__ == "__main__":
    main()
