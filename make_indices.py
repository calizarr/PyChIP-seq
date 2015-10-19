#!/usr/bin/python
from __future__ import print_function
import argparse
import subprocess


def options():
    parser = argparse.ArgumentParser(description="Generate Indices for STAR")
    parser.add_argument("nThreads", help="Number of threads to use for STAR")
    parser.add_argument("read_length", help="Input max(ReadLength)-1")
    args = parser.parse_args()
    return args


def main():
    args = options()
    # masked_genome = "/shares/tmockler_share/clizarraga/bdaccessions/Ref/Brachypodium_distachyon.mainGenome.masked.fasta"
    masked_genome = "/shares/tmockler_share/clizarraga/Brachypodium_distachyon/Phytozome/v3.1/assembly/Bdistachyon_314_v3.0.hardmasked.fa"
    star_indices = "/shares/tmockler_share/clizarraga/Brachypodium_distachyon/Phytozome/v3.1/assembly/STAR_indices_masked"
    # gtf_file = "/shares/tmockler_share/clizarraga/Brachypodium_distachyon/Phytozome/v3.1/annotation/Bdistachyon_314_v3.1.gene_exons.gtf"
    gff_file = "/shares/tmockler_share/clizarraga/Brachypodium_distachyon/Phytozome/v3.1/annotation/Bdistachyon_314_v3.1.gene_exons.gff3"
    # Generate Indices
    cmd = "STAR --runThreadN {0} --runMode genomeGenerate --genomeDir {1} --genomeFastaFiles {2} --sjdbGTFfile {3} --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang {4}" \
          .format(args.nThreads, star_indices, masked_genome, gff_file, args.read_length)
    print("Running cmd: ")
    print(cmd)
    subprocess.call(cmd, shell=True)


if __name__ == "__main__":
    main()
