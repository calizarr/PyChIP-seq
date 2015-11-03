#!/usr/bin/python
from __future__ import print_function
import argparse
import subprocess


def options():
    parser = argparse.ArgumentParser(description="Generate Indices for STAR",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("nThreads", help="Number of threads to use for STAR")
    parser.add_argument("read_length", help="Input max(ReadLength)-1")
    parser.add_argument("-st",
                        metavar="Star Indices",
                        default = "/shares/tmockler_share/clizarraga/Brachypodium_distachyon/Phytozome/v3.1/assembly/STAR_indices_masked",
                        help="Directory for STAR indices")
    parser.add_argument("-mg",
                        metavar = "Genome (masked)",
                        default = "/shares/tmockler_share/clizarraga/Brachypodium_distachyon/Phytozome/v3.1/assembly/Bdistachyon_314_v3.0.hardmasked.fa",
                        help="Masked genome fasta file")
    parser.add_argument("-gff",
                        default = "/shares/tmockler_share/clizarraga/Brachypodium_distachyon/Phytozome/v3.1/annotation/Bdistachyon_314_v3.1.gene_exons.gff3",
                        help="GFF file for STAR")
    parser.add_argument("-usegtf", action="store_true", help="Use GTF with STAR instead of GFF (pass GTF to -gff)")
    args = parser.parse_args()
    return args


def main():
    args = options()
    # Generate Indices
    if args.usegtf:
        cmd = "STAR --runThreadN {0} --runMode genomeGenerate --genomeDir {1} --genomeFastaFiles {2} --sjdbGTFfile {3} --sjdbOverhang {4}" \
          .format(args.nThreads, args.st, args.mg, args.gff, args.read_length)
    else:
        cmd = "STAR --runThreadN {0} --runMode genomeGenerate --genomeDir {1} --genomeFastaFiles {2} --sjdbGTFfile {3} --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang {4}" \
          .format(args.nThreads, args.st, args.mg, args.gff, args.read_length)
    print("Running cmd: ")
    print(cmd)
    subprocess.call(cmd, shell=True)


if __name__ == "__main__":
    main()
