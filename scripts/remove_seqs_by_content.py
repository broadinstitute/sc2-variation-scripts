#!/usr/bin/env python

import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO

parser = argparse.ArgumentParser(
            description="""This cleans sequences dumped from GISAID
                            It replaces spaces in sequence ids/descriptions with underscores and makes all sequences upper-case.
                            It additionally provides options to include/exclude sequences by pattern in the id/description,
                            and by inclusion in a date range."""
            )

parser.add_argument('fasta_in', type=argparse.FileType('r'), help='Input fasta.')
parser.add_argument('--out_fasta', dest='fasta_out', type=str, default=None, help="Output fasta. If not provided, a second fasta will be created in the same directory as the input.")
parser.add_argument('--gap_drop_threshold', dest='gap_threshold', type=float, nargs="?", default=None, const=0, help="Threshold at which sequences will be dropped if great than this fraction is made of gap characters (-). (Default: %s)")
parser.add_argument('--ambig_drop_threshold', dest='ambig_threshold', type=float, nargs="?", default=None, const=0, help="Threshold at which sequences will be dropped if great than this fraction is made of ambig characters (N). (Default: %s)")

def write_and_drop_seqs(fasta_in, fasta_out, gap_threshold=None, ambig_threshold=None):
    print("ambig_threshold",ambig_threshold)
    print("gap_threshold",gap_threshold)

    with open(fasta_out, "w") as handle:
        fasta_out = FastaIO.FastaWriter(handle, wrap=80) # wrap=None
        fasta_out.write_header()
        for record in SeqIO.parse(fasta_in.name, "fasta"):
            ambig_fraction = record.seq.count("N")/float(len(record))
            gap_fraction   = record.seq.count("-")/float(len(record))

            

            if (ambig_threshold!=None and ambig_fraction>ambig_threshold) or (gap_threshold!=None and gap_fraction>gap_threshold):
                print("omitting",record.id,"ambig:",ambig_fraction,"gap:",gap_fraction)
                continue
            else:
                #print("writing",record.id,"ambig:",ambig_fraction,"gap:",gap_fraction,"gap_chrs:",record.seq.count("-"))
                fasta_out.write_record(record)

if __name__ == "__main__":
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    print("Reformatting",args.fasta_in.name)
    write_and_drop_seqs(args.fasta_in, args.fasta_out, args.gap_threshold, args.ambig_threshold)