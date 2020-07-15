#!/usr/bin/env python
import sys,os
from Bio import SeqIO
from Bio.SeqIO import FastaIO

def single_line_records(fasta_in):
    in_fasta_basename = os.path.splitext(os.path.basename(fasta_in))[0]
    out_basedir = os.path.realpath(os.path.dirname(fasta_in))
    out_filepath = os.path.join(out_basedir,in_fasta_basename+"_single_lines.fasta")

    if os.path.exists(out_filepath):
        raise IOError("%s already exists; skipping..." % out_filepath)

    with open(out_filepath, "w") as handle:
        fasta_out = FastaIO.FastaWriter(handle, wrap=None)
        fasta_out.write_header()
        for record in SeqIO.parse(fasta_in, "fasta"): 
            fasta_out.write_record(record)

if __name__ == "__main__":
    if len(sys.argv) == 2:
        print("Reformatting",sys.argv[1])
        single_line_records(sys.argv[1])
    else:
        print("Usage: {} file_to_format.fasta".format(sys.argv[0]))
        exit(1)