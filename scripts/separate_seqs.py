#!/usr/bin/env python
import sys,os
from Bio import SeqIO
from Bio.SeqIO import FastaIO

def split_records(fasta_in):
    for record in SeqIO.parse(fasta_in, "fasta"):
        in_fasta_basename = os.path.splitext(os.path.basename(fasta_in))[0]
        out_basedir = os.path.realpath(os.path.join(os.path.dirname(fasta_in),in_fasta_basename))
        if not os.path.isdir(out_basedir):
            os.makedirs(out_basedir, exist_ok=True)
        out_filepath = os.path.join(out_basedir,record.id+".fasta")
        print("%s %i -> %s" % (record.id, len(record), out_filepath))
        if not os.path.exists(out_filepath):
            with open(out_filepath, "w") as handle:
                fasta_out = FastaIO.FastaWriter(handle, wrap=None)
                fasta_out.write_header()
                fasta_out.write_record(record)
                #SeqIO.write(record, handle, "fasta")
        else:
            #raise IOError("%s already exists; skipping..." % out_filepath)
            print("%s already exists; skipping..." % out_filepath)

if __name__ == "__main__":
    if len(sys.argv) == 2:
        print("Splitting",sys.argv[1])
        split_records(sys.argv[1])
    else:
        print("Usage: {} file_to_split.fasta".format(sys.argv[0]))
        exit(1)