#!/usr/bin/env python
import sys,os
from Bio import SeqIO
from Bio.SeqIO import FastaIO

def subset_to_ids_not_in_file(fasta_in, ids_file, fasta_out):
    in_fasta_basename = os.path.splitext(os.path.basename(fasta_in))[0]
    out_basedir = os.path.realpath(os.path.dirname(fasta_in))
    #out_filepath = os.path.join(out_basedir,in_fasta_basename+"_subset.fasta")
    out_filepath = fasta_out

    ids_to_include = set()

    with open(ids_file) as ids_file:
        for line in ids_file:
            ids_to_include.add(line.rstrip().replace("\n",""))

    #if os.path.exists(out_filepath):
    #    raise IOError("%s already exists; skipping..." % out_filepath)

    with open(out_filepath, "w") as handle:
        fasta_out = FastaIO.FastaWriter(handle, wrap=None)
        fasta_out.write_header()
        for record in SeqIO.parse(fasta_in, "fasta"):
            if record.id not in ids_to_include:
                fasta_out.write_record(record)

if __name__ == "__main__":
    if len(sys.argv) == 4:
        print("Subsetting",sys.argv[1],"to include IDs in",sys.argv[2])
        print("Output:", sys.argv[3])
        subset_to_ids_not_in_file(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        print("Usage: {} in.fasta file_containing_ids out.fasta".format(sys.argv[0]))
        exit(1)