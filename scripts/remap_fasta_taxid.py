#!/usr/bin/env python
import sys,os
from Bio import SeqIO
from Bio.SeqIO import FastaIO

def remap_tax_id(fasta_in, remap_table):
    in_fasta_basename = os.path.splitext(os.path.basename(fasta_in))[0]
    out_basedir = os.path.realpath(os.path.dirname(fasta_in))
    out_filepath = os.path.join(out_basedir,in_fasta_basename+"_remapped_taxids.fasta")

    if os.path.exists(out_filepath):
        raise IOError("%s already exists; skipping..." % out_filepath)

    id_map = dict()
    with open(remap_table, "r") as map_handle:
        for line in map_handle:
            old,new = line.split("\t")
            if old in id_map:
                raise LookupError("%s already found in map" % old)
            id_map[old] = new

    with open(out_filepath, "w") as handle:
        fasta_out = FastaIO.FastaWriter(handle, wrap=None)
        fasta_out.write_header()
        for record in SeqIO.parse(fasta_in, "fasta"):
            record.id=id_map[record.id]
            fasta_out.write_record(record)

if __name__ == "__main__":
    if len(sys.argv) == 3:
        print("Reformatting",sys.argv[1])
        remap_tax_id(sys.argv[1], sys.argv[2])
    else:
        print("Usage: {} file_to_format.fasta old_new.tsv".format(sys.argv[0]))
        print("   Where old_new.tsv has old->new mappings, tab-delimited in two columns: (old, new)")
        exit(1)