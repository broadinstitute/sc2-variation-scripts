#!/usr/bin/env python
import sys,os,csv
from Bio import SeqIO
from Bio.SeqIO import FastaIO

def remap_tax_id(fasta_in, remap_table):
    in_fasta_basename = os.path.splitext(os.path.basename(fasta_in))[0]
    out_basedir = os.path.realpath(os.path.dirname(fasta_in))
    out_filepath = os.path.join(out_basedir,in_fasta_basename+"_annotated_for_beast.fasta")

    if os.path.exists(out_filepath):
        raise IOError("%s already exists; skipping..." % out_filepath)

    id_map = dict()
    with open(remap_table, "r") as map_handle:
        for line in map_handle:
            seqid,*other_fields = line.split("\t")
            if seqid in id_map:
                raise LookupError("%s already found in map" % seqid)
            if seqid=="taxa":
                # if this is a figtree-formatted annotation file, skip the header row
                # that has nothing but column labels
                continue
            id_map[seqid] = other_fields

    with open(out_filepath, "w") as handle:
        fasta_out = FastaIO.FastaWriter(handle, wrap=None)
        fasta_out.write_header()
        for record in SeqIO.parse(fasta_in, "fasta"):
            if record.id in id_map:
                if sum([len(x) for x in id_map[record.id]])>0:
                    new_description="|".join(id_map[record.id]).replace("||","|?|").replace("||","|?|")
                    record.description=new_description
                    fasta_out.write_record(record)
            else:
                print("Warning: '{}' not found in {}".format(record.id,os.path.basename(remap_table)))

if __name__ == "__main__":
    if len(sys.argv) == 3:
        print("Reformatting",sys.argv[1])
        remap_tax_id(sys.argv[1], sys.argv[2])
    else:
        print("Usage: {} file_to_format.fasta field_mappings.tsv".format(sys.argv[0]))
        print("   Where field_mappings.tsv has id->field mappings, tab-delimited in multiple columns: (id, field1, field2, ...)")
        print("   Note: fields become pipe-delimited.")
        print("   As an example:")
        print("     '>seq1 original description'")
        print("   annoted with this matching line from a TSV file:")
        print("     'seq1	case	second	11-Jul-2014'")
        print("   would result in an annotated fasta file containing:")
        print("     '>seq1 case|second|11-Jul-2014'")
        print("   ")
        exit(1)