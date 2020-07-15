#!/usr/bin/env python
import sys,os,hashlib
from Bio import SeqIO
from Bio.SeqIO import FastaIO

hashes_seen_before=set()


def first_occurrences_only(fasta_in):
    in_fasta_basename = os.path.splitext(os.path.basename(fasta_in))[0]
    out_basedir = os.path.realpath(os.path.dirname(fasta_in))
    out_filepath = os.path.join(out_basedir,in_fasta_basename+"_first_occurrences_only.fasta")

    total_seq_count=0

    if os.path.exists(out_filepath):
        raise IOError("%s already exists; skipping..." % out_filepath)

    with open(out_filepath, "w") as handle:
        fasta_out = FastaIO.FastaWriter(handle, wrap=None)
        fasta_out.write_header()
        for record in SeqIO.parse(fasta_in, "fasta"):
            total_seq_count+=1
            record_hash = hashlib.sha256(str(record.seq).encode("UTF-8")).hexdigest()
            if record_hash not in hashes_seen_before:
                hashes_seen_before.add(record_hash)
                fasta_out.write_record(record)
            else:
                pass
                #print("{} is identical to a sequence earlier in the input file; skipping...".format(record.id))

    print("{} seqs seen".format(total_seq_count))
    print("{} unique seqs found".format(len(hashes_seen_before)))
    print("{} identical duplicates removed".format(total_seq_count-len(hashes_seen_before)))

if __name__ == "__main__":
    if len(sys.argv) == 2:
        print("Reformatting",sys.argv[1])
        first_occurrences_only(sys.argv[1])
        print("Done.")
    else:
        print("Usage: {} file_to_format.fasta".format(sys.argv[0]))
        exit(1)