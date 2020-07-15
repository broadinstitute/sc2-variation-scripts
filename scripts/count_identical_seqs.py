#!/usr/bin/env python
import sys,os,hashlib
from Bio import SeqIO
from Bio.SeqIO import FastaIO
from collections import defaultdict, OrderedDict

def identical_seq_counts(fasta_in, id_lookup_files=None, tsv_out=None):
    unique_seq_ids_by_hash=defaultdict(set)
    unique_seqs_by_hash={}
    unique_seq_id_lookup_matched_counts = defaultdict(OrderedDict)
    unique_seq_id_lookup_matched_ids    = defaultdict(OrderedDict)

    in_fasta_basename = os.path.splitext(os.path.basename(fasta_in))[0]
    out_basedir = os.path.realpath(os.path.dirname(fasta_in))
    out_filepath = os.path.join(out_basedir,in_fasta_basename+"_identical_seq_counts.tsv") #tsv_out

    lookup_ids=OrderedDict()

    if id_lookup_files:
        for id_filepath in id_lookup_files:
            with open(id_filepath,"r") as inf:
                basename=os.path.splitext(os.path.basename(id_filepath))[0]
                lookup_ids[basename] = inf.read().splitlines()

    if os.path.exists(out_filepath):
        raise IOError("%s already exists; skipping..." % out_filepath)

    total_seq_count=0
    for record in SeqIO.parse(fasta_in, "fasta"):
        total_seq_count+=1
        record_hash = str(hashlib.sha256(str(record.seq).encode("UTF-8")).hexdigest())
        unique_seq_ids_by_hash[record_hash].add(str(record.id))
        #unique_seqs_by_hash[record_hash]=str(record.seq)

        # id_lookup_counts = OrderedDict()
        # id_lookup_ids    = OrderedDict()
        # for lookup_name in lookup_ids.keys():
        #     id_lookup_counts[lookup_name] = 0
        #     id_lookup_ids[lookup_name] = set()
        # id_lookup_counts["other"]=0
        # id_lookup_ids["other"]=set()

        if id_lookup_files:
            for lookup_name in list(lookup_ids.keys())+["other"]:
                if lookup_name not in unique_seq_id_lookup_matched_counts[record_hash]:
                    unique_seq_id_lookup_matched_counts[record_hash][lookup_name] = 0
                if lookup_name not in unique_seq_id_lookup_matched_ids[record_hash]:
                    unique_seq_id_lookup_matched_ids[record_hash][lookup_name] = set()

            match_found = False
            for lookup_name,lookup_ids_list in lookup_ids.items():
                if str(record.id) in lookup_ids_list:
                    match_found = True
                    unique_seq_id_lookup_matched_counts[record_hash][lookup_name] += 1
                    unique_seq_id_lookup_matched_ids[record_hash][lookup_name].add(str(record.id))
            if not match_found:
                unique_seq_id_lookup_matched_counts[record_hash]["other"] += 1
                unique_seq_id_lookup_matched_ids[record_hash]["other"].add(str(record.id))

        #unique_seq_id_lookup_matched_counts[record_hash] = id_lookup_counts
        #unique_seq_id_lookup_matched_ids[record_hash]    = id_lookup_ids

    with open(out_filepath, "w") as outf:
        header=["seq_sha256","num_seqs","seq_IDs"]
        if id_lookup_files:
            for e in list(lookup_ids.keys())+["other"]:
                header.extend([e,e+"_ct"])
        #header=["seq_sha256","num_seqs","seq_IDs","sequence"]
        outf.write("\t".join(header)+"\n")
        for seqhash,seq_ids in unique_seq_ids_by_hash.items():
            #outf.write("\t".join([seqhash,str(len(seq_ids)),",".join(seq_ids),unique_seqs_by_hash[seqhash]])+"\n")

            matches_str_elem = []
            for m,ids_to_write,count_to_write in zip(unique_seq_id_lookup_matched_ids[seqhash].keys(), unique_seq_id_lookup_matched_ids[seqhash].values(),unique_seq_id_lookup_matched_counts[seqhash].values()):
                matches_str_elem.extend([",".join(ids_to_write),str(count_to_write)])
            matches_str="\t".join(matches_str_elem)
            
            outf.write("\t".join([seqhash,str(len(seq_ids)),",".join(seq_ids),matches_str])+"\n")
            print("\t".join([seqhash,str(len(seq_ids)),",".join(seq_ids),matches_str])+"\n")
    

    print("{} seqs seen".format(total_seq_count))
    print("{} unique seqs found".format(len(unique_seq_ids_by_hash.keys())))
    print("{} duplicate sequences".format(total_seq_count-len(unique_seq_ids_by_hash.keys())))

if __name__ == "__main__":
    if len(sys.argv) >= 2:
        print("Processing",sys.argv[1])
        identical_seq_counts(sys.argv[1],sys.argv[2:])
        print("Done.")
    else:
        print("Usage: {} infasta.fasta [id_lookup_file1,id_lookup_file2, ...]".format(sys.argv[0]))
        exit(1)