#!/usr/bin/env python
import sys,os,hashlib,csv
from collections import defaultdict, OrderedDict
import argparse
#import shelve

from Bio import SeqIO
from Bio.SeqIO import FastaIO

# requirements: biopython

class HaplotypeSequence(object):
    def __init__(self,seq):
        self.seq=str(seq.seq)
        self.attribute_counts=defaultdict(int)
        self.seqid=seq.id
        self.sha256sum=str(hashlib.sha256(str(self.seq).encode("UTF-8")).hexdigest())

    @property
    def sha256(self):
        return self.sha256sum

    def increment_attribute(self, attribute_name):
        self.attribute_counts[attribute_name]+=1

    def decrement_attribute(self, attribute_name):
        self.attribute_counts[attribute_name]-=1

    def set_attribute_count(self,attribute_name,count):
        self.attribute_counts[attribute_name]=count

    def set_sequence_id(self, seq_id):
        self.seqid = seq_id

    def __str__(self):
        return "\t".join(["ID:"+self.seqid,"sha256:"+self.sha256sum,",".join("%s:%s"%(k,v) for k,v in self.attribute_counts.items())])


def main():
    seqdb={}

    parser = argparse.ArgumentParser(
            description="""This tabulates numbers of identical sequences and creates a PopART-compatible
                            NEXUS file, optionally adding trait attribute counts."""
        )
    parser.add_argument('in_fasta', 
                        type=str, 
                        help='Input fasta.')
    parser.add_argument('metadata_tsv_path', 
                        type=str, 
                        help="""tab-delimited metadata file; 
                                must include one column with values
                                matching seq id values. Extra metadata ignored.""")
    parser.add_argument('out_nexus', 
                        type=str, 
                        help="Output nexus file.")
    parser.add_argument('--seqIdColumnName', 
                        dest='seq_id_column_name', 
                        type=str, 
                        help="""name of the column in the metadata file to match against sequence IDs. 
                                If not provided the first column will be used.""")
    parser.add_argument('--traitColumnNamesToInclude', 
                        dest='trait_column_names_to_include', 
                        type=str, 
                        nargs='+', 
                        help="names of columns in metadata file to include (header values)")
    parser.add_argument('--booleanTraitColumnNames', 
                        dest='boolean_trait_column_names', 
                        type=str, 
                        nargs='+', 
                        help="names of columns in metadata file (header values) to interpret as boolean values")
    parser.add_argument('--traitValuesToUseAsBooleanTrue', 
                        dest='boolean_true_values', 
                        type=str, 
                        nargs='+', 
                        default=["TRUE","YES"], 
                        help="""values to interpret as boolean true values. 
                            If these values are seen, the column name is used as the value for tabulation (default: %(default)s)""")
    parser.add_argument('--includeMissing', 
                        dest='include_missing',
                        action='store_true',  
                        help="""If specified, sequences with missing metadata will be included in the output""")
    parser.add_argument('--includeUnknown', 
                        dest='include_unknown',
                        action='store_true',  
                        help="""If specified, an "UNKNOWN" trait will be included in the output for entries without any counts""")

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    metadata = read_in_metadata(args.metadata_tsv_path, 
                                boolean_trait_column_names    = args.boolean_trait_column_names, 
                                seq_id_column_name            = args.seq_id_column_name, 
                                boolean_true_values           = args.boolean_true_values, 
                                trait_column_names_to_include = args.trait_column_names_to_include)
    for record in read_in_seqs(args.in_fasta):
        if record.id not in metadata and not args.include_missing:
            print("MISSING FROM METADATA",record.id)
            continue
        
        s=HaplotypeSequence(record)

        if s.seqid in seqdb:
            s=seqdb[s.seqid]

        if record.id in metadata:
            record_metadata=metadata[record.id]
            for key,val in record_metadata.items():
                if key in args.boolean_trait_column_names:
                    if val==True:
                        s.increment_attribute(key)
                    else:
                        if key not in s.attribute_counts:
                            s.set_attribute_count(key,0)
                else:
                    s.increment_attribute(val)
        elif args.include_missing:
            pass
    
        seqdb[s.seqid] = s

    #for key in seqdb:
    #    print(seqdb[key])

    write_popart_nexus(seqdb, args.out_nexus, args.in_fasta, " ".join(sys.argv),args.include_unknown)

def write_popart_nexus(seqdb, out_nexus, in_fasta, arg_string,include_unknown):
    pass
    with open(out_nexus, "w") as fout:            

        fout.write("#NEXUS\n[File created using msa_fasta_to_popart_nexus.py using %s and %s]\n\nBEGIN TAXA;\n" % (in_fasta, arg_string))

        seq_ids=[]
        #for key in seqdb:
        #    print(seqdb[key])
        fout.write("DIMENSIONS NTAX=%d;\n\nTAXLABELS\n%s\n;\n\nEND;\n\n" % (len(seqdb), '\n'.join(seqdb.keys())))
        
        seq_length=0
        for key,value in seqdb.items():
            seq_length = len(value.seq)
            break # we only need one row
        fout.write("BEGIN CHARACTERS;\n")
        fout.write("DIMENSIONS NCHAR=%d;\n"%seq_length)
        fout.write("FORMAT DATATYPE=DNA MISSING=N GAP=-;\n")
        fout.write("MATRIX")
        for key in seqdb:
            fout.write("\n%s %s\n" % (key,seqdb[key].seq)) # maybe use "_"join(seqdb[key].sequence_ids)
        fout.write("\n;\n")
        fout.write("\nEND;\n\n")
        
        trait_keys=set()
        for key,value in seqdb.items():
            trait_keys |= set(value.attribute_counts.keys())
            #break # we only need one row to key the key list; missing entries return 0 due to defaultdict()
        trait_keys=list(trait_keys)
        fout.write("BEGIN TRAITS;\n")
        trait_nums=len(trait_keys)
        if include_unknown:
            trait_nums+=1
        fout.write("Dimensions NTRAITS=%d;\n" % trait_nums)
        fout.write("Format labels=yes missing=? separator=Comma;\n")
        fout.write("TraitLabels %s;\n" % " ".join(trait_keys+([] if not include_unknown else ["UNKNOWN"])))
        fout.write("Matrix\n\n")

        for key,value in seqdb.items():
            seq_obj = seqdb[key]
            trait_states_for_seq = [value.attribute_counts[trait] for trait in trait_keys]
            if include_unknown and sum(trait_states_for_seq)==0:
                trait_states_for_seq.append(1)
            else:
                trait_states_for_seq.append(0)
            fout.write("%s %s\n" % (key, ",".join([str(trait_state) for trait_state in trait_states_for_seq])))
            
        fout.write(";\n\nEND;\n\n")    
        

def read_in_metadata(metadata_tsv_path, 
                        boolean_trait_column_names    = None, 
                        seq_id_column_name            = None, 
                        boolean_true_values           = None, 
                        trait_column_names_to_include = None):
    metadata=OrderedDict()
    seq_id_column_name            = seq_id_column_name or None
    boolean_true_values           = boolean_true_values or []
    boolean_trait_column_names    = boolean_trait_column_names or []
    trait_column_names_to_include = trait_column_names_to_include or []

    with open(metadata_tsv_path,'rt') as metadata_file:
        reader = csv.DictReader(metadata_file, delimiter='\t')
        for row in reader:
            # use the first column as the seq IDs if it is not specified
            if seq_id_column_name is None:
                seq_id_column_name=list(row.keys())[0]
            # if the user has not specified the traits to include
            # include all of them except the seq ID column
            if len(trait_column_names_to_include)==0:
                trait_column_names_to_include=[k for k in row.keys() if k not in [seq_id_column_name]]
            
            # only store the desired trait values with seq ID as key
            dict_to_store={}
            for key,val in row.items():
                if key in trait_column_names_to_include:
                    if key in boolean_trait_column_names:
                        if val in boolean_true_values:
                            dict_to_store[key] = True
                        else:
                            dict_to_store[key] = False
                    else:
                        dict_to_store[key] = val

            metadata[row[seq_id_column_name]] = dict_to_store
    return metadata

def read_in_seqs(in_fasta):
    for record in SeqIO.parse(in_fasta, "fasta"):
        yield record

if __name__ == "__main__":
    main()
