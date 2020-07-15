#!/usr/bin/env python
import argparse
import sys,os,copy,re,datetime
from Bio import SeqIO
from Bio.Seq import MutableSeq, Seq
from Bio.Alphabet import DNAAlphabet
from Bio.SeqIO import FastaIO

CHARACTER_TO_USE="_"

# TODO: incorporate bits from here:
# https://github.com/pvanheus/ncov/blob/master/combine_gisaid.py

def bp_range(s):
    try:
        start,end = map(int, s.split('-'))
        return start,end
    except:
        raise argparse.ArgumentTypeError("Coordinates must be in the format 'start-end'")

parser = argparse.ArgumentParser(
            description="""This cleans sequences dumped from GISAID
                            It replaces spaces in sequence ids/descriptions with underscores and makes all sequences upper-case.
                            It additionally provides options to include/exclude sequences by pattern in the id/description,
                            and by inclusion in a date range."""
            )

parser.add_argument('fasta_in', type=argparse.FileType('r'), help='Input fasta.')
parser.add_argument('--out_fasta', dest='fasta_out', type=str, default=None, help="Output fasta. If not provided, a second fasta will be created in the same directory as the input.")
parser.add_argument('--include_regex', dest='filter_include_expression', type=str, default=None, help="String with regular expression to match on; if supplied only sequence IDs/descriptions matching this expression will be INCLUDED. (Note: This happens before spaces are replaced.)")
parser.add_argument('--exclude_regex', dest='filter_exclude_expression', type=str, default=None, help="String with regular expression to match on; if supplied sequence IDs/descriptions matching this expression will be EXCLUDED (Note: This happens after inclusion matching, and before spaces are replaced.)")
parser.add_argument('--include_coordinates', dest='bp_ranges', type=bp_range, nargs='+', help="Ranges (inclusive; one-indexed) of coordinates to be included in the output (concatenated)")
parser.add_argument('--start_date', dest='start_date', type=lambda s: datetime.datetime.strptime(s, '%Y-%m-%d'), default=None, help="Date (not time) in ISO8601 format (e.g. 2012-01-01). The time to use as a lower bound for inclusion (inclusive). Note: Sequences without date-like patterns will be excluded.")
parser.add_argument('--end_date', dest='end_date', type=lambda s: datetime.datetime.strptime(s, '%Y-%m-%d'), default=None, help="Date (not time) in ISO8601 format (e.g. 2012-01-01). The time to use as a lower bound for inclusion (inclusive). Note: Sequences without date-like patterns will be excluded.")
parser.add_argument('--ungap', dest='ungap', type=str, nargs="?", const="-", default=None, help="If set, any gaps will be removed following all other processing steps (i.e. gaps are removed after pulling out specified regions). This can be given a value to change the gap character used (Default: '-').")


def clean_seqs(fasta_in,fasta_out=None,filter_include_expression=None,filter_exclude_expression=None,bp_ranges=None,start_date=None,end_date=None,ungap=None):
    iso_date_re = re.compile(r'(\d{4}-\d{2}-\d{2})')

    bp_ranges = bp_ranges or []
    
    bp_range_str          = "_".join([str(t[0])+"-"+str(t[1])+"bp" for t in bp_ranges])
    start_date_str        = "" if not start_date else "starting_"+start_date.strftime("%Y-%m-%d")
    end_date_str          = "" if not end_date else "ending_"+end_date.strftime("%Y-%m-%d")
    filter_include_str    = "" if not filter_include_expression else "only_subset_by_filter"
    filter_exclude_str    = "" if not filter_exclude_expression else "excluding_some_by_filter"
    output_summary_string = "_".join(s for s in [bp_range_str,start_date_str,end_date_str,filter_include_str,filter_exclude_str] if len(s)>0)

    if len(output_summary_string)>0:
        output_summary_string="_"+output_summary_string

    in_fasta_basename = os.path.splitext(os.path.basename(fasta_in.name))[0]
    out_basedir       = os.path.realpath(os.path.dirname(fasta_in.name))

    out_filepath = fasta_out or os.path.join(out_basedir,in_fasta_basename+"_cleaned"+output_summary_string+".fasta")

    if os.path.exists(out_filepath):
        raise IOError("%s already exists; skipping..." % out_filepath)

    if filter_include_expression:
        filter_include_re = re.compile(filter_include_expression)
    if filter_exclude_expression:
        filter_exclude_re = re.compile(filter_exclude_expression)

    with open(out_filepath, "w") as handle:
        fasta_out = FastaIO.FastaWriter(handle, wrap=80) # wrap=None
        fasta_out.write_header()
        for record in SeqIO.parse(fasta_in.name, "fasta"):
            should_output=True
            if filter_include_expression:
                should_output=False
                if filter_include_re.search(record.id) or filter_include_re.search(record.description):
                    should_output=True
            
            if filter_exclude_expression and (filter_exclude_re.search(record.id) or filter_exclude_re.search(record.description)):
                should_output=False

            if start_date:
                for field in [record.description,record.id]:
                    match = iso_date_re.search(field)
                    if match:
                        seq_date = datetime.datetime.strptime(match.group(0), "%Y-%m-%d")
                        if seq_date<start_date:
                            should_output=False

            if end_date:
                for field in [record.description,record.id]:
                    match = iso_date_re.search(field)
                    if match:
                        seq_date = datetime.datetime.strptime(match.group(0), "%Y-%m-%d")
                        if seq_date>end_date:
                            should_output=False
            
            if should_output:
                
                if len(bp_ranges)==0:
                    record.seq=MutableSeq(str(record.seq).upper(), DNAAlphabet())
                else:
                    output_seq=MutableSeq("", DNAAlphabet())
                    for start,end in bp_ranges:
                        start-=1 # remove one since biopython seqs are zero-indexed
                        # end-=1 # remove one since biopython seqs are zero-indexed; not needed because slice upper is exclusive
                        start=max(start,0) # bound to limit of sequence
                        end=min(end,len(record)) # bound to limit of sequence
                        output_seq+=record.seq[start:end]
                    record.seq=Seq(str(output_seq).upper(),DNAAlphabet())

                if ungap!=None:
                    record.seq=Seq(str(record.seq).upper(),DNAAlphabet()).ungap(ungap)

                #record.id=copy.deepcopy(record.id).replace(" ",CHARACTER_TO_USE)
                #record.description=copy.deepcopy(record.description).replace(" ",CHARACTER_TO_USE)
                # set the id to the description, which is the ID in the case of GISAID
                # and remove the description. 
                record.id=copy.deepcopy(record.description).replace(" ",CHARACTER_TO_USE)
                record.description=""
                fasta_out.write_record(record)


if __name__ == "__main__":
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    print("Reformatting",args.fasta_in.name)
    clean_seqs(args.fasta_in, args.fasta_out, args.filter_include_expression, args.filter_exclude_expression, 
                args.bp_ranges, args.start_date, args.end_date, args.ungap)
