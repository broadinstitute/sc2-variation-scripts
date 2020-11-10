#!/usr/bin/env python

import os, sys, re, copy, glob
import csv
import statistics
import re
import datetime
import dateutil.parser
import hashlib
from collections import defaultdict, OrderedDict

from Bio import SeqIO
from Bio.SeqIO import FastaIO

import numpy as np
import pandas as pd
from pandas.tseries.offsets import *
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
from matplotlib.backends.backend_pdf import PdfPages

placeholder_date_vals = datetime.datetime.strptime('01/15/01', '%m/%d/%y')

def read_metadata(metadata_filepath, geolevel=None):
    geolevel = geolevel or "region"

    regions = set()

    sample_metadata = {}

    with open(metadata_filepath) as inf:
        metadata_reader = csv.DictReader(inf, dialect='excel-tab')

        for idx,row in enumerate(metadata_reader):
            sample_metadata[row["strain"]] = row
            if len(row[geolevel])>0:
                regions.add(row[geolevel])

    return (regions,sample_metadata)


def read_seq_sub_freqs(fasta_in, gene_range, aa_sub, regions, sample_metadata, out_filepath=None, exclude_older_than=None, geolevel=None):
    """ gene and aa sub positions should be one-indexed """

    geolevel = geolevel or "region"

    freq_by_date = {} 
    cumulative_freq_by_date = {} 

    matching_codons_seen = defaultdict(int)
    other_codons_seen = defaultdict(int)

    numerator_per_date = defaultdict(int)
    denominator_per_date = defaultdict(int)

    numerator_regions_by_date = defaultdict(lambda: OrderedDict({region:0 for region in list(regions)}))
    denominator_regions_by_date = defaultdict(lambda: OrderedDict({region:0 for region in list(regions)}))
    frequency_in_regions_by_date = defaultdict(lambda: OrderedDict({region:0 for region in list(regions)}))

    gene_start,gene_end=gene_range # one-indexed

    aa_sub_re = re.compile(r"^(?P<anc>[ARNDCQEGHILKMFPSTWYV])?(?P<aa_pos>\d+)(?P<aa_sub>[ARNDCQEGHILKMFPSTWYV])$")
    m = re.match(aa_sub_re, aa_sub)
    assert m.group("aa_pos") and m.group("aa_sub")
    aa_expected = m.group("aa_sub")
    aa_expected_pos = int(m.group("aa_pos"))

    out_filepath = out_filepath or "data/output"
    matching_fasta = out_filepath+"/"+aa_sub+"_matching_seqs.fasta"

    with open(matching_fasta, "w") as handle:
        fasta_out = FastaIO.FastaWriter(handle, wrap=None)
        fasta_out.write_header()

        for idx,record in enumerate(SeqIO.parse(fasta_in, "fasta")):
            if record.id in sample_metadata:
                sample_date=None
                sample_region=None
                try:
                    sample_date = datetime.datetime.strptime(sample_metadata[record.id]["date"], "%Y-%m-%d")
                    # parse date, defaulting to set values for missing values
                    # sample_date = dateutil.parser.parse(sample_metadata[record.id]["date"].replace("-XX",""), default=placeholder_date_vals)
                    sample_region = sample_metadata[record.id][geolevel]
                except:
                    # continue to next seq in the event date could not be parsed
                    #print("could not parse date for",record.id)
                    continue

                if (exclude_older_than is not None and sample_date < datetime.datetime.strptime(exclude_older_than, "%Y-%m-%d")) or (sample_date > datetime.datetime.now()):
                    continue

                gene_seq=record.seq[gene_start-1:gene_end]#.ungap("-")

                denominator_regions_by_date[sample_date][sample_region]+=1

                # to examine degeneracy
                codon=gene_seq[((aa_expected_pos-1)*3):((aa_expected_pos-1)*3+3)]
                aa_seen=""
                aa_seen = codon.translate()

                if aa_seen == aa_expected:
                    numerator_per_date[sample_date]+=1

                    record.id=copy.deepcopy(record.description)+"|"+sample_date.strftime('%Y-%m-%d')
                    record.description=""
                    fasta_out.write_record(record)

                    numerator_regions_by_date[sample_date][sample_region]+=1

                    matching_codons_seen[str(codon.upper())]+=1
                else:
                    other_codons_seen[str(codon.upper())]+=1

                denominator_per_date[sample_date]+=1

    print("matching_codons_seen seen",matching_codons_seen)
    print("other_codons_seen seen",other_codons_seen)

    cumulative_numerator=0
    cumulative_denominator=0
    for date,denominator in sorted(denominator_per_date.items()):
        freq_by_date[date] = float(float(numerator_per_date.get(date,0))/float(denominator))
        cumulative_numerator+=numerator_per_date.get(date,0)
        cumulative_denominator+=denominator
        cumulative_freq_by_date[date] = float(cumulative_numerator)/float(cumulative_denominator)

        for region in regions:
            if denominator_regions_by_date[date][region]==0:
                frequency_in_regions_by_date[date][region]="NaN"
            else:    
                frequency_in_regions_by_date[date][region]=float(numerator_regions_by_date[date][region])/float(denominator_regions_by_date[date][region])

    last_N_dates_with_records=7
    print("mean frequency of the last {} dates with records:".format(last_N_dates_with_records), statistics.mean([e[1] for e in sorted(freq_by_date.items())[len(freq_by_date)-last_N_dates_with_records:]]))
    print("number in the last {} dates with records:".format(last_N_dates_with_records), sum([e[1] for e in sorted(denominator_per_date.items())[len(denominator_per_date)-last_N_dates_with_records:]]))

    denominator_regions_by_date_df = pd.DataFrame(denominator_regions_by_date)
    denominator_regions_by_date_df.transpose()[exclude_older_than:].sort_index().to_csv(out_filepath+"/{geolevel}_by_date.tsv".format(geolevel=geolevel), sep="\t", header=True, index=True)

    numerator_regions_by_date_df = pd.DataFrame(numerator_regions_by_date)
    numerator_regions_by_date_df.transpose()[exclude_older_than:].sort_index().to_csv(out_filepath+"/"+aa_sub+"_per_{geolevel}_by_date.tsv".format(geolevel=geolevel), sep="\t", header=True, index=True)

    frequency_in_regions_by_date_df = pd.DataFrame(frequency_in_regions_by_date)
    frequency_in_regions_by_date_df.transpose()[exclude_older_than:].sort_index().to_csv(out_filepath+"/"+aa_sub+"_{geolevel}_frequency_by_date.tsv".format(geolevel=geolevel), sep="\t", header=True, index=True)

    freq_by_date_df = pd.DataFrame.from_dict(freq_by_date, orient='index')
    freq_by_date_df.index = pd.to_datetime(freq_by_date_df.index)
    freq_by_date_df=freq_by_date_df.sort_index()

    time_period_ago_date=(datetime.datetime.now()-datetime.timedelta(days=last_N_dates_with_records)).strftime('%Y-%m-%d')
    denominator_per_date_df = pd.DataFrame.from_dict(denominator_per_date,orient='index')
    denominator_per_date_df.index = pd.to_datetime(denominator_per_date_df.index)
    denominator_per_date_df=denominator_per_date_df.sort_index()
    #denominator_per_date_df=denominator_per_date_df.transpose().sort_index()
    denominator_per_date_df.to_csv(out_filepath+"/total_sequences_by_date.tsv", sep="\t", header=True, index=True)
    print(freq_by_date_df)
    print(freq_by_date_df.loc[time_period_ago_date:,])
    print("mean frequency of the last {} days:{}".format(last_N_dates_with_records, freq_by_date_df.loc[time_period_ago_date:].mean()))
    print("number in the last {} days: {}".format(last_N_dates_with_records,denominator_per_date_df.loc[time_period_ago_date:].count()))

    print(aa_sub,"seen",cumulative_numerator,"/",cumulative_denominator,"=",float(cumulative_numerator)/cumulative_denominator)

    return (freq_by_date, cumulative_freq_by_date)

def write_freqs(tsv_out, freq_by_date, cumulative_freq_by_date, exclude_older_than=None):
    with open(tsv_out,"w") as outf:
        fieldnames = ['date', 'frequency', 'cumulative_freq']
        writer = csv.DictWriter(outf, fieldnames=fieldnames, dialect="excel-tab")

        writer.writeheader()
        
        for date,freq in freq_by_date.items():
            if (exclude_older_than is None or date > datetime.datetime.strptime(exclude_older_than, "%Y-%m-%d")) and (date < datetime.datetime.now()):
                writer.writerow({'date': date.strftime('%Y-%m-%d'), 'frequency': freq, "cumulative_freq": cumulative_freq_by_date[date]})

def plot_stacked_bar(tsv_in, avg_period=1, pdf_out="regions.pdf", tsv_avg_out=None, ylabel="number", title="collections by region over time", normalize=True, stacked=True, monochrome=True,denominator_tsv=None,vertical_line_tsv=None):
    series = pd.read_csv(tsv_in, delimiter="\t", header=0, index_col=0, parse_dates=True, squeeze=True)

    denominator_regions_by_date = None
    if denominator_tsv is not None:
        denominator_regions_by_date = pd.read_csv(denominator_tsv, delimiter="\t", header=0, index_col=0, parse_dates=True, squeeze=True)

    vertical_line_values = None
    if vertical_line_tsv is not None:
        vertical_line_values = pd.read_csv(vertical_line_tsv, delimiter="\t", header=None, index_col=0, parse_dates=True, squeeze=True)

    df = pd.DataFrame(series)
    df.dropna(axis=1,how="all",inplace=True) # remove columns with all-null values
    if denominator_regions_by_date is not None:
        denominator_regions_by_date.dropna(axis=1,how="all",inplace=True) # remove columns with all-null values
        denominator_regions_by_date = denominator_regions_by_date.loc[:, (denominator_regions_by_date != 0).any(axis=0)] # remove columns with all zeroes
    # remove columns with only one item
    # for col in df.columns:
    #     if len(df[col].unique()) == 1:
    #         print("dropping column",col)
    #         df.drop(col,inplace=True,axis=1)

    if monochrome:
        plt.style.use('grayscale')
    else:
        plt.style.use('seaborn-colorblind')
    
    # interpolate frequencies for missing dates
    df = df.reindex(pd.date_range(df.index.min(), df.index.max()), fill_value="NaN")
    df_orig = df.copy(deep=True)
    df = df.astype(float)
    df_end_interpolated = df.copy(deep=True)
    df_end_interpolated = df_end_interpolated.astype(float)
    df = df.interpolate(method='linear', limit_area="inside", axis=0)

    df_end_interpolated = df_end_interpolated.interpolate(method='linear', limit_area="inside", axis=0).ffill() #limit=avg_period*5,

    avg_series = df.rolling(avg_period).mean()

    if denominator_regions_by_date is not None:
        last_valid_indices = denominator_regions_by_date.apply(lambda x: x[x != 0].index[-1] if x.count()>0 else pd.nan)
    else:
        last_valid_indices = avg_series.apply(lambda x: x[x != 0].index[-1] if x.count()>0 else pd.nan)

    for col in df.columns:
        df_end_interpolated.loc[:,col] = np.nan
        
        last_valid_index = last_valid_indices.loc[col]
        if not pd.isna(last_valid_index):
            # set values up to the last valid index as null for the df representing the end of the graph            
            df_end_interpolated.loc[:last_valid_index,col] = np.nan
            # set the df from the end index to the interpolated values from the main dataframe
            df_end_interpolated.loc[last_valid_index:,col] = avg_series.loc[last_valid_index,col]#.copy(deep=True)
            # set the values at the end of the main dataframe to null in the region covered by df_end_interpolated
            avg_series.loc[last_valid_index+DateOffset(days=1):,col] = np.nan 

    if tsv_avg_out is not None:
        avg_series.to_csv(tsv_avg_out, sep="\t", header=True, index=True)

    print(" Dataset covers %s days" % int((df.index.max()-df.index.min()).days) )

    if normalize:
        avg_series = avg_series.div(avg_series.sum(axis=1), axis=0)

    # === to zero out dates without values ===
    #for col in df.columns:
    #    for date in df.index:
    #        if pd.isna(df.loc[date,col]) or df.loc[date,col]==None or df.loc[date,col]=="":
    #            avg_series[col][date]="NaN"


    with PdfPages(pdf_out) as pdf:
        colors=None
        if not stacked and monochrome:
            colors=["#a6d5e4"]*len(avg_series.columns)
            avg_series=avg_series.loc[:,(avg_series >= 0.75).idxmax().sort_values(ascending=True, na_position="first").index]

        avg_series.plot(stacked=stacked, 
                        subplots=not stacked, 
                        color="#4045b9" if not stacked else None, 
                        alpha=0.7, 
                        sharex=False, 
                        sharey=False, 
                        grid=True, 
                        linewidth=1 if not stacked else 0)
        
        fig=plt.gcf() 

        for idx,fig_ax in enumerate(plt.gcf().axes):
            print("Plotting: ",fig_ax.legendlabels[0])

            if not stacked:
                if len(fig_ax.legendlabels)==1:
                    if vertical_line_values is not None and fig_ax.legendlabels[0] in vertical_line_values:
                        fig_ax.axvline(x=vertical_line_values[fig_ax.legendlabels[0]],ymin=0,ymax=1,linestyle="--",color="#8d8d8d")

                    if denominator_regions_by_date is not None:
                        denominator_records_by_day=denominator_regions_by_date[fig_ax.legendlabels[0]]

                        denominator_records_by_day = denominator_records_by_day.reindex(pd.date_range(df.index.min(), df.index.max()), fill_value="NaN")
                        denominator_records_by_day = denominator_records_by_day.astype(float)

                        fig_ax.set(title="{} (N={})".format(fig_ax.legendlabels[0],int(denominator_records_by_day.sum())))
                        
                        denominator_records_by_day.plot(ax=fig_ax,color="#e93224",alpha=0.0,linewidth=1.5,secondary_y=True) # needed to set up the right axis for the bar plot
                        fig_ax.right_ax.bar(denominator_regions_by_date.loc[:,fig_ax.legendlabels[0]].index,denominator_regions_by_date.loc[:,fig_ax.legendlabels[0]].values,color="#e93224",alpha=0.8,linewidth=1.5,width=-0.8,align="edge")

                        fig_ax.set_ylim(0,1)
                        #fig_ax.right_ax.set_yscale('log')
                        fig_ax.set_ylabel(ylabel="frequency")
                        fig_ax.right_ax.set_ylabel(ylabel="collections",color="#e93224")
                        if denominator_records_by_day.sum()>0:
                            fig_ax.right_ax.set_ylim(0,denominator_regions_by_date[fig_ax.legendlabels[0]].max())
                        else:
                            fig_ax.right_ax.set_ylim(0,10)

                        fig_ax.set_zorder(2)
                        fig_ax.right_ax.set_zorder(fig_ax.get_zorder()+1) # -1 for bars behind
                        fig_ax.patch.set_visible(False) # ensure first axis does not hide second axis
                    else:
                        fig_ax.set(title="{}".format(fig_ax.legendlabels[0]))
                    fig_ax.get_legend().remove()

            if normalize:
                fig_ax.set(xlim = (df.index.min(), df.index.max()), ylim=(0,1))
            else:
                fig_ax.set(xlim = (df.index.min(), df.index.max()))


            fig_ax.fill_between(avg_series.loc[:,fig_ax.legendlabels[0]].index,avg_series.loc[:,fig_ax.legendlabels[0]].values,0, color="#a6d5e4",linewidth=0,alpha=0.75,facecolor="#a6d5e4")

            if not df_end_interpolated.loc[:,fig_ax.legendlabels[0]].isnull().values.all():
                fig_ax.fill_between(df_end_interpolated.loc[:,fig_ax.legendlabels[0]].index,df_end_interpolated.loc[:,fig_ax.legendlabels[0]].values,0, color="#a6d5e4",alpha=0.35,linewidth=0,facecolor="#a6d5e4")

            fig_ax.set_axisbelow(True)
            fig_ax.minorticks_on()
            # Customize the major grid
            fig_ax.grid(which='major', linestyle='-', linewidth='0.5', color='#5a5a5a')
            # Customize the minor grid
            fig_ax.grid(which='minor', linestyle=':', linewidth='0.5', color='#adadad', axis="y")

        if not stacked:
            size = fig.get_size_inches()
            fig.set_size_inches(size[0]*1, size[1]*1.5*(len(avg_series.columns)//6), forward=True)

        plt.tight_layout()
        print("saving", pdf_out, "...")
        pdf.savefig()

def plot_from_tsv(tsv_in, avg_period=1, pdf_out="freq.pdf", tsv_avg_out=None, ylabel="frequency", title="frequency",num_items=None):
    series = pd.read_csv(tsv_in, delimiter="\t", header=0, index_col=0, parse_dates=True, squeeze=True)

    df = pd.DataFrame(series)

    plt.style.use('grayscale')

    df = df.reindex(pd.date_range(df.index.min(), df.index.max()), fill_value="NaN")
    df = df.astype(float)
    df = df.interpolate(method='linear', axis=0)

    avg_series = df.rolling(avg_period).mean()

    for col in df.columns:
        for date in df.index:
            if pd.isna(df.loc[date,col]) or df.loc[date,col]==None or df.loc[date,col]=="":
                avg_series[col][date]="NaN"

    if tsv_avg_out is not None:
        avg_series.to_csv(tsv_avg_out, sep="\t", header=True, index=True)

    print(" df.index.max()-df.index.min()", (df.index.max()-df.index.min()).days )

    color_dict = {"frequency":"#4045b9","cumulative_freq":"#222222"} #d6ff00

    avg_series.plot(color=[color_dict.get(x, '#333333') for x in avg_series.columns],linewidth=1)
    ax = plt.gca()
    ax.set(xlim = (df.index.min(), df.index.max()), ylim=(0,1))
    ax.set(ylabel=ylabel)
    ax.get_legend().remove()
    ax.set_axisbelow(True)
    ax.minorticks_on()
    # Customize the major grid
    ax.grid(which='major', linestyle='-', linewidth='0.5', color='black')
    # Customize the minor grid
    ax.grid(which='minor', linestyle=':', linewidth='0.5', color='gray', axis="y")


    if (num_items is not None and num_items>0):
        plt.title(title+" (N=%s)" % int(num_items))
    else:
        plt.title(title)
    plt.fill_between(avg_series.index,avg_series["frequency"].values,0, color="#a6d5e4",linewidth=0,facecolor="#a6d5e4")
    plt.tight_layout()
    print("saving", pdf_out, "...")
    plt.savefig(pdf_out)

if __name__ == "__main__":
    # =============================
    # Spike D614G
    gene_range=(21563,25384)
    aa_sub="D614G"
    #aa_sub="Y453F"

    # ORF1b P314L
    #gene_range=(13468,21555)
    #aa_sub="P314L"

    # ORF1a L3606F
    #gene_range=(266,13468)
    #aa_sub="L3606F"

    # N gene A222V
    # see: https://www.medrxiv.org/content/10.1101/2020.10.25.20219063v1.full.pdf
    #gene_range=(28274,29533)
    #aa_sub="A220V"
    # =============================

    msa_fasta_pattern = "data/*_intermediate_output/msa.ends-masked.sites-masked.fasta"
    msa_fastas = glob.glob(msa_fasta_pattern)
    input_msa = ""
    if len(msa_fastas) == 1:
        input_msa = msa_fastas[0] # sys.argv[1]
    else:
        raise FileNotFoundError("%s not found as one distinct file" % msa_fasta_pattern)



    gisaid_metadata_pattern = "data/metadata*.tsv"
    metadata_tsvs = glob.glob(gisaid_metadata_pattern)
    metadata_tsv = ""
    if len(msa_fastas) == 1:
        metadata_tsv = metadata_tsvs[0] # sys.argv[1]
    else:
        raise FileNotFoundError("%s not found as one distinct file" % gisaid_metadata_pattern)

    travel_restriction_dates = "data/travel_restriction_dates/border_restriction_dates.tsv" # sys.argv[3]
    out_filepath="data/output"

    re_compute_freqs=True

    for geolevel in ["region"]:
    #for geolevel in ["country","region"]:
        
        regions=True
        # set bool to re-compute frequencies from input data, set False to only run plotting
        if re_compute_freqs:
            regions,sample_metadata        = read_metadata(metadata_tsv,geolevel=geolevel)
            freq_by_date, cumulative_freq_by_date = read_seq_sub_freqs(input_msa, gene_range, aa_sub, regions, sample_metadata, out_filepath=out_filepath, geolevel=geolevel, exclude_older_than="2019-09-01")
            write_freqs(out_filepath+"/{aa_sub}_freqs.tsv".format(aa_sub=aa_sub), 
                        freq_by_date, 
                        cumulative_freq_by_date, 
                        exclude_older_than="2019-09-01")
        
        avg_period=7

        seqs_by_day = pd.read_csv(out_filepath+"/total_sequences_by_date.tsv", delimiter="\t", header=0, index_col=0, parse_dates=True, squeeze=True)
        total_num_seqs = seqs_by_day.sum()
        print("total_num_seqs",total_num_seqs)

        plot_from_tsv( out_filepath+"/{aa_sub}_freqs.tsv".format(aa_sub=aa_sub),
                        avg_period=avg_period,
                        ylabel="frequency of {aa_sub}".format(aa_sub=aa_sub),
                        title="global frequency of {aa_sub} over time".format(aa_sub=aa_sub),
                        pdf_out=out_filepath+"/{aa_sub}_freq_global_{avg}-day_avg.pdf".format(aa_sub=aa_sub,avg=7),
                        tsv_avg_out=out_filepath+"/{aa_sub}_freq_global_{avg}-day_avg_interpolated.tsv".format(aa_sub=aa_sub,avg=7),
                        num_items=total_num_seqs)
        plot_stacked_bar(out_filepath+"/{aa_sub}_{geolevel}_frequency_by_date.tsv".format(aa_sub=aa_sub,geolevel=geolevel),
                        avg_period=avg_period,
                        ylabel="frequency",
                        title="{aa_sub} frequency by {geolevel} over time".format(aa_sub=aa_sub,geolevel=geolevel),
                        pdf_out=out_filepath+"/{aa_sub}_frequency_by_{geolevel}_{avg}-day_avg.pdf".format(aa_sub=aa_sub,avg=avg_period,geolevel=geolevel),
                        tsv_avg_out=out_filepath+"/{aa_sub}_frequency_by_{geolevel}_{avg}-day_avg_interpolated.tsv".format(aa_sub=aa_sub,avg=avg_period,geolevel=geolevel),
                        normalize=False,
                        stacked=False,
                        denominator_tsv=out_filepath+"/{geolevel}_by_date.tsv".format(geolevel=geolevel),
                        vertical_line_tsv=travel_restriction_dates)
        
        
