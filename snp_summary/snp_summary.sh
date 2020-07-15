#!/bin/bash

# way to get the absolute path to this script that should
# work regardless of whether or not this script has been sourced
# Find original directory of bash script, resovling symlinks
# http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in/246128#246128
function absolute_path() {
    local SOURCE="$1"
    while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
        DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
        if [[ "$OSTYPE" == "darwin"* ]]; then
            SOURCE="$(readlink "$SOURCE")"
        else
            SOURCE="$(readlink -f "$SOURCE")"
        fi
        [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
    done
    echo "$SOURCE"
}
SOURCE="${BASH_SOURCE[0]}"
SCRIPT=$(absolute_path "$SOURCE")
SCRIPT_DIRNAME="$(dirname "$SOURCE")"
SCRIPTPATH="$(cd -P "$(echo $SCRIPT_DIRNAME)" &> /dev/null && pwd)"
SCRIPT="$SCRIPTPATH/$(basename "$SCRIPT")"

pushd "$SCRIPTPATH" &> /dev/null

intermediate_dir="$(dirname $(ls -1 data/*intermediate_output/msa.ends-masked.sites-masked_fewer_ambig.fasta))"

# replace ambig codes with N
sed '/^[^>]/s/[R|Y|W|S|M|K|H|B|V|D]/N/g' ${intermediate_dir}/msa.ends-masked.sites-masked_fewer_ambig.fasta > ${intermediate_dir}/msa.ends-masked.sites-masked_fewer_ambig_ATCGN-only.fasta

# write SNPs based on MSA
snp-sites -v -o ${intermediate_dir}/snps_only_ATCGN.vcf ${intermediate_dir}/msa.ends-masked.sites-masked_fewer_ambig_ATCGN-only.fasta

# split multiallelic rows into multiple rows
bcftools norm ${intermediate_dir}/snps_only_ATCGN.vcf -o ${intermediate_dir}/snps_only_ATCGN_norm.vcf -m -snps
# count number of ACTG SNPs
#tail -n +8 snps_only_ATCGN_norm.vcf | cut -f5 | grep -E "[ACTG]" | wc -l

# count distinct positions
#tail -n +8 snps_only_ATCGN_norm.vcf | cut -f1-5 | grep -E '\d+\s\d+\s\.\s.\s[ACTG]' | cut -f2 | sort | uniq | wc -l


#rename chromosome column:
#cat ${intermediate_dir}/snps_only_ATCGN.vcf | sed "s/^1/NC_045512.2/" > ${intermediate_dir}/snps_only_ATCGN_updated_chr.vcf
cat ${intermediate_dir}/snps_only_ATCGN_norm.vcf | sed "s/^1/NC_045512.2/" > ${intermediate_dir}/snps_only_ATCGN_norm_updated_chr.vcf

# remove ambig from vcf:
awk '$1 ~ /^#/ {print $0;next} {if ($4 ~ /A|C|T|G/ && $5 ~ /A|C|T|G/) print $0}' ${intermediate_dir}/snps_only_ATCGN_norm_updated_chr.vcf > ${intermediate_dir}/snps_only_ATCGN_norm_updated_chr_no_ambig_snps.vcf

mkdir -p data/output/snpEffOutput

#annotate via snpEff:
snpEff ann -v -csvStats data/output/snpEffOutput/snp_stats.csv -stats data/output/snpEffOutput/snp_stats.html NC_045512.2 ${intermediate_dir}/snps_only_ATCGN_norm_updated_chr_no_ambig_snps.vcf > data/output/snps.ann.vcf

# extract a few positions of interest
#bcftools view -r "NC_045512.2:23403,NC_045512.2:14408,NC_045512.2:241,NC_045512.2:3037" snpEffOutput/snps.ann.vcf.gz > ./snpEffOutput/snps_of_interest.ann.vcf

popd