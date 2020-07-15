#!/bin/bash

#set -e -o pipefail

MASK_BEGINNING_BP_NUM=100 #https://github.com/nextstrain/ncov/blob/master/config/config.yaml 
#MASK_BEGINNING_BP_NUM=55 # http://virological.org/t/issues-with-sars-cov-2-sequencing-data/473
MASK_END_BP_NUM=50 # https://github.com/nextstrain/ncov/blob/master/config/config.yaml
#MASK_END_BP_NUM=99 #http://virological.org/t/issues-with-sars-cov-2-sequencing-data/473
# Specific sites to mask in the reference genome's coordinates.
# These are 1-indexed coordinates of sites that have been identified as prone to sequencing errors.

# from ncov config: https://github.com/nextstrain/ncov/blob/master/config/config.yaml
MASK_SITES="13402 24389 24390" 
# from virological: http://virological.org/t/issues-with-sars-cov-2-sequencing-data/473
# likely recurrent sequencing errors:
#MASK_SITES="$MASK_SITES 187 1059 2094 3037 3130 6990 8022 10323 10741 11074 13408 14786 19684 20148 21137 24034 24378 25563 26144 26461 26681 28077 28826 28854 29700" 
# masking any homoplasic positions that are exclusive to a single sequencing lab or geographic location
#MASK_SITES="$MASK_SITES 4050 13402"
# despite having strong phylogenetic signal, are also strongly homoplasic
#MASK_SITES="$MASK_SITES 11083 15324 21575"


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

if [ $# -eq 0 ]; then
    echo "Usage: $0 context_seqs.fasta our_seqs.fasta reference.fasta output.fasta --isAligned"
    echo "--isAligned"
    exit 1
fi

OTHER_FASTA="$1"
OUR_FASTA="$2"
REFERENCE_FASTA="$3"
OUTPUT_FASTA="$4"
INPUT_IS_ALIGNED=$5

INTERMEDIATE_OUTPUT_DIR="$(dirname $OUTPUT_FASTA)/$(date +'%Y-%m-%d_%H-%M-%S')_intermediate_output"

mkdir -p "${INTERMEDIATE_OUTPUT_DIR}"

# ========= DOWNLOAD ====================
# (download from GISAID manually)
# =======================================

# ========= NORMALIZE ===================
# normalize fasta seq names via Nextstrain script
"${SCRIPTPATH}/nextstrain_ncov/normalize_gisaid_fasta.sh" "${OTHER_FASTA}" "${INTERMEDIATE_OUTPUT_DIR}/$(basename ${OTHER_FASTA}).normalized" 28000
"${SCRIPTPATH}/nextstrain_ncov/normalize_gisaid_fasta.sh" "${OUR_FASTA}" "${INTERMEDIATE_OUTPUT_DIR}/$(basename ${OUR_FASTA}).normalized" 28000
# =======================================

OTHER_FASTA_TO_MERGE=""
OUR_FASTA_TO_MERGE=""
if [[ "$INPUT_IS_ALIGNED" == "--isAligned" ]]; then
    # ========= UNGAP INPUT =================
    # ungap other seqs
    "${SCRIPTPATH}/clean_gisaid.py" --ungap --out_fasta "${INTERMEDIATE_OUTPUT_DIR}/$(basename ${OTHER_FASTA}).normalized.ungapped.fasta" "${INTERMEDIATE_OUTPUT_DIR}/$(basename ${OTHER_FASTA}).normalized"
    # ungap our seqs
    "${SCRIPTPATH}/clean_gisaid.py" --ungap --out_fasta "${INTERMEDIATE_OUTPUT_DIR}/$(basename ${OUR_FASTA}).ungapped.fasta" "${INTERMEDIATE_OUTPUT_DIR}/$(basename ${OUR_FASTA}).normalized"
    # =======================================
    OTHER_FASTA_TO_MERGE="${INTERMEDIATE_OUTPUT_DIR}/$(basename ${OTHER_FASTA}).normalized.ungapped.fasta"
    OUR_FASTA_TO_MERGE="${INTERMEDIATE_OUTPUT_DIR}/$(basename ${OUR_FASTA}).ungapped.fasta"
else
    # otherwise use directly
    OTHER_FASTA_TO_MERGE="${INTERMEDIATE_OUTPUT_DIR}/$(basename ${OTHER_FASTA}).normalized"
    OUR_FASTA_TO_MERGE="${INTERMEDIATE_OUTPUT_DIR}/$(basename ${OUR_FASTA}).normalized"
fi

# ========= MERGE IN OUR DATA ===========
# merge our data
cat "${OTHER_FASTA_TO_MERGE}" "${OUR_FASTA_TO_MERGE}" > "${INTERMEDIATE_OUTPUT_DIR}/merged.fasta"
# =======================================

# ========= CREATE ALIGNMENT TO REF =====
# align all sequences to the reference
mafft --preservecase --auto --thread -1 --keeplength --addfragments "${INTERMEDIATE_OUTPUT_DIR}/merged.fasta" "${REFERENCE_FASTA}" > "${INTERMEDIATE_OUTPUT_DIR}/msa.fasta"
# =======================================

# ========= MASK ENDS AND SITES =========
# mask the ends
# see: http://virological.org/t/issues-with-sars-cov-2-sequencing-data/473
"${SCRIPTPATH}/nextstrain_ncov/mask-alignment.py" --alignment "${INTERMEDIATE_OUTPUT_DIR}/msa.fasta" --output "${INTERMEDIATE_OUTPUT_DIR}/msa.ends-masked.fasta" --mask-from-beginning $MASK_BEGINNING_BP_NUM --mask-from-end $MASK_END_BP_NUM

# mask problematic sites 
"${SCRIPTPATH}/nextstrain_ncov/mask-alignment.py" --alignment "${INTERMEDIATE_OUTPUT_DIR}/msa.ends-masked.fasta" --output "${INTERMEDIATE_OUTPUT_DIR}/msa.ends-masked.sites-masked.fasta" --mask-from-beginning 0 --mask-sites ${MASK_SITES}
# =======================================

# ========= REMOVE AMBIG SEQS ===========
# remove any sequences with >2% N (~600bp)
#"${SCRIPTPATH}/remove_seqs_by_content.py" "${INTERMEDIATE_OUTPUT_DIR}/${OTHER_FASTA}.normalized" --out_fasta "${INTERMEDIATE_OUTPUT_DIR}/${OTHER_FASTA}.normalized.fewer_ambig" --gap_drop_threshold 0.0 --ambig_drop_threshold 0.02
"${SCRIPTPATH}/remove_seqs_by_content.py" "${INTERMEDIATE_OUTPUT_DIR}/msa.ends-masked.sites-masked.fasta" --out_fasta "${INTERMEDIATE_OUTPUT_DIR}/msa.ends-masked.sites-masked_fewer_ambig.fasta" --gap_drop_threshold 0.0 --ambig_drop_threshold 0.02
# =======================================

# ========= TRIMAL ======================
# run trimal?
#trimal -fasta -automated1 -in "${INTERMEDIATE_OUTPUT_DIR}/msa.ends-masked.sites-masked.fasta" -out "${INTERMEDIATE_OUTPUT_DIR}/msa.ends-masked.sites-masked.trimal-trimmed.fasta"
# =======================================

# ========= UNGAP (unalign) =============
#"${SCRIPTPATH}/clean_gisaid.py" --ungap --out_fasta "${OUTPUT_FASTA}" "${INTERMEDIATE_OUTPUT_DIR}/msa.ends-masked.sites-masked.trimal-trimmed.fasta"
#"${SCRIPTPATH}/clean_gisaid.py" --ungap --out_fasta "${INTERMEDIATE_OUTPUT_DIR}/seqs.ends-masked.sites-masked.ungapped.fasta" "${INTERMEDIATE_OUTPUT_DIR}/msa.ends-masked.sites-masked.fasta"
"${SCRIPTPATH}/clean_gisaid.py" --ungap --out_fasta "${OUTPUT_FASTA}" "${INTERMEDIATE_OUTPUT_DIR}/msa.ends-masked.sites-masked.fasta"
# =======================================
