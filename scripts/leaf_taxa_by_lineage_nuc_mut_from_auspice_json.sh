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

if [ $# -eq 0 ]; then
    echo "Usage: $(basename $0) input_file_for_auspice.json nucleotide_mutation_ [PATTERN_TO_INCLUDE]"
    echo "        ex. $(basename $0) input_file_for_auspice.json C1059T [MA_MGH]"
    exit 1
fi

JSON_FILE="${1}"
NUT_MUT_STR="${2}"
PATTERN_TO_MATCH="${3:-"MA_MGH"}"

# echo "JSON_FILE $JSON_FILE"
# echo "NUT_MUT_STR $NUT_MUT_STR"
# echo "PATTERN_TO_MATCH $PATTERN_TO_MATCH"

#echo ""
#jq -r --arg NUT_MUT_MATCH_STR "${NUT_MUT_STR}" '.tree | recurse(.children[]?) | select( .branch_attrs.mutations.nuc[]? | contains($NUT_MUT_MATCH_STR)) | recurse(.children[]?) | select(.name | contains("NODE")| not) | .name' ${JSON_FILE} | grep "${PATTERN_TO_MATCH}" | sort | tr '\n' ','
jq -r --arg NUT_MUT_MATCH_STR "${NUT_MUT_STR}" '.tree | recurse(.children[]?) | select( .branch_attrs.mutations.nuc[]? | contains($NUT_MUT_MATCH_STR)) | recurse(.children[]?) | select(.name | contains("NODE")| not) | .name' ${JSON_FILE} | grep "${PATTERN_TO_MATCH}" | sort | paste -sd "," -
#printf "\n"