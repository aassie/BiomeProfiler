#!/bin/bash
# ------------------------------------------------------------------
# [Adrien Assié] WormBiom: Estimated profiles of C. elegans related
#               Microbiome based on genome and 16S rRNA Amplicons
# ------------------------------------------------------------------

# Last update 14/01/25 -- V1.0

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

VERSION=1.0
SUBJECT=WormBiom
USAGE="wormbiom.sh -a ASV-TABLE -m ASV-METADATA -s ASV-FASTA -b BACTERIALIST -w WB-DATABASE -o OUTFOLDER -s1 CATEGORY -t TABLETYPE [-v] [-st]"
set -e

# Define script usage function first
function script_usage() {
    cat << EOF
    __        __                   ____  _
    \ \      / /__  _ __ _ __ ___ | __ )(_) ___  _ __ ___   ___
     \ \ /\ / / _ \| '__| '_ \` _ \|  _ \| |/ _ \| '_ \` _ \ / _ \ 
      \ V  V / (_) | |  | | | | | | |_) | | (_) | | | | | |  __/
       \_/\_/ \___/|_|  |_| |_| |_|____/|_|\___/|_| |_| |_|\___|
    =====================================================================
    ==                           v1.0                                   ==
    =====================================================================
    WormBiom - Pipeline to generate estimated metabolic pathways from 16S
    rRNA biom tables.
    Usage: $USAGE
    Options:
       -h|--help                Display this help

       -a|--ASVtable            Required. Biom/tsv/csv File

       -m|--ASVmeta            Required. Count table related Metadata

       -s|--ASVsequence        Required. related ASV sequences

       -b|--BacteriaList       Required. Collection Name or list of bacteria
                               name used in the experiment. One per line

       -w|--wbdb              Required. Path to WormBiome genome database

       -o|--OutFolder          Output folder
                               Default: ./OUT/

       -s1|--Selector1         Name of column to use as category used
                               for graph plot categorisation and basic
                               statistics

       -t|--TableType          Input table type (BIOM/TSV/CSV)

       -c|--CountThreshold     Minimum read count threshold for samples
                               Default: 1000

       -v|--Verbose            Verbose mode

       -st|--Stats             Perform statistical analysis using the s1 column option
EOF
}

#Setting up variable
## General
GMeta=$DIR/Data/Bacteria_Phylogeny.txt
SSUCount=$DIR/Data/16S.Counts.txt
SSUSeq=$DIR/Data/16S.sequences.fasta

MetacycDB=$DIR/Data/Metacyc.database.tsv
MetacycMeta=$DIR/Data/Metacyc.metadata.tsv

#WB=$DIR/Data/WB.tsv
KeggMeta=$DIR/Data/Kegg.metadata.tsv

RunKegg=1
RunMeta=1
RunGFF=1
RunHyp=1
STATS=0

OUT="./WormBiom.Result/"
DEBUG=0

CountThreshold=1000
Tableopt1="BIOM"  # Default table type
VERBOSE=0

## Options

WBoption2=NULL
BNAME=$(date +%Y%m%d)

## Options parsing
while [ "$1" != "" ]; do
    case $1 in
        -h|--help)
         script_usage
         exit 0
         ;;
        -a|--ASVtable)
          shift
          ASVTable=$1
          ;;
        -m|--ASVmeta)
          shift
          ASVMeta=$1
          ;;
        -s|--ASVsequence)
          shift
          ASVsequence=$1
          ;;
        -b|--BacteriaList)
          shift
          BACList=$1
          ;;
        -w|--wbdb)
            shift
            WB=$1
            ;;
        -o|--OutFolder)
          shift
          OUT=$1
          ;;
        -s1|--Selector1)
          shift
          Sel1=$1
          ;;
        -t|--TableType)
          shift
          Tableopt1=$1
          ;;
        -c|--CountThreshold)
          shift
          CountThreshold=$1
          ;;
        -v|--Verbose)
          VERBOSE=1
          ;;
        -st|--Stats)
          STATS=1
          ;;
         *)
          echo "Invalid command: no parameter included with argument $1"
          exit 1
          ;;
    esac
    shift
done

# AFTER options parsing, check required arguments
REQUIRED_ARGS=( 
    "ASVTable:(-a|--ASVtable)" 
    "ASVMeta:(-m|--ASVmeta)" 
    "ASVsequence:(-s|--ASVsequence)" 
    "BACList:(-b|--BacteriaList)"
    "WB:(-w|--wbdb)"
)

# Function to check required arguments
check_required_args() {
    local missing=()
    for arg in "${REQUIRED_ARGS[@]}"; do
        local var_name="${arg%%:*}"
        local opt_name="${arg#*:}"
        if [ -z "${!var_name}" ]; then
            missing+=("$opt_name")
        fi
    done

    if [ ${#missing[@]} -ne 0 ]; then
        echo "Error: Missing required arguments:"
        printf '%s\n' "${missing[@]}"
        echo
        echo "use -h option for help"
        exit 1
    fi
}

# Now check arguments after they've been parsed
check_required_args

# Validate input files exist
for file in "$ASVTable" "$ASVMeta" "$ASVsequence"; do
    if [ ! -f "$file" ]; then
        echo "Error: File does not exist: $file"
        exit 1
    fi
done

# Validate table type
if [[ ! "$Tableopt1" =~ ^(BIOM|TSV|CSV)$ ]]; then
    echo "Error: Invalid table type '$Tableopt1'. Must be one of: BIOM, TSV, CSV"
    exit 1
fi

# If we get here, all required arguments are present and valid
echo "All required arguments provided and validated"
echo "----"

date 

#Creating out folder
mkdir -p $OUT

#Checking for Blast
echo "Looking for Blast"
echo
MKBLASDB=$(which makeblastdb)
TBLASTN=$(which blastn)
if [ -z ${MKBLASDB+x} ]; then echo "makeblastdb has not been found, stopping" && exit 1; else echo "makeblastdb has been found"; fi
if [ -z ${TBLASTN+x} ]; then echo "blastn has not been found, stopping" && exit 1; else echo  "blastn has been found"; fi
echo
echo "----"
#Checking Bacteria list
echo "Checking Collection"
echo
echo $DIR/Data/Collection/$BACList\.fasta

if ls $DIR/Data/Collection/$BACList\.fasta 1> /dev/null 2>&1;
then echo "Found Bacteria $BACList Collection";
  BNAME=$BACList;
else echo "Creating Bacteria Collection";
  REFsequence=$DIR/Data/Collection/$BNAME\.fasta
  [ -e "$REFMissing" ] && rm -f $REFMissing;
  [ -e "$REFsequence" ] && rm -f $REFsequence;
  for i in  $(cat $BACList);
   do a=$(cat $SSUSeq | grep -w -A 1 $i || true);
    if (( $(grep -c . <<<"$a") > 1 ));
      then echo $a >> $REFsequence;
    else echo $i >> $REFMissing;
    fi;
  done;
  echo "New Collection Created"
  report1=$(grep -c "^" $BACList)
  if ls $REFMissing 1> /dev/null 2>&1; then report2=$(grep -c "^" $REFMissing); else report2=0; fi
  echo $report2"/"$report1" were missing in the reference"
  gsed -i "s/ /\n/g" $REFsequence
fi
echo

echo "----"
#Checking Blast Database
echo "Checking Blast Database"
echo
if ls $DIR/BlastDB/$BNAME* 1> /dev/null 2>&1;
  then echo "Found Blast database, using it." ;
else echo "Blast database not found, creating...";
    makeblastdb -dbtype nucl -in $REFsequence -title $BNAME -out $DIR/BlastDB/$BNAME;
fi
echo
echo "----"
#Blasting
echo "Running blast to map ASV to Genomes"
BlastResult=$OUT/Blast.Result.txt
blastn -db $DIR/BlastDB/$BNAME -query $ASVsequence -outfmt 6 -perc_identity 99 -qcov_hsp_perc 100 -out $BlastResult
echo "Blast done"
echo
echo "----"
#Running R script
echo "Running WormBiom.R"
echo $VERSION
echo

SourceDir=$DIR/R/WB.script.R
echo here: $Tableopt1 $SourceDir $DEBUG $STATS
echo
#line for debug
echo "Rscript $DIR/R/WB.R $ASVTable $ASVMeta $BlastResult $OUT $GMeta $KeggMeta $WB $Sel1 $SSUCount $DEBUG $STATS $CountThreshold $Tableopt1 $SourceDir"

Rscript $DIR/R/WB.R \
    "$ASVTable" \
    "$ASVMeta" \
    "$BlastResult" \
    "$OUT" \
    "$GMeta" \
    "$KeggMeta" \
    "$WB" \
    "$Sel1" \
    "$SSUCount" \
    "$DEBUG" \
    "$STATS" \
    "$CountThreshold" \
    "$Tableopt1" \
    "$SourceDir"
echo
echo "----"
#Creating report
#TBA

