#!/bin/bash
#$ -N dl_meth
#$ -q fat002.q
#$ -cwd
#$ -j y
#$ -V
#$ -S /bin/bash
source ~/.bashrc
conda activate r_env

set -u

# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders (PD x IBD)
# Data download on the server
#
# Usage:
#   qsub 03_download_data.sh              # default: download up to tier 2
#   TIER=1 qsub 03_download_data.sh       # only the minimal set needed to run the pipeline
#   TIER=3 qsub 03_download_data.sh       # plus the large GSE145361 matrix
#   INCLUDE_RAW=1 TIER=4 qsub ...         # including RAW.tar (about 32 GB, not needed by this design)
#
# Tiers:
#   1 = IBD primary analysis set + PD confounding adjudication set, enough to run the whole pipeline
#   2 = plus the PD replication set, the paediatric IBD cohort and the negative-control disease
#   3 = the 6.8 GB beta matrix of GSE145361 (largest sample size but no age field)
#   4 = all RAW.tar files (IDAT)
# ==============================================================

TIER=${TIER:-2}
INCLUDE_RAW=${INCLUDE_RAW:-0}

BASE=data/methylation
mkdir -p "$BASE"

# Established approach: segmented, resumable download with aria2c + an outer until loop that survives connection resets
dl(){
  local url="$1" outdir="$2" outname="$3"; local tries=0
  if [ -s "$outdir/$outname" ] && [ ! -f "$outdir/$outname.aria2" ]; then
    echo "  [skip] $outname is complete"
    return 0
  fi
  echo "[$(date '+%F %T')] >>> $outname"
  until aria2c -x16 -s16 --max-connection-per-server=16 --min-split-size=20M -c \
        --file-allocation=none --retry-wait=15 --dir="$outdir" -o "$outname" "$url"; do
    tries=$((tries+1)); echo "  retry #$tries in 15s ..."; sleep 15
    [ "$tries" -ge 50 ] && { echo "  !! GIVE UP $outname"; break; }
  done
}

geo_url(){
  local gse="$1" f="$2"
  local stub="${gse%???}nnn"
  echo "https://ftp.ncbi.nlm.nih.gov/geo/series/${stub}/${gse}/suppl/${f}"
}

get(){
  local tier="$1" gse="$2" f="$3"
  [ "$tier" -le "$TIER" ] || return 0
  case "$f" in
    *RAW.tar) [ "$INCLUDE_RAW" -eq 1 ] || return 0 ;;
  esac
  local d="$BASE/$gse"; mkdir -p "$d"
  dl "$(geo_url "$gse" "$f")" "$d" "$f"
}

echo "========== TIER=$TIER  INCLUDE_RAW=$INCLUDE_RAW =========="

# ---- tier 1: needed to run the pipeline ----
get 1 GSE87648  GSE87648_processedMethWb.txt.gz
get 1 GSE72774  GSE72774_datBetaNormalized.csv.gz
get 1 GSE72774  GSE72774_Explanation_of_characteristic_variables.pdf

# ---- tier 2: remaining cohorts of the primary analysis ----
get 2 GSE111629 GSE111629_PEGblood_450kMethylationDataBackgroundNormalized.txt.gz
get 2 GSE112611 GSE112611_beta_values.txt.gz
get 2 GSE152027 GSE152027_IOP_processed_signals.csv.gz

# ---- tier 3: the largest PD cohort ----
get 3 GSE145361 GSE145361_Vallerga2020_NCOMMS_AvgBeta_Matrix-file.txt.gz

# ---- tier 4: raw IDAT archives ----
get 4 GSE87648  GSE87648_RAW.tar
get 4 GSE72774  GSE72774_RAW.tar
get 4 GSE111629 GSE111629_RAW.tar
get 4 GSE112611 GSE112611_RAW.tar
get 4 GSE152027 GSE152027_RAW.tar
get 4 GSE145361 GSE145361_RAW.tar

# ==============================================================
# Verification
# ==============================================================
echo
echo "========== Verification =========="
fail=0
while IFS= read -r p; do
  sz=$(stat -c%s "$p")
  printf "  %-70s %12s bytes" "$p" "$sz"
  case "$p" in
    *.gz)
      if gzip -t "$p" 2>/dev/null; then echo "   gzip OK"
      else echo "   gzip corrupt; delete the file and re-run"; fail=1; fi ;;
    *.tar)
      if tar -tf "$p" >/dev/null 2>&1; then echo "   tar OK"
      else echo "   tar corrupt; delete the file and re-run"; fail=1; fi ;;
    *) echo "   (not verified)" ;;
  esac
done < <(find "$BASE" -type f ! -name '*.aria2' | sort)

# Markers of unfinished segments
if find "$BASE" -name '*.aria2' | grep -q .; then
  echo
  echo "The following files are incomplete (.aria2 control files present); re-running this script resumes them:"
  find "$BASE" -name '*.aria2' | sed 's/\.aria2$//' | sed 's/^/  /'
  fail=1
fi

echo
echo "Size of each data set:"
du -sh "$BASE"/* 2>/dev/null

echo
if [ "$fail" -eq 0 ]; then echo "[$(date '+%F %T')] ALL DONE, all files verified"
else echo "[$(date '+%F %T')] some files are incomplete or failed verification; re-run this script to resume"; fi
