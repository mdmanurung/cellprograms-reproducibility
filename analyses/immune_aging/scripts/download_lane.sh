#!/bin/bash
# Download 7 query-cohort h5ads (2 streams), pseudobulk + composition counts
# each as it lands, then delete the h5ad to manage disk.
set -u
BASE="https://artyomovlab.wustl.edu/scn-m/var/datasets/curated/maksim/GZMB_GZMK"
DL=/workspace/aging_pbmc
PB=/workspace/aging_pbmc/pseudobulk
mkdir -p $DL/{abf300,aida,cima,kcl,soundlife,spac,ucsf} $PB

process_one () {  # dir h5ad_path cohort_name
  local dir=$1 h5=$2 nm=$3
  echo "[$(date +%H:%M)] pseudobulking $nm ..."
  python3 /workspace/repo2/analyses/immune_aging/scripts/pseudobulk_h5ad.py \
      "$h5" "$PB/$dir" "$nm" --min-cells 10 > "$PB/$dir/pseudobulk_log.txt" 2>&1 \
    && echo "[$(date +%H:%M)] $nm pseudobulk OK" || { echo "[$(date +%H:%M)] $nm PSEUDOBULK FAILED"; return 1; }
  python3 /workspace/composition_counts.py "$h5" "$nm" \
      > /workspace/aging_pbmc/composition/${dir}_ct_counts.csv 2>/workspace/aging_pbmc/composition/${dir}_counts_err.log \
    && echo "[$(date +%H:%M)] $nm composition counts OK" || echo "[$(date +%H:%M)] $nm composition counts FAILED (non-fatal)"
  rm -f "$h5"
  echo "[$(date +%H:%M)] $nm COMPLETE (h5ad deleted)"
}

dl_one () {  # dir filename
  local dir=$1 f=$2
  curl -s --retry 3 --retry-delay 10 -C - -o "$DL/$dir/$f" "$BASE/$dir/files/$f" \
    && echo "[$(date +%H:%M)] downloaded $f" || { echo "[$(date +%H:%M)] DOWNLOAD FAILED: $f"; return 1; }
}

stream_A () {  # big files
  process_one cima    CIMA_reannotated.h5ad CIMA
  process_one aida    AIDA_reannotated.h5ad AIDA
  echo "STREAM_A DONE"
}

stream_B () {
  dl_one soundlife SoundLife_reannotated.h5ad && process_one soundlife $DL/soundlife/SoundLife_reannotated.h5ad SoundLife
  dl_one kcl       KCL_reannotated.h5ad       && process_one kcl       $DL/kcl/KCL_reannotated.h5ad KCL
  dl_one abf300    ABF300_reannotated.h5ad    && process_one abf300    $DL/abf300/ABF300_reannotated.h5ad ABF300
  dl_one ucsf      UCSF_reannotated.h5ad      && process_one ucsf      $DL/ucsf/UCSF_reannotated.h5ad UCSF
  dl_one spac      SPAC_filtered_donors_reannotated.h5ad && process_one spac $DL/spac/SPAC_filtered_donors_reannotated.h5ad SPAC
  echo "STREAM_B DONE"
}

stream_A & stream_B & wait
echo "ALL DOWNLOADS+PROCESSING DONE"
