#!/bin/bash
# Tier-2 (interacting) campaign driver  -- WP-L.
#
# For each ensemble in the priority list: run rhmc with resume-from-checkpoint
# (the committed trajectory counter is read back through rmeas -ckptInfo, which
# validates the whole manifest; the run is skipped once the target is reached),
# then rmeas per saved configuration (per-config TSVs are skipped when they
# already exist, so this is restartable too), then rmeas -analyze.
#
# Usage:
#   campaign/t2.sh                  # everything, in priority order
#   campaign/t2.sh L1g15m00 pureL1  # only the named ensembles
#   DRYRUN=1 campaign/t2.sh         # print the commands without running
#   RADIAL_T2_OUT=/new/path campaign/t2.sh  # explicit fresh campaign root
#
# This standard-overlap campaign MUST NOT reuse the legacy output/radial/t2
# checkpoints or measurements, which used additive masses.  Changing a mass,
# rational window, or other frozen manifest field likewise requires a new
# RADIAL_T2_OUT; the old directory is preserved as an immutable campaign.
#
# Sizing (single core, measured 2026-08-21; see doc/06 WP-L):
#   L=1, nt=60, N_f=2 HMC: 41.5 s/trajectory alone (68-96 s with 3-4 runs
#   sharing the machine; m>0 is cheaper).  L=2: 424-479 s/trajectory measured.
#   Measurement per L=1 config (m=0, thermalized): currents+disc+ward 50-115 s,
#   scalars 30-50 s (now including the volume estimators), cond 140 s (m=0;
#   several x cheaper at m>=0.1), gluon 0.03 s, wspec 3.9 s.
#
# Trajectory targets below are the preliminary (laptop) profile; the
# production values are in the comments of each ensemble.  The legacy
# additive-mass campaign in output/radial/t2 was never completed under this
# script and is not resumable with the current binaries.

set -euo pipefail

WORKTREE="$(cd "$(dirname "$0")"/../../../.. && pwd)"
BUILD="$WORKTREE/build_mac"
BIN="$BUILD/bin"
OUT="${RADIAL_T2_OUT:-$WORKTREE/output/radial/t2-standard-overlap}"
export SDKROOT="${SDKROOT:-/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk}"
DRYRUN="${DRYRUN:-0}"
NNOISE="${NNOISE:-8}"         # condensate noise vectors (the stored legacy files used 8)

run() {
  echo "+ $*"
  if [ "$DRYRUN" != 1 ]; then "$@"; fi
}

# side effects (directories, logs) only outside a dry run
prep() {
  if [ "$DRYRUN" != 1 ]; then mkdir -p "$1"; fi
}

logto() {
  if [ "$DRYRUN" != 1 ]; then tee -a "$1"; else cat; fi
}

build() {
  ( cd "$BUILD" && run make experimental/radial/rhmc \
    && run make experimental/radial/rmeas )
}

# committed trajectory counter of an ensemble (0 if no checkpoint); the read
# goes through loadCheckpoint, so a manifest mismatch aborts loudly here
# (rmeas exits 1 and prints "checkpoint mismatch: <field>").
ckpt_traj() {
  local dir=$1; shift
  if [ ! -f "$dir/ckpt" ]; then echo 0; return; fi
  local out
  out=$("$BIN/rmeas" "$@" -dir:"$dir" -ckptInfo:true | sed -n 's/^ckptinfo: traj: //p')
  if [ -z "$out" ]; then
    echo "ckpt_traj: no counter reported for $dir (manifest mismatch?)" >&2
    exit 1
  fi
  echo "$out"
}

# run_hmc <dir> <target_traj> <manifest+hmc args...>
run_hmc() {
  local dir=$1 target=$2; shift 2
  prep "$dir"
  local have remaining
  have=$(ckpt_traj "$dir" "$@")
  remaining=$(( target - have ))
  if [ "$remaining" -le 0 ]; then
    echo "== $dir: $have/$target trajectories -- HMC complete, skipping"
    return
  fi
  echo "== $dir: $have/$target trajectories -- running $remaining more"
  run "$BIN/rhmc" "$@" -ntraj:"$remaining" \
      -ckpt:"$dir/ckpt" -cfg:"$dir/cfg" 2>&1 | logto "$dir/hmc.log"
}

# run_meas <dir> <tmin> <obs> <manifest args...> [extra rmeas args appended]
run_meas() {
  local dir=$1 tmin=$2 obsv=$3; shift 3
  run "$BIN/rmeas" "$@" -dir:"$dir" -obs:"$obsv" -tmin:"$tmin" \
      2>&1 | logto "$dir/meas.log"
}

run_analyze() {
  local dir=$1 tmin=$2; shift 2
  run "$BIN/rmeas" "$@" -dir:"$dir" -obs: -analyze:true -tmin:"$tmin" \
      2>&1 | logto "$dir/analyze.log"
}

# ---------------------------------------------------------------------------
# ensembles.  MANIFEST args are shared verbatim by rhmc and rmeas -- that is
# what makes the checkpoint validation effective.  a_t = 0.2, M = 1,
# rationals 31/11 on [0.3, 12.5] (widened to [0.15, 14] for g2R = 3, where
# WP-H expects sigma_min well below the g2R=1 value 0.58; maxRelErr(11) there
# is ~4.5e-4 -- an MD-force error only, corrected by the order-31 Metropolis
# test at some acceptance cost).  Masses are STANDARD overlap masses,
# D(m)=(1-m/2)D_ov+m (rho=1), in units of R (R = 1): the deck's mR is m here.
# The legacy campaign's -masses were additive mu; mu = 0.1 corresponds to
# m = 0.0952 (doc/02 section 4).  Hasenbusch ladder [m_sea, 0.5], tau = 1,
# steps 4 x innerSteps 5 -- exactly WP-H's demo settings.
# ---------------------------------------------------------------------------

ens_L1g15m00() {  # T2.4 T2.6 T2.9 T2.1 T2.7 T2.8 (+T2.5 via disc); anchor point
  local dir="$OUT/L1g15m00"
  local args=(-lev:1 -nt:60 -at:0.2 -g2R:1.5 -nf:2 -masses:0.0,0.5 -seed:1001)
  # preliminary profile: 80 trajectories = 14 measured configs (production: 170+)
  run_hmc "$dir" 80 "${args[@]}" -warmup:10 -measEvery:5 -ckptFreq:10 -windowEvery:10
  run_meas "$dir" 11 currents,scalars,gluon,wspec "${args[@]}" -disc:true
  run_analyze "$dir" 11 "${args[@]}" -disc:true
}

cond_scan() {      # T2.3: condensate at the sea mass, published slope points
  local m=$1 seed=$2
  local tag="L1g15m0${m/0./}"
  local ntraj=60   # preliminary profile (production: 300+)
  local dir="$OUT/$tag"
  local args=(-lev:1 -nt:60 -at:0.2 -g2R:1.5 -nf:2 -masses:"$m",0.5 -seed:"$seed")
  run_hmc "$dir" "$ntraj" "${args[@]}" -warmup:10 -measEvery:5 -ckptFreq:10 -windowEvery:10
  run_meas "$dir" 11 cond "${args[@]}" -nnoise:"$NNOISE"
  run_analyze "$dir" 11 "${args[@]}" -nnoise:"$NNOISE"
}

ens_L1g05m00() {   # the g2R trend, weak end (preliminary: 60 traj)
  local dir="$OUT/L1g05m00"
  local args=(-lev:1 -nt:60 -at:0.2 -g2R:0.5 -nf:2 -masses:0.0,0.5 -seed:1006)
  run_hmc "$dir" 60 "${args[@]}" -warmup:10 -measEvery:5 -ckptFreq:10 -windowEvery:10
  run_meas "$dir" 11 currents,scalars,gluon,wspec "${args[@]}" -disc:true
  run_analyze "$dir" 11 "${args[@]}" -disc:true
}

ens_L1g10m00() {
  local dir="$OUT/L1g10m00"
  local args=(-lev:1 -nt:60 -at:0.2 -g2R:1.0 -nf:2 -masses:0.0,0.5 -seed:1007)
  run_hmc "$dir" 100 "${args[@]}" -warmup:10 -measEvery:5 -ckptFreq:10 -windowEvery:10
  run_meas "$dir" 11 currents,scalars,gluon,wspec "${args[@]}" -disc:true
  run_analyze "$dir" 11 "${args[@]}" -disc:true
}

ens_L1g10nf4() {   # the N_f trend (slide 12/14); ~79 s/traj expected (2 pf copies)
  local dir="$OUT/L1g10nf4"
  local args=(-lev:1 -nt:60 -at:0.2 -g2R:1.0 -nf:4 -masses:0.0,0.5 -seed:1008)
  # 100 traj ~ 2.2 h; 18 measured configs -- thin statistics, trend-level only
  run_hmc "$dir" 100 "${args[@]}" -warmup:10 -measEvery:5 -ckptFreq:10 -windowEvery:10
  run_meas "$dir" 11 currents,scalars,gluon "${args[@]}" -disc:true
  run_analyze "$dir" 11 "${args[@]}" -disc:true
}

ens_L1g10nf6() {   # ~116 s/traj expected (3 pf copies)
  local dir="$OUT/L1g10nf6"
  local args=(-lev:1 -nt:60 -at:0.2 -g2R:1.0 -nf:6 -masses:0.0,0.5 -seed:1009)
  # 70 traj ~ 2.3 h; 12 measured configs -- trend-level only
  run_hmc "$dir" 70 "${args[@]}" -warmup:10 -measEvery:5 -ckptFreq:10 -windowEvery:10
  run_meas "$dir" 11 currents,scalars,gluon "${args[@]}" -disc:true
  run_analyze "$dir" 11 "${args[@]}" -disc:true
}

ens_L2g30m00() {   # ONE L=2 demonstration point (T2.6's l=3 splitting, (a/R)^2)
  local dir="$OUT/L2g30m00"
  # WIDENED WINDOW [0.15, 14]: sigma_min at g2R=3 is unknown until windowCheck
  # measures it.  If rhmc still hard-stops, widen only in a NEW
  # RADIAL_T2_OUT (for example .../t2-standard-overlap-window2); never resume
  # the stopped checkpoint with a changed rational.
  local args=(-lev:2 -nt:60 -at:0.2 -g2R:3.0 -nf:2 -masses:0.0,0.5 -seed:1010
              -ratLo:0.15 -ratHi:14.0)
  # preliminary profile: 48 traj (~6 h at 424-479 s/traj) -> 10 measured
  # configs; the deck-level L=2 statistics (the 3% l=3 splitting) are a cluster job
  run_hmc "$dir" 48 "${args[@]}" -warmup:8 -measEvery:4 -ckptFreq:5 -windowEvery:5
  run_meas "$dir" 9 currents,scalars,gluon "${args[@]}" -disc:true
  # dense wspec at L=2 (dim 5040) costs minutes per configuration -- measure
  # every 3rd config (3-4 of 10), which is also all the deck used (T2.1)
  run_meas "$dir" 9 wspec "${args[@]}" -tstride:3
  run_analyze "$dir" 9 "${args[@]}" -disc:true
}

pure_ens() {       # exact-heatbath pure gauge: slide-9 free baselines (T2.2)
  local L=$1 nconf=$2 seed=$3
  local dir="$OUT/pureL$L"
  local g2R
  g2R=$(echo "$L" | awk '{printf "%.1f", $1*1.0}')   # g2a = 1.0 row
  local args=(-lev:"$L" -nt:60 -at:0.2 -g2R:"$g2R" -nf:0 -masses:0.0 -seed:"$seed")
  prep "$dir"
  run_meas "$dir" 0 gluon "${args[@]}" -pure:true -nconf:"$nconf"
  run_analyze "$dir" 0 "${args[@]}" -pure:true
}

# ---------------------------------------------------------------------------

# the default list is the N_f = 2 preliminary profile; the N_f = 4, 6 ensembles
# (L1g10nf4, L1g10nf6: the deck's N_f trend, below laptop statistics) and
# pureL4 are run by naming them explicitly
all=(L1g15m00 L1g15m01 L1g15m02 L1g15m03 L1g15m04
     L1g05m00 L1g10m00 L2g30m00 pureL1 pureL2)
want=("${@:-${all[@]}}")

build

for e in "${want[@]}"; do
  echo ""
  echo "############ $e ############"
  case "$e" in
    L1g15m00) ens_L1g15m00 ;;
    L1g15m01) cond_scan 0.1 1002 ;;
    L1g15m02) cond_scan 0.2 1003 ;;
    L1g15m03) cond_scan 0.3 1004 ;;
    L1g15m04) cond_scan 0.4 1005 ;;
    L1g05m00) ens_L1g05m00 ;;
    L1g10m00) ens_L1g10m00 ;;
    L1g10nf4) ens_L1g10nf4 ;;
    L1g10nf6) ens_L1g10nf6 ;;
    L2g30m00) ens_L2g30m00 ;;
    pureL1)   pure_ens 1 256 2001 ;;
    pureL2)   pure_ens 2 256 2002 ;;
    pureL4)   pure_ens 4 128 2003 ;;
    *) echo "unknown ensemble: $e" >&2; exit 1 ;;
  esac
done

echo ""
echo "campaign pass complete.  Per-ensemble summaries: $OUT/*/analysis/summary.tsv"
