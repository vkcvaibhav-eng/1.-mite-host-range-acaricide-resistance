#!/bin/bash
# Runs the BEAST 2 chains on a Google Cloud VM. The VM start-up script (pasted once when the VM is created) is:
#   #!/bin/bash
#   curl -fsSL -o /root/mite_start.sh https://raw.githubusercontent.com/vkcvaibhav-eng/1.-mite-host-range-acaricide-resistance/main/analysis/beast/gcloud_startup.sh
#   bash /root/mite_start.sh
# First boot: installs BEAST, downloads the repository, starts the chains.
# Later boots: chains continue from their saved state (-resume) until CHAIN_LENGTH in run_config.env is reached.
# When all chains are finished: convergence summary, results packed in /opt/mite/results.tar.gz, machine off.
# Progress is written every 15 minutes to /opt/mite/status.log and to the serial console (lines "[mite-beast]").
set -uo pipefail
WORK=${MITE_WORK:-/opt/mite}
REPO_URL=https://github.com/vkcvaibhav-eng/1.-mite-host-range-acaricide-resistance.git
mkdir -p "$WORK"; cd "$WORK"
say() { echo "[mite-beast $(date '+%F %T')] $*" | tee -a "$WORK/status.log"; [ -w /dev/console ] && echo "[mite-beast] $*" > /dev/console 2>/dev/null; true; }

if [ -z "${MITE_LOCAL_REPO:-}" ]; then
  export DEBIAN_FRONTEND=noninteractive
  apt-get update -qq >/dev/null 2>&1; apt-get install -y -qq git curl python3 >/dev/null 2>&1
  if [ -d repo/.git ]; then git -C repo fetch -q origin main && git -C repo reset -q --hard origin/main
  else git clone -q --depth 1 "$REPO_URL" repo; fi
  REPO="$WORK/repo"
else
  REPO="$MITE_LOCAL_REPO"          # test mode: use a local copy, no apt, no clone
fi
# shellcheck disable=SC1091
source "$REPO/analysis/beast/run_config.env"
[ -n "${MITE_CHAIN_LENGTH:-}" ] && CHAIN_LENGTH=$MITE_CHAIN_LENGTH
[ -n "${MITE_SHUTDOWN:-}" ] && SHUTDOWN_WHEN_DONE=$MITE_SHUTDOWN
[ -n "${MITE_LOG_EVERY:-}" ] && LOG_EVERY=$MITE_LOG_EVERY
[ -n "${MITE_SEEDS:-}" ] && SEEDS=$MITE_SEEDS
say "commit $(git -C "$REPO" rev-parse --short HEAD 2>/dev/null); chain length $CHAIN_LENGTH; seeds $SEEDS"

BEAST=${MITE_BEAST:-$WORK/beast/bin}
if [ ! -x "$BEAST/beast" ]; then
  curl -fsSL -o beast.tgz "$BEAST_URL" || { say "ERROR: BEAST download failed"; exit 1; }
  echo "$BEAST_SHA256  beast.tgz" | sha256sum -c - >/dev/null || { say "ERROR: BEAST checksum does not match"; exit 1; }
  tar -xzf beast.tgz && rm beast.tgz
fi

mkdir -p runs
pids=()
# states completed by a chain; also cuts its logs back to the saved state (see chain_state.py)
done_states() { python3 "$REPO/analysis/beast/chain_state.py" "runs/chain_$1" "$ANALYSIS" "$1" "$LOG_EVERY"; }
for s in $SEEDS; do
  d="runs/chain_$s"; mkdir -p "$d"
  have=$(done_states "$s"); have=${have:-0}
  left=$((CHAIN_LENGTH - have))
  if [ "$left" -le 0 ]; then say "chain $s already has $have states"; continue; fi
  # BEAST -resume runs chainLength further states, so the input file is written with what is left
  python3 "$REPO/scripts/05e_beast_xml.py" "$ANALYSIS" --chain "$left" --log-every "$LOG_EVERY" >/dev/null || { say "ERROR: could not write the input file"; exit 1; }
  cp "$REPO/analysis/beast/$ANALYSIS.xml" "$d/$ANALYSIS.xml"
  if [ "$have" -gt 0 ]; then mode=-resume; else mode=-overwrite; rm -f "$d/$ANALYSIS.xml.state"; fi
  echo "$have" > "$d/run_offset.txt"
  ( cd "$d" && "$BEAST/beast" -seed "$s" -threads "$THREADS_PER_CHAIN" -statefile "$ANALYSIS.xml.state" $mode "$ANALYSIS.xml" >> "beast_$s.out" 2>&1 ) &
  pids+=($!)
  say "chain $s started ($mode, $have states done, $left to go)"
done

progress() {
  line=""
  for s in $SEEDS; do
    last=$(tail -n 1 "runs/chain_$s/$ANALYSIS.$s.log" 2>/dev/null | cut -f1)
    line="$line  $s: ${last:-0}"
  done
  say "states done of $CHAIN_LENGTH -$line"
}
( while sleep "${MITE_PROGRESS_SECONDS:-900}"; do progress; done ) &
watcher=$!
[ ${#pids[@]} -gt 0 ] && wait "${pids[@]}"
pkill -P $watcher 2>/dev/null; kill $watcher 2>/dev/null
progress

failed=0
for s in $SEEDS; do
  have=$(done_states "$s"); have=${have:-0}
  if [ "$have" -lt "$CHAIN_LENGTH" ]; then failed=1; say "ERROR: chain $s stopped at $have states; see runs/chain_$s/beast_$s.out"; tail -n 5 "runs/chain_$s/beast_$s.out" | tee -a status.log; fi
done

# convergence summary and one file to download
mkdir -p results
for s in $SEEDS; do
  "$BEAST/loganalyser" -b "$BURNIN_PERCENT" "runs/chain_$s/$ANALYSIS.$s.log" > "results/ess_chain_$s.txt" 2>&1
  cp "runs/chain_$s/$ANALYSIS.$s.log" results/
  "$BEAST/logcombiner" -log "runs/chain_$s/$ANALYSIS.$s.trees" -o "results/$ANALYSIS.$s.thinned.trees" -b "$BURNIN_PERCENT" -resample 100000 > "results/logcombiner_$s.out" 2>&1
  low=$(awk '$1=="item"{h=1;next} h && $9 ~ /^[0-9.]+$/ && $9 < 200 {n++} END {print n+0}' "results/ess_chain_$s.txt")
  min=$(awk '$1=="item"{h=1;next} h && $9 ~ /^[0-9.]+$/ {if (m=="" || $9+0 < m+0) {m=$9; it=$1}} END {print it " " m}' "results/ess_chain_$s.txt")
  post=$(awk '$1=="posterior" {print $9}' "results/ess_chain_$s.txt")
  say "chain $s: ESS of posterior $post; lowest ESS $min; $low logged quantities below 200 (burn-in $BURNIN_PERCENT%)"
done
cp status.log "runs/chain_${SEEDS%% *}/$ANALYSIS.xml" "$REPO/analysis/beast/run_config.env" results/
tar -czf results.tar.gz results && say "results packed: $WORK/results.tar.gz ($(du -h results.tar.gz | cut -f1))"

# switch off only if chains ran during this boot; after a restart to download results the machine stays on
if [ "$failed" = 0 ] && [ "$SHUTDOWN_WHEN_DONE" = yes ] && [ ${#pids[@]} -gt 0 ]; then say "all chains finished; switching the machine off"; shutdown -h now
elif [ "$failed" = 0 ]; then say "all chains finished; machine left on so the results can be downloaded - stop it afterwards"; fi
