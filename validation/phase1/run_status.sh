#!/usr/bin/env bash
# =============================================================================
# Which validation runs are alive, and what are they doing?
# -----------------------------------------------------------------------------
#   bash run_status.sh          # every run this directory knows about
#   bash run_status.sh v3       # only runs whose tag starts with "v3"
#
# WHY THIS EXISTS -- AND WHY IT DOES NOT USE pgrep.
#
# The obvious way to ask "is the run still going?" is to poll the process table:
#
#     until ! pgrep -f "run_v3_bkmr.R" >/dev/null; do sleep 10; done   # BROKEN
#
# One such loop, on its own, works. Two of them never exit.
#
# `pgrep -f` matches against FULL COMMAND LINES, and each waiting shell's command
# line contains the text `run_v3_bkmr.R`. pgrep does not match its own ancestors,
# so a lone waiter does not find itself -- which is why this bug hides during
# testing. But two waiters are SIBLINGS, not ancestors: each one's pgrep finds
# the OTHER, both conditions stay true forever, and both spin until killed.
# (Verified on this machine: a lone `pgrep -f X` from a shell whose command line
# contains X returns 0 matches; two such sibling shells each return 2.)
#
# The failure is silent and looks exactly like a job that is still running --
# the worst possible disguise for a job that finished an hour ago.
#
# Every runner here already writes `logs/<tag>_<stamp>.pid`, so the question can
# be answered without matching anything:
#
#     kill -0 "$PID"      # true iff a process with that PID exists
#
# `kill -0` sends no signal; it only tests existence and permission. No pattern,
# no sibling collision, no false positive from a grep of the process table.
#
# If you genuinely need a pattern match, bracket one character so that the
# pattern does not appear literally in the command lines doing the matching:
#
#     pgrep -f "[r]un_v3_bkmr.R"
#
# That works because `[r]un_v3_bkmr.R` is a regex matching the text
# `run_v3_bkmr.R`, while the watcher's own command line contains `[r]un...` and
# so is not matched -- by itself or by its siblings. Note the limit: bracketing
# protects against watchers, not against genuinely unrelated processes that
# mention the same filename. Prefer the PID file: PIDs are exact, patterns are
# approximate.
#
# CAVEAT worth knowing: PIDs are reused by the OS. A stale .pid file whose number
# has been recycled will report "alive" for an unrelated process. The log's
# modification time is printed alongside precisely so that a run claiming to be
# alive while its log has not moved in hours is visibly suspicious.
# =============================================================================
set -uo pipefail
cd "$(cd "$(dirname "$0")" && pwd)"

FILTER="${1:-}"

shopt -s nullglob 2>/dev/null || setopt NULL_GLOB 2>/dev/null || true
PIDFILES=(logs/*.pid)

if [[ ${#PIDFILES[@]} -eq 0 ]]; then
  echo "No .pid files in logs/ -- no run has been started from this directory."
  exit 0
fi

printf "%-26s %8s  %-7s  %-19s  %s\n" "RUN" "PID" "STATE" "LOG LAST WRITTEN" "LOG"
printf "%-26s %8s  %-7s  %-19s  %s\n" "---" "---" "-----" "----------------" "---"

alive_count=0
for pf in "${PIDFILES[@]}"; do
  base="$(basename "$pf" .pid)"
  [[ -n "$FILTER" && "$base" != "$FILTER"* ]] && continue

  pid="$(tr -d '[:space:]' < "$pf" 2>/dev/null)"
  log="logs/${base}.log"

  if [[ -n "$pid" ]] && kill -0 "$pid" 2>/dev/null; then
    state="ALIVE"; alive_count=$((alive_count + 1))
  else
    state="done"
  fi

  if [[ -f "$log" ]]; then
    mtime="$(date -r "$log" '+%Y-%m-%d %H:%M:%S' 2>/dev/null || echo '?')"
  else
    mtime="(no log)"; log="-"
  fi

  printf "%-26s %8s  %-7s  %-19s  %s\n" "$base" "${pid:-?}" "$state" "$mtime" "$log"
done

echo
if [[ $alive_count -eq 0 ]]; then
  echo "Nothing running."
else
  echo "$alive_count run(s) alive.  Follow one with:  tail -f logs/<run>.log"
  echo "Stop one with:  kill \$(cat logs/<run>.pid)"
  echo
  echo "If a run is ALIVE but its log has not been written to recently, suspect a"
  echo "recycled PID rather than a working job -- check the log's last lines."
fi
