#!/usr/bin/env python3
"""Number of states a chain has completed, after making its log files fit the saved state.

BEAST 2 saves the chain state every few log intervals. If a run is stopped from outside (machine switched off),
the trace and tree logs can run past the saved state, or end at different samples, and BEAST then refuses to
resume. This script cuts both logs back to the last sample covered by the saved state and prints that state
number. The state file counts states within one run only, so the number completed before the run is kept in
run_offset.txt (written by gcloud_startup.sh when a run is started).

Usage: python3 chain_state.py <chain directory> <analysis> <seed> <log_every>
"""
import os
import re
import sys
from pathlib import Path

d, name, seed, every = Path(sys.argv[1]), sys.argv[2], sys.argv[3], int(sys.argv[4])
state, offset_file = d / f"{name}.xml.state", d / "run_offset.txt"
if not state.exists():
    print(0)
    sys.exit()
offset = int(offset_file.read_text()) if offset_file.exists() else 0
m = re.search(r"sample='(\d+)'", state.read_text()[:200])
in_run = int(m.group(1)) if m else 0
if offset_file.exists() and os.path.getmtime(state) <= os.path.getmtime(offset_file):
    done = offset                      # stopped before the first save of the last run
else:
    done = offset + in_run
    if (in_run + 1) % every == 0:      # periodic saves are written one state before the log line
        done += 1
done = done // every * every

log = d / f"{name}.{seed}.log"
if log.exists():
    keep = []
    for line in log.read_text().splitlines(keepends=True):
        first = line.split("\t", 1)[0]
        if first.isdigit() and int(first) > done:
            continue
        keep.append(line)
    if keep and not keep[-1].endswith("\n"):
        keep.pop()                     # half-written last line
    log.write_text("".join(keep))
trees = d / f"{name}.{seed}.trees"
if trees.exists():
    keep = []
    for line in trees.read_text().splitlines(keepends=True):
        t = re.match(r"tree STATE_(\d+)", line)
        if t and (int(t.group(1)) > done or not line.rstrip().endswith(";")):
            continue
        keep.append(line)
    trees.write_text("".join(keep))
print(done)
