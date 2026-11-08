#!/bin/bash
# SPDX-License-Identifier: MIT
# randomized battery with --check; every program (gen, tw2approx, verify_td) runs under
# run_limited.sh (32GB / 10 min).  Exit code: 0 if no case failed, 1 otherwise.
cd "$(dirname "$0")/.."
LIM=./run_limited.sh
D=$(mktemp -d)
trap 'rm -rf "$D"' EXIT
fail=0; total=0; skip=0
failed() { echo "FAIL $*"; fail=$((fail+1)); }
graph() { # gen args... -> $D/g.gr; a failing generator counts as one failed case
  if ! $LIM ./bin/gen "$@" > $D/g.gr; then total=$((total+1)); failed "gen $*"; return 1; fi
}
check_run() { # gr K extra_args...   (K: a known upper bound on tw; K = tw for the grids)
  local gr=$1 K=$2; shift 2
  local out rc v vrc w vw
  total=$((total+1))
  out=$($LIM ./bin/tw2approx --check "$@" -o $D/o.td $gr 2>&1); rc=$?
  if [ $rc -eq 3 ] && grep -q "memory limit" <<< "$out"; then skip=$((skip+1)); return; fi
  if [ $rc -ne 0 ]; then failed "rc=$rc $gr $*"; tail -3 <<< "$out"; return; fi
  v=$($LIM ./bin/verify_td $gr $D/o.td); vrc=$?
  if [ $vrc -ne 0 ] || [[ "$v" != "OK width="* ]]; then failed "invalid TD (verify_td rc=$vrc) $gr $*: $v"; return; fi
  vw=${v#OK width=}; vw=${vw%% *}   # width measured by the verifier
  # width, k and lower bound reported by tw2approx: exactly one result line with three integers
  w=$(sed -n 's/^c result: width=\([0-9]*\)  k=\([0-9]*\)  (=> width <= 2k+1, tw >= \([0-9]*\)).*/\1 \2 \3/p' <<< "$out")
  if ! [[ "$vw" =~ ^[0-9]+$ && "$w" =~ ^[0-9]+\ [0-9]+\ [0-9]+$ ]]; then
    failed "cannot read the widths (verify_td: '$v', tw2approx: '$w') $gr $*"; return
  fi
  set -- $w
  if [ $vw -ne $1 ]; then failed "verified width $vw != reported width $1 on $gr"; return; fi
  if [ $3 -gt $K ]; then failed "lower bound $3 > K=$K >= tw on $gr"; return; fi
  if [ $vw -gt $((2*$2+1)) ]; then failed "width $vw > 2k+1 (k=$2) on $gr"; return; fi
  if [ $vw -gt $((2*K+1)) ]; then failed "width $vw > 2K+1 (K=$K >= tw) on $gr"; return; fi
}
for s in $(seq 1 50); do
  K=$((1 + s % 3)); N=$((12 + (s*7) % 20)); P=0.$((5 + s % 5))
  graph ktree $N $K $P $s && check_run $D/g.gr $K --init random --seed $s --no-lb
done
for rc in "2 5" "3 3" "3 6" "4 4" "2 9"; do
  set -- $rc; T=$(( $1 < $2 ? $1 : $2 ))
  graph grid $1 $2 && for s in 1 2 3; do check_run $D/g.gr $T --init random --seed $s --no-lb; done
done
for s in $(seq 1 10); do
  # unknown tw: K = 99 is no bound, so only validity and width <= 2k+1 are checked
  graph gnm 14 $((14 + s*2)) $s && check_run $D/g.gr 99 --init random --seed $s --no-lb
done
echo "battery: total=$total passed=$((total-fail-skip)) skipped(memory guard)=$skip failed=$fail"
[ $fail -eq 0 ]
