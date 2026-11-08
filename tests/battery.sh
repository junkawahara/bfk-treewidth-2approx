#!/bin/bash
# randomized battery with --check; every run limited to 32GB / 10 min by run_limited.sh
cd "$(dirname "$0")/.."
D=$(mktemp -d)
fail=0; total=0; skip=0
check_run() { # gr known_tw_upper extra_args...
  local gr=$1 K=$2; shift 2
  total=$((total+1))
  out=$(./run_limited.sh ./bin/tw2approx --check "$@" -o $D/o.td $gr 2>&1); rc=$?
  if [ $rc -eq 3 ] && echo "$out" | grep -q "memory limit"; then skip=$((skip+1)); return; fi
  if [ $rc -ne 0 ]; then echo "FAIL rc=$rc $gr $*"; echo "$out" | tail -3; fail=$((fail+1)); return; fi
  v=$(./bin/verify_td $gr $D/o.td)
  w=$(echo "$out" | sed -n 's/.*result: width=\([0-9-]*\).*k=\([0-9-]*\).*tw >= \([0-9]*\).*/\1 \2 \3/p')
  set -- $w
  if [[ "$v" != OK* ]]; then echo "FAIL invalid TD $gr: $v"; fail=$((fail+1)); return; fi
  if [ $3 -gt $K ]; then echo "FAIL wrong lower bound $3 > tw<=$K on $gr"; fail=$((fail+1)); return; fi
  if [ $1 -gt $((2*$2+1)) ]; then echo "FAIL width $1 > 2k+1 on $gr"; fail=$((fail+1)); return; fi
  if [ $1 -gt $((2*K+1)) ]; then echo "FAIL width $1 > 2tw+1 on $gr"; fail=$((fail+1)); return; fi
}
for s in $(seq 1 50); do
  K=$((1 + s % 3)); N=$((12 + (s*7) % 20)); P=0.$((5 + s % 5))
  ./bin/gen ktree $N $K $P $s > $D/g.gr
  check_run $D/g.gr $K --init random --seed $s --no-lb
done
for rc in "2 5" "3 3" "3 6" "4 4" "2 9"; do
  set -- $rc; ./bin/gen grid $1 $2 > $D/g.gr; T=$(( $1 < $2 ? $1 : $2 ))
  for s in 1 2 3; do check_run $D/g.gr $T --init random --seed $s --no-lb; done
done
for s in $(seq 1 10); do
  ./bin/gen gnm 14 $((14 + s*2)) $s > $D/g.gr   # unknown tw: only validity / 2k+1 checked
  check_run $D/g.gr 99 --init random --seed $s --no-lb
done
echo "battery: total=$total passed=$((total-fail-skip)) skipped(memory guard)=$skip failed=$fail"
rm -rf $D
