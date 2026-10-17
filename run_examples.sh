#!/usr/bin/env bash
#
# Run the redistributable examples and compare the final log-likelihood
# with the reference values in examples/reference_values.tsv.
#
#   bash run_examples.sh                       # uses ./local/bin/amlet
#   AMLET=$HOME/amlet/bin/amlet bash run_examples.sh
#   RUN=/tmp/amlet-run bash run_examples.sh    # where outputs are written
#
# The Swissmetro data are not shipped: put swissmetro.dat (from
# https://transp-or.epfl.ch/data/swissmetro.dat) in examples/swissmetro/,
# or set SWISSMETRO_DAT=/path/to/swissmetro.dat; the script tries to
# download it otherwise, and skips the Swissmetro examples if it cannot.

set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AMLET="${AMLET:-$ROOT/local/bin/amlet}"
RUN="${RUN:-$ROOT/run}"
SWISSMETRO_URL="https://transp-or.epfl.ch/data/swissmetro.dat"

[ -x "$AMLET" ] || { echo "ERROR: AMLET executable not found at $AMLET (run build.sh, or set AMLET=...)" >&2; exit 2; }
mkdir -p "$RUN"
RUN="$(cd "$RUN" && pwd)"

prepare_swissmetro() {
  local dir="$RUN/swissmetro" dat="${SWISSMETRO_DAT:-$ROOT/examples/swissmetro/swissmetro.dat}"
  [ -s "$dir/swissmetro_train96" ] && return 0
  if [ ! -s "$dat" ]; then
    dat="$dir/swissmetro.dat"
    echo "Downloading Swissmetro data from $SWISSMETRO_URL"
    curl -fsSL -m 120 -o "$dat" "$SWISSMETRO_URL" 2>/dev/null \
      || wget -q -T 120 -O "$dat" "$SWISSMETRO_URL" 2>/dev/null \
      || { rm -f "$dat"; return 1; }
  fi
  python3 "$ROOT/examples/swissmetro/swissmetro_to_amlet.py" "$dat" "$dir/swissmetro_train96"
}

pass=0; fail=0; skip=0
for ex in australia simulated swissmetro; do
  mkdir -p "$RUN/$ex"
  cp -p "$ROOT/examples/$ex/"* "$RUN/$ex/"
done

if ! prepare_swissmetro; then
  echo "Swissmetro data unavailable: skipping the Swissmetro examples."
  skip_swiss=1
else
  skip_swiss=0
fi

printf '\n%-11s %-13s %18s %18s  %s\n' example driver expected obtained status
while IFS=$'\t' read -r ex driver output expected tol source; do
  case "$ex" in ''|\#*) continue ;; esac
  if [ "$ex" = swissmetro ] && [ "$skip_swiss" = 1 ]; then
    printf '%-11s %-13s %18s %18s  %s\n' "$ex" "$driver" "$expected" "-" SKIPPED
    skip=$((skip+1)); continue
  fi
  (cd "$RUN/$ex" && "$AMLET" "$driver" > "$driver.log" 2>&1)
  rc=$?
  got="$(awk '/^Optimal/ {print $2; exit}' "$RUN/$ex/$output" 2>/dev/null)"
  status="$(python3 -c '
import sys
exp, got, tol = sys.argv[1], sys.argv[2], float(sys.argv[3])
try:
    ok = abs(float(got) - float(exp)) <= tol * abs(float(exp))
except ValueError:
    ok = False
print("OK" if ok else "FAILED")' "$expected" "${got:-nan}" "$tol")"
  [ "$rc" -ne 0 ] && status="FAILED (exit code $rc)"
  printf '%-11s %-13s %18s %18s  %s\n' "$ex" "$driver" "$expected" "${got:--}" "$status"
  case "$status" in OK) pass=$((pass+1)) ;; *) fail=$((fail+1)) ;; esac
done < "$ROOT/examples/reference_values.tsv"

printf '\n%d passed, %d failed, %d skipped. Outputs in %s\n' "$pass" "$fail" "$skip" "$RUN"
[ "$fail" -eq 0 ]
