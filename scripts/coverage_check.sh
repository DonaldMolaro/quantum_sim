#!/bin/sh
set -e

if ! command -v gcov >/dev/null 2>&1; then
  echo "coverage_check: gcov not found in PATH"
  exit 2
fi

find . -name "*.gcov" -delete

objs=$(find . -name "*.o" -print)
if [ -z "$objs" ]; then
  echo "coverage_check: no object files found"
  exit 2
fi

gcov -b -c $objs >/tmp/gcov.out 2>/dev/null || true

root=$(pwd)
gcov_files=$(find . -name "*.gcov" -print)
if [ -z "$gcov_files" ]; then
  echo "coverage_check: no gcov files produced"
  exit 2
fi

filtered=""
exclude_shor=0
if [ -z "$QSIM_SLOW_TESTS" ] || [ "$QSIM_SLOW_TESTS" = "0" ]; then
  exclude_shor=1
fi

for f in $gcov_files; do
  src=$(awk -F: '/Source:/ {print $4; exit}' "$f")
  [ -z "$src" ] && continue
  resolved="$src"
  case "$resolved" in
    /*) ;;
    *) resolved="$root/$resolved" ;;
  esac
  case "$resolved" in
    "$root"/*)
      if [ -f "$resolved" ]; then
        case "$resolved" in
          */tests/*|*/cli/*|*/scripts/*|*/dist/*)
            ;;
          */all_tests.cc|*/main.cc)
            ;;
          *.cc)
            if [ "$exclude_shor" -eq 1 ] && [ "$resolved" = "$root/algorithms/shor_classical.cc" ]; then
              :
            elif [ "$exclude_shor" -eq 1 ] && [ "$resolved" = "$root/algorithms/shor_quantum.cc" ]; then
              :
            elif [ "$exclude_shor" -eq 1 ] && [ "$resolved" = "$root/demos/shor_demo.cc" ]; then
              :
            elif [ "$exclude_shor" -eq 1 ] && [ "$resolved" = "$root/demos/latin_demo.cc" ]; then
              :
            elif [ "$exclude_shor" -eq 1 ] && [ "$resolved" = "$root/demos/grover_demo.cc" ]; then
              :
            else
              filtered="$filtered $f"
            fi
            ;;
        esac
      fi
      ;;
  esac
 done

if [ -z "$filtered" ]; then
  echo "coverage_check: no project sources matched coverage filter"
  exit 2
fi

summary=$(awk -F: '
  $1 ~ /^[[:space:]]*-/ { next }
  $1 ~ /^[[:space:]]*#/ { total++; next }
  $1 ~ /^[[:space:]]*[0-9]/ { total++; hit++; next }
  END {
    if (total == 0) {
      printf("0 0\n");
    } else {
      printf("%d %d\n", hit, total);
    }
  }
' $filtered)

hit=$(echo "$summary" | awk '{print $1}')
total=$(echo "$summary" | awk '{print $2}')

if [ "$total" -eq 0 ]; then
  echo "coverage_check: no executable lines detected"
  exit 2
fi

pct=$((hit * 100 / total))
printf "coverage_check: %d/%d lines executed (%d%%)\n" "$hit" "$total" "$pct"

if [ "$pct" -ne 100 ]; then
  echo "coverage_check: coverage below 100%"
  exit 1
fi
