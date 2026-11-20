#!/bin/sh
# Run every example and render browser visualizations into viewer/exports/.
#
#   scripts/render_examples.sh                 # everything
#   scripts/render_examples.sh --only grover   # just examples whose name contains "grover"
#   scripts/render_examples.sh --out /tmp/viz  # write somewhere else
#
# Then open viewer/exports/index.html. See scripts/render_examples.py for details.
set -eu

ROOT_DIR="$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)"
cd "$ROOT_DIR"

if ! command -v python3 >/dev/null 2>&1; then
  echo "error: python3 is required to build the HTML pages" >&2
  exit 1
fi

# Always run make: it is a no-op when up to date, and it rebuilds a stale binary
# (an old ./quantum_sim would silently produce fewer visualizations).
echo "Building quantum_sim..."
make quantum_sim >/dev/null

exec python3 scripts/render_examples.py "$@"
