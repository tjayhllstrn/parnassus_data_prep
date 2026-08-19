#!/usr/bin/env bash
# Entry point for the root-inspection skill.
#
#   skills/root-inspection/scripts/root_inspect.sh <subcommand> [args...]
#
# Subcommands: info | count | select | dump | plot | plot-all
# Run with -h (or a subcommand + -h) for full usage.
#
# This wrapper only activates the ROOT environment and forwards to
# root_tools.py; all logic (and all safety checks) live there.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
# Note: parnassus-core/setup.sh also exports REPO_ROOT (pointing at
# parnassus-core), so root_tools.py keys off PARNASSUS_DATA_PREP_ROOT instead.
REPO_ROOT="${PARNASSUS_DATA_PREP_ROOT:-$(cd "$SCRIPT_DIR/../../.." && pwd -P)}"
export REPO_ROOT

# shellcheck source=_env.sh
source "$SCRIPT_DIR/_env.sh"

exec python "$SCRIPT_DIR/root_tools.py" "$@"
