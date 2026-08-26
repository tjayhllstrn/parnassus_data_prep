#!/usr/bin/env bash
# Sourced helper: make ROOT + PyROOT available for the root-inspection skill.
#
# ROOT is NOT on the bare cluster PATH. It comes from the parnassus-core conda
# environment, which is activated by sourcing:
#     ../parnassus-core/setup.sh
# (override with PARNASSUS_CORE_SETUP=/path/to/setup.sh)
#
# This helper is idempotent: if `root` and PyROOT already work, it does nothing.

_ri_have_root() {
  command -v root >/dev/null 2>&1 && python -c 'import ROOT' >/dev/null 2>&1
}

ensure_root_env() {
  if _ri_have_root; then
    return 0
  fi

  local setup="${PARNASSUS_CORE_SETUP:-}"
  if [[ -z "$setup" ]]; then
    setup="$(cd "${REPO_ROOT:-.}/.." 2>/dev/null && pwd)/parnassus-core/setup.sh"
  fi

  if [[ ! -f "$setup" ]]; then
    echo "ERROR: cannot find the ROOT setup script: $setup" >&2
    echo "       Set PARNASSUS_CORE_SETUP to the parnassus-core setup.sh path." >&2
    return 1
  fi

  # conda/module init scripts are not -u/-e safe; relax the shell while sourcing.
  local had_u=0 had_e=0
  [[ $- == *u* ]] && had_u=1
  [[ $- == *e* ]] && had_e=1
  set +u
  set +e
  # shellcheck disable=SC1090
  source "$setup" >/dev/null 2>&1
  [[ $had_u == 1 ]] && set -u
  [[ $had_e == 1 ]] && set -e

  if ! _ri_have_root; then
    echo "ERROR: ROOT/PyROOT still unavailable after sourcing $setup" >&2
    echo "       Try in an interactive shell:  source $setup && root-config --version" >&2
    return 1
  fi
  return 0
}

ensure_root_env || exit 1
