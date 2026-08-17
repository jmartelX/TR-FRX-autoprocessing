# ============================================================
# setup_env.sh — prepare the environment for trfrx_pipeline.py
#
# SOURCE this file (do not execute it), so the environment stays
# active in your current shell:
#
#     source setup_env.sh
#
# What it does:
#   1. loads every external tool used by the TR-FRX scripts as a module
#      (phenix, CCP4, PyMOL, XDS, autoPROC, H5ToXds);
#   2. if the Python env already exists  -> just activates it;
#      otherwise                         -> creates it, installs all
#                                            dependencies, and activates it.
#
# After sourcing, run any of the linked scripts in the same shell, e.g.:
#     python /path/to/trfrx_pipeline.py --dry-run
#
# Edit MODULES below if a tool has a different name on a given machine.
# ============================================================

VENV="$HOME/.venv/trfrx"
PKGS="numpy pandas gemmi scipy matplotlib dask seaborn h5py"
MODULES="phenix ccp4 pymol xds autoPROC H5ToXds"

# --- make sure we were sourced (otherwise activation would be lost) ---
_sourced=0
if [ -n "${ZSH_VERSION:-}" ]; then
    case "${ZSH_EVAL_CONTEXT:-}" in *:file*) _sourced=1 ;; esac
elif [ -n "${BASH_VERSION:-}" ]; then
    (return 0 2>/dev/null) && _sourced=1
fi
if [ "$_sourced" -eq 0 ]; then
    echo "Please SOURCE this script so the env stays active:"
    echo "    source ${BASH_SOURCE[0]:-setup_env.sh}"
    exit 1
fi

# --- 1. external tools (each loaded independently; missing ones just warn) ---
if type module >/dev/null 2>&1; then
    for _m in $MODULES; do
        if module load "$_m" 2>/dev/null; then
            echo "[trfrx] loaded module: $_m"
        else
            echo "[trfrx] note: module '$_m' not available here (skipped)"
        fi
    done
else
    echo "[trfrx] note: no 'module' command on this machine — skipping module loads."
fi
unset PYTHONHOME PYTHONPATH 2>/dev/null || true

# --- 2. create-or-activate the env ---
if [ -f "$VENV/bin/activate" ]; then
    echo "[trfrx] environment found — activating $VENV"
    source "$VENV/bin/activate"
    # ensure h5py is present in a pre-existing env (added later for Eiger header reading; idempotent)
    if ! python -c "import h5py" 2>/dev/null; then
        echo "[trfrx] h5py missing in existing env — installing ..."
        python -m pip install h5py >/dev/null 2>&1 \
            && echo "[trfrx] h5py installed." \
            || echo "[trfrx] warning: could not install h5py (Eiger header reading will fall back)."
    fi
else
    echo "[trfrx] no environment yet — creating $VENV and installing dependencies ..."
    _base="$(command -v phenix.python || command -v python3 || true)"
    if [ -z "$_base" ]; then
        echo "[trfrx] ERROR: no python found. Load a phenix/python module first." >&2
        return 1
    fi
    "$_base" -m venv "$VENV" || { echo "[trfrx] ERROR: could not create venv." >&2; return 1; }
    source "$VENV/bin/activate"
    python -m pip install --upgrade pip
    python -m pip install $PKGS || { echo "[trfrx] ERROR: pip install failed." >&2; return 1; }
    echo "[trfrx] environment ready."
fi

echo "[trfrx] using python: $(command -v python)"
