#!/bin/bash
# ============================================================
# setup_svd_env.sh — Create and populate the SVD Python environment
#
# Run once:
#   bash setup_svd_env.sh
#
# To reactivate later in a new session:
#   source ~/.venv/svd_pipeline/bin/activate
#   # or, if you use phenix:
#   module load phenix
#   source ~/.venv/svd_pipeline/bin/activate
#
# To run the SVD script inside the environment:
#   source ~/.venv/svd_pipeline/bin/activate
#   ./run_series_diffmaps.py
# ============================================================

set -euo pipefail

VENV_DIR="$HOME/.venv/svd_pipeline"

# ── Check Python ─────────────────────────────
# Prefer phenix.python if available (has gemmi, numpy, pandas already)
if command -v phenix.python &>/dev/null; then
    BASE_PYTHON=$(command -v phenix.python)
    echo "Using phenix.python as base: $BASE_PYTHON"
else
    BASE_PYTHON=$(command -v python3)
    echo "phenix.python not found, using system python3: $BASE_PYTHON"
fi

# ── Create venv ──────────────────────────────
if [ -d "$VENV_DIR" ]; then
    echo "Virtual environment already exists at $VENV_DIR"
    echo "Delete it first to recreate: rm -rf $VENV_DIR"
else
    echo "Creating virtual environment at $VENV_DIR ..."
    "$BASE_PYTHON" -m venv "$VENV_DIR" --system-site-packages
    echo "Done."
fi

# ── Activate ─────────────────────────────────
source "$VENV_DIR/bin/activate"
echo "Environment activated: $VENV_DIR"

# ── Upgrade pip ──────────────────────────────
pip install --quiet --upgrade pip

# ── Install packages ─────────────────────────
echo ""
echo "Installing packages..."

pip install --quiet \
    numpy \
    pandas \
    gemmi \
    "dask[array]" \
    scipy \
    matplotlib \
    seaborn \
    wxmplot

echo ""
echo "Installed packages:"
pip list | grep -E "numpy|pandas|gemmi|dask|scipy|matplotlib|seaborn|wxmplot"

# ── Print reactivation instructions ──────────
echo ""
echo "============================================================"
echo "  Setup complete!"
echo ""
echo "  To reactivate this environment in a new session:"
echo ""
echo "    module load phenix          # if needed for mtzdmp / cad / phenix"
echo "    module load ccp4            # if needed"
echo "    source $VENV_DIR/bin/activate"
echo ""
echo "  Then run your script:"
echo "    ./run_series_diffmaps.py"
echo ""
echo "  To deactivate:"
echo "    deactivate"
echo "============================================================"