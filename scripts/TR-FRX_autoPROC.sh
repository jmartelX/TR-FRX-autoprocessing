#!/bin/bash
#
# TR-FRX autoPROC chunked processing + consolidated reports.
#
# Handles BOTH input types automatically:
#   * CBF   : a template with #### placeholders, e.g. SAMPLE_1_####.cbf
#   * HDF5  : an Eiger master file,               e.g. SAMPLE_001_master.h5
# The mode is detected from the file extension (override via the 1st argument).
#
# After all chunks finish, a final SLURM job regroups the final statistics into
#   autoproc_chunks/reports/
#     staraniso_report.pdf   staraniso_statistics.csv   (STARANISO route)
#     truncate_report.pdf    truncate_statistics.csv    (classical autoPROC)
#     parsing_diagnostics.txt
# Each PDF holds radiation-damage plots vs image number (high-resolution limit,
# R-factors, I/sigma, CC(1/2), Wilson B, unit-cell edges & volume, images used)
# plus per-chunk statistics tables, and records which image ranges were used.
#
# Usage:
#   # Old way (run from the directory containing the images / master file):
#   ./TR-FRX_autoPROC.sh                           # use IMAGE_TEMPLATE set below
#   ./TR-FRX_autoPROC.sh SAMPLE_001_master.h5   # override the input file
#
#   # Read images from one place, write results to another:
#   ./TR-FRX_autoPROC.sh <input_dir> <output_dir>
#   ./TR-FRX_autoPROC.sh <input_dir> <output_dir> SAMPLE_001_master.h5
#
#   # Reports ONLY — regenerate the statistics from chunks already processed,
#   # without re-running autoPROC. Asks for an image limit (Enter = all):
#   ./TR-FRX_autoPROC.sh --stats                   # autoproc_chunks in ./
#   ./TR-FRX_autoPROC.sh --stats <dir>             # dir = output dir OR the
#                                                  #   autoproc_chunks folder
#   ./TR-FRX_autoPROC.sh --stats --max-image 2000  # skip the prompt
#
# --max-image N keeps only chunks that END at or before N (a chunk straddling N
# is excluded, never truncated): N=2000 keeps ...1801-2000 and drops 2001-2200.
# Limited reports are written with an _upto<N> suffix, so the full-series report
# is never overwritten.

set -euo pipefail

# =========================== CONFIGURATION ===========================
RUN_ID=006

# CBF template (with ####) OR Eiger master .h5 file. Auto-detected below.
IMAGE_TEMPLATE="SAMPLE_1_####.cbf"

FIRST_IMG=1
LAST_IMG=3000

REF_FIRST_IMG=1
REF_LAST_IMG=300
CHUNK_SIZE=300

# Space group and unit cell — EDIT these to match your crystal (example values).
SYMM="I41"
CELL="114 114 118 90 90 90"

# --- Real per-dataset acquisition time -------------------------------------
# The report computes a true time window [t_start, t_end] for every dataset and
# writes reports/time_windows.csv. Leave FRAME_TIME_MS / OSC_PER_IMAGE blank to
# read them from the image header (CBF "Exposure_period" / Eiger "frame_time");
# set them only to override when the header is missing or wrong.
#
# Time origin: image T0_IMAGE happens at T0_SECONDS, so t(i)=T0_SECONDS+(i-T0_IMAGE)*frame_time.
# Leave T0_IMAGE BLANK (default) => t=0 at the END of the reference chunk, i.e. the
# reference dataset is 0 s and every later chunk's time is the time at its LAST image
# ((image_last - ref_last_image) * frame_time). Set T0_IMAGE/T0_SECONDS explicitly only
# for a pump/mix experiment whose trigger is a known image at a known clock time.
FRAME_TIME_MS=""      # per-image frame period in ms (incl. dead time)
OSC_PER_IMAGE=""      # oscillation per image in degrees
T0_IMAGE=""           # blank = end of reference chunk (t0=0 there); or a trigger image index
T0_SECONDS=0          # absolute seconds assigned to T0_IMAGE

# SLURM partition — EDIT to match your cluster.
SLURM_PARTITION="nice"
SLURM_CPUS=24
SLURM_MEM=24000
# =====================================================================

# --------------------------- Command-line arguments ---------------------------
# Four ways to launch the script:
#   ./TR-FRX_autoPROC.sh                                     # legacy mode:
#                                                            #   images AND results
#                                                            #   in the current dir
#   ./TR-FRX_autoPROC.sh master.h5                           # + override input file
#   ./TR-FRX_autoPROC.sh <input_dir> <output_dir>           # read here, write there
#   ./TR-FRX_autoPROC.sh <input_dir> <output_dir> master.h5 # + override input file
# -----------------------------------------------------------------------------
IMAGES_DIR="$(pwd)"
OUTPUT_DIR="$(pwd)"

# --- Options: --stats (reports only) and --max-image N (image limit) ----------
STATS_ONLY=0
MAX_IMAGE=""
POSITIONAL=()
while [ "$#" -gt 0 ]; do
    case "$1" in
        --stats|--stats-only|--reports-only)
            STATS_ONLY=1; shift ;;
        --max-image|--max-images)
            MAX_IMAGE="${2:-}"
            if [ -z "$MAX_IMAGE" ]; then
                echo "ERROR: --max-image requires a value (e.g. --max-image 2000)" >&2
                exit 1
            fi
            shift 2 ;;
        --max-image=*|--max-images=*)
            MAX_IMAGE="${1#*=}"; shift ;;
        # --- Timing / image-header overrides (work in --stats too) -------------
        --frame-time-ms)   FRAME_TIME_MS="${2:-}"; shift 2 ;;
        --frame-time-ms=*) FRAME_TIME_MS="${1#*=}"; shift ;;
        --osc-per-image)   OSC_PER_IMAGE="${2:-}"; shift 2 ;;
        --osc-per-image=*) OSC_PER_IMAGE="${1#*=}"; shift ;;
        --images-dir)      IMAGES_DIR="${2:-}"; shift 2 ;;
        --images-dir=*)    IMAGES_DIR="${1#*=}"; shift ;;
        --image-template)   IMAGE_TEMPLATE="${2:-}"; shift 2 ;;
        --image-template=*) IMAGE_TEMPLATE="${1#*=}"; shift ;;
        --t0-image)   T0_IMAGE="${2:-}"; shift 2 ;;
        --t0-image=*) T0_IMAGE="${1#*=}"; shift ;;
        --t0-seconds)   T0_SECONDS="${2:-}"; shift 2 ;;
        --t0-seconds=*) T0_SECONDS="${1#*=}"; shift ;;
        -h|--help)
            # Script header up to the first non-comment line
            # (no hard-coded line number: stays correct if the help text grows).
            awk 'NR==1 {next} /^#/ {sub(/^# ?/, ""); print; next} {exit}' \
                "${BASH_SOURCE[0]}"
            exit 0 ;;
        -*)
            echo "ERROR: unknown option '$1' (see --help)" >&2
            exit 1 ;;
        *)
            POSITIONAL+=("$1"); shift ;;
    esac
done
set -- ${POSITIONAL[@]+"${POSITIONAL[@]}"}

if [ -n "$MAX_IMAGE" ] && ! printf '%s' "$MAX_IMAGE" | grep -Eq '^[0-9]+$'; then
    echo "ERROR: --max-image must be a positive integer (got '$MAX_IMAGE')" >&2
    exit 1
fi

if [ "$STATS_ONLY" -eq 1 ]; then
    # Reports mode: a single optional argument = output directory OR the
    # autoproc_chunks folder itself. Statistics need no images; time_windows.csv,
    # however, needs the frame cadence: either --frame-time-ms, or
    # --images-dir/--image-template (header reading), or the values from the
    # CONFIGURATION block at the top of the script.
    case "$#" in
        0) : ;;
        1) OUTPUT_DIR="$1" ;;
        *)
            echo "Usage: $0 --stats [<dir>] [--max-image N] \\" >&2
            echo "          [--frame-time-ms MS] [--osc-per-image DEG] \\" >&2
            echo "          [--images-dir DIR] [--image-template TPL] \\" >&2
            echo "          [--t0-image N] [--t0-seconds S]" >&2
            exit 1 ;;
    esac
else
    case "$#" in
        0) : ;;                                       # all defaults (current dir)
        1) IMAGE_TEMPLATE="$1" ;;                      # override input file only
        2) IMAGES_DIR="$1"; OUTPUT_DIR="$2" ;;         # read from $1, write to $2
        3) IMAGES_DIR="$1"; OUTPUT_DIR="$2"; IMAGE_TEMPLATE="$3" ;;
        *)
            echo "Usage: $0 [<input_dir> <output_dir>] [image_file]" >&2
            exit 1
            ;;
    esac
fi

# --- Absolute paths (SLURM compute nodes must be able to see these dirs) ---
# In --stats mode no image is read, so we only check the output directory
# (which must already exist; we don't create it, to avoid typos).
if [ "$STATS_ONLY" -eq 0 ]; then
    if [ ! -d "$IMAGES_DIR" ]; then
        echo "ERROR: input directory not found: $IMAGES_DIR" >&2
        exit 1
    fi
    IMAGES_DIR="$(cd "$IMAGES_DIR" && pwd)"
    mkdir -p "$OUTPUT_DIR"
else
    if [ ! -d "$OUTPUT_DIR" ]; then
        echo "ERROR: directory not found: $OUTPUT_DIR" >&2
        exit 1
    fi
    # In --stats, if an images directory is provided (config or --images-dir),
    # make it absolute so the header reading (time_windows.csv) can find it;
    # otherwise continue (the user may pass --frame-time-ms instead).
    if [ -d "$IMAGES_DIR" ]; then
        IMAGES_DIR="$(cd "$IMAGES_DIR" && pwd)"
    fi
fi
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- Detect mode from the input file extension (processing only) ---
if [ "$STATS_ONLY" -eq 0 ]; then
    case "$IMAGE_TEMPLATE" in
        *.h5|*.H5|*.nxs|*.NXS) INPUT_MODE="HDF5" ;;
        *.cbf|*.CBF)           INPUT_MODE="CBF" ;;
        *)
            echo "ERROR: unrecognized input type for '$IMAGE_TEMPLATE' (expected .cbf or .h5)" >&2
            exit 1
            ;;
    esac
    echo "Input mode: $INPUT_MODE ($IMAGE_TEMPLATE)"
    module load autoPROC
fi

# --- All processing runs go into a sub-folder of the output directory ---
# In --stats mode, the argument may point either at a folder CONTAINING the
# autoPROC_<first>_<last> chunks (whatever its name: autoproc_chunks, autoproc_copy,
# a timepoint folder...), or at its parent folder. Detection is by content,
# not by name, so all directory layouts are accepted.
PROCESS_DIR="${OUTPUT_DIR}/autoproc_chunks"

if [ "$STATS_ONLY" -eq 1 ]; then
    # ---------------- Reports-only mode: no autoPROC ----------------
    # Two directory layouts are accepted (detection by content, not by name):
    #   flat   : <dir>/autoPROC_1_150           (TR-FRX_autoPROC.sh)
    #   nested : <dir>/1_150/autoPROC_1_150      (copy from trfrx_full_pipeline)
    if ls -d "$OUTPUT_DIR"/autoPROC_*_* >/dev/null 2>&1; then
        PROCESS_DIR="$OUTPUT_DIR"                        # chunks directly here
    elif [ -d "$PROCESS_DIR" ]; then
        :                                                # <dir>/autoproc_chunks
    elif ls -d "$OUTPUT_DIR"/*/autoPROC_*_* >/dev/null 2>&1; then
        PROCESS_DIR="$OUTPUT_DIR"                        # nested layout
    fi
    if [ ! -d "$PROCESS_DIR" ]; then
        echo "ERROR: chunk directory not found: $PROCESS_DIR" >&2
        echo "         (run the script from the output directory, or pass it" >&2
        echo "          as an argument: $0 --stats /path/to/output)" >&2
        exit 1
    fi
    # Available chunks (autoPROC_<first>_<last>), sorted by first image.
    # maxdepth 2: covers both the flat layout AND the nested layout
    # (<dir>/1_150/autoPROC_1_150) produced by trfrx_full_pipeline.
    # -L: follows symlinks (chunk folders are sometimes linked).
    CHUNK_DIRS=$(find -L "$PROCESS_DIR" -maxdepth 2 -type d -name 'autoPROC_*_*' \
                 -exec basename {} \; 2>/dev/null \
                 | sed -E 's/^autoPROC_([0-9]+)_([0-9]+)$/\1 \2/' \
                 | grep -E '^[0-9]+ [0-9]+$' | sort -n -u)
    if [ -z "$CHUNK_DIRS" ]; then
        echo "ERROR: no autoPROC_<first>_<last> chunk in $PROCESS_DIR" >&2
        exit 1
    fi
    N_CHUNKS=$(printf '%s\n' "$CHUNK_DIRS" | wc -l | tr -d ' ')
    FIRST_AVAIL=$(printf '%s\n' "$CHUNK_DIRS" | awk '{print $1}' | sort -n | head -1)
    LAST_AVAIL=$(printf '%s\n' "$CHUNK_DIRS" | awk '{print $2}' | sort -n | tail -1)
    echo ""
    echo "=========== Reports only (no autoPROC reprocessing) ==========="
    echo "  Chunk directory  : $PROCESS_DIR"
    echo "  Chunks available : $N_CHUNKS   (images $FIRST_AVAIL -> $LAST_AVAIL)"
    printf '%s\n' "$CHUNK_DIRS" | awk '{printf "      images %s-%s\n", $1, $2}'
    echo "==============================================================="
    # Interactive limit prompt (Enter = all). A chunk that exceeds the limit
    # is EXCLUDED (never truncated): limit 2000 -> ...1801-2000 kept,
    # 2001-2200 dropped.
    if [ -z "$MAX_IMAGE" ]; then
        while true; do
            printf "Maximum image for statistics [Enter = all (%s)]: " "$LAST_AVAIL"
            if ! read -r ANSWER; then ANSWER=""; echo; fi
            ANSWER="$(printf '%s' "$ANSWER" | tr -d '[:space:]')"
            if [ -z "$ANSWER" ]; then
                MAX_IMAGE=""            # no limit
                break
            fi
            if printf '%s' "$ANSWER" | grep -Eq '^[0-9]+$'; then
                KEPT=$(printf '%s\n' "$CHUNK_DIRS" \
                       | awk -v lim="$ANSWER" '$2 <= lim' | wc -l | tr -d ' ')
                if [ "$KEPT" -eq 0 ]; then
                    FIRST_END=$(printf '%s\n' "$CHUNK_DIRS" | head -1 | awk '{print $2}')
                    echo "  $ANSWER excludes all chunks (the first one ends at $FIRST_END). Try again."
                    continue
                fi
                LAST_KEPT=$(printf '%s\n' "$CHUNK_DIRS" \
                            | awk -v lim="$ANSWER" '$2 <= lim' | tail -1)
                MAX_IMAGE="$ANSWER"
                echo "  -> $KEPT chunk(s) kept; last = images ${LAST_KEPT% *}-${LAST_KEPT#* }"
                break
            fi
            echo "  '$ANSWER' is not an integer — type a number (e.g. 2000) or Enter."
        done
    fi
    if [ -n "$MAX_IMAGE" ]; then
        echo "  Limit applied: images <= $MAX_IMAGE  (files suffixed _upto${MAX_IMAGE})"
    else
        echo "  No limit: all chunks are included."
    fi
    echo ""
else
    mkdir -p "$PROCESS_DIR"
fi

# --- Processing job IDs (for the final report) ---
CHUNK_JOBIDS=()
CHUNK_LABELS=()

FIRST_OUTDIR=""
FIRST_JOBID=""

if [ "$STATS_ONLY" -eq 0 ]; then

# --- Process the reference dataset ---
REF_JOB_BASENAME="${REF_FIRST_IMG}-${REF_LAST_IMG}_autoproc"
REF_JOB_SCRIPT="${PROCESS_DIR}/${REF_JOB_BASENAME}.sh"
REF_OUTDIR="${PROCESS_DIR}/autoPROC_${REF_FIRST_IMG}_${REF_LAST_IMG}"

cat > "$REF_JOB_SCRIPT" <<EOF
#!/bin/bash
set -euo pipefail
module load autoPROC
cd "$PROCESS_DIR"
process -Id "$RUN_ID,$IMAGES_DIR,$IMAGE_TEMPLATE,$REF_FIRST_IMG,$REF_LAST_IMG" symm="$SYMM" cell="$CELL" AutoProcScale_RunStaraniso=yes -d "$REF_OUTDIR"
EOF

sbatch_output=$(sbatch -p "$SLURM_PARTITION" -n 1 -c "$SLURM_CPUS" --mem="$SLURM_MEM" "$REF_JOB_SCRIPT")
FIRST_JOBID=$(echo "$sbatch_output" | awk '{print $4}')
FIRST_OUTDIR="$REF_OUTDIR"
CHUNK_JOBIDS+=("$FIRST_JOBID")
CHUNK_LABELS+=("${REF_FIRST_IMG}-${REF_LAST_IMG} (reference)")
echo "  [submitted] job $FIRST_JOBID: images ${REF_FIRST_IMG}-${REF_LAST_IMG} (reference)"

# --- Process the following datasets in chunks ---
i=$((REF_LAST_IMG + 1))

while [ "$i" -le "$LAST_IMG" ]; do
    j=$((i + CHUNK_SIZE - 1))
    if [ "$j" -gt "$LAST_IMG" ]; then
        j="$LAST_IMG"
    fi

    JOB_BASENAME="${i}-${j}_autoproc"
    JOB_SCRIPT="${PROCESS_DIR}/${JOB_BASENAME}.sh"
    OUTDIR="${PROCESS_DIR}/autoPROC_${i}_${j}"

    cat > "$JOB_SCRIPT" <<EOF
#!/bin/bash
set -euo pipefail
module load autoPROC
cd "$PROCESS_DIR"
process -Id "$RUN_ID,$IMAGES_DIR,$IMAGE_TEMPLATE,$i,$j" -ref "$FIRST_OUTDIR/truncate-unique.mtz" AutoProcScale_RunStaraniso=yes -d "$OUTDIR"
EOF

    chunk_output=$(sbatch -p "$SLURM_PARTITION" -n 1 -c "$SLURM_CPUS" --mem="$SLURM_MEM" --dependency=afterok:"$FIRST_JOBID" "$JOB_SCRIPT")
    CHUNK_JOBID=$(echo "$chunk_output" | awk '{print $4}')
    CHUNK_JOBIDS+=("$CHUNK_JOBID")
    CHUNK_LABELS+=("${i}-${j}")
    echo "  [submitted] job $CHUNK_JOBID: images ${i}-${j} (depends on $FIRST_JOBID)"

    i=$((j + 1))
done

fi   # end of the autoPROC submission block (skipped in --stats mode)

# =====================================================================
# --- Report generator (written to PROCESS_DIR, then launched) ---------
# =====================================================================
cat > "$PROCESS_DIR/trfrx_damage_report.py" <<'PYEOF'
#!/usr/bin/env python3
"""TR-FRX radiation-damage report from autoPROC chunk statistics.

Standalone tool (development version — not yet wired into TR-FRX_autoPROC.sh).
Run it after processing, pointed at an existing ``autoproc_chunks`` directory:

    python trfrx_damage_report.py --process-dir /path/to/autoproc_chunks
    python trfrx_damage_report.py --process-dir ... --dataset SAMPLE

It needs only libraries already provided by setup_env.sh (matplotlib, gemmi,
numpy for the plots / Wilson-B; the CSV and scraped stats work with the
Python standard library alone). Every optional dependency is imported
defensively, so a missing one degrades to a warning rather than a crash.

For a given ``autoproc_chunks`` directory (the PROCESS_DIR created by
TR-FRX_autoPROC.sh), this script scans every ``autoPROC_<first>_<last>`` chunk
sub-directory, extracts the final merging statistics for both processing routes
(STARANISO and classical TRUNCATE) and writes, into a ``reports/`` sub-folder:

    reports/staraniso_statistics.csv   reports/staraniso_report.pdf
    reports/truncate_statistics.csv    reports/truncate_report.pdf
    reports/parsing_diagnostics.txt

Because TR-FRX processes successive image ranges ("chunks") that accumulate
dose, each chunk is a dose point. The reports therefore plot a set of
*radiation-damage* indicators as a function of the mean image number (a dose
proxy):

    * High-resolution limit vs dose               (rising = damage)
    * R_merge / R_meas / R_pim vs dose            (rising = damage)
    * Mean I/sigma(I) vs dose                     (falling = damage)
    * ISa, asymptotic I/sigma vs dose             (falling = damage;
      Diederichs 2010, Acta Cryst. D66, 733-740)
    * CC(1/2) vs dose                             (falling = damage)
    * CC(anom), anomalous half-set CC vs dose     (falling = anomalous
      signal decay, an early/sensitive damage reporter)
    * Completeness (spherical) vs dose            (falling = damage/rejections)
    * Multiplicity vs dose
    * Wilson B-factor vs dose                     (rising = loss of
      high-resolution order; the global B-factor "scaling"/decay curve)
    * Unit-cell edges (true A lengths) vs dose    (expansion = damage)
    * Images actually used per chunk vs dose

Scalar statistics come from the AIMLESS "Summary data" tables and the per-route
``*.table1`` files (CC(ano), completeness, multiplicity included). ISa is scraped
from the XDS ``CORRECT.LP`` error model. Unit cells come from the same files. The
Wilson B-factor is scraped when autoPROC prints it and otherwise estimated from
the merged MTZ via ``gemmi`` (relative B: the offset from omitting the atomic
scattering term cancels in dB across chunks).

Any statistic that cannot be found is reported as ``N/A`` rather than causing a
failure; parsing_diagnostics.txt lists every file scanned and value picked so
the parser can be tightened for non-standard filenames.
"""

import argparse
import csv
import datetime
import math
import os
import re
import sys


# --------------------------------------------------------------------------
# Statistics extracted from the AIMLESS/table1 "Overall InnerShell OuterShell"
# rows. Each entry: (csv_column, [keywords that must appear in the label],
#                    [keywords that must NOT appear], want_outer_shell)
# --------------------------------------------------------------------------
METRICS = [
    ("resolution_low",        ["low resolution"],              [],            False),
    ("resolution_high",       ["high resolution"],             [],            True),
    ("Rmerge",                ["rmerge", "all"],               [],            True),
    ("Rmeas",                 ["rmeas", "all"],                [],            True),
    ("Rpim",                  ["rpim", "all"],                 [],            True),
    ("Mean_I_over_sigma",     ["mean", "sd(i)"],               [],            True),
    ("CC_half",               ["cc(1/2)"],                     [],            True),
    ("CC_ano",                ["cc(ano)"],                     [],            True),
    ("Completeness",              ["completeness"],                           ["anomalous", "ellipsoidal"], True),
    ("Completeness_ellipsoidal",  ["completeness", "ellipsoidal"],            ["anomalous"],                True),
    ("Multiplicity",              ["multiplicity"],                           ["anomalous"],                True),
    ("Anom_completeness",             ["anomalous", "completeness"],              ["ellipsoidal"],          True),
    ("Anom_completeness_ellipsoidal", ["anomalous", "completeness", "ellipsoidal"], [],                    True),
    ("Anom_multiplicity",             ["anomalous", "multiplicity"],              [],                       True),
    ("N_observations",        ["total", "observations"],       [],            False),
    ("N_unique",              ["total", "unique"],             [],            False),
]

NUM = r"[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?"
ROW_RE = re.compile(r"^\s*(.*?)\s+(%s)\s+(%s)\s+(%s)\s*$" % (NUM, NUM, NUM))
CHUNK_RE = re.compile(r"^autoPROC_(\d+)_(\d+)$")

# Wilson B-factor — several autoPROC/CTRUNCATE/STARANISO phrasings.
WILSON_RES = [
    re.compile(r"wilson\s*b[- ]?factor[^0-9\-+]*(%s)" % NUM, re.I),
    re.compile(r"\bb[- ]?wilson\b[^0-9\-+]*(%s)" % NUM, re.I),
    re.compile(r"wilson\s*plot[^0-9]*?b[- ]?factor[^0-9\-+]*(%s)" % NUM, re.I),
    re.compile(r"estimated\s+b[- ]?factor[^0-9\-+]*(%s)" % NUM, re.I),
]
# Six cell parameters after a "unit cell" / "average unit cell" marker.
CELL6 = r"(%s)[ \t]+(%s)[ \t]+(%s)[ \t]+(%s)[ \t]+(%s)[ \t]+(%s)" % ((NUM,) * 6)
CELL_AVG_RE = re.compile(r"average\s+unit\s+cell[^\n]*?" + CELL6, re.I)
CELL_ANY_RE = re.compile(r"unit[- ]cell[^\n]*?" + CELL6, re.I)

# Real number of images that survived into scaling (autoPROC/XDS reject some).
# AIMLESS batch count is the primary source (1 batch ~ 1 image); XDS logs back
# it up. Every match is recorded in the diagnostics so the source is auditable.
IMAGES_USED_RES = [
    ("aimless-batches", re.compile(r"number\s+of\s+batches\s*[:=]?\s*(\d+)", re.I)),
    ("images-used",     re.compile(r"number\s+of\s+images\s+used\s*[:=]?\s*(\d+)", re.I)),
    ("images-accepted", re.compile(r"images\s+(?:used|accepted)\s*[:=]?\s*(\d+)", re.I)),
    ("xds-of-images",   re.compile(r"of\s+(\d+)\s+images", re.I)),
]

# ISa — asymptotic I/sigma(I) = 1/sqrt(a*b) from the XDS CORRECT.LP error model
# (Diederichs, 2010, Acta Cryst. D66, 733-740). CORRECT.LP prints a header row
# "a  b  ISa" followed by the three values; ISa is the third. Higher = better;
# a fall across the series indicates radiation damage / degrading measurability.
ISA_RE = re.compile(r"\ba\s+b\s+ISa\b[^\n]*\n\s*[-+.\deE]+\s+[-+.\deE]+\s+([\d.]+)",
                    re.I)


def to_float(value):
    """Best-effort conversion of a scraped statistic to float, else None."""
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    text = str(value).strip().replace(",", "")
    if text == "" or text.upper() == "N/A" or set(text) <= set("-"):
        return None
    try:
        return float(text)
    except ValueError:
        m = re.search(NUM, text)
        return float(m.group()) if m else None


def cell_volume(cell):
    """Triclinic unit-cell volume from (a, b, c, alpha, beta, gamma) in deg."""
    if not cell or any(v is None for v in cell):
        return None
    a, b, c, al, be, ga = cell
    ca, cb, cg = (math.cos(math.radians(x)) for x in (al, be, ga))
    fac = 1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg
    if fac <= 0:
        return None
    return a * b * c * math.sqrt(fac)


def find_summary_blocks(text):
    """Return [{'title', 'metrics'}] for every AIMLESS Summary data table in
    *text*. metrics maps label -> (overall, inner, outer)."""
    lines = text.splitlines()
    n = len(lines)
    blocks = []
    i = 0
    while i < n:
        if "Summary data for" in lines[i]:
            title = lines[i].strip()
            header = None
            for j in range(i + 1, min(i + 12, n)):
                if "Overall" in lines[j] and ("Shell" in lines[j] or "Low" in lines[j]):
                    header = j
                    break
            if header is not None:
                metrics = {}
                k = header + 1
                while k < n:
                    line = lines[k]
                    stripped = line.strip()
                    if (stripped.startswith("$$")
                            or "SUMMARY_END" in stripped
                            or stripped.startswith("====")
                            or stripped.startswith("Estimates of resolution")
                            or "Summary data for" in line):
                        break
                    m = ROW_RE.match(line)
                    if m:
                        label = re.sub(r"\s+", " ", m.group(1).strip())
                        if label:
                            metrics[label] = (m.group(2), m.group(3), m.group(4))
                    k += 1
                blocks.append({"title": title, "metrics": metrics})
                i = k
                continue
        i += 1
    return blocks


def is_staraniso(path, title):
    return "staraniso" in path.lower() or "staraniso" in title.lower()


def parse_table1(text):
    """Parse a ``*.table1`` (label Overall InnerShell OuterShell) file."""
    lines = text.splitlines()
    header = None
    for idx, line in enumerate(lines):
        if "Overall" in line and "Shell" in line:
            header = idx
            break
    if header is None:
        return {}
    metrics = {}
    for line in lines[header + 1:]:
        stripped = line.strip()
        if not stripped or set(stripped) <= set("-"):
            continue
        m = ROW_RE.match(line)
        if m:
            label = re.sub(r"\s+", " ", m.group(1).strip())
            if label:
                metrics[label] = (m.group(2), m.group(3), m.group(4))
    return metrics


def scrape_wilson_b(text):
    """Return the last Wilson B-factor mentioned in *text*, or None."""
    best = None
    for rex in WILSON_RES:
        found = rex.findall(text)
        if found:
            best = to_float(found[-1])
    return best


def scrape_staraniso_wilson_b(text):
    """STARANISO Popov-Bourenkov Wilson B from staraniso_alldata.log.

    STARANISO prints a 'Wilson plots using as expected intensity' section whose
    columns are 'intercept, gradient, WilsonB, Goodness-of-Fit'; the 'PB' row
    (Popov & Bourenkov 2003 expected intensities) gives WilsonB (3rd column).
    The last such block in the file is the converged one.
    """
    heads = list(re.finditer(r"intercept\s*,\s*gradient\s*,\s*WilsonB", text, re.I))
    if not heads:
        return None
    tail = text[heads[-1].end():]
    m = re.search(r"^\s*PB\s+(%s)\s+(%s)\s+(%s)\s+(%s)" % (NUM, NUM, NUM, NUM),
                  tail, re.M)
    return to_float(m.group(3)) if m else None


def scrape_cell(text):
    """Return (a, b, c, al, be, ga) preferring an 'Average unit cell' line."""
    m = CELL_AVG_RE.search(text) or CELL_ANY_RE.search(text)
    if not m:
        return None
    cell = tuple(to_float(g) for g in m.groups())
    if any(v is None or v <= 0 for v in cell):
        return None
    return cell


def scrape_isa(text):
    """Return the last ISa (asymptotic I/sigma) value in *text*, or None.
    Only XDS CORRECT.LP carries it, so this is a no-op on other files."""
    found = ISA_RE.findall(text)
    return to_float(found[-1]) if found else None


def scrape_images_used(text):
    """Return (n_images_used, source_tag) or (None, None)."""
    for tag, rex in IMAGES_USED_RES:
        m = rex.search(text)
        if m:
            return int(m.group(1)), tag
    return None, None


def source_priority(fname):
    """Higher = more authoritative source for the merging statistics."""
    f = fname.lower()
    if f in ("truncate-unique.table1", "staraniso_alldata-unique.table1"):
        return 100
    if f.endswith("-unique.table1"):
        return 60
    if f.endswith(".table1"):
        return 40
    if f == "aimless.log":
        return 20
    if f == "aimless_alldata.log":
        return 15
    if f.startswith("aimless"):
        return 10
    return 5


def collect_chunk_blocks(chunk_dir, diagnostics):
    """Scan a chunk directory and return the best (truncate, staraniso) data:
    merging metrics, unit cell, merged-MTZ path and images-used count."""
    best = {"truncate": None, "staraniso": None}
    best_rank = {"truncate": (-1, -1), "staraniso": (-1, -1)}
    extra = {r: {"cell": None, "mtz": None, "wilson_b": None, "images_used": None,
                 "wb_source": None, "isa": None,
                 "cell_rank": -1, "mtz_rank": -1, "wb_rank": -1, "img_rank": -1,
                 "isa_rank": -1}
             for r in ("truncate", "staraniso")}

    def consider(route, metrics, fname, rel):
        if not metrics:
            return
        diagnostics.append("    [%s] %s (%d metrics) <- %s"
                           % (route, os.path.basename(chunk_dir), len(metrics), rel))
        rank = (source_priority(fname), len(metrics))
        if rank > best_rank[route]:
            best_rank[route] = rank
            best[route] = {"metrics": metrics, "source": rel}

    def consider_extra(route, fname, rel, cell, wilson, images, isa=None):
        prio = source_priority(fname)
        e = extra[route]
        if cell is not None and prio > e["cell_rank"]:
            e["cell"], e["cell_rank"] = cell, prio
        # ISa is route-agnostic (one XDS integration feeds both routes).
        if isa is not None and prio > e["isa_rank"]:
            e["isa"], e["isa_rank"] = isa, prio
            diagnostics.append("    [%s] ISa = %s <- %s" % (route, isa, rel))
        # Wilson B: penalise the _early/_late radiation-damage half-sets so the
        # full-data value (e.g. truncate.log) wins over truncate_late.log.
        wb_prio = prio - 3 if ("_late" in fname.lower()
                               or "_early" in fname.lower()) else prio
        if wilson is not None and wb_prio > e["wb_rank"]:
            e["wilson_b"], e["wb_rank"] = wilson, wb_prio
            method = "Popov-Bourenkov" if route == "staraniso" else "CTRUNCATE"
            e["wb_source"] = "%s (%s)" % (os.path.basename(fname), method)
            diagnostics.append("    [%s] Wilson B = %s <- %s (%s)"
                               % (route, wilson, rel, method))
        if images is not None:
            n, tag = images
            # Prefer AIMLESS batch count; keep the first plausible hit.
            weight = {"aimless-batches": 3, "images-used": 2,
                      "images-accepted": 2, "xds-of-images": 1}.get(tag, 0)
            if weight > e["img_rank"]:
                e["images_used"], e["img_rank"] = n, weight
                diagnostics.append("    [%s] images used = %d <- %s (%s)"
                                   % (route, n, rel, tag))

    def consider_mtz(route, fname, fpath):
        prio = 100 if fname.lower().endswith("-unique.mtz") else 50
        e = extra[route]
        if prio > e["mtz_rank"]:
            e["mtz"], e["mtz_rank"] = fpath, prio

    def handle_text(fname, rel, text, path_hint=""):
        """Scrape one text file (loose on disk or read from summary.tar.gz)."""
        low = fname.lower()
        route_of_file = "staraniso" if "staraniso" in low else "truncate"
        cell = scrape_cell(text)
        wilson = (scrape_staraniso_wilson_b(text) if route_of_file == "staraniso"
                  else scrape_wilson_b(text))
        images = scrape_images_used(text)
        images = images if images[0] is not None else None
        isa = scrape_isa(text)
        # Unit cell / images used / ISa are route-agnostic — record for both so the
        # STARANISO route is not left blank when only a truncate log has them.
        for r in ("truncate", "staraniso"):
            consider_extra(r, fname, rel, cell,
                           wilson if r == route_of_file else None, images, isa)
        if low.endswith(".table1"):
            consider(route_of_file, parse_table1(text), fname, rel)
        elif low.endswith(".log") or low.endswith(".lp"):
            for block in find_summary_blocks(text):
                route = ("staraniso" if is_staraniso(path_hint or fname,
                                                     block["title"])
                         else "truncate")
                consider(route, block["metrics"], fname, rel)

    tarballs = []
    for root, _dirs, files in os.walk(chunk_dir):
        for fname in files:
            low = fname.lower()
            fpath = os.path.join(root, fname)
            rel = os.path.relpath(fpath, chunk_dir)
            if low.endswith(".mtz"):
                route = "staraniso" if "staraniso" in low else "truncate"
                consider_mtz(route, fname, fpath)
                continue
            if low.endswith(".tar.gz") or low.endswith(".tgz"):
                tarballs.append(fpath)
                continue
            if not (low.endswith(".table1") or low.endswith(".log")
                    or low.endswith(".lp")):
                continue
            try:
                with open(fpath, "r", errors="replace") as fh:
                    text = fh.read(2_000_000)     # cap: XDS .LP can be huge
            except (IOError, OSError):
                continue
            handle_text(fname, rel, text, fpath)

    # Fallback for pipeline-copied trees (trfrx_full_pipeline copies only a few
    # files per chunk): the loose *.table1 / aimless logs are absent, but
    # summary.tar.gz still carries them. Only opened when a route has no metrics
    # yet, and only its small text members are read.
    _need_tar = (best["truncate"] is None or best["staraniso"] is None
                 or extra["truncate"]["isa"] is None)
    if tarballs and _need_tar:
        import tarfile
        for tpath in tarballs:
            try:
                tf = tarfile.open(tpath, "r:*")
            except (tarfile.TarError, IOError, OSError) as exc:
                diagnostics.append("    could not open %s (%s)"
                                   % (os.path.basename(tpath), exc))
                continue
            try:
                for member in tf:
                    if not member.isfile():
                        continue
                    mlow = member.name.lower()
                    if not (mlow.endswith(".table1") or mlow.endswith(".log")
                            or mlow.endswith(".lp")):     # .lp -> CORRECT.LP (ISa)
                        continue
                    if member.size > 2_000_000:           # skip huge INTEGRATE.LP
                        continue
                    fh = tf.extractfile(member)
                    if fh is None:
                        continue
                    try:
                        text = fh.read().decode("utf-8", "replace")
                    finally:
                        fh.close()
                    base = os.path.basename(member.name)
                    # No early exit: every candidate must be seen so the usual
                    # source_priority ranking picks the SAME file the flat layout
                    # would use (truncate-unique.table1 beats aimless_alldata-*).
                    handle_text(base, "%s:%s" % (os.path.basename(tpath), base),
                                text, member.name)
            except (tarfile.TarError, IOError, OSError, EOFError) as exc:
                diagnostics.append("    error reading %s (%s)"
                                   % (os.path.basename(tpath), exc))
            finally:
                tf.close()

    result = {}
    for route in ("truncate", "staraniso"):
        data = best[route] or {"metrics": {}, "source": "N/A"}
        e = extra[route]
        data.update({"cell": e["cell"], "mtz": e["mtz"],
                     "wilson_b": e["wilson_b"], "wb_source": e["wb_source"],
                     "images_used": e["images_used"], "isa": e["isa"]})
        has_any = (best[route] or e["cell"] or e["mtz"]
                   or e["wilson_b"] is not None or e["images_used"] is not None
                   or e["isa"] is not None)
        result[route] = data if has_any else None
    return result


def match_metric(metrics, must, mustnot):
    for label, values in metrics.items():
        low = label.lower()
        if all(x in low for x in must) and not any(x in low for x in mustnot):
            return values
    return None


def build_row(chunk_name, first, last, block, wilson_b, images_used,
              wilson_source="N/A", t_start=None, t_end=None, t_mid=None):
    """Return an ordered list of (column, value) for one chunk."""
    def _ts(v):
        return "%.3f" % v if v is not None else "N/A"
    row = [
        ("chunk", chunk_name),
        ("image_first", first),
        ("image_last", last),
        ("n_images", last - first + 1),
        ("images_used", images_used if images_used is not None else "N/A"),
        ("t_start_s", _ts(t_start)),
        ("t_end_s", _ts(t_end)),
        ("t_mid_s", _ts(t_mid)),
    ]
    metrics = block["metrics"] if block else {}
    for col, must, mustnot, want_outer in METRICS:
        values = match_metric(metrics, must, mustnot)
        row.append((col, values[0] if values else "N/A"))
        if want_outer:
            row.append((col + "_outer", values[2] if values else "N/A"))
    cell = block.get("cell") if block else None
    vol = cell_volume(cell) if cell else None
    isa = block.get("isa") if block else None
    row.append(("Isa", "%.2f" % isa if isa is not None else "N/A"))
    row.append(("Wilson_B", "%.2f" % wilson_b if wilson_b is not None else "N/A"))
    row.append(("Wilson_B_source", wilson_source))
    names = ["cell_a", "cell_b", "cell_c", "cell_al", "cell_be", "cell_ga"]
    for idx, nm in enumerate(names):
        row.append((nm, "%.3f" % cell[idx] if cell else "N/A"))
    row.append(("cell_volume", "%.1f" % vol if vol is not None else "N/A"))
    row.append(("source_file", block["source"] if block else "N/A"))
    return row


def write_csv(path, rows):
    if not rows:
        headers = ["chunk", "image_first", "image_last", "n_images"]
        rows_out = []
    else:
        headers = [c for c, _ in rows[0]]
        rows_out = [[v for _, v in r] for r in rows]
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(headers)
        writer.writerows(rows_out)


# --------------------------------------------------------------------------
# Real per-dataset time: read the frame period + oscillation from ONE image
# header (they are constant across a rotation series), then map every chunk's
# image window to a time window. CBF (miniCBF ASCII) and Eiger HDF5 masters are
# both supported; every read is defensive and degrades to a warning, never a
# crash. Explicit --frame-time-ms / --osc-per-image overrides always win.
# --------------------------------------------------------------------------
def _resolve_image_path(images_dir, template):
    """Best-effort path to a single image whose header we can read.

    Eiger: the template IS the master .h5 (one file). CBF: the template carries
    a #### placeholder — glob the images dir for the first matching frame."""
    if not template:
        return None
    base = os.path.join(images_dir, template) if images_dir else template
    if os.path.isfile(base):
        return base
    if "#" in template:
        import glob
        pat = os.path.join(images_dir or ".", re.sub(r"#+", "*", template))
        hits = sorted(glob.glob(pat))
        if hits:
            return hits[0]
    return base   # may not exist; caller handles a missing file


def _read_eiger_header(path):
    """(dt_s, osc_deg, meta) from an Eiger HDF5 master. Multi-key: Dectris NeXus
    layouts vary (frame_time vs count_time; omega_increment vs omega_range_average
    vs a per-frame omega array)."""
    meta = {"reader": "eiger"}
    try:
        import h5py
    except Exception as e:                       # h5py not installed in this env
        meta["note"] = "h5py unavailable (%s)" % e
        return None, None, meta
    try:
        with h5py.File(path, "r") as f:
            def g(p):
                try:
                    d = f[p]
                    return d[()] if getattr(d, "shape", None) == () else d[...]
                except Exception:
                    return None
            dt = None
            for k in ("/entry/instrument/detector/frame_time",
                      "/entry/instrument/detector/count_time",
                      "/entry/instrument/detector/detectorSpecific/frame_time"):
                v = g(k)
                if v is not None:
                    dt = float(v); meta["dt_key"] = k; break
            osc = None
            for k in ("/entry/sample/goniometer/omega_increment",
                      "/entry/sample/goniometer/omega_range_average",
                      "/entry/instrument/detector/goniometer/omega_increment",
                      "/entry/instrument/detector/goniometer/omega_range_average"):
                v = g(k)
                if v is not None:
                    osc = abs(float(v)); meta["osc_key"] = k; break
            if osc is None:                       # last resort: diff the omega array
                arr = g("/entry/sample/goniometer/omega")
                try:
                    if arr is not None and len(arr) >= 2:
                        osc = abs(float(arr[1]) - float(arr[0]))
                        meta["osc_key"] = "omega[]diff"
                except TypeError:
                    pass
            nimg = g("/entry/instrument/detector/detectorSpecific/nimages")
            if nimg is not None:
                meta["nimages"] = int(nimg)
            wl = g("/entry/instrument/beam/incident_wavelength")
            if wl is not None:
                meta["wavelength"] = float(wl)
        return dt, osc, meta
    except Exception as e:
        meta["note"] = "read failed (%s)" % e
        return None, None, meta


def _read_cbf_header(path):
    """(dt_s, osc_deg, meta) from a PILATUS/Eiger miniCBF ASCII header. Handles
    both "# Key value unit" (Dectris) and "Key = value" forms."""
    meta = {"reader": "cbf"}
    try:
        with open(path, "rb") as fh:
            text = fh.read(8192).decode("latin-1", "replace")
    except Exception as e:
        meta["note"] = "open failed (%s)" % e
        return None, None, meta

    def find(key):
        m = re.search(r"(?mi)^\s*#?\s*%s\s*[=:\s]\s*(%s)" % (re.escape(key), NUM),
                      text)
        return float(m.group(1)) if m else None

    dt = find("Exposure_period")                 # frame-to-frame time (incl. dead time)
    osc = find("Angle_increment")
    et = find("Exposure_time")
    sa = find("Start_angle")
    if et is not None:
        meta["exposure_time_s"] = et
    if sa is not None:
        meta["start_angle_deg"] = sa
    return dt, (abs(osc) if osc is not None else None), meta


def frame_period_from_header(images_dir, template, frame_time_ms=None,
                             osc_per_image=None):
    """(dt_s, osc_deg, meta): per-image frame period (s) and oscillation (deg).
    Overrides win; otherwise auto-detect CBF vs Eiger by extension. dt_s is None
    only when neither a header nor an override is available (source=UNRESOLVED)."""
    meta = {"source": ""}
    dt_s = (frame_time_ms / 1000.0) if frame_time_ms else None
    osc = osc_per_image if osc_per_image else None
    if dt_s is not None and osc is not None:
        meta["source"] = "override"
        return dt_s, osc, meta
    path = _resolve_image_path(images_dir, template)
    meta["image"] = path
    hdt = hosc = None
    if path:
        ext = os.path.splitext(path)[1].lower()
        if ext in (".h5", ".hdf5", ".nxs"):
            hdt, hosc, hmeta = _read_eiger_header(path)
        elif ext == ".cbf":
            hdt, hosc, hmeta = _read_cbf_header(path)
        else:
            hmeta = {"reader": "unknown-ext"}
        meta.update(hmeta)
        if not (path and os.path.isfile(path)):
            meta.setdefault("note", "image not found: %s" % path)
    else:
        meta["note"] = "no image template/dir given"
    if dt_s is None:
        dt_s = hdt
    if osc is None:
        osc = hosc
    if dt_s is not None:
        meta["source"] = meta.get("reader", "header")
    else:
        meta["source"] = "UNRESOLVED"
    return dt_s, osc, meta


def time_window(image_first, image_last, dt_s, t0_image=1, t0_s=0.0):
    """(t_start, t_end, t_mid, duration) in seconds for an image window, using the
    linear model t(i) = t0_s + (i - t0_image) * dt_s. All None when dt_s is None."""
    if dt_s is None:
        return None, None, None, None
    t_start = t0_s + (image_first - t0_image) * dt_s
    t_end = t0_s + (image_last - t0_image) * dt_s
    t_mid = t0_s + (0.5 * (image_first + image_last) - t0_image) * dt_s
    duration = (image_last - image_first + 1) * dt_s
    return t_start, t_end, t_mid, duration


def write_time_windows_csv(path, chunk_list, dt_s, osc_deg, t0_image, t0_s,
                           source, images_used_by=None):
    """Canonical per-dataset time table. One row per chunk (dataset). Written on
    every run; blank times + source=UNRESOLVED when no header/override was found."""
    images_used_by = images_used_by or {}
    cols = ["dataset_id", "image_first", "image_last", "n_images", "images_used",
            "frame_time_s", "osc_deg", "t_start_s", "t_end_s", "t_mid_s",
            "t0_image", "t0_s", "source"]

    def f(v, nd=4):
        return "" if v is None else ("%.*f" % (nd, v))

    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for name, first, last in chunk_list:
            ts, te, tm, _dur = time_window(first, last, dt_s, t0_image, t0_s)
            iu = images_used_by.get(name)
            w.writerow([name.replace("autoPROC_", ""), first, last,
                        last - first + 1, ("" if iu is None else iu),
                        f(dt_s, 6), f(osc_deg, 4), f(ts), f(te), f(tm),
                        t0_image, f(t0_s), source])


# --------------------------------------------------------------------------
# gemmi-based quantity: Wilson B (from merged intensities).
# --------------------------------------------------------------------------
def _read_mtz_intensities(mtz_path):
    """Return (miller_dict {(h,k,l): I}, inv_d2_array, I_array) or (None,..)."""
    try:
        import gemmi
        import numpy as np
    except ImportError:
        return None, None, None
    try:
        mtz = gemmi.read_mtz_file(mtz_path)
    except (RuntimeError, IOError, OSError, ValueError):
        return None, None, None
    icol = next((c for c in mtz.columns if c.type == "J"), None)
    if icol is None:
        return None, None, None
    try:
        hkl = mtz.make_miller_array()
        inv_d2 = mtz.make_1_d2_array()
        ivals = np.asarray(icol.array, dtype=float)
    except (RuntimeError, ValueError):
        return None, None, None
    good = np.isfinite(ivals)
    data = {}
    for (h, k, l), inten in zip(hkl[good], ivals[good]):
        data[(int(h), int(k), int(l))] = float(inten)
    return data, inv_d2[good], ivals[good]


def wilson_b_from_arrays(inv_d2, ivals, d_cut=4.0, nbins=20):
    """Relative Wilson B from a straight-line fit of ln<I> vs (sin th/lambda)^2.

    Uses reflections with d < d_cut (the ~linear high-resolution part of the
    Wilson plot). The atomic-scattering term is omitted, so the absolute value
    carries an offset that is identical across chunks (same composition and
    binning) and therefore cancels in dB. Returns B (A^2) or None."""
    try:
        import numpy as np
    except ImportError:
        return None
    if inv_d2 is None or len(inv_d2) < 200:
        return None
    sel = inv_d2 >= (1.0 / (d_cut * d_cut))
    if sel.sum() < 200:
        sel = np.ones(len(inv_d2), dtype=bool)
    x = inv_d2[sel] / 4.0                    # (sin theta / lambda)^2
    y = ivals[sel]
    edges = np.linspace(x.min(), x.max(), nbins + 1)
    idx = np.clip(np.digitize(x, edges) - 1, 0, nbins - 1)
    xs, ys = [], []
    for b in range(nbins):
        m = idx == b
        if m.sum() < 5:
            continue
        mi = y[m].mean()
        if mi <= 0:
            continue
        xs.append(0.5 * (edges[b] + edges[b + 1]))
        ys.append(math.log(mi))
    if len(xs) < 5:
        return None
    slope, _ = np.polyfit(np.array(xs), np.array(ys), 1)
    b = -slope / 2.0
    return b if math.isfinite(b) else None


def gemmi_available():
    try:
        import gemmi  # noqa: F401
        import numpy   # noqa: F401
        return True
    except ImportError:
        return False


def wilson_b_by_chunk(rows_meta, diagnostics, route):
    """Estimate Wilson B from each chunk's merged MTZ (only where not scraped).

    rows_meta: list of dicts with keys chunk, mtz, wilson_scraped.
    Returns {chunk: B}."""
    result = {}
    if not gemmi_available():
        diagnostics.append("    [%s] gemmi unavailable in this Python — Wilson B "
                           "(MTZ fallback) skipped; run with the trfrx venv to "
                           "enable it" % route)
        return result
    for meta in rows_meta:
        mtz = meta["mtz"]
        if not mtz or meta["wilson_scraped"] is not None:
            continue
        _data, inv_d2, ivals = _read_mtz_intensities(mtz)
        if inv_d2 is None:
            diagnostics.append("    [%s] could not read MTZ %s" % (route, mtz))
            continue
        b = wilson_b_from_arrays(inv_d2, ivals)
        if b is not None:
            result[meta["chunk"]] = b
            diagnostics.append("    [%s] Wilson B = %.2f <- %s (from merged I)"
                               % (route, b, os.path.basename(mtz)))
    return result


# --------------------------------------------------------------------------
# PDF report.
# --------------------------------------------------------------------------
# Colour palette (Tol "bright", teal as the 3rd colour). Applied per panel by
# series order; dashes are used ONLY where two lines would coincide (cell edge b
# on a, images used on requested). Everything else is a solid continuous step.
_PALETTE = ["#4477AA", "#EE6677", "#44AA99"]
_SOLID   = "-"
_DASH    = (0, (4, 2.2))
_INK, _MUT, _GRID, _SPINE = "#1a2530", "#5c6b76", "#eceff2", "#c4ced4"

# Linear image->seconds transform for the secondary (top) x-axis. Set once in
# main() to (dt_s, t0_image, t0_s) when a frame period is known; left None (no
# seconds axis) when time is UNRESOLVED. Read by _draw_panel.
_XTIME = None


def _step_xy(rows, col, pct=False):
    """(xs, ys) tracing a continuous STEP line: the value is held flat across each
    chunk's image window [image_first, image_last] and steps at every boundary, so
    a value is never placed at a single ambiguous x (it spans the window it was
    measured over). Drops N/A. *pct* expresses the series as % change from the
    first window."""
    vals = []
    for r in rows:
        d = dict(r)
        y = to_float(d.get(col))
        a = to_float(d.get("image_first"))
        b = to_float(d.get("image_last"))
        if y is None or a is None or b is None:
            continue
        vals.append((a, b, y))
    if pct:
        if not vals or vals[0][2] == 0:
            return [], []
        base = vals[0][2]
        vals = [(a, b, 100.0 * (y - base) / base) for a, b, y in vals]
    xs, ys = [], []
    for a, b, y in vals:
        xs += [a, b]
        ys += [y, y]
    return xs, ys


# Damage-report panels. Each: (title, y-label, invert-y, series, slug).
# series item: (legend, column, pct, colour_idx, dash). Colour is by series order
# from _PALETTE; dashes only separate coincident lines (a=b, requested=used).
# The slug is the filename stem for the per-panel publication figures.
_DAMAGE_PANELS = [
    ("High-resolution limit", "d_high (Å)", True,
        [("", "resolution_high", False, 0, _SOLID)], "resolution"),
    ("R-factors", "R", False,
        [("Rmerge", "Rmerge", False, 0, _SOLID),
         ("Rmeas", "Rmeas", False, 1, _SOLID),
         ("Rpim", "Rpim", False, 2, _SOLID)], "R-factors"),
    ("Mean I / σ(I)", "I/σ", False,
        [("overall", "Mean_I_over_sigma", False, 0, _SOLID),
         ("outer", "Mean_I_over_sigma_outer", False, 1, _SOLID)], "I-over-sigma"),
    ("ISa  (asymptotic I/σ)", "ISa", False,
        [("", "Isa", False, 0, _SOLID)], "ISa"),
    ("CC½", "CC½", False,
        [("overall", "CC_half", False, 0, _SOLID),
         ("outer", "CC_half_outer", False, 1, _SOLID)], "CC-half"),
    ("CC(anom)", "CC(ano)", False,
        [("overall", "CC_ano", False, 0, _SOLID),
         ("outer", "CC_ano_outer", False, 1, _SOLID)], "CC-ano"),
    # Spherical AND ellipsoidal completeness (overall). For STARANISO the
    # spherical value is low by design (anisotropic cut-off) while the ellipsoidal
    # value — completeness within the resolution surface actually used — stays
    # high; showing both makes that explicit. TRUNCATE has no ellipsoidal value,
    # so only the spherical line is drawn there. Outer-shell values are in the
    # per-chunk statistics table.
    ("Completeness", "completeness (%)", False,
        [("spherical", "Completeness", False, 0, _SOLID),
         ("ellipsoidal", "Completeness_ellipsoidal", False, 1, _SOLID)],
        "completeness"),
    ("Multiplicity", "multiplicity", False,
        [("overall", "Multiplicity", False, 0, _SOLID),
         ("outer", "Multiplicity_outer", False, 1, _SOLID)], "multiplicity"),
    ("Wilson B-factor", "B (Å²)", False,
        [("", "Wilson_B", False, 0, _SOLID)], "Wilson-B"),
    # Unit-cell edges as TRUE Ångström lengths (not % change). a = b by symmetry,
    # so b is dashed to stay visible where it coincides with a.
    ("Unit-cell edges", "edge (Å)", False,
        [("a", "cell_a", False, 0, _SOLID),
         ("b", "cell_b", False, 1, _DASH),
         ("c", "cell_c", False, 2, _SOLID)], "cell-edges"),
    ("Images per chunk", "N images", False,
        [("requested", "n_images", False, 0, _SOLID),
         ("used", "images_used", False, 1, _DASH)], "images"),
]


def _x_range(rows):
    """(x0, x1, xpad): tight image-number range straight from the windows."""
    xf = [to_float(dict(r).get("image_first")) for r in rows]
    xl = [to_float(dict(r).get("image_last")) for r in rows]
    xf = [v for v in xf if v is not None]
    xl = [v for v in xl if v is not None]
    x0 = min(xf) if xf else 0.0
    x1 = max(xl) if xl else 1.0
    return x0, x1, 0.012 * ((x1 - x0) or 1.0)


def _draw_panel(ax, panel, rows, xr):
    """Draw one damage panel (step lines + shared styling) onto *ax*.

    Shared by the combined PDF page and the standalone publication figures, so
    the two never drift apart. Returns True if any series was drawn."""
    ptitle, ylab, invert, specs, _slug = panel
    x0, x1, xpad = xr
    C = _PALETTE
    drew = False
    for lab, col, pct, ci, dash in specs:
        xs, ys = _step_xy(rows, col, pct)
        if not xs:
            continue
        ax.plot(xs, ys, color=C[ci % len(C)], lw=1.8, ls=dash,
                label=(lab or None), solid_capstyle="round")
        drew = True
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(_SPINE)
        ax.spines[sp].set_linewidth(0.9)
    ax.grid(axis="y", color=_GRID, lw=0.9)
    ax.set_axisbelow(True)
    ax.tick_params(length=0, labelsize=7.5, colors=_MUT)
    ax.set_title(ptitle, loc="left", fontsize=10.5, fontweight="bold",
                 color=_INK, pad=7)
    ax.set_ylabel(ylab, fontsize=8, color=_MUT)
    ax.set_xlabel("image number", fontsize=8, color=_MUT)
    ax.set_xlim(x0 - xpad, x1 + xpad)
    if _XTIME and _XTIME[0]:
        dt, ti, t0 = _XTIME
        secax = ax.secondary_xaxis(
            "top",
            functions=(lambda x, dt=dt, ti=ti, t0=t0: t0 + (x - ti) * dt,
                       lambda t, dt=dt, ti=ti, t0=t0: ti + (t - t0) / dt))
        secax.set_xlabel("time (s)", fontsize=8, color=_MUT)
        secax.tick_params(length=0, labelsize=7.5, colors=_MUT)
        for sp in secax.spines.values():
            sp.set_visible(False)
    if invert and drew:
        ax.invert_yaxis()
    if not drew:
        ax.text(0.5, 0.5, "no data", ha="center", va="center",
                fontsize=8, color=_MUT, transform=ax.transAxes)
    if drew and sum(1 for s in specs if s[0]) > 1:
        ax.legend(fontsize=7.5, frameon=False, loc="best", handlelength=2.0)
    return drew


def _damage_plots(pdf, plt, rows, dataset, title, wilson_method=""):
    """Combined report page: one continuous step line per metric, flat across
    each image window (see _step_xy). Panels: high-res limit, R-factors,
    I/sigma, CC(1/2), Wilson B, cell edges, images."""
    nrows = int(math.ceil(len(_DAMAGE_PANELS) / 2.0))
    # Page height scales with the number of rows so panels keep their aspect as
    # metrics are added/removed (2.35*nrows + 1.9 == 11.3 in at 4 rows).
    fig = plt.figure(figsize=(8.4, 2.35 * nrows + 1.9))
    fig.patch.set_facecolor("white")
    htop = 1.0 - 0.75 / (2.35 * nrows + 1.9)          # header just under the top
    fig.text(0.07, htop, dataset, fontsize=17, fontweight="bold", color=_INK)
    fig.text(0.07, htop - 0.017,
             "%s route  ·  radiation-damage indicators vs image number" % title,
             fontsize=9, color=_MUT)
    xr = _x_range(rows)
    for idx, panel in enumerate(_DAMAGE_PANELS, 1):
        _draw_panel(fig.add_subplot(nrows, 2, idx), panel, rows, xr)
    fig.subplots_adjust(left=0.09, right=0.955, top=htop - 0.035, bottom=0.04,
                        hspace=0.50, wspace=0.30)
    pdf.savefig(fig)
    plt.close(fig)


def _damage_figs(plt, rows, out_dir, route, suffix=""):
    """One standalone publication figure per panel — a 600-dpi PNG and a vector
    SVG — written to <out_dir>/plots_<route><suffix>/<slug>.{png,svg}."""
    sub = os.path.join(out_dir, "plots_%s%s" % (route, suffix))
    os.makedirs(sub, exist_ok=True)
    xr = _x_range(rows)
    for panel in _DAMAGE_PANELS:
        slug = panel[4]
        fig = plt.figure(figsize=(6.0, 3.6))
        fig.patch.set_facecolor("white")
        _draw_panel(fig.add_subplot(1, 1, 1), panel, rows, xr)
        fig.tight_layout()
        fig.savefig(os.path.join(sub, "%s.png" % slug), dpi=600,
                    bbox_inches="tight", facecolor="white")
        fig.savefig(os.path.join(sub, "%s.svg" % slug),
                    bbox_inches="tight", facecolor="white")   # vector, DPI-free
        plt.close(fig)
    return sub


# Display rows for the statistics table: (label, kind, *cols). kind is
# 'range' (low-high), 'ovouter' (overall (outer)), 'cell' (a b c) or 'plain'.
TABLE_SPEC = [
    ("Resolution (A)",          "range",   "resolution_low", "resolution_high"),
    ("Rmerge",                  "ovouter", "Rmerge"),
    ("Rmeas",                   "ovouter", "Rmeas"),
    ("Rpim",                    "ovouter", "Rpim"),
    ("I / sigma(I)",            "ovouter", "Mean_I_over_sigma"),
    ("ISa",                     "plain",   "Isa"),
    ("CC(1/2)",                 "ovouter", "CC_half"),
    ("CC(ano)",                 "ovouter", "CC_ano"),
    ("Completeness (%)",        "ovouter", "Completeness"),
    ("Completeness ellip. (%)", "ovouter", "Completeness_ellipsoidal"),
    ("Multiplicity",            "ovouter", "Multiplicity"),
    ("Wilson B (A^2)",          "plain",   "Wilson_B"),
    ("Unit cell a b c (A)",     "cell",    "cell_a", "cell_b", "cell_c"),
    ("Cell volume (A^3)",       "plain",   "cell_volume"),
    ("N observations",          "plain",   "N_observations"),
    ("N unique",                "plain",   "N_unique"),
    ("Images requested",        "plain",   "n_images"),
    ("Images used",             "plain",   "images_used"),
]


def _cell_text(spec, d):
    label, kind = spec[0], spec[1]
    if kind == "range":
        lo, hi = d.get(spec[2], "N/A"), d.get(spec[3], "N/A")
        return "%s - %s" % (lo, hi)
    if kind == "ovouter":
        ov = d.get(spec[2], "N/A")
        out = d.get(spec[2] + "_outer", "N/A")
        return "%s (%s)" % (ov, out)
    if kind == "cell":
        # Stack the three edges on their own lines so full precision still fits
        # inside one narrow column (the row is given extra height below).
        return "\n".join(str(d.get(c, "N/A")) for c in spec[2:])
    return str(d.get(spec[2], "N/A"))


def _stats_tables(pdf, plt, rows, dataset, title):
    # Balance chunks across pages so the last page is never a lone column
    # (e.g. 16 chunks -> 4+4+4+4, not 5+5+5+1).
    max_per_page = 5
    n = len(rows)
    n_pages = max(1, math.ceil(n / max_per_page))
    per_page = math.ceil(n / n_pages)
    first_w = 0.30
    # Fixed per-column width (based on a full page) so short pages keep the same
    # column size and are left-justified rather than stretched across the sheet.
    fixed_w = (1.0 - first_w) / max_per_page

    for start in range(0, n, per_page):
        page_rows = rows[start:start + per_page]
        dicts = [dict(r) for r in page_rows]
        col_labels = ["Statistic"] + [d["chunk"].replace("autoPROC_", "")
                                      for d in dicts]
        table_data = [[spec[0]] + [_cell_text(spec, d) for d in dicts]
                      for spec in TABLE_SPEC]

        fig = plt.figure(figsize=(8.27, 11.69))
        fig.suptitle("%s  —  %s\nstatistics per chunk (overall, outer shell in "
                     "parentheses)   [%d-%d]"
                     % (dataset, title, start + 1, start + len(page_rows)),
                     fontsize=11)
        ax = fig.add_axes([0.02, 0.03, 0.96, 0.9])
        ax.axis("off")
        col_widths = [first_w] + [fixed_w] * len(page_rows)
        tbl = ax.table(cellText=table_data, colLabels=col_labels,
                       colWidths=col_widths, loc="center", cellLoc="center")
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(7)
        tbl.scale(1, 1.5)
        # The multi-line unit-cell row (a/b/c stacked) needs extra height.
        cell_row = next((i for i, s in enumerate(TABLE_SPEC)
                         if s[1] == "cell"), None)
        cell_row = cell_row + 1 if cell_row is not None else None   # +1 for header
        # Left-align the statistic-name column, bold the header, grow the
        # stacked unit-cell row.
        for (r, c), cell in tbl.get_celld().items():
            if c == 0 and r > 0:
                cell.set_text_props(ha="left")
                cell.PAD = 0.03
            if r == 0:
                cell.set_text_props(weight="bold")
            if r == cell_row:
                cell.set_height(cell.get_height() * 2.8)
        pdf.savefig(fig)
        plt.close(fig)


def write_pdf(path, dataset, title, rows, image_summary, wilson_method="",
              limit_note="", route="", suffix=""):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
    except ImportError:
        sys.stderr.write("WARNING: matplotlib not available; skipping %s\n" % path)
        return False

    stamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    with PdfPages(path) as pdf:
        # ---- Cover / image-usage page ----
        fig = plt.figure(figsize=(8.27, 11.69))
        fig.text(0.5, 0.93, dataset, ha="center", fontsize=22, weight="bold")
        fig.text(0.5, 0.895, "%s route — radiation-damage report" % title,
                 ha="center", fontsize=13)
        ax = fig.add_axes([0.08, 0.06, 0.86, 0.80])
        ax.axis("off")
        lines = [
            "Generated: %s" % stamp,
            "Chunks (image ranges): %d" % len(image_summary),
            "Total images requested: %d" % sum(n for _, _, n in image_summary),
        ]
        if limit_note:
            lines += ["", limit_note]
        lines += ["", "Image ranges:"]
        for name, rng, n in image_summary:
            lines.append("    %-24s images %-14s (%d)" % (name, rng, n))
        ax.text(0.0, 1.0, "\n".join(lines), va="top", ha="left",
                family="monospace", fontsize=9)
        pdf.savefig(fig)
        plt.close(fig)

        if rows:
            _damage_plots(pdf, plt, rows, dataset, title, wilson_method)
            _stats_tables(pdf, plt, rows, dataset, title)

    # Publication figures: one 600-dpi PNG + one SVG per panel, per route.
    if rows and route:
        try:
            sub = _damage_figs(plt, rows, os.path.dirname(path), route, suffix)
            print("  per-panel figures (PNG+SVG) -> %s" % sub)
        except Exception as exc:                      # never fail the PDF over PNGs
            sys.stderr.write("WARNING: per-panel figures failed: %s\n" % exc)
    return True


def detect_dataset(process_dir, override):
    if override:
        return override
    p = process_dir.rstrip(os.sep)
    base = os.path.basename(p)
    if base == "autoproc_chunks":
        parent = os.path.basename(os.path.dirname(p))
        return parent or base
    return base


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--process-dir", required=True,
                        help="Path to the autoproc_chunks directory.")
    parser.add_argument("--out", default=None,
                        help="Output directory (default: <process-dir>/reports).")
    parser.add_argument("--dataset", default=None,
                        help="Dataset name shown at the top of the reports "
                             "(default: parent folder of autoproc_chunks).")
    parser.add_argument("--max-image", type=int, default=None, metavar="N",
                        help="Stop the statistics at image N: keep only chunks "
                             "that END at or before N (a chunk straddling N is "
                             "excluded, never truncated). E.g. --max-image 2000 "
                             "keeps ...1801-2000 and drops 2001-2200. Output "
                             "files get an _upto<N> suffix so the full report "
                             "is never overwritten.")
    # --- Real-time (acquisition time) options -----------------------------
    parser.add_argument("--image-template", default=None,
                        help="CBF template with #### OR an Eiger master .h5, used "
                             "to read the frame period / oscillation for "
                             "time_windows.csv.")
    parser.add_argument("--images-dir", default=None,
                        help="Directory holding the raw images (for --image-template).")
    parser.add_argument("--frame-time-ms", type=float, default=None, metavar="MS",
                        help="Override per-image frame period in ms (else read from header).")
    parser.add_argument("--osc-per-image", type=float, default=None, metavar="DEG",
                        help="Override oscillation per image in degrees (else from header).")
    parser.add_argument("--t0-image", type=int, default=None, metavar="N",
                        help="Image index taken as t=0. Default: the reference chunk's "
                             "LAST image, so the reference dataset is 0 s and each later "
                             "chunk's time is the time at its last image.")
    parser.add_argument("--t0-seconds", type=float, default=0.0, metavar="S",
                        help="Absolute seconds assigned to --t0-image (default 0).")
    args = parser.parse_args()

    process_dir = os.path.abspath(args.process_dir)
    out_dir = args.out or os.path.join(process_dir, "reports")
    os.makedirs(out_dir, exist_ok=True)
    dataset = detect_dataset(process_dir, args.dataset)

    # Read the frame period + oscillation ONCE (constant across a rotation series).
    # The time origin (t0_image) is resolved later, once the chunks are known.
    # Failures degrade to a warning; the CSV is still written with blank times
    # (source=UNRESOLVED).
    global _XTIME
    dt_s, osc_deg, tmeta = frame_period_from_header(
        args.images_dir, args.image_template, args.frame_time_ms, args.osc_per_image)
    t0_s = args.t0_seconds
    tw_source = tmeta.get("source", "UNRESOLVED")
    if dt_s is None:
        sys.stderr.write("WARNING: no frame period (%s) — time_windows.csv will have "
                         "blank times; pass --frame-time-ms/--osc-per-image or check "
                         "--image-template/--images-dir.\n" % tmeta.get("note", "no header"))

    def _scan_chunks(root):
        """(name, first, last, path) for autoPROC_<first>_<last> dirs in *root*."""
        found = []
        try:
            names = sorted(os.listdir(root))
        except OSError:
            return found
        for name in names:
            full = os.path.join(root, name)
            m = CHUNK_RE.match(name)
            if m and os.path.isdir(full):
                found.append((name, int(m.group(1)), int(m.group(2)), full))
        return found

    # Two layouts are supported:
    #   flat   : <process_dir>/autoPROC_1_150            (TR-FRX_autoPROC.sh)
    #   nested : <process_dir>/1_150/autoPROC_1_150      (trfrx_full_pipeline copy)
    entries = _scan_chunks(process_dir)
    if not entries:
        for sub in sorted(os.listdir(process_dir)):
            subdir = os.path.join(process_dir, sub)
            if os.path.isdir(subdir):
                entries.extend(_scan_chunks(subdir))
    entries.sort(key=lambda c: c[1])
    chunks = [(n, a, b) for n, a, b, _p in entries]
    all_chunks = list(chunks)   # full dataset list for time_windows.csv (pre --max-image)

    # Time origin: default t0 = the reference chunk's LAST image (all_chunks is sorted
    # by image_first, so all_chunks[0] is the reference). Then the reference dataset
    # ends at t=0 and each later chunk's t_end = (image_last - ref_last_image)*dt_s.
    # An explicit --t0-image (pump/mix trigger) overrides this default.
    if args.t0_image is not None:
        t0_image = args.t0_image
    else:
        t0_image = all_chunks[0][2] if all_chunks else 1
    _XTIME = (dt_s, t0_image, t0_s) if dt_s else None
    chunk_path = dict((n, p) for n, _a, _b, p in entries)
    if not chunks:
        sys.stderr.write("ERROR: no autoPROC_<first>_<last> chunk directory found "
                         "in %s (searched one level deep).\n" % process_dir)
        return 1

    # --max-image: keep whole chunks only. A chunk is kept when its LAST image is
    # <= N, so the cut lands on a chunk boundary and no partial window is ever
    # plotted (limit 2000 keeps 1801-2000, drops 2001-2200).
    limit_note = ""
    suffix = ""
    dropped = []
    if args.max_image is not None:
        keep = [c for c in chunks if c[2] <= args.max_image]
        dropped = [c for c in chunks if c[2] > args.max_image]
        if not keep:
            sys.stderr.write(
                "ERROR: --max-image %d excludes every chunk (the first one ends "
                "at %s).\n" % (args.max_image,
                               chunks[0][2] if chunks else "n/a"))
            return 1
        chunks = keep
        suffix = "_upto%d" % args.max_image
        limit_note = ("Limited to images <= %d  (last chunk kept: %s; %d chunk(s) "
                      "excluded)" % (args.max_image, chunks[-1][0], len(dropped)))
        print("Limit --max-image %d: keeping %d chunk(s) up to image %d; "
              "dropping %d." % (args.max_image, len(chunks), chunks[-1][2],
                                len(dropped)))
        for nm, a, b in dropped:
            print("    excluded: %s (images %d-%d)" % (nm, a, b))

    diagnostics = ["TR-FRX damage report - %s"
                   % datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
                   "Dataset: %s" % dataset,
                   "Process dir: %s" % process_dir,
                   "Chunks found: %d" % len(chunks),
                   "Image limit: %s" % (limit_note or "none (all chunks)"),
                   "Per-chunk sources:"]
    for nm, a, b in dropped:
        diagnostics.append("    EXCLUDED by --max-image: %s (images %d-%d)"
                           % (nm, a, b))

    best_by_chunk = {}
    image_summary = []
    images_used_by_chunk = {}
    meta = {"truncate": [], "staraniso": []}
    for name, first, last in chunks:
        image_summary.append((name, "%d-%d" % (first, last), last - first + 1))
        best = collect_chunk_blocks(chunk_path[name], diagnostics)
        best_by_chunk[name] = (first, last, best)
        iu = None
        for route in ("staraniso", "truncate"):
            bb = best[route] if best[route] else {}
            if bb.get("images_used") is not None:
                iu = bb.get("images_used"); break
        images_used_by_chunk[name] = iu
        mid = 0.5 * (first + last)
        for route in ("truncate", "staraniso"):
            b = best[route] if best[route] else {}
            meta[route].append({"chunk": name, "mid": mid,
                                "mtz": b.get("mtz"),
                                "wilson_scraped": b.get("wilson_b")})

    diagnostics.append("MTZ-derived Wilson B (fallback when not scraped):")
    wilson_mtz = {}
    for route in ("truncate", "staraniso"):
        wilson_mtz[route] = wilson_b_by_chunk(meta[route], diagnostics, route)

    def route_wilson(route, name):
        b = best_by_chunk[name][2][route]
        scraped = b.get("wilson_b") if b else None
        if scraped is not None:
            return scraped, (b.get("wb_source") or ("%s log" % route))
        mtzb = wilson_mtz[route].get(name)
        if mtzb is not None:
            return mtzb, "MTZ fit (FALLBACK - different method!)"
        return None, "N/A"

    def rows_for(route):
        rows = []
        for name, first, last in chunks:
            b = best_by_chunk[name][2][route]
            # Each route uses ONLY its own Wilson B (truncate -> CTRUNCATE,
            # staraniso -> STARANISO Popov-Bourenkov). No cross-route borrowing:
            # the two are computed with different normalisations, so mixing them
            # would be meaningless. MTZ fit is a clearly-flagged last resort only.
            wilson, wsource = route_wilson(route, name)
            images = b.get("images_used") if b else None
            ts, te, tm, _dur = time_window(first, last, dt_s, t0_image, t0_s)
            rows.append(build_row(name, first, last, b, wilson, images, wsource,
                                  t_start=ts, t_end=te, t_mid=tm))
        return rows

    wmethods = {
        "truncate":  "CTRUNCATE Wilson plot (average-atom, composition from cell)",
        "staraniso": "STARANISO Wilson B (Popov & Bourenkov 2003 normalisation)",
    }
    for route, pretty in [("truncate", "Classical autoPROC (TRUNCATE)"),
                          ("staraniso", "STARANISO")]:
        rows = rows_for(route)
        write_csv(os.path.join(out_dir, "%s_statistics%s.csv" % (route, suffix)),
                  rows)
        write_pdf(os.path.join(out_dir, "%s_report%s.pdf" % (route, suffix)),
                  dataset, pretty, rows, image_summary, wmethods[route],
                  limit_note, route=route, suffix=suffix)

    with open(os.path.join(out_dir,
                           "parsing_diagnostics%s.txt" % suffix), "w") as fh:
        fh.write("\n".join(diagnostics) + "\n")

    # Canonical per-dataset time table — always written, over the FULL chunk set
    # (independent of --max-image), and consumed by trfrx_full_pipeline.py et al.
    tw_path = os.path.join(out_dir, "time_windows.csv")
    write_time_windows_csv(tw_path, all_chunks, dt_s, osc_deg, t0_image, t0_s,
                           tw_source, images_used_by_chunk)
    total_s = (len(all_chunks) and dt_s
               and sum((b - a + 1) for _n, a, b in all_chunks) * dt_s)
    _t0kind = ("end of reference (auto)" if args.t0_image is None else "explicit")
    print("time_windows.csv: dt=%s s, osc=%s deg, %d datasets, total=%s s, source=%s"
          % ("%.6g" % dt_s if dt_s else "N/A",
             "%.4g" % osc_deg if osc_deg else "N/A",
             len(all_chunks),
             "%.4g" % total_s if total_s else "N/A", tw_source))
    print("  t0 = image %s (%s), t0_s=%s; timepoint = END of each window (t_end_s)"
          % (t0_image, _t0kind, ("%.4g" % t0_s)))

    print("Reports written to %s" % out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
PYEOF

# Optional image limit, passed through to the report generator.
REPORT_ARGS=""
if [ -n "$MAX_IMAGE" ]; then
    REPORT_ARGS="--max-image $MAX_IMAGE"
fi

# Image context + timing for computing the time windows (time_windows.csv).
# The quotes are preserved as-is in generate_reports.sh (unquoted EOF heredoc)
# then re-interpreted at run time — safe even if a path contains a space.
REPORT_ARGS="$REPORT_ARGS --image-template \"$IMAGE_TEMPLATE\" --images-dir \"$IMAGES_DIR\" --t0-seconds \"$T0_SECONDS\""
# --t0-image only when set; blank => report defaults it to the reference chunk's last image.
if [ -n "$T0_IMAGE" ]; then
    REPORT_ARGS="$REPORT_ARGS --t0-image \"$T0_IMAGE\""
fi
if [ -n "$FRAME_TIME_MS" ]; then
    REPORT_ARGS="$REPORT_ARGS --frame-time-ms \"$FRAME_TIME_MS\""
fi
if [ -n "$OSC_PER_IMAGE" ]; then
    REPORT_ARGS="$REPORT_ARGS --osc-per-image \"$OSC_PER_IMAGE\""
fi

REPORT_JOB_SCRIPT="${PROCESS_DIR}/generate_reports.sh"
cat > "$REPORT_JOB_SCRIPT" <<EOF
#!/bin/bash
set -euo pipefail
# autoPROC is NOT required for the reports (python + matplotlib + gemmi are
# enough): load the module if it exists, without letting the script fail.
if type module >/dev/null 2>&1; then module load autoPROC 2>/dev/null || true; fi
cd "$PROCESS_DIR"
# Prefer the trfrx venv (created by setup_env.sh): it carries matplotlib and
# gemmi, needed for the PDF plots and the MTZ-derived Wilson B. Fall back to any
# python3/python on PATH (the CSVs and scraped stats still work without them).
if [ -x "\$HOME/.venv/trfrx/bin/python" ]; then
    PY="\$HOME/.venv/trfrx/bin/python"
elif command -v python3 >/dev/null 2>&1; then
    PY=python3
else
    PY=python
fi
"\$PY" "$PROCESS_DIR/trfrx_damage_report.py" --process-dir "$PROCESS_DIR" $REPORT_ARGS
EOF
chmod +x "$REPORT_JOB_SCRIPT"

if [ "$STATS_ONLY" -eq 1 ]; then
    # Without a limit, the files take their default names: warn before
    # overwriting existing reports (e.g. those copied by trfrx_full_pipeline).
    if [ -z "$MAX_IMAGE" ]; then
        # NB: no `ls` here — a glob with no match would return a non-zero code
        # and `set -e`/`pipefail` would kill the script silently.
        EXISTING=0
        for _f in "$PROCESS_DIR"/reports/*_report.pdf \
                  "$PROCESS_DIR"/reports/*_statistics.csv; do
            if [ -e "$_f" ]; then EXISTING=$((EXISTING + 1)); fi
        done
        if [ "$EXISTING" -gt 0 ]; then
            echo "WARNING: $EXISTING report file(s) already exist in"
            echo "            $PROCESS_DIR/reports and will be OVERWRITTEN."
            echo "            (a --max-image limit adds a _upto<N> suffix and"
            echo "             keeps the existing files)"
            printf "Continue and overwrite? [y/N]: "
            if ! read -r CONFIRM; then CONFIRM=""; echo; fi
            case "$(printf '%s' "$CONFIRM" | tr '[:upper:]' '[:lower:]')" in
                y|yes) : ;;
                *) echo "Aborted (no file modified)."; exit 0 ;;
            esac
        fi
    fi
    # Reports only: nothing to wait for, generate right away (a few seconds).
    echo "Generating reports..."
    bash "$REPORT_JOB_SCRIPT"
    echo ""
    echo "Reports in: $PROCESS_DIR/reports"
    if [ -n "$MAX_IMAGE" ]; then
        echo "  (files suffixed _upto${MAX_IMAGE}; the full report is kept)"
    fi
    exit 0
fi

# NB: afterany (not afterok) => the report is ALWAYS generated once all chunks
# have finished, even if some failed (in that case we simply get N/A for those
# chunks instead of no report at all).
DEP=$(IFS=:; echo "${CHUNK_JOBIDS[*]}")
report_output=$(sbatch -p "$SLURM_PARTITION" -n 1 -c 1 --mem=8000 --dependency=afterany:"$DEP" "$REPORT_JOB_SCRIPT")
REPORT_JOBID=$(echo "$report_output" | awk '{print $4}')
echo "  [submitted] job $REPORT_JOBID: report generation (depends on $DEP)"

# --- Summary of all submitted jobs ---
echo ""
echo "==================== Submitted jobs ===================="
printf "  %-12s %s\n" "JOB ID" "IMAGES"
for idx in "${!CHUNK_JOBIDS[@]}"; do
    printf "  %-12s %s\n" "${CHUNK_JOBIDS[$idx]}" "${CHUNK_LABELS[$idx]}"
done
printf "  %-12s %s\n" "$REPORT_JOBID" "final reports"
echo "======================================================="
echo "Total: $(( ${#CHUNK_JOBIDS[@]} + 1 )) jobs (${#CHUNK_JOBIDS[@]} processing + 1 report)"
echo ""
echo "Track with: squeue -j $(IFS=,; echo "${CHUNK_JOBIDS[*]},$REPORT_JOBID")"
echo ""
echo "Processing submitted. Results in: $PROCESS_DIR"
echo "Final reports (generated after processing) in: $PROCESS_DIR/reports"
echo ""
echo "Tip: to regenerate only the reports + time_windows.csv (without reprocessing):"
echo "         $0 --stats \"$OUTPUT_DIR\""
echo "         (for real times, add --frame-time-ms MS or --images-dir/--image-template)"
