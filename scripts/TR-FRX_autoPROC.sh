#!/bin/bash
#
# TR-FRX autoPROC chunked processing + consolidated reports.
#
# Handles BOTH input types automatically:
#   * CBF   : a template with #### placeholders, e.g. PfuGRHPR_006_1_####.cbf
#   * HDF5  : an Eiger master file,               e.g. PfuGR_003_001_master.h5
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
#   ./TR-FRX_autoPROC.sh PfuGR_003_001_master.h5   # override the input file
#
#   # Read images from one place, write results to another:
#   ./TR-FRX_autoPROC.sh <input_dir> <output_dir>
#   ./TR-FRX_autoPROC.sh <input_dir> <output_dir> PfuGR_003_001_master.h5

set -euo pipefail

# =========================== CONFIGURATION ===========================
RUN_ID=006

# CBF template (with ####) OR Eiger master .h5 file. Auto-detected below.
IMAGE_TEMPLATE="PfuGRHPR_006_1_####.cbf"

FIRST_IMG=1
LAST_IMG=3000

REF_FIRST_IMG=1
REF_LAST_IMG=300
CHUNK_SIZE=300

SYMM="I41"
CELL="114 114 118 90 90 90"

SLURM_PARTITION="nice"
SLURM_CPUS=24
SLURM_MEM=24000
# =====================================================================

# ----------------------- Arguments en ligne de commande -----------------------
# Quatre façons de lancer le script :
#   ./TR-FRX_autoPROC.sh                                     # ancienne méthode :
#                                                            #   images ET résultats
#                                                            #   dans le dossier courant
#   ./TR-FRX_autoPROC.sh master.h5                           # + override du fichier
#   ./TR-FRX_autoPROC.sh <input_dir> <output_dir>           # lire ici, écrire là
#   ./TR-FRX_autoPROC.sh <input_dir> <output_dir> master.h5 # + override du fichier
# -----------------------------------------------------------------------------
IMAGES_DIR="$(pwd)"
OUTPUT_DIR="$(pwd)"

case "$#" in
    0) : ;;                                       # tout par défaut (dossier courant)
    1) IMAGE_TEMPLATE="$1" ;;                      # override du fichier uniquement
    2) IMAGES_DIR="$1"; OUTPUT_DIR="$2" ;;         # lire dans $1, écrire dans $2
    3) IMAGES_DIR="$1"; OUTPUT_DIR="$2"; IMAGE_TEMPLATE="$3" ;;
    *)
        echo "Usage : $0 [<input_dir> <output_dir>] [image_file]" >&2
        exit 1
        ;;
esac

# --- Chemins absolus (les nœuds de calcul SLURM doivent voir ces dossiers) ---
if [ ! -d "$IMAGES_DIR" ]; then
    echo "ERREUR : dossier d'entrée introuvable : $IMAGES_DIR" >&2
    exit 1
fi
IMAGES_DIR="$(cd "$IMAGES_DIR" && pwd)"
mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- Détection du mode d'après l'extension du fichier d'entrée ---
case "$IMAGE_TEMPLATE" in
    *.h5|*.H5|*.nxs|*.NXS) INPUT_MODE="HDF5" ;;
    *.cbf|*.CBF)           INPUT_MODE="CBF" ;;
    *)
        echo "ERREUR : type d'entrée non reconnu pour '$IMAGE_TEMPLATE' (attendu .cbf ou .h5)" >&2
        exit 1
        ;;
esac
echo "Mode d'entrée : $INPUT_MODE ($IMAGE_TEMPLATE)"

module load autoPROC

# --- Tous les traitements seront créés dans un sous-dossier du dossier de sortie ---
PROCESS_DIR="${OUTPUT_DIR}/autoproc_chunks"
mkdir -p "$PROCESS_DIR"

# --- Identifiants des jobs de traitement (pour le rapport final) ---
CHUNK_JOBIDS=()
CHUNK_LABELS=()

FIRST_OUTDIR=""
FIRST_JOBID=""

# --- Traitement du jeu de référence ---
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
CHUNK_LABELS+=("${REF_FIRST_IMG}-${REF_LAST_IMG} (référence)")
echo "  [soumis] job $FIRST_JOBID : images ${REF_FIRST_IMG}-${REF_LAST_IMG} (référence)"

# --- Traitement des jeux suivants par blocs ---
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
    echo "  [soumis] job $CHUNK_JOBID : images ${i}-${j} (dépend de $FIRST_JOBID)"

    i=$((j + 1))
done

# =====================================================================
# --- Générateur de rapports (écrit dans PROCESS_DIR, puis lancé) ------
# =====================================================================
cat > "$PROCESS_DIR/trfrx_damage_report.py" <<'PYEOF'
#!/usr/bin/env python3
"""TR-FRX radiation-damage report from autoPROC chunk statistics.

Standalone tool (development version — not yet wired into TR-FRX_autoPROC.sh).
Run it after processing, pointed at an existing ``autoproc_chunks`` directory:

    python trfrx_damage_report.py --process-dir /path/to/autoproc_chunks
    python trfrx_damage_report.py --process-dir ... --dataset CaMDH_012

It needs only libraries already provided by setup_env.sh (matplotlib, gemmi,
numpy for the plots / Wilson-B / R_d; the CSV and scraped stats work with the
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
    * CC(1/2) vs dose                             (falling = damage)
    * Wilson B-factor and its increase dB vs dose (rising = loss of
      high-resolution order; the global B-factor "scaling"/decay curve)
    * Unit-cell edges and volume, % change vs dose (expansion = damage)
    * Cross-sweep R_d vs dose                      (rising = damage)
    * Images actually used per chunk vs dose

Scalar statistics come from the AIMLESS "Summary data" tables and the per-route
``*.table1`` files. Unit cells are scraped from the same files. The Wilson
B-factor is scraped when autoPROC prints it and otherwise estimated directly
from the merged MTZ (relative B: the offset from omitting the atomic scattering
term cancels in dB across chunks). The cross-sweep R_d and the Wilson B both use
``gemmi`` to read the merged MTZs.

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
                 "wb_source": None,
                 "cell_rank": -1, "mtz_rank": -1, "wb_rank": -1, "img_rank": -1}
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

    def consider_extra(route, fname, rel, cell, wilson, images):
        prio = source_priority(fname)
        e = extra[route]
        if cell is not None and prio > e["cell_rank"]:
            e["cell"], e["cell_rank"] = cell, prio
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

    for root, _dirs, files in os.walk(chunk_dir):
        for fname in files:
            low = fname.lower()
            fpath = os.path.join(root, fname)
            rel = os.path.relpath(fpath, chunk_dir)
            if low.endswith(".mtz"):
                route = "staraniso" if "staraniso" in low else "truncate"
                consider_mtz(route, fname, fpath)
                continue
            if not (low.endswith(".table1") or low.endswith(".log")
                    or low.endswith(".lp")):
                continue
            try:
                with open(fpath, "r", errors="replace") as fh:
                    text = fh.read(2_000_000)     # cap: XDS .LP can be huge
            except (IOError, OSError):
                continue
            route_of_file = "staraniso" if "staraniso" in low else "truncate"
            cell = scrape_cell(text)
            wilson = (scrape_staraniso_wilson_b(text) if route_of_file == "staraniso"
                      else scrape_wilson_b(text))
            images = scrape_images_used(text)
            images = images if images[0] is not None else None
            # Unit cell / images used are route-agnostic — record for both so the
            # STARANISO route is not left blank when only a truncate log has them.
            for r in ("truncate", "staraniso"):
                consider_extra(r, fname, rel, cell,
                               wilson if r == route_of_file else None, images)
            if low.endswith(".table1"):
                consider(route_of_file, parse_table1(text), fname, rel)
            elif low.endswith(".log"):
                for block in find_summary_blocks(text):
                    route = ("staraniso" if is_staraniso(fpath, block["title"])
                             else "truncate")
                    consider(route, block["metrics"], fname, rel)

    result = {}
    for route in ("truncate", "staraniso"):
        data = best[route] or {"metrics": {}, "source": "N/A"}
        e = extra[route]
        data.update({"cell": e["cell"], "mtz": e["mtz"],
                     "wilson_b": e["wilson_b"], "wb_source": e["wb_source"],
                     "images_used": e["images_used"]})
        has_any = (best[route] or e["cell"] or e["mtz"]
                   or e["wilson_b"] is not None or e["images_used"] is not None)
        result[route] = data if has_any else None
    return result


def match_metric(metrics, must, mustnot):
    for label, values in metrics.items():
        low = label.lower()
        if all(x in low for x in must) and not any(x in low for x in mustnot):
            return values
    return None


def build_row(chunk_name, first, last, block, wilson_b, images_used,
              wilson_source="N/A"):
    """Return an ordered list of (column, value) for one chunk."""
    row = [
        ("chunk", chunk_name),
        ("image_first", first),
        ("image_last", last),
        ("n_images", last - first + 1),
        ("images_used", images_used if images_used is not None else "N/A"),
    ]
    metrics = block["metrics"] if block else {}
    for col, must, mustnot, want_outer in METRICS:
        values = match_metric(metrics, must, mustnot)
        row.append((col, values[0] if values else "N/A"))
        if want_outer:
            row.append((col + "_outer", values[2] if values else "N/A"))
    cell = block.get("cell") if block else None
    vol = cell_volume(cell) if cell else None
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
# gemmi-based quantities: Wilson B (from merged intensities) and cross-sweep R_d.
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
def _series(rows, col):
    """(x, y) of mean image number vs numeric metric, dropping N/A."""
    xs, ys = [], []
    for r in rows:
        d = dict(r)
        y = to_float(d.get(col))
        first = to_float(d.get("image_first"))
        last = to_float(d.get("image_last"))
        if y is None or first is None or last is None:
            continue
        xs.append(0.5 * (first + last))
        ys.append(y)
    return xs, ys


def _pct_change_series(rows, col):
    xs, ys = _series(rows, col)
    if not ys or ys[0] == 0:
        return [], []
    base = ys[0]
    return xs, [100.0 * (y - base) / base for y in ys]


def _damage_plots(pdf, plt, rows, dataset, title, wilson_method=""):
    """Radiation-damage plots (metrics vs image number)."""
    fig = plt.figure(figsize=(8.27, 11.69))
    fig.suptitle("%s  —  %s\nradiation-damage indicators vs image number"
                 % (dataset, title), fontsize=13, y=0.98)
    xlabel = "Image number"

    def ax_at(idx):
        ax = fig.add_subplot(4, 2, idx)
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=7)
        ax.set_xlabel(xlabel, fontsize=7)
        return ax

    # 1) High-resolution limit vs dose. Y-axis inverted so that better
    #    resolution (smaller d) is at the top: the curve then falls as damage
    #    pushes the limit to higher d.
    ax = ax_at(1)
    xs, ys = _series(rows, "resolution_high")
    if xs:
        ax.plot(xs, ys, "o-", ms=3, color="C0")
        ax.invert_yaxis()
    ax.set_title("High-resolution limit (better = up)", fontsize=8)
    ax.set_ylabel("d_high (A)", fontsize=7)

    # 2) R-factors vs dose (rising = damage).
    ax = ax_at(2)
    for col, lab in [("Rmerge", "Rmerge"), ("Rmeas", "Rmeas"), ("Rpim", "Rpim")]:
        xs, ys = _series(rows, col)
        if xs:
            ax.plot(xs, ys, "o-", ms=3, label=lab)
    ax.set_title("R-factors (overall)", fontsize=8)
    ax.set_ylabel("R", fontsize=7)
    ax.legend(fontsize=6)

    # 3) Mean I/sigma(I) vs dose (falling = damage).
    ax = ax_at(3)
    for col, lab in [("Mean_I_over_sigma", "overall"),
                     ("Mean_I_over_sigma_outer", "outer shell")]:
        xs, ys = _series(rows, col)
        if xs:
            ax.plot(xs, ys, "o-", ms=3, label=lab)
    ax.set_title("Mean I / sigma(I)", fontsize=8)
    ax.set_ylabel("I/sigma(I)", fontsize=7)
    ax.legend(fontsize=6)

    # 4) CC(1/2) vs dose (falling = damage).
    ax = ax_at(4)
    for col, lab in [("CC_half", "overall"), ("CC_half_outer", "outer shell")]:
        xs, ys = _series(rows, col)
        if xs:
            ax.plot(xs, ys, "o-", ms=3, label=lab)
    ax.set_title("CC(1/2)", fontsize=8)
    ax.set_ylabel("CC(1/2)", fontsize=7)
    ax.legend(fontsize=6)

    # 5) Wilson B (global decay / B-scaling curve).
    ax = ax_at(5)
    xs, ys = _series(rows, "Wilson_B")
    if xs:
        ax.plot(xs, ys, "o-", ms=3, color="C3")
    else:
        ax.text(0.5, 0.5, "Wilson B unavailable", ha="center", va="center",
                fontsize=7, transform=ax.transAxes)
    ax.set_title(("Wilson B-factor\n(%s)" % wilson_method) if wilson_method
                 else "Wilson B-factor", fontsize=6.5)
    ax.set_ylabel("Wilson B (A^2)", fontsize=7)

    # 6) Unit-cell edges, % change vs dose. Always draw all three edges; distinct
    #    line styles/markers (and some transparency) keep 'a' visible even when a
    #    symmetry constraint makes it coincide with another edge.
    ax = ax_at(6)
    for k, style in [("a", "o-"), ("b", "s--"), ("c", "^:")]:
        xk, yk = _pct_change_series(rows, "cell_" + k)
        if xk:
            ax.plot(xk, yk, style, ms=3, alpha=0.8, label=k)
    ax.set_title("Unit-cell edges (% change)", fontsize=8)
    ax.set_ylabel("delta edge (%)", fontsize=7)
    ax.legend(fontsize=6)

    # 7) Unit-cell volume, % change vs dose.
    ax = ax_at(7)
    xs, ys = _pct_change_series(rows, "cell_volume")
    if xs:
        ax.plot(xs, ys, "o-", ms=3, color="C2")
    ax.set_title("Unit-cell volume (% change)", fontsize=8)
    ax.set_ylabel("delta V (%)", fontsize=7)

    # 8) Images actually used per chunk vs dose.
    ax = ax_at(8)
    xs, ys = _series(rows, "images_used")
    xr, yr = _series(rows, "n_images")
    if xr:
        ax.plot(xr, yr, "o--", ms=3, color="0.6", label="requested")
    if xs:
        ax.plot(xs, ys, "o-", ms=3, color="C1", label="used")
    ax.set_title("Images per chunk", fontsize=8)
    ax.set_ylabel("N images", fontsize=7)
    ax.legend(fontsize=6)

    fig.subplots_adjust(left=0.09, right=0.9, top=0.90, bottom=0.06,
                        hspace=0.5, wspace=0.45)
    pdf.savefig(fig)
    plt.close(fig)


# Display rows for the statistics table: (label, kind, *cols). kind is
# 'range' (low-high), 'ovouter' (overall (outer)), 'cell' (a b c) or 'plain'.
TABLE_SPEC = [
    ("Resolution (A)",          "range",   "resolution_low", "resolution_high"),
    ("Rmerge",                  "ovouter", "Rmerge"),
    ("Rmeas",                   "ovouter", "Rmeas"),
    ("Rpim",                    "ovouter", "Rpim"),
    ("I / sigma(I)",            "ovouter", "Mean_I_over_sigma"),
    ("CC(1/2)",                 "ovouter", "CC_half"),
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


def write_pdf(path, dataset, title, rows, image_summary, wilson_method=""):
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
            "",
            "Image ranges:",
        ]
        for name, rng, n in image_summary:
            lines.append("    %-24s images %-14s (%d)" % (name, rng, n))
        ax.text(0.0, 1.0, "\n".join(lines), va="top", ha="left",
                family="monospace", fontsize=9)
        pdf.savefig(fig)
        plt.close(fig)

        if rows:
            _damage_plots(pdf, plt, rows, dataset, title, wilson_method)
            _stats_tables(pdf, plt, rows, dataset, title)
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
    args = parser.parse_args()

    process_dir = os.path.abspath(args.process_dir)
    out_dir = args.out or os.path.join(process_dir, "reports")
    os.makedirs(out_dir, exist_ok=True)
    dataset = detect_dataset(process_dir, args.dataset)

    chunks = []
    for name in sorted(os.listdir(process_dir)):
        m = CHUNK_RE.match(name)
        full = os.path.join(process_dir, name)
        if m and os.path.isdir(full):
            chunks.append((name, int(m.group(1)), int(m.group(2))))
    chunks.sort(key=lambda c: c[1])

    diagnostics = ["TR-FRX damage report - %s"
                   % datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
                   "Dataset: %s" % dataset,
                   "Process dir: %s" % process_dir,
                   "Chunks found: %d" % len(chunks),
                   "Per-chunk sources:"]

    best_by_chunk = {}
    image_summary = []
    meta = {"truncate": [], "staraniso": []}
    for name, first, last in chunks:
        image_summary.append((name, "%d-%d" % (first, last), last - first + 1))
        best = collect_chunk_blocks(os.path.join(process_dir, name), diagnostics)
        best_by_chunk[name] = (first, last, best)
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
            rows.append(build_row(name, first, last, b, wilson, images, wsource))
        return rows

    wmethods = {
        "truncate":  "CTRUNCATE Wilson plot (average-atom, composition from cell)",
        "staraniso": "STARANISO Wilson B (Popov & Bourenkov 2003 normalisation)",
    }
    for route, pretty in [("truncate", "Classical autoPROC (TRUNCATE)"),
                          ("staraniso", "STARANISO")]:
        rows = rows_for(route)
        write_csv(os.path.join(out_dir, "%s_statistics.csv" % route), rows)
        write_pdf(os.path.join(out_dir, "%s_report.pdf" % route),
                  dataset, pretty, rows, image_summary, wmethods[route])

    with open(os.path.join(out_dir, "parsing_diagnostics.txt"), "w") as fh:
        fh.write("\n".join(diagnostics) + "\n")

    print("Reports written to %s" % out_dir)


if __name__ == "__main__":
    main()
PYEOF

REPORT_JOB_SCRIPT="${PROCESS_DIR}/generate_reports.sh"
cat > "$REPORT_JOB_SCRIPT" <<EOF
#!/bin/bash
set -euo pipefail
module load autoPROC
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
"\$PY" "$PROCESS_DIR/trfrx_damage_report.py" --process-dir "$PROCESS_DIR"
EOF

# NB : afterany (et non afterok) => le rapport est TOUJOURS généré une fois que
# tous les chunks sont terminés, même si certains ont échoué (on aura alors
# simplement des N/A pour ces chunks au lieu d'aucun rapport du tout).
DEP=$(IFS=:; echo "${CHUNK_JOBIDS[*]}")
report_output=$(sbatch -p "$SLURM_PARTITION" -n 1 -c 1 --mem=8000 --dependency=afterany:"$DEP" "$REPORT_JOB_SCRIPT")
REPORT_JOBID=$(echo "$report_output" | awk '{print $4}')
echo "  [soumis] job $REPORT_JOBID : génération des rapports (dépend de $DEP)"

# --- Récapitulatif de tous les jobs soumis ---
echo ""
echo "==================== Jobs soumis ===================="
printf "  %-12s %s\n" "JOB ID" "IMAGES"
for idx in "${!CHUNK_JOBIDS[@]}"; do
    printf "  %-12s %s\n" "${CHUNK_JOBIDS[$idx]}" "${CHUNK_LABELS[$idx]}"
done
printf "  %-12s %s\n" "$REPORT_JOBID" "rapports finaux"
echo "====================================================="
echo "Total : $(( ${#CHUNK_JOBIDS[@]} + 1 )) jobs (${#CHUNK_JOBIDS[@]} traitements + 1 rapport)"
echo ""
echo "Suivi : squeue -j $(IFS=,; echo "${CHUNK_JOBIDS[*]},$REPORT_JOBID")"
echo ""
echo "Traitements soumis. Résultats dans : $PROCESS_DIR"
echo "Rapports finaux (générés après traitement) dans : $PROCESS_DIR/reports"
