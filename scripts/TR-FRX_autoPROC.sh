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
# including how many images were used and which image ranges.
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
cat > "$PROCESS_DIR/autoproc_reports.py" <<'PYEOF'
#!/usr/bin/env python3
"""Consolidate autoPROC chunk statistics into per-route reports.

For a given ``autoproc_chunks`` directory (the PROCESS_DIR created by
TR-FRX_autoPROC.sh), this script scans every ``autoPROC_<first>_<last>`` chunk
sub-directory, extracts the final merging statistics for both processing routes
and writes, into a ``reports/`` sub-folder:

    reports/staraniso_statistics.csv   reports/staraniso_report.pdf
    reports/truncate_statistics.csv    reports/truncate_report.pdf
    reports/parsing_diagnostics.txt

* "staraniso" = the STARANISO (anisotropy-corrected) route.
* "truncate"  = the classical autoPROC route (truncate-unique.mtz).

Each report also records how many images were used and which image ranges
(derived from the chunk directory names) went into the merge.

The statistics are read from the AIMLESS-style "Summary data" table that
autoPROC writes into its log files. Any statistic that cannot be found is
reported as ``N/A`` rather than causing a failure; the diagnostics file lists
every candidate file that was scanned so the parser can be tightened if a site
uses non-standard filenames.

Usage:
    python autoproc_reports.py --process-dir /path/to/autoproc_chunks
"""

import argparse
import csv
import datetime
import os
import re
import sys


# --------------------------------------------------------------------------
# Statistics extracted for the consolidated "Table 1" style report.
# Each entry: (csv_column, [keywords that must appear in the AIMLESS label],
#              [keywords that must NOT appear], want_outer_shell)
# Labels are matched case-insensitively against the AIMLESS row label.
# --------------------------------------------------------------------------
METRICS = [
    ("resolution_low",        ["low resolution"],              [],            False),
    ("resolution_high",       ["high resolution"],             [],            True),
    ("Rmerge",                ["rmerge", "all"],               [],            True),
    ("Rmeas",                 ["rmeas", "all"],                [],            True),
    ("Rpim",                  ["rpim", "all"],                 [],            True),
    ("Mean_I_over_sigma",     ["mean", "sd(i)"],               [],            True),
    ("CC_half",               ["cc(1/2)"],                     [],            True),
    # STARANISO reports both spherical and ellipsoidal completeness; the
    # classical route reports a single "Completeness" (matched by the spherical
    # rows too, since it carries neither keyword). "ellipsoidal" is therefore
    # N/A for the classical route.
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


def find_summary_blocks(text):
    """Return a list of {'title', 'metrics'} for every AIMLESS Summary data
    table found in *text*. metrics maps label -> (overall, inner, outer)."""
    lines = text.splitlines()
    n = len(lines)
    blocks = []
    i = 0
    while i < n:
        if "Summary data for" in lines[i]:
            title = lines[i].strip()
            # Locate the "Overall  InnerShell  OuterShell" header nearby.
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
                    # End of the summary table. AIMLESS closes it with a
                    # "$$ <!--SUMMARY_END-->" line; the free-text sections that
                    # follow ("Estimates of resolution limits", "Average unit
                    # cell", ...) contain stray numbers we must NOT parse.
                    if (stripped.startswith("$$")
                            or "SUMMARY_END" in stripped
                            or stripped.startswith("====")
                            or stripped.startswith("Estimates of resolution")
                            or "Summary data for" in line):
                        break
                    # Rows like "Rmerge in top intensity bin  0.255  -  -" carry
                    # dashes instead of 3 numbers: skip them WITHOUT ending the
                    # table (the old code stopped here and lost CC(1/2),
                    # completeness, multiplicity, anomalous stats, ...).
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
    """Parse an autoPROC/Global-Phasing ``*.table1`` summary file.

    These files hold the official per-route "Table 1" statistics as a simple
    ``label  Overall  InnerShell  OuterShell`` table under an
    "Overall InnerShell OuterShell" header. Returns a metrics dict
    label -> (overall, inner, outer), or {} if no table is found.
    """
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
        if not stripped or set(stripped) <= set("-"):   # blank or "----" rule
            continue
        m = ROW_RE.match(line)
        if m:
            label = re.sub(r"\s+", " ", m.group(1).strip())
            if label:
                metrics[label] = (m.group(2), m.group(3), m.group(4))
    return metrics


def source_priority(fname):
    """Higher = more authoritative source for the merging statistics.

    autoPROC's per-route ``*-unique.table1`` files are the official Table-1
    summaries and are preferred over scraping AIMLESS logs; the plain
    ``aimless.log`` is preferred over the ``_early``/``_late`` (radiation-damage
    half-sets) and ``_alldata`` variants."""
    f = fname.lower()
    if f in ("truncate-unique.table1", "staraniso_alldata-unique.table1"):
        return 100
    if f.endswith("-unique.table1"):        # e.g. aimless_alldata-unique.table1
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
    """Scan a chunk directory and return the best (truncate, staraniso) metrics.

    Both ``*.table1`` files (preferred) and ``*.log`` AIMLESS summaries are
    considered; the route is taken from the filename ("staraniso" vs anything
    else) and the highest-priority source wins (ties break on metric count)."""
    best = {"truncate": None, "staraniso": None}
    best_rank = {"truncate": (-1, -1), "staraniso": (-1, -1)}

    def consider(route, metrics, fname, rel):
        if not metrics:
            return
        diagnostics.append("    [%s] %s (%d metrics) <- %s"
                           % (route, os.path.basename(chunk_dir),
                              len(metrics), rel))
        rank = (source_priority(fname), len(metrics))
        if rank > best_rank[route]:
            best_rank[route] = rank
            best[route] = {"metrics": metrics, "source": rel}

    for root, _dirs, files in os.walk(chunk_dir):
        for fname in files:
            low = fname.lower()
            if not (low.endswith(".table1") or low.endswith(".log")):
                continue
            fpath = os.path.join(root, fname)
            rel = os.path.relpath(fpath, chunk_dir)
            try:
                with open(fpath, "r", errors="replace") as fh:
                    text = fh.read()
            except (IOError, OSError):
                continue
            if low.endswith(".table1"):
                route = "staraniso" if "staraniso" in low else "truncate"
                consider(route, parse_table1(text), fname, rel)
            else:  # .log
                for block in find_summary_blocks(text):
                    route = ("staraniso" if is_staraniso(fpath, block["title"])
                             else "truncate")
                    consider(route, block["metrics"], fname, rel)
    return best


def match_metric(metrics, must, mustnot):
    for label, values in metrics.items():
        low = label.lower()
        if all(x in low for x in must) and not any(x in low for x in mustnot):
            return values
    return None


def build_row(chunk_name, first, last, block):
    """Return an ordered list of (column, value) for one chunk."""
    row = [
        ("chunk", chunk_name),
        ("image_first", first),
        ("image_last", last),
        ("n_images", last - first + 1),
    ]
    metrics = block["metrics"] if block else {}
    for col, must, mustnot, want_outer in METRICS:
        values = match_metric(metrics, must, mustnot)
        overall = values[0] if values else "N/A"
        row.append((col, overall))
        if want_outer:
            outer = values[2] if values else "N/A"
            row.append((col + "_outer", outer))
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


def write_pdf(path, title, rows, image_summary):
    """Render a multi-page PDF: an image-usage page plus statistics tables.
    Falls back gracefully (no PDF, only a warning) if matplotlib is absent."""
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
        fig = plt.figure(figsize=(8.27, 11.69))  # A4 portrait
        fig.suptitle(title, fontsize=16, y=0.97)
        ax = fig.add_axes([0.06, 0.06, 0.88, 0.84])
        ax.axis("off")
        lines = [
            "Generated: %s" % stamp,
            "Chunks merged: %d" % len(image_summary),
            "Total images used: %d" % sum(n for _, _, n in image_summary),
            "",
            "Image ranges used:",
        ]
        for name, rng, n in image_summary:
            lines.append("    %-24s images %-14s (%d)" % (name, rng, n))
        ax.text(0.0, 1.0, "\n".join(lines), va="top", ha="left",
                family="monospace", fontsize=9)
        pdf.savefig(fig)
        plt.close(fig)

        # ---- Statistics tables (chunks as columns, metrics as rows) ----
        if rows:
            metric_cols = [c for c, _ in rows[0]
                           if c not in ("chunk", "image_first", "image_last", "source_file")]
            chunks_per_page = 6
            for start in range(0, len(rows), chunks_per_page):
                page_rows = rows[start:start + chunks_per_page]
                col_labels = [""] + [dict(r)["chunk"].replace("autoPROC_", "")
                                     for r in page_rows]
                table_data = []
                for metric in metric_cols:
                    line = [metric] + [str(dict(r).get(metric, "N/A")) for r in page_rows]
                    table_data.append(line)
                fig = plt.figure(figsize=(8.27, 11.69))
                fig.suptitle("%s - statistics (%d-%d)"
                             % (title, start + 1, start + len(page_rows)), fontsize=12)
                ax = fig.add_axes([0.03, 0.03, 0.94, 0.9])
                ax.axis("off")
                tbl = ax.table(cellText=table_data, colLabels=col_labels,
                               loc="center", cellLoc="center")
                tbl.auto_set_font_size(False)
                tbl.set_fontsize(7)
                tbl.scale(1, 1.4)
                pdf.savefig(fig)
                plt.close(fig)
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--process-dir", required=True,
                        help="Path to the autoproc_chunks directory.")
    parser.add_argument("--out", default=None,
                        help="Output directory (default: <process-dir>/reports).")
    args = parser.parse_args()

    process_dir = os.path.abspath(args.process_dir)
    out_dir = args.out or os.path.join(process_dir, "reports")
    os.makedirs(out_dir, exist_ok=True)

    chunks = []
    for name in sorted(os.listdir(process_dir)):
        m = CHUNK_RE.match(name)
        full = os.path.join(process_dir, name)
        if m and os.path.isdir(full):
            chunks.append((name, int(m.group(1)), int(m.group(2))))
    chunks.sort(key=lambda c: c[1])

    diagnostics = ["autoPROC report generation - %s"
                   % datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
                   "Process dir: %s" % process_dir,
                   "Chunks found: %d" % len(chunks),
                   "Summary-data blocks located:"]

    truncate_rows, staraniso_rows = [], []
    image_summary = []
    for name, first, last in chunks:
        image_summary.append((name, "%d-%d" % (first, last), last - first + 1))
        best = collect_chunk_blocks(os.path.join(process_dir, name), diagnostics)
        truncate_rows.append(build_row(name, first, last, best["truncate"]))
        staraniso_rows.append(build_row(name, first, last, best["staraniso"]))

    write_csv(os.path.join(out_dir, "truncate_statistics.csv"), truncate_rows)
    write_csv(os.path.join(out_dir, "staraniso_statistics.csv"), staraniso_rows)
    write_pdf(os.path.join(out_dir, "truncate_report.pdf"),
              "Classical autoPROC (TRUNCATE)", truncate_rows, image_summary)
    write_pdf(os.path.join(out_dir, "staraniso_report.pdf"),
              "STARANISO", staraniso_rows, image_summary)

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
if command -v python3 >/dev/null 2>&1; then PY=python3; else PY=python; fi
"\$PY" "$PROCESS_DIR/autoproc_reports.py" --process-dir "$PROCESS_DIR"
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
