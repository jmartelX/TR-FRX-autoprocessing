#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
trfrx_full_pipeline.py — ONE command for a whole TR-FRX series:

    Stage 1  copy the autoPROC_* chunks into <output>/<chunks_name>/   (copy_files.py)
    Stage 2  CAD-clean the chosen MTZ per timepoint into <output>/dfo/ (fetch_clean_mtz.py)
    Stage 3  dimple the first timepoint (MODEL + chosen MTZ); final.pdb -> <output>/dfo/
    Stage 4  Fo-Fo diff maps + SVD + peak analysis, run inside <output>/dfo/

It merges three scripts into one self-contained tool:
    copy_files.py  +  fetch_clean_mtz.py  +  the diffmap/SVD/peak pipeline
(which itself already merged diffmaps.py + SVD_all_in_one.py). Their conventions
are reused verbatim: pathlib.Path throughout, argparse CLI, --dry-run support,
the MTZ trailing-index regexes, and the label / resolution auto-detection helpers.

At startup you are asked whether to use the truncate (default) or the staraniso
MTZ; that choice drives BOTH dimple and the cleaning step (use --file to skip
the prompt).

Difference maps go to   ./output_dfo/
SVD output goes to      ./output_svd/
Peak tables / images /  ./output_dfo/reports/  and  ./output_svd/reports/
per-map PDFs and the    series recap live under the two output folders.

Supported MTZ trailing-index patterns (from the source scripts):
    _N.mtz          e.g. CaMDH_073_137_2.mtz
    _start-end.mtz  e.g. CaMDH_073_275_1-250.mtz
    _start_end.mtz  e.g. CaMDH_073_275_1_300.mtz

Typical use (non-expert friendly — auto-detects ref MTZ, labels, res):
    ./trfrx_full_pipeline.py MODEL.pdb  /path/to/autoPROC_chunks  /path/to/output
    ./trfrx_full_pipeline.py MODEL.pdb  ../chunks  ./out  --file truncate
    ./trfrx_full_pipeline.py MODEL.pdb  ../chunks  ./out  --dry-run

Absolute and relative paths both work. Everything else is auto-detected.

The original single-stage pipeline CLI is preserved for advanced / cluster use
and is selected automatically whenever --dir, --only or --svd-only is passed:
    ./trfrx_full_pipeline.py --dir /path/to/dfo --sigma 3.5 --cluster
"""

from __future__ import annotations

import argparse
import logging
import math
import re
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path

# ═════════════════════════════════════════════════════════════════════════
# Defaults (see the module docstring / README for the radius rationale)
# ═════════════════════════════════════════════════════════════════════════
SIGMA_DEFAULT   = "auto"  # peak cutoff: "auto" (noise-based) or a fixed number of sigma
SIGMA_FLOOR     = 3.0    # auto cutoff never drops below this (community standard)
AUTO_EXPECTED_FALSE = 3.0  # auto cutoff admits ~this many noise peaks per map
                           # (1 = only the single strongest noise peak = strictest)
RADIUS_DEFAULT  = 4.0    # "nearest residue" base search radius, A (H-bond+vdW; scales up)
N_PEAKS_DEFAULT = 20     # peaks IMAGED per map (the table keeps all above the cutoff)
MAX_TABLE_PEAKS = 200    # hard safety cap on tabulated peaks per map
LOW_RES_DEFAULT = 10.0   # low-resolution cutoff for the difference maps

# Merge radius is deliberately MUCH smaller than the residue radius: its only
# job is to collapse the several grid points of a single peak into one row.
# Auto = max(MERGE_RADIUS_MIN, grid_spacing). Override with --merge-radius.
MERGE_RADIUS_MIN = 1.5   # A

# PyMOL images — wwPDB-validation-report look (grey model, tight mesh cage)
IMG_W, IMG_H = 1000, 800
DISPLAY_SIGMA = 3.0     # contour the mesh in the figures at this |sigma| (fixed)
DISPLAY_NEAR  = 5.0     # only residues within this many A of the peak are shown
MESH_CARVE    = 1.8     # A — tight cage hugging the shown atoms (validation look)
# "lone" atoms of interest — waters + common ions — shown as spheres and labelled
# (they have no bonds, so `show sticks` alone would leave them invisible).
_LONE_SEL = "(resn HOH or elem Na+K+Mg+Ca+Mn+Fe+Zn+Cu+Ni+Co+Cd+Cl+Br+I)"
# Ions / single atoms of interest (waters EXCLUDED) — these get a name label.
_ION_SEL  = "(elem Na+K+Mg+Ca+Mn+Fe+Zn+Cu+Ni+Co+Cd+Cl+Br+I)"

# Signal-to-noise: noise floor is estimated from density further than
# SOLVENT_DIST from any atom (bulk solvent); signal from the top-N near peaks.
SOLVENT_DIST = 5.0   # A
SNR_TOP_N    = 5
# Persistence: how strong the σ-vs-time monotonic trend must be to call a site
# "rising"/"decaying" rather than "flickering" (|Spearman rho| threshold).
TREND_RHO = 0.6

# CaMDH active-site residues used by the existing dFo PyMOL scripts — only used
# as a convenience default if you later enable the active-site flag.
DEFAULT_ACTIVE_SITE = [82, 88, 148, 151, 175]


# ═════════════════════════════════════════════════════════════════════════
# § A. Dependency checks  (extends run_series_diffmaps.check_dependencies)
# ═════════════════════════════════════════════════════════════════════════
def check_dependencies(need_images: bool = True) -> bool:
    ok = True

    programs = {
        "mtzdmp":                     "module load ccp4",
        "phenix.fobs_minus_fobs_map": "module load phenix",
    }
    for prog, fix in programs.items():
        if not shutil.which(prog):
            print(f"  MISSING program : {prog}  ->  {fix}")
            ok = False

    # Now REQUIRED for peak finding / merging / reports.
    packages = {
        "numpy":      "pip install numpy",
        "pandas":     "pip install pandas",
        "gemmi":      "pip install gemmi   (or: conda install -c conda-forge gemmi)",
        "scipy":      "pip install scipy   (peak finding + KDTree merge)",
        "matplotlib": "pip install matplotlib   (peak tables + PDF reports)",
    }
    for pkg, fix in packages.items():
        try:
            __import__(pkg)
        except ImportError:
            print(f"  MISSING package  : {pkg}  ->  {fix}")
            ok = False

    if need_images and not check_pymol(verbose=True):
        # Not fatal on its own — user can run with --skip-images. We warn here
        # and let main() decide, but a clear message is printed either way.
        pass

    optional = {
        "dask":    "pip install dask[array]  (SVD uses numpy fallback if absent)",
        "seaborn": "pip install seaborn      (nicer SVD plots)",
    }
    for pkg, note in optional.items():
        try:
            __import__(pkg)
        except ImportError:
            print(f"  optional missing : {pkg}  ->  {note}")

    return ok


def check_pymol(verbose: bool = False) -> bool:
    """True if a usable PyMOL is available (executable or importable module)."""
    if shutil.which("pymol"):
        return True
    try:
        import pymol  # noqa: F401  (module build, e.g. `module load pymol`)
        return True
    except Exception:
        if verbose:
            print("  MISSING program : pymol  ->  module load pymol  "
                  "(or run with --skip-images)")
        return False


# ═════════════════════════════════════════════════════════════════════════
# § B. Reused verbatim from diffmaps.py / SVD_all_in_one.py
# ═════════════════════════════════════════════════════════════════════════
_INDEX_PATTERNS = [
    re.compile(r"_(\d+)-(\d+)\.mtz$"),   # _start-end.mtz
    re.compile(r"_(\d+)_(\d+)\.mtz$"),   # _start_end.mtz
    re.compile(r"_(\d+)\.mtz$"),          # _N.mtz
]


def mtz_sort_index(mtz: Path) -> int | None:
    for pat in _INDEX_PATTERNS:
        m = pat.search(mtz.name)
        if m:
            return int(m.group(1))
    return None


def _run_silent(cmd: list[str]) -> str:
    try:
        cp = subprocess.run(cmd, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return cp.stdout + cp.stderr
    except FileNotFoundError:
        return ""


def detect_fsig_labels(mtz: Path) -> str:
    """Return 'FP,SIGFP' or 'F,SIGF' by inspecting the MTZ header."""
    txt = _run_silent(["mtzdmp", str(mtz)]) if shutil.which("mtzdmp") \
          else _run_silent(["phenix.mtz.dump", str(mtz)])
    words = set(txt.split())
    if "FP" in words and "SIGFP" in words:
        return "FP,SIGFP"
    if "F" in words and "SIGF" in words:
        return "F,SIGF"
    raise RuntimeError(
        f"Could not find FP/SIGFP or F/SIGF in {mtz.name}.\n"
        f"Tokens found (first 60): {' '.join(list(words)[:60])}"
    )


def detect_resolution_limit(mtz: Path) -> float:
    """High-resolution limit — gemmi first, mtzdmp fallback (from the sources)."""
    try:
        import gemmi
        res = gemmi.read_mtz_file(str(mtz)).resolution_high()
        if res and res > 0:
            return round(res, 4)
    except Exception:
        pass

    txt   = _run_silent(["mtzdmp", str(mtz)])
    lines = txt.splitlines()
    for i, line in enumerate(lines):
        if "Resolution Range" in line or "Resolution range" in line:
            for bline in lines[i:i + 4]:
                nums   = re.findall(r"\d+\.\d+", bline)
                floats = [float(x) for x in nums if 1.0 <= float(x) <= 999.0]
                if len(floats) >= 2:
                    return min(floats)
    raise RuntimeError(
        f"Could not detect resolution limit from {mtz.name}.\n"
        f"mtzdmp output (first 40 lines):\n" + "\n".join(lines[:40])
    )


def find_pdb(directory: Path) -> Path:
    pdbs = sorted(directory.glob("*.pdb"))
    if not pdbs:
        raise RuntimeError("No .pdb file found in directory.")
    if len(pdbs) == 1:
        return pdbs[0]
    model = directory / "model.pdb"
    if model.exists():
        return model
    raise RuntimeError(
        f"Multiple PDB files found but model.pdb does not exist.\n"
        f"PDBs present: {', '.join(p.name for p in pdbs)}"
    )


def collect_source_mtz(directory: Path) -> list[Path]:
    """Source MTZ files sorted by trailing index. Excludes dFo_* files."""
    indexed = []
    for mtz in directory.glob("*.mtz"):
        if mtz.name.startswith("dFo_"):
            continue
        idx = mtz_sort_index(mtz)
        if idx is not None:
            indexed.append((idx, mtz))
    if not indexed:
        raise RuntimeError("No MTZ files with a recognised trailing index found.")
    indexed.sort(key=lambda t: t[0])
    return [mtz for _, mtz in indexed]


def compute_diffmap(target, ref, model, labels, reslim, dry_run, out_dir, log_dir=None):
    """phenix.fobs_minus_fobs_map — reused verbatim from SVD_all_in_one.py.

    The diffmap .log is written to *log_dir* (default out_dir) so all logs can be
    gathered in one place; the .eff and the dFo_*.mtz stay in out_dir.
    """
    stem     = target.stem
    out_dir.mkdir(parents=True, exist_ok=True)
    log_dir  = log_dir or out_dir
    log_dir.mkdir(parents=True, exist_ok=True)
    eff_path = out_dir / f"dFo-Fo_{stem}.eff"
    log_path = log_dir / f"{stem}-{ref.stem}_diffmap.log"
    prefix   = out_dir / f"dFo_{stem}"

    eff_content = (
        f"f_obs_1_file_name = {target}\n"
        f"f_obs_1_label = {labels}\n"
        f"f_obs_2_file_name = {ref}\n"
        f"f_obs_2_label = {labels}\n"
        f"high_resolution = {reslim}\n"
        f"low_resolution = {LOW_RES_DEFAULT}\n"
        f"sigma_cutoff = 3.0\n"
        f"phase_source = {model}\n"
        f"ignore_non_isomorphous_unit_cells = True\n"
        f"advanced {{\n"
        f"  multiscale = True\n"
        f"}}\n"
    )

    if dry_run:
        print(f"  [dry-run] would write {eff_path.name} and run phenix")
        return

    eff_path.write_text(eff_content)
    with log_path.open("w") as log:
        try:
            subprocess.run(
                ["phenix.fobs_minus_fobs_map",
                 str(target), str(ref), str(model), str(eff_path),
                 f"file_name_prefix={prefix}", "job_id=1"],
                stdout=log, stderr=log, check=True,
            )
        except subprocess.CalledProcessError:
            raise RuntimeError(f"phenix failed — see {log_path.name}")
        except FileNotFoundError:
            raise RuntimeError("phenix.fobs_minus_fobs_map not found in PATH.")

    print(f"  MTZ: {prefix}_1.mtz")
    print(f"  LOG: {log_path.name}")


# ═════════════════════════════════════════════════════════════════════════
# § C. run_svd — reused from SVD_all_in_one.py, with a return value added
#      (singular values) so the series recap can report SV weights.
#      Existing outputs (maps, rSV.csv, plot) are unchanged.
# ═════════════════════════════════════════════════════════════════════════
def run_svd(dfo_dir: Path, out_dir: Path, time_step_ms: float | None = None):
    """Scan dfo_dir for dFo_*.mtz, run SVD, write maps/csv/plot. Returns [S]."""
    import numpy as np
    import pandas as pd
    import gemmi

    try:
        import matplotlib
        matplotlib.use("Agg")            # headless — never touch X11 (SLURM/SSH safe)
        import matplotlib.pyplot as plt
        HAVE_PLOT = True
    except Exception as _e:
        print(f"  SVD: matplotlib unavailable ({_e}) — no plot.")
        HAVE_PLOT = False
    # seaborn is OPTIONAL (nicer palette only) — its absence must NOT skip the plot
    try:
        import seaborn as sns
        HAVE_SNS = True
    except Exception:
        HAVE_SNS = False

    try:
        from scipy.optimize import curve_fit
        HAVE_SCIPY = True
    except Exception:
        HAVE_SCIPY = False

    PLOT_SV_LIMIT = 7
    DO_FIT_SV0    = True
    FIT_MODEL     = "decay"

    def dataset_number(name):
        ints = re.findall(r"\d+", name)
        return int(ints[-2]) if len(ints) >= 2 else None

    _grid_size = None

    def mtz_to_matrix(path):
        nonlocal _grid_size
        mtz = gemmi.read_mtz_file(str(path))
        if _grid_size is None:
            grid = mtz.transform_f_phi_to_map("FoFo", "PHFc")
            arr  = np.array(grid, copy=False)
            _grid_size = list(arr.shape)
            print(f"  SVD: grid size set to {_grid_size} from {path.name}")
            return arr
        grid = mtz.transform_f_phi_to_map("FoFo", "PHFc", exact_size=_grid_size)
        return np.array(grid, copy=False)

    def matrix_to_mtz(matrix, template, out):
        mtz  = gemmi.read_mtz_file(str(template))
        rlim = mtz.resolution_high()
        eg   = mtz.get_f_phi_on_grid("FoFo", "PHFc", _grid_size)
        rs   = gemmi.transform_f_phi_grid_to_map(eg)
        for z in range(rs.nw):
            for y in range(rs.nv):
                for x in range(rs.nu):
                    rs.set_value(x, y, z, float(matrix[x, y, z]))
        recigrid = gemmi.transform_map_to_f_phi(rs, half_l=True)
        data     = recigrid.prepare_asu_data(dmin=rlim)
        mtz.set_data(data)
        mtz.write_to_file(str(out))

    entries = []
    for p in sorted(dfo_dir.glob("dFo_*.mtz")):
        n = dataset_number(p.stem)
        if n is None:
            print(f"  SVD: skipping {p.name} (cannot parse dataset number)")
            continue
        entries.append((n, p))

    if not entries:
        print("  SVD: no dFo_*.mtz files found — skipping.")
        return None

    step    = time_step_ms
    x_label = f"Time (ms), step={step} ms" if step is not None else "Dataset"
    print(f"  SVD: X axis = {x_label}, {len(entries)} maps")

    df = pd.DataFrame({
        "n":       [e[0] for e in entries],
        "x_num":   [float(e[0] - 1) if step is None else float((e[0] - 1) * step)
                    for e in entries],
        "x_label": [e[1].stem if step is None else f"{(e[0]-1)*step} ms"
                    for e in entries],
        "path":    [e[1] for e in entries],
    }).sort_values("n").reset_index(drop=True)

    out_dir.mkdir(parents=True, exist_ok=True)

    first_map = mtz_to_matrix(df.path[0])
    shape, m  = first_map.shape, int(np.prod(first_map.shape))
    n_files   = len(df)
    A         = np.zeros((m, n_files), dtype=np.float32)
    A[:, 0]   = first_map.reshape(-1)
    for j in range(1, n_files):
        mat = mtz_to_matrix(df.path[j])
        if mat.shape != shape:
            raise ValueError(f"Shape mismatch: {df.path[j].name}")
        A[:, j] = mat.reshape(-1)

    try:
        import dask.array as da
        U, S, VT = da.linalg.svd(da.array(A))
        U = np.array(U); S = np.array(S); VT = np.array(VT)
        print("  SVD: using dask backend")
    except ImportError:
        print("  SVD: dask not found, using numpy.linalg.svd")
        U, S, VT = np.linalg.svd(A, full_matrices=False)

    sigmat = np.zeros((n_files, n_files), dtype=np.float32)
    for i in range(len(S)):
        sigmat[i, i] = S[i]
    US = np.array(np.matrix(U) @ np.matrix(sigmat))

    template = df.path[0]
    for rank in range(n_files):
        matrix_to_mtz(U[:, rank].reshape(shape),  template, out_dir / f"lSV_{rank}.mtz")
        matrix_to_mtz(US[:, rank].reshape(shape), template, out_dir / f"SVscaled_lSV_{rank}.mtz")

    if n_files >= 2:
        Ar1 = np.array(US @ VT.T[1, :])
        matrix_to_mtz(Ar1.reshape(shape), df.path[1], out_dir / f"reformed_{df.path[1].stem}.mtz")

    t_index   = np.array(df.x_num, dtype=float)
    x_labels  = list(df.x_label)
    scaled_tf = (sigmat @ VT)[:n_files, :]
    index_label = "time_ms" if step is not None else "dataset"
    pd.DataFrame({f"rSV{i}": scaled_tf[i] for i in range(n_files)},
                 index=x_labels).to_csv(out_dir / "rSV.csv", index_label=index_label)

    # Singular values -> file (additive; nothing existing depends on it) + return
    sv_list = [float(s) for s in S]
    s_total = float(S.sum()) or 1.0
    pd.DataFrame({
        "component": list(range(len(sv_list))),
        "singular_value": sv_list,
        "fraction": [s / s_total for s in sv_list],
    }).to_csv(out_dir / "singular_values.csv", index=False)

    print("  SVD: singular values (relative weight):")
    for i, s in enumerate(S):
        bar = "#" * int(30 * s / s_total)
        print(f"    SV{i}: {s:.4g}  ({100*s/s_total:.1f}%)  {bar}")
    print(f"  SVD: wrote lSV/SVscaled maps and rSV.csv -> {out_dir}")

    if not HAVE_PLOT:
        print("  SVD: matplotlib not available — skipping plots.")
        return sv_list

    def monoexp(t, A_, tau, C):
        return A_ * np.exp(-t / tau) + C if FIT_MODEL == "decay" \
               else A_ * (1 - np.exp(-t / tau)) + C

    svlim   = min(PLOT_SV_LIMIT, n_files)
    palette = ([sns.color_palette("bright", n_colors=svlim)[i] for i in range(svlim)]
               if HAVE_SNS else [plt.cm.tab10(i / 10) for i in range(svlim)])

    order    = np.argsort(t_index)
    t        = t_index[order]
    t_labels = [x_labels[i] for i in order]

    fig, ax = plt.subplots(figsize=(10, 6))
    for i in range(svlim):
        yi = np.array(scaled_tf[i], dtype=float)[order]
        ax.plot(range(len(t)), yi, "o--", color=palette[i],
                markersize=8, linewidth=2, alpha=0.7, label=f"SV{i}")
        if i == 0 and DO_FIT_SV0 and HAVE_SCIPY and step is not None and len(t) >= 4:
            try:
                C0   = float(np.mean(yi[-max(1, len(yi) // 10):]))
                A0   = float(yi[0] - C0) if FIT_MODEL == "decay" else float(yi[-1] - C0)
                tau0 = float(max((t.max() - t.min()) / 3.0, 1.0))
                popt, _ = curve_fit(monoexp, t, yi, p0=[A0, tau0, C0],
                                    bounds=([-np.inf, 1e-6, -np.inf],
                                            [np.inf, np.inf, np.inf]), maxfev=20000)
                tfit = np.linspace(float(t.min()), float(t.max()), 300)
                xfit = np.interp(tfit, t, range(len(t)))
                ax.plot(xfit, monoexp(tfit, *popt), "-", color=palette[0],
                        linewidth=2.5, alpha=0.95, label=f"SV0 fit (tau={popt[1]:.2f} ms)")
                with open(out_dir / "SV0_monoexp_fit.txt", "w") as f:
                    f.write(f"model {FIT_MODEL}\nA {popt[0]}\ntau_ms {popt[1]}\nC {popt[2]}\n")
                print(f"  SVD: SV0 fit -> A={popt[0]:.4g}, tau={popt[1]:.4g} ms, C={popt[2]:.4g}")
            except Exception as e:
                print(f"  SVD: SV0 fit failed ({e})")

    ax.set_xticks(range(len(t_labels)))
    ax.set_xticklabels(t_labels, rotation=45, ha="right", fontsize=8)
    ax.set_xlabel(x_label); ax.set_ylabel("Magnitude")
    ax.set_title("Right Singular Vectors"); ax.legend(); ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / "rSV_plot.pdf")            # PDF (print quality)
    fig.savefig(out_dir / "rSV_plot.png", dpi=150)   # PNG (for the HTML report)
    plt.close(fig)
    print(f"  SVD: plot saved -> {out_dir / 'rSV_plot.pdf'} (+ .png)")
    return sv_list


# ═════════════════════════════════════════════════════════════════════════
# § D. Feature 1 — interactive resolution handling
# ═════════════════════════════════════════════════════════════════════════
def resolution_table(mtz_files: list[Path]) -> dict[str, float]:
    res_values: dict[str, float] = {}
    for mtz in mtz_files:
        try:
            res_values[mtz.name] = detect_resolution_limit(mtz)
        except RuntimeError as e:
            print(f"WARNING: could not detect resolution for {mtz.name}: {e}")
    return res_values


def print_resolution_table(res_values: dict[str, float]) -> None:
    print("\nHigh-resolution limit per MTZ (best -> worst):")
    print(f"  {'res (A)':>8}   file")
    for name, r in sorted(res_values.items(), key=lambda x: x[1]):
        print(f"  {r:>8.3f}   {name}")


def prompt_resolution_choice(mtz_files: list[Path],
                             res_values: dict[str, float],
                             interactive: bool) -> tuple[float, list[Path]]:
    """
    Return (reslim, kept_files).

    Default (Enter / non-interactive / cluster): worst common resolution across
    all files — the conservative choice required for SVD grid consistency.
    Two prompts at most: (1) cutoff to use, (2) optionally drop worse files.
    """
    worst = max(res_values.values())
    worst_name = max(res_values, key=res_values.get)

    if not interactive:
        print(f"Resolution (worst common, auto): {worst:.3f} A  <- {worst_name}")
        return worst, mtz_files

    print_resolution_table(res_values)
    print(f"\nWorst (common) resolution = {worst:.3f} A  ({worst_name})")

    # Prompt 1 — cutoff
    try:
        ans = input(f"High-resolution cutoff for ALL maps in A "
                    f"[Enter = {worst:.3f}]: ").strip()
    except EOFError:
        ans = ""
    reslim = float(ans) if ans else worst

    # Prompt 2 — drop worse files
    droppable = {n: r for n, r in res_values.items() if r > reslim + 1e-6}
    kept = list(mtz_files)
    if droppable:
        print("These files are worse than the chosen cutoff:")
        for n, r in sorted(droppable.items(), key=lambda x: -x[1]):
            print(f"    {r:.3f} A  {n}")
        try:
            drop = input("Drop them before proceeding? [y/N]: ").strip().lower()
        except EOFError:
            drop = ""
        if drop.startswith("y"):
            kept = [m for m in mtz_files if m.name not in droppable]
            print(f"  Dropped {len(droppable)} file(s); {len(kept)} remain.")
    return reslim, kept


# ═════════════════════════════════════════════════════════════════════════
# § E. Feature 2 — peak detection, merging, residue assignment
# ═════════════════════════════════════════════════════════════════════════
@dataclass
class Peak:
    x: float
    y: float
    z: float
    sigma: float                 # signed (positive or negative lobe)
    volume: float = 0.0          # blob volume above the contour, A^3
    resname: str = "?"
    chain: str = "?"
    seqid: int = 0
    dist: float = float("nan")   # symmetry-correct nearest-atom distance, A
    short_label: str = "?"
    atom: str = "?"
    near: bool = False           # within the base chemical radius of the model
    search_radius: float = 0.0   # volume-scaled radius used to find the residue
    ncs_count: int = 1           # in how many NCS molecules this site is seen
    ncs_of: int = 1              # number of NCS copies of that molecule
    ncs_chains: str = ""         # which chains show it, e.g. "A,B,C"

    @property
    def ncs_unique(self) -> bool:
        """Event seen in only one molecule out of several NCS copies."""
        return self.ncs_of > 1 and self.ncs_count == 1

    @property
    def sign(self) -> str:
        return "+" if self.sigma >= 0 else "-"

    @property
    def blob_radius(self) -> float:
        """Effective radius (A) of a sphere with this peak's volume."""
        return (3.0 * self.volume / (4.0 * math.pi)) ** (1.0 / 3.0) if self.volume > 0 else 0.0


def peak_search_radius(volume: float, base_radius: float) -> float:
    """Volume-scaled nearest-residue radius: blob surface + chemical margin."""
    r_blob = (3.0 * volume / (4.0 * math.pi)) ** (1.0 / 3.0) if volume > 0 else 0.0
    return min(max(r_blob + base_radius, base_radius), 3.0 * base_radius)


def peak_image_carve(volume: float) -> float:
    """Volume-scaled isomesh carve / zoom radius for the figure (A)."""
    r_blob = (3.0 * volume / (4.0 * math.pi)) ** (1.0 / 3.0) if volume > 0 else 0.0
    return min(max(r_blob + 1.5, 2.5), 8.0)


# Finer grid than the default: smoother mesh + better-localised peaks.
_MAP_SAMPLE_RATE = 3.0


def load_grid_from_mtz(mtz_path: Path):
    """Return (grid, rms, voxel_volume_A3) for a FoFo/PHFc MTZ (dFo or SVD map)."""
    import numpy as np
    import gemmi
    mtz  = gemmi.read_mtz_file(str(mtz_path))
    grid = mtz.transform_f_phi_to_map("FoFo", "PHFc", sample_rate=_MAP_SAMPLE_RATE)
    arr  = np.array(grid, copy=False)
    rms  = float(arr.std()) or 1.0
    vox  = grid.unit_cell.volume / (grid.nu * grid.nv * grid.nw)
    return grid, rms, vox


def write_ccp4(grid, out_path: Path) -> None:
    """Dump a gemmi grid to a CCP4 .ccp4 map for PyMOL isomesh."""
    import gemmi
    ccp4 = gemmi.Ccp4Map()
    ccp4.grid = grid
    ccp4.update_ccp4_header()
    ccp4.write_ccp4_map(str(out_path))


def _grid_position(grid, u: int, v: int, w: int):
    import gemmi
    frac = gemmi.Fractional(u / grid.nu, v / grid.nv, w / grid.nw)
    return grid.unit_cell.orthogonalize(frac)


def find_peaks(grid, rms: float, sigma: float, voxel_vol: float = 1.0) -> list[Peak]:
    """
    One peak per connected blob beyond the contour (|value| >= sigma*rms), placed
    at the blob's extremum and tagged with the blob volume (A^3). Using connected
    components — rather than a local-maximum filter — means a single feature yields
    a single row even when its top is broad/flat.
    """
    import numpy as np
    arr = np.asarray(grid, dtype=np.float32)
    thr = sigma * rms

    try:
        from scipy.ndimage import label, maximum_position, minimum_position
    except Exception:
        # crude fallback: threshold voxels only (scipy is a hard dependency though)
        peaks = []
        for mask, _ in ((arr >= thr, True), (arr <= -thr, False)):
            us, vs, ws = np.nonzero(mask)
            for u, v, w in zip(us.tolist(), vs.tolist(), ws.tolist()):
                p = _grid_position(grid, u, v, w)
                peaks.append(Peak(p.x, p.y, p.z, float(arr[u, v, w]) / rms))
        return peaks

    peaks: list[Peak] = []
    for thr_mask, want_pos in ((arr >= thr, True), (arr <= -thr, False)):
        comp, ncomp = label(thr_mask)
        if ncomp == 0:
            continue
        sizes = np.bincount(comp.ravel())
        labels = list(range(1, ncomp + 1))
        extrema = (maximum_position if want_pos else minimum_position)(arr, comp, labels)
        for lbl_id, (u, v, w) in zip(labels, extrema):
            p = _grid_position(grid, u, v, w)
            vol = float(sizes[lbl_id] * voxel_vol)
            peaks.append(Peak(p.x, p.y, p.z, float(arr[u, v, w]) / rms, volume=vol))
    return peaks


def _canonical_asu_xyz(peaks: list[Peak], st):
    """
    Fold each peak into a canonical position in the asymmetric unit by applying
    all space-group operations and picking the lexicographically-smallest wrapped
    fractional coordinate. Symmetry copies of one peak collapse to the same point,
    so they can be deduplicated. Returns a list of orthogonal (x,y,z) tuples.
    """
    import gemmi
    sg  = st.find_spacegroup()
    ops = list(sg.operations()) if sg else [gemmi.Op("x,y,z")]
    out = []
    for p in peaks:
        fr  = st.cell.fractionalize(gemmi.Position(p.x, p.y, p.z))
        best = None
        for op in ops:
            t = op.apply_to_xyz([fr.x, fr.y, fr.z])
            key = tuple(round(c % 1.0, 4) for c in t)
            if best is None or key < best:
                best = key
        orth = st.cell.orthogonalize(gemmi.Fractional(*best))
        out.append((orth.x, orth.y, orth.z))
    return out


def merge_peaks(peaks: list[Peak], merge_radius: float, st=None) -> list[Peak]:
    """
    Deduplicate peaks — WITHIN sign only (a +/- lobe pair is a real dipole and is
    never merged). When *st* is given, peaks are first folded into the asymmetric
    unit so that the many symmetry copies of a single feature collapse to one row
    (crucial for high-symmetry space groups). Keeps the highest-|sigma| member of
    each cluster. Direct KDTree union-find, no peakmax.
    """
    import numpy as np
    if not peaks:
        return []

    # Position used for clustering: ASU-folded when we have symmetry, else raw.
    if st is not None:
        canon = _canonical_asu_xyz(peaks, st)
    else:
        canon = [(p.x, p.y, p.z) for p in peaks]

    out: list[Peak] = []
    for want_pos in (True, False):
        idx = [i for i, p in enumerate(peaks) if (p.sigma >= 0) == want_pos]
        if not idx:
            continue
        coords = np.array([canon[i] for i in idx])
        for members in _cluster_indices(coords, merge_radius):
            grp = [peaks[idx[m]] for m in members]
            out.append(max(grp, key=lambda p: abs(p.sigma)))
    out.sort(key=lambda p: abs(p.sigma), reverse=True)
    return out


def one_letter(resname: str) -> str:
    import gemmi
    try:
        info = gemmi.find_tabulated_residue(resname)
        if info and info.is_amino_acid():
            return info.one_letter_code.upper()
    except Exception:
        pass
    return resname


def load_structure(model_path: Path):
    import gemmi
    st = gemmi.read_structure(str(model_path))
    st.setup_entities()
    return st


WATER_RESNAMES = {"HOH", "WAT", "DOD", "H2O", "SOL", "TIP", "TIP3"}


def assign_nearest_residue(peaks: list[Peak], st, base_radius: float,
                           ignore_hydrogens: bool = True,
                           ignore_waters: bool = True) -> None:
    """
    Fill each peak's nearest-residue fields via gemmi.NeighborSearch.

    Fixes vs. the first version:
      * distance is measured to the actual symmetry image (mark.pos), NOT the
        asymmetric-unit copy (which gave nonsense 50-120 A distances);
      * the search radius is volume-scaled per peak (bigger peaks reach further);
      * hydrogens are ignored (report the nearest heavy atom);
      * waters are ignored (assign to the nearest protein residue / ligand / ion,
        not to a water) — waters are still shown in the figures for context;
      * a `near` flag marks peaks within the base chemical radius of the model.
    """
    import gemmi
    if not peaks:
        return
    model = st[0]
    max_r = 3.0 * base_radius                       # NeighborSearch upper bound
    ns = gemmi.NeighborSearch(model, st.cell, max_r).populate()
    for pk in peaks:
        pk.search_radius = peak_search_radius(pk.volume, base_radius)
        pos = gemmi.Position(pk.x, pk.y, pk.z)
        best, best_d = None, 1e9
        for mark in ns.find_atoms(pos, "\0", radius=pk.search_radius):
            cra = mark.to_cra(model)
            if ignore_hydrogens and cra.atom.is_hydrogen():
                continue
            if ignore_waters and cra.residue.name in WATER_RESNAMES:
                continue
            d = pos.dist(mark.pos)                   # symmetry-image distance
            if d < best_d:
                best_d, best = d, cra
        if best is not None:
            pk.resname = best.residue.name
            pk.chain   = best.chain.name
            pk.seqid   = best.residue.seqid.num
            pk.atom    = best.atom.name
            pk.dist    = round(best_d, 2)
            pk.near    = best_d <= base_radius
            pk.short_label = f"{one_letter(best.residue.name)}{best.residue.seqid.num}"
        else:
            pk.resname = "(none)"
            pk.short_label = "-"
            pk.near = False


NCS_RMSD_MAX = 3.0   # A — chains superposing below this are treated as NCS copies


def _chain_ca(chain):
    """{seqid_num: CA Position} for a polymer chain."""
    d = {}
    for res in chain.get_polymer():
        a = res.find_atom("CA", "*")
        if a is not None:
            d[res.seqid.num] = a.pos
    return d


def detect_ncs(st) -> dict:
    """
    Group the model's chains into NCS families by superposition and return
    {chain_name: (group_id, transform_to_ref, group_size)}. The transform maps
    that chain's coordinates into its family's reference-chain frame, so peaks in
    different molecules can be compared in a common frame. Chains that don't
    superpose on anything (single copies, ligand-only chains) get group_size 1.
    """
    import gemmi
    model = st[0]
    cas = {}
    for ch in model:
        d = _chain_ca(ch)
        if len(d) >= 5:                       # need a real polymer
            cas[ch.name] = d
    names = list(cas)
    out = {n: (i, gemmi.Transform(), 1) for i, n in enumerate(names)}  # identity default
    assigned = set()
    gid = 0
    for i, ref in enumerate(names):
        if ref in assigned:
            continue
        family = [ref]
        assigned.add(ref)
        out[ref] = (gid, gemmi.Transform(), 1)          # ref -> identity
        for other in names[i + 1:]:
            if other in assigned:
                continue
            common = sorted(set(cas[ref]) & set(cas[other]))
            if len(common) < 5:
                continue
            try:
                sup = gemmi.superpose_positions([cas[ref][s] for s in common],
                                                [cas[other][s] for s in common])
            except Exception:
                continue
            if sup.rmsd <= NCS_RMSD_MAX:
                family.append(other)
                assigned.add(other)
                out[other] = (gid, sup.transform, 0)    # other -> ref frame
        size = len(family)
        for n in family:
            g, t, _ = out[n]
            out[n] = (g, t, size)
        gid += 1
    return out


def annotate_ncs(peaks: list[Peak], st, radius: float) -> None:
    """
    Tag each peak with how many NCS molecules show the same event. Peaks are
    mapped into their family's reference frame and clustered there (within sign),
    so the same site across copies groups together. Sets ncs_count / ncs_of /
    ncs_chains; peaks seen in a single copy of a multi-copy molecule are flagged
    via Peak.ncs_unique.
    """
    import numpy as np
    import gemmi
    ncs = detect_ncs(st)
    if not any(g[2] > 1 for g in ncs.values()):
        return                                          # no NCS — nothing to do

    ref_xyz = []
    for pk in peaks:
        pos = gemmi.Position(pk.x, pk.y, pk.z)
        info = ncs.get(pk.chain)
        rp = info[1].apply(pos) if info else pos
        ref_xyz.append([rp.x, rp.y, rp.z])

    for want_pos in (True, False):
        idx = [i for i, p in enumerate(peaks) if (p.sigma >= 0) == want_pos]
        if not idx:
            continue
        coords = np.array([ref_xyz[i] for i in idx])
        for members in _cluster_indices(coords, radius):
            grp = [peaks[idx[m]] for m in members]
            chains = sorted({p.chain for p in grp if p.chain not in ("?", "")})
            size = max((ncs.get(p.chain, (0, None, 1))[2] for p in grp), default=1)
            for p in grp:
                p.ncs_count = len(chains) if chains else 1
                p.ncs_of = size
                p.ncs_chains = ",".join(chains)


def auto_sigma_threshold(grid, rms: float, st, floor: float = SIGMA_FLOOR,
                         expected_false: float = AUTO_EXPECTED_FALSE,
                         solvent_dist: float = SOLVENT_DIST) -> float:
    """
    Data-driven peak cutoff from extreme-value statistics of the noise.

    The signal-free bulk-solvent region (> solvent_dist from any atom) gives the
    noise standard deviation. We set the cutoff so that only ~expected_false noise
    peaks are expected in the whole map: sigma * sqrt(2 ln(N_independent /
    expected_false)). expected_false = 1 gives the single strongest noise peak
    (strictest, ~5 sigma); a few (default 3) relaxes it by ~0.3-0.5 sigma so
    genuine mid-sigma features survive. Robust (uses the noise std, not one hot
    voxel) and adapts per map. Never below `floor`.
    """
    import numpy as np
    import gemmi
    arr  = np.asarray(grid, dtype=np.float32)
    mask = gemmi.FloatGrid(grid.nu, grid.nv, grid.nw)
    mask.set_unit_cell(grid.unit_cell)
    mask.spacegroup = grid.spacegroup
    for cra in st[0].all():
        mask.set_points_around(cra.atom.pos, radius=solvent_dist, value=1.0)
    marr   = np.asarray(mask, dtype=np.float32)
    region = arr[marr == 0]
    if region.size == 0 or rms <= 0:
        return floor
    noise_std = float(region.std())
    # independent points ~ voxels / oversampling^3 (grid is ~sample_rate x Nyquist)
    n_ind   = max(region.size / (_MAP_SAMPLE_RATE ** 3), 10.0)
    z       = math.sqrt(2.0 * math.log(max(n_ind / expected_false, math.e)))
    return round(max(floor, (noise_std / rms) * z), 2)


def resolve_sigma_cutoff(sigma_arg, grid, rms: float, st) -> float:
    """Turn the --sigma argument ('auto' or a number) into a concrete cutoff."""
    if isinstance(sigma_arg, str) and sigma_arg.strip().lower() == "auto":
        return auto_sigma_threshold(grid, rms, st)
    return float(sigma_arg)


# ── Map-quality checklist helpers (diffmap-stage) ────────────────────────
def cell_isomorphism(mtz_path: Path, ref_path: Path) -> dict:
    """Unit-cell comparison of a dataset vs the reference (isomorphism check).

    Returns the deviation AND the actual cell values (A1), so the recap can show
    real numbers, not just a worst-case percentage.
    """
    import gemmi
    c1 = gemmi.read_mtz_file(str(mtz_path)).cell
    c0 = gemmi.read_mtz_file(str(ref_path)).cell
    edge = max(abs(c1.a - c0.a) / c0.a, abs(c1.b - c0.b) / c0.b,
               abs(c1.c - c0.c) / c0.c) * 100.0
    ang  = max(abs(c1.alpha - c0.alpha), abs(c1.beta - c0.beta),
               abs(c1.gamma - c0.gamma))
    fmt = lambda c: (round(c.a, 2), round(c.b, 2), round(c.c, 2),
                     round(c.alpha, 2), round(c.beta, 2), round(c.gamma, 2))
    return {
        "max_edge_pct": round(edge, 3), "max_angle_deg": round(ang, 3),
        "cell": fmt(c1), "ref_cell": fmt(c0),
    }


def cell_str(cell) -> str:
    """'a, b, c, α, β, γ' for a (a,b,c,al,be,ga) tuple, or '-' if missing."""
    if not cell:
        return "-"
    a, b, c, al, be, ga = cell
    return f"{a:.1f} {b:.1f} {c:.1f} / {al:.1f} {be:.1f} {ga:.1f}"


def isomorphism_status(iso: dict) -> str:
    e, a = iso["max_edge_pct"], iso["max_angle_deg"]
    if e < 0.5 and a < 0.2:
        return "OK"
    if e < 1.0 and a < 0.5:
        return "WARN"
    return "FAIL"


_FOBS_BIN_RE = re.compile(
    r"^\s*\d+:\s+[\d.]+\s*-\s*[\d.]+\s+[\d.]+\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s*$")


def parse_fobs_agreement(log_path: Path):
    """Dataset agreement (scaling quality) between the two Fobs of a difference map.

    Parses the 'Fobs1_vs_Fobs2 statistics' table from the phenix diffmap log and
    returns reflection-weighted CC and R over all bins, plus the high-resolution
    shell (worst ~10% of reflections, where agreement degrades). This replaces the
    old regex that mistakenly grabbed the model r_work.
    """
    if not log_path.is_file():
        return None
    lines = log_path.read_text(errors="ignore").splitlines()
    try:
        start = next(i for i, ln in enumerate(lines)
                     if ln.strip().startswith("Fobs1_vs_Fobs2"))
    except StopIteration:
        return None

    rows = []  # (n_refl, CC, R) low-res -> high-res
    for ln in lines[start:]:
        m = _FOBS_BIN_RE.match(ln)
        if m:
            rows.append((int(m.group(1)), float(m.group(2)), float(m.group(3))))
    if not rows:
        return None

    def wmean(sel, idx):
        n = sum(r[0] for r in sel) or 1
        return sum(r[0] * r[idx] for r in sel) / n

    ntot = sum(r[0] for r in rows)
    # High-res shell: accumulate the last bins until ~10% of reflections.
    target, acc, hi = 0.1 * ntot, 0, []
    for r in reversed(rows):
        hi.append(r)
        acc += r[0]
        if acc >= target:
            break
    return {
        "cc_overall": round(wmean(rows, 1), 3), "r_overall": round(wmean(rows, 2), 4),
        "cc_high":    round(wmean(hi, 1), 3),   "r_high":    round(wmean(hi, 2), 4),
        "n_refl": ntot, "n_bins": len(rows),
    }


def write_peak_table(peaks: list[Peak], csv_path: Path, txt_path: Path,
                     map_name: str) -> None:
    import pandas as pd
    rows = []
    for rank, p in enumerate(peaks, 1):
        rows.append({
            "rank": rank, "sigma": round(p.sigma, 2), "sign": p.sign,
            "volume_A3": round(p.volume, 1),
            "x": round(p.x, 2), "y": round(p.y, 2), "z": round(p.z, 2),
            "nearest_res": p.resname, "chain": p.chain, "seqid": p.seqid,
            "atom": p.atom, "short": p.short_label,
            "dist_A": p.dist, "near": p.near,
            "ncs": f"{p.ncs_count}/{p.ncs_of}", "ncs_chains": p.ncs_chains,
            "ncs_unique": p.ncs_unique,
            "search_radius_A": round(p.search_radius, 1),
        })
    df = pd.DataFrame(rows)
    df.to_csv(csv_path, index=False)

    lines = [f"Peak table — {map_name}", "=" * 76,
             f"{'#':>3} {'sigma':>7} {'vol(A3)':>8} {'short':>7} "
             f"{'atom':>5} {'dist':>6} {'near':>5} {'ncs':>6}"]
    for r in rows:
        note = " *unique*" if r["ncs_unique"] else ""
        lines.append(f"{r['rank']:>3} {r['sigma']:>7.2f} {r['volume_A3']:>8.1f} "
                     f"{r['short']:>7} {r['atom']:>5} "
                     f"{('-' if r['dist_A'] != r['dist_A'] else format(r['dist_A'],'.2f')):>6} "
                     f"{'yes' if r['near'] else 'no':>5} {r['ncs']:>6}{note}")
    txt_path.write_text("\n".join(lines) + "\n")


# ═════════════════════════════════════════════════════════════════════════
# § E-bis. Extra "is this map interesting?" analyses
#   1) signal-to-noise heuristic   2) active-site proximity flag
#   (persistence across the series is computed later, over all maps)
# ═════════════════════════════════════════════════════════════════════════
def compute_snr(grid, rms: float, st, peaks: list[Peak],
                solvent_dist: float = SOLVENT_DIST, top_n: int = SNR_TOP_N) -> dict:
    """
    Signal-to-noise heuristic for one map.

    noise floor = RMS of density further than *solvent_dist* from any atom
    (bulk-solvent region), expressed in units of the map sigma.
    signal      = mean of the top-N near-model peak |sigma|.
    snr         = signal / noise floor.  High signal + low far-from-model
                  density  ->  the map is worth a look.
    """
    import numpy as np
    import gemmi
    arr  = np.asarray(grid, dtype=np.float32)
    mask = gemmi.FloatGrid(grid.nu, grid.nv, grid.nw)
    mask.set_unit_cell(grid.unit_cell)
    mask.spacegroup = grid.spacegroup
    for cra in st[0].all():
        mask.set_points_around(cra.atom.pos, radius=solvent_dist, value=1.0)
    marr = np.asarray(mask, dtype=np.float32)

    far = arr[marr == 0]
    far_rms = float(far.std()) if far.size else float("nan")
    noise = (far_rms / rms) if (rms and far_rms == far_rms) else float("nan")

    near = sorted((abs(p.sigma) for p in peaks if p.dist == p.dist), reverse=True)
    signal = float(np.mean(near[:top_n])) if near else 0.0
    snr = (signal / noise) if (noise and noise == noise) else float("nan")
    return {"snr": snr, "noise_floor": noise, "signal": signal}


def active_site_hits(peaks: list[Peak], st, active_site: list[int] | None,
                     radius: float, sigma_cut: float) -> list[tuple[Peak, float]]:
    """Peaks with |sigma| >= sigma_cut within *radius* of any active-site residue."""
    if not active_site:
        return []
    import gemmi
    aset = set(active_site)
    positions = [cra.atom.pos for cra in st[0].all()
                 if cra.residue.seqid.num in aset]
    if not positions:
        return []
    hits = []
    for p in peaks:
        if abs(p.sigma) < sigma_cut:
            continue
        pos = gemmi.Position(p.x, p.y, p.z)
        d = min(pos.dist(q) for q in positions)
        if d <= radius:
            hits.append((p, round(d, 2)))
    return hits


def _cluster_indices(coords, radius: float) -> list[list[int]]:
    """Union-find spatial clustering of points within *radius*. Returns groups."""
    import numpy as np
    n = len(coords)
    parent = list(range(n))
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    try:
        from scipy.spatial import cKDTree
        pairs = cKDTree(coords).query_pairs(radius)
    except Exception:
        pairs = set()
        for i in range(n):
            for j in range(i + 1, n):
                if np.linalg.norm(coords[i] - coords[j]) <= radius:
                    pairs.add((i, j))
    for i, j in pairs:
        parent[find(i)] = find(j)
    clusters: dict[int, list[int]] = {}
    for i in range(n):
        clusters.setdefault(find(i), []).append(i)
    return list(clusters.values())


def _classify_trend(orders: list[int], seq: list[float]) -> str:
    """Rising / decaying / flickering from the σ-vs-time monotonic trend."""
    import numpy as np
    if len(seq) < 3:
        return "-"
    try:
        from scipy.stats import spearmanr
        rho = spearmanr(orders, seq).correlation
    except Exception:
        diffs = np.diff(seq)
        up, dn = int((diffs > 0).sum()), int((diffs < 0).sum())
        rho = (up - dn) / max(len(diffs), 1)
    if rho is None or rho != rho:
        return "flat"
    if rho >= TREND_RHO:
        return f"rising (rho={rho:+.2f})"
    if rho <= -TREND_RHO:
        return f"decaying (rho={rho:+.2f})"
    return f"flickering (rho={rho:+.2f})"


def find_persistent_sites(dfo_results: list[dict], radius: float) -> list[dict]:
    """
    Cluster every map's peaks across the (time-ordered) series and score each
    site: in how many timepoints it appears, its persistence fraction, and the
    monotonic trend of its |sigma| over time. Monotonic sites = real kinetics;
    flickering = noise.  dfo_results must be in timepoint order.
    """
    import numpy as np
    total = len(dfo_results)
    pts, ref = [], []
    for order, r in enumerate(dfo_results):
        for p in r["peaks"]:
            pts.append([p.x, p.y, p.z])
            ref.append((order, p))
    if not pts:
        return []

    sites = []
    for members in _cluster_indices(np.array(pts), radius):
        by_order: dict[int, list[float]] = {}
        for idx in members:
            order, p = ref[idx]
            by_order.setdefault(order, []).append(abs(p.sigma))
        orders = sorted(by_order)
        n_tp   = len(orders)
        if n_tp < 2:
            continue          # seen in only one map -> not a "persistent" site
        seq    = [max(by_order[o]) for o in orders]
        best_i = max(members, key=lambda i: abs(ref[i][1].sigma))
        sites.append({
            "short":       ref[best_i][1].short_label,
            "n_tp":        n_tp,
            "persistence": n_tp / total if total else 0.0,
            "trend":       _classify_trend(orders, seq),
            "max_sigma":   max(abs(ref[i][1].sigma) for i in members),
        })
    sites.sort(key=lambda s: (-s["n_tp"], -s["max_sigma"]))
    return sites


# ═════════════════════════════════════════════════════════════════════════
# § F. Feature 3 — PyMOL peak images + per-map PDF + series recap
# ═════════════════════════════════════════════════════════════════════════
def build_pymol_script(model_path: Path, ccp4_path: Path, peaks: list[Peak],
                       rms: float, display_sigma: float, png_paths: list[Path],
                       display_near: float = DISPLAY_NEAR) -> str:
    """
    Generate a single PyMOL script that renders one image per peak.

    Style = wwPDB validation report: model as sticks with GREY carbons and
    heteroatoms coloured by element (O red, N blue, S yellow, P orange), a thin
    mesh cage carved tightly around the atoms, ray-traced on white with depth
    cueing. Only the map colours are ours (flashy_blue + / flashy_orange -).

      * the mesh is contoured at a FIXED |display_sigma| (default 3.5), decoupled
        from each peak's detected height;
      * carved tightly (~MESH_CARVE A) around the shown atoms + peak for the
        validation "cage" look, not a large ball;
      * only residues within DISPLAY_NEAR A of the peak are shown (uncluttered);
      * waters are shown as red spheres and labelled if they're a contact;
      * symexp builds symmetry neighbours; hydrogens removed; centred on the peak.
    """
    # IMPORTANT: PyMOL normalises CCP4 maps to sigma by default, which made the
    # old code (level = display_sigma * rms) contour at ~0.1 sigma -> a solid
    # blob. We turn normalisation OFF and contour in the map's ABSOLUTE units,
    # where the correct level is display_sigma * rms (rms = grid std, from us).
    body = [
        "from pymol import cmd, util",
        "cmd.set('normalize_ccp4_maps', 0)",     # contour in absolute units (see note)
        "cmd.bg_color('white')",
        "cmd.set('ray_opaque_background', 1)",
        "cmd.set_color('flashy_blue',   [0.10, 0.45, 1.0])",
        "cmd.set_color('flashy_orange', [1.0, 0.40, 0.0])",
        f"cmd.load({str(model_path)!r}, 'enzyme')",
        f"cmd.load({str(ccp4_path)!r}, 'diffmap')",
        "cmd.hide('everything')",
        "cmd.remove('hydro')",
        "cmd.set('ray_shadows', 0)",
        "cmd.set('antialias', 2)",
        "cmd.set('depth_cue', 0)",               # off: difference cages read flatter/clearer
        "cmd.set('two_sided_lighting', 1)",
        "cmd.set('mesh_width', 0.5)",
        "cmd.set('mesh_quality', 3)",
        "cmd.set('stick_radius', 0.16)",
        "cmd.set('valence', 0)",
        "cmd.set('nonbonded_size', 0.5)",
        "cmd.set('sphere_scale', 0.30)",
        "cmd.set('label_size', 16)",
        "cmd.set('label_color', 'black')",
        "cmd.set('label_font_id', 7)",           # bold sans-serif
        "cmd.set('label_outline_color', 'white')",
        "cmd.set('float_labels', 1)",
    ]

    for i, pk in enumerate(peaks):
        level    = display_sigma * rms                   # absolute contour level
        level_i  = level if pk.sigma >= 0 else -level
        color    = "flashy_blue" if pk.sigma >= 0 else "flashy_orange"
        carve_r  = peak_image_carve(pk.volume)           # cage radius around the peak
        sym_cut  = display_near + carve_r + 2.0
        png      = png_paths[i]
        near_sel = f"byres ((enzyme or sym*) within {display_near:.1f} of pk)"

        body += [
            "cmd.delete('sym*'); cmd.delete('pk'); cmd.delete('mesh'); cmd.delete('near')",
            f"cmd.pseudoatom('pk', pos=[{pk.x:.3f}, {pk.y:.3f}, {pk.z:.3f}])",
            "cmd.hide('everything')",
            # symmetry neighbours around the peak (needs CRYST1 in the PDB)
            f"try:\n    cmd.symexp('sym', 'enzyme', 'pk', {sym_cut:.1f})\nexcept Exception:\n    pass",
            f"cmd.select('near', {near_sel!r})",
            f"cmd.select('lone', 'near and {_LONE_SEL}')",    # waters + ions (single atoms)
            "cmd.show('sticks', 'near')",                     # bonded atoms (residues, ligands)
            "cmd.show('nb_spheres', 'lone')",                 # lone atoms as spheres
            "cmd.color('grey70', 'near')",                    # grey carbons ...
            "cmd.util.cnc('near')",                            # ... heteroatoms by element
            "cmd.color('red', 'near and resn HOH')",          # waters red
            # difference-density cage: a sphere of radius carve_r AROUND THE PEAK,
            # so the blob is always shown (not clipped to nearby atoms).
            f"cmd.isomesh('mesh', 'diffmap', {level_i:.4f}, 'pk', carve={carve_r:.2f})",
            f"cmd.color({color!r}, 'mesh')",
            # a small marker at the peak centre so the eye finds it immediately
            "cmd.show('nonbonded', 'pk')",
            "cmd.color('yellow', 'pk')",
            # fresh labels each frame: residue names (protein Cα, one per residue)
            # and ion / single-atom names — never waters.
            "cmd.label('all', '')",
            "cmd.label('near and polymer and name CA', '\"%s%s\" % (resn, resi)')",
            f"cmd.label('near and {_ION_SEL}', '\"%s%s\" % (resn, resi)')",
            # consistent framing: orient on the local scene, then fit around the peak
            "cmd.orient('near or pk')",
            f"cmd.zoom('pk', {carve_r + 2.5:.2f})",
            f"cmd.png({str(png)!r}, width={IMG_W}, height={IMG_H}, dpi=150, ray=1)",
        ]
    return "\n".join(body) + "\n"


def render_peak_images(model_path: Path, ccp4_path: Path, peaks: list[Peak],
                       rms: float, display_sigma: float, img_dir: Path,
                       display_near: float = DISPLAY_NEAR) -> list[Path]:
    """Run PyMOL once to render every peak image. Returns existing PNG paths."""
    img_dir.mkdir(parents=True, exist_ok=True)
    png_paths = [img_dir / f"peak_{i+1:02d}.png" for i in range(len(peaks))]
    if not peaks:
        return []
    script = build_pymol_script(model_path, ccp4_path, peaks, rms, display_sigma,
                                png_paths, display_near=display_near)
    script_path = img_dir / "_render_peaks.py"   # .py so multiline python runs
    script_path.write_text(script)

    exe = shutil.which("pymol")
    cmd = [exe or "pymol", "-cq", str(script_path)]
    try:
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                       check=True)
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"    peak images: PyMOL render failed ({e}). "
              f"Try 'module load pymol' or --skip-images.")
        return []
    return [p for p in png_paths if p.exists()]


def build_map_pdf(map_name: str, peaks: list[Peak], images: list[Path],
                  pdf_path: Path, display_sigma: float = DISPLAY_SIGMA) -> None:
    """One PDF per map: peak table first, then peak images (2-up grid)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.image as mpimg

    ROWS_PER_PAGE = 28
    with PdfPages(pdf_path) as pdf:
        # Table pages — the full list (paginated), not just the imaged peaks.
        if not peaks:
            fig, ax = plt.subplots(figsize=(8.27, 11.69))
            ax.axis("off")
            ax.set_title(f"Peak table — {map_name}", fontsize=13, weight="bold", pad=14)
            ax.text(0.5, 0.5, "No peaks above cutoff.", ha="center")
            pdf.savefig(fig); plt.close(fig)
        for start in range(0, len(peaks), ROWS_PER_PAGE):
            chunk = peaks[start:start + ROWS_PER_PAGE]
            fig, ax = plt.subplots(figsize=(8.27, 11.69))  # A4 portrait
            ax.axis("off")
            suffix = (f"  (rows {start+1}-{start+len(chunk)} of {len(peaks)})"
                      if len(peaks) > ROWS_PER_PAGE else f"  ({len(peaks)} peaks)")
            ax.set_title(f"Peak table — {map_name}{suffix}",
                         fontsize=13, weight="bold", pad=14)
            cell_text = [[f"{start+i+1}", f"{p.sigma:+.2f}", f"{p.volume:.0f}",
                          p.short_label, f"{p.resname} {p.seqid}", p.atom,
                          "-" if p.dist != p.dist else f"{p.dist:.2f}",
                          "yes" if p.near else "no",
                          f"{p.ncs_count}/{p.ncs_of}" + ("*" if p.ncs_unique else "")]
                         for i, p in enumerate(chunk)]
            table = ax.table(cellText=cell_text,
                             colLabels=["#", "sigma", "vol A3", "short",
                                        "residue", "atom", "dist A", "near", "ncs"],
                             loc="upper center", cellLoc="center")
            table.auto_set_font_size(False); table.set_fontsize(9)
            table.scale(1, 1.4)
            pdf.savefig(fig); plt.close(fig)

        # Following pages — images, 2 columns x 3 rows per page
        per_page = 6
        for start in range(0, len(images), per_page):
            chunk = images[start:start + per_page]
            fig, axes = plt.subplots(3, 2, figsize=(8.27, 11.69))
            axes = axes.ravel()
            for ax in axes:
                ax.axis("off")
            for k, img in enumerate(chunk):
                idx = start + k
                axes[k].imshow(mpimg.imread(img))
                pk = peaks[idx]
                dtxt = "-" if pk.dist != pk.dist else f"{pk.dist:.1f} A"
                tag  = "" if pk.near else " [far]"
                axes[k].set_title(
                    f"#{idx+1}  {pk.short_label} ({dtxt}){tag}\n"
                    f"detected {pk.sigma:+.1f}$\\sigma$ · shown at "
                    f"{display_sigma:.1f}$\\sigma$", fontsize=9)
            pdf.savefig(fig); plt.close(fig)
    print(f"    report: {pdf_path}")


def build_series_recap(dfo_results: list[dict], sv_values, recap_path: Path,
                       radius: float) -> None:
    """
    One higher-level file across the whole series so you don't open every PDF.
    Contains: best peak + S/N per timepoint, SV weights, persistent-site table
    (recurrence + rising/decaying/flickering trend), and active-site hits.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    import pandas as pd

    def _r(x, n):
        """Round for display; None/NaN -> None (blank in CSV)."""
        try:
            v = float(x)
            return round(v, n) if v == v else None
        except (TypeError, ValueError):
            return None

    # Per-timepoint summary CSV — plain-language headers, rounded (A3/A4)
    csv_path = recap_path.with_suffix(".csv")
    pd.DataFrame([{
        "timepoint": r["name"],
        "peaks found": r["n_peaks"],
        "strongest peak (sigma)": _r(r["best_sigma"], 2),
        "strongest peak near": r["best_res"],
        "signal / noise": _r(r.get("snr"), 2),
        "noise level (sigma)": _r(r.get("noise_floor"), 3),
        "active-site hit?": ("yes" if r.get("active_flag") else "no"),
        "active-site peaks": ";".join(r.get("active_hits", [])),
    } for r in dfo_results]).to_csv(csv_path, index=False)

    # Persistence across the series (time-ordered) + CSV
    sites = find_persistent_sites(dfo_results, radius)
    if sites:
        pd.DataFrame(sites).to_csv(
            recap_path.with_name(recap_path.stem + "_persistence.csv"), index=False)

    # Map-quality checklist CSV — cell values (A1), dataset agreement CC/R (A2)
    have_quality = any("iso" in r for r in dfo_results)
    if have_quality:
        pd.DataFrame([{
            "timepoint": r["name"],
            "resolution (A)": _r(r.get("resolution"), 2),
            "unit cell (a b c / al be ga)": cell_str(r.get("iso", {}).get("cell")),
            "unit-cell change (%)": _r(r.get("iso", {}).get("max_edge_pct"), 2),
            "angle change (deg)": _r(r.get("iso", {}).get("max_angle_deg"), 2),
            "same crystal form?": r.get("iso_status"),
            "dataset agreement CC": _r((r.get("scaling") or {}).get("cc_overall"), 3),
            "dataset agreement R": _r((r.get("scaling") or {}).get("r_overall"), 3),
            "high-res agreement R": _r((r.get("scaling") or {}).get("r_high"), 3),
            "peak threshold (sigma)": _r(r.get("cutoff"), 2),
        } for r in dfo_results]).to_csv(
            recap_path.with_name("map_quality.csv"), index=False)

    any_active = any(r.get("active_flag") for r in dfo_results)

    with PdfPages(recap_path) as pdf:
        # ── Page 1: map-quality checklist ────────────────────────────────
        if have_quality:
            fig, ax = plt.subplots(figsize=(8.27, 11.69))
            ax.axis("off")
            ax.set_title("Map-quality checklist", fontsize=15, weight="bold", pad=16)
            cell_text = []
            for r in dfo_results:
                iso = r.get("iso", {})
                sc  = r.get("scaling") or {}
                cell_text.append([
                    r["name"][:20],
                    f"{r['resolution']:.2f}" if r.get("resolution") else "-",
                    cell_str(iso.get("cell")),
                    f"{iso.get('max_edge_pct', float('nan')):.2f}",
                    r.get("iso_status", "-"),
                    f"{sc['cc_overall']:.3f}" if sc.get("cc_overall") is not None else "-",
                    f"{sc['r_overall']:.3f}" if sc.get("r_overall") is not None else "-",
                    f"{r.get('cutoff', float('nan')):.2f}",
                ])
            tbl = ax.table(cellText=cell_text,
                           colLabels=["timepoint", "res (Å)", "unit cell (a b c/α β γ)",
                                      "cell Δ%", "same form?", "agree CC", "agree R",
                                      "peak σ"],
                           loc="upper center", cellLoc="center")
            tbl.auto_set_font_size(False); tbl.set_fontsize(6.5); tbl.scale(1, 1.4)
            # colour the "same form?" cell by status
            status_col = 4
            for ridx, r in enumerate(dfo_results, start=1):
                s = r.get("iso_status", "-")
                col = {"OK": "#c8f7c5", "WARN": "#ffe9a8", "FAIL": "#f7bcbc"}.get(s)
                if col:
                    tbl[(ridx, status_col)].set_facecolor(col)
            ax.text(0.05, 0.02,
                    "cell Δ% = largest unit-cell edge change vs the reference dataset;  "
                    "same crystal form: OK <0.5% / WARN <1% / FAIL otherwise.  "
                    "agree CC/R = how well the two datasets of each difference map agree.",
                    fontsize=6.5, style="italic", transform=ax.transAxes)
            pdf.savefig(fig); plt.close(fig)

        # ── Page 2: peak / SVD summary ───────────────────────────────────
        fig, ax = plt.subplots(figsize=(8.27, 11.69))
        ax.axis("off")
        ax.set_title("TR-FRX series recap", fontsize=15, weight="bold", pad=16)

        y = 0.94
        ax.text(0.05, y, "Best peak per timepoint  (map: best-sigma @ res | S/N)",
                fontsize=12, weight="bold"); y -= 0.028
        for r in dfo_results:
            mark = " *AS*" if r.get("active_flag") else ""
            snr = r.get("snr", float("nan"))
            ax.text(0.06, y, f"{r['name'][:34]:<36} {r['best_sigma']:+5.1f}s @ "
                             f"{r['best_res']:<6} | S/N {snr:4.1f} "
                             f"({r['n_peaks']}pk){mark}",
                    fontsize=7.5, family="monospace"); y -= 0.020

        y -= 0.02
        ax.text(0.05, y, "Singular-value weights", fontsize=12, weight="bold")
        y -= 0.028
        if sv_values:
            tot = sum(sv_values) or 1.0
            for i, s in enumerate(sv_values[:8]):
                bar = "#" * int(30 * s / tot)
                ax.text(0.06, y, f"SV{i}: {100*s/tot:5.1f}%  {bar}",
                        fontsize=8, family="monospace"); y -= 0.020

        y -= 0.02
        ax.text(0.05, y, "Persistent sites  (site | timepoints | trend | max|sigma|)",
                fontsize=12, weight="bold"); y -= 0.028
        for s in sites[:12]:
            ax.text(0.06, y, f"{s['short']:<7} {s['n_tp']:>2}/{len(dfo_results)} "
                             f"({100*s['persistence']:3.0f}%)  "
                             f"{s['trend']:<20} max={s['max_sigma']:.1f}",
                    fontsize=7.5, family="monospace"); y -= 0.020

        if any_active:
            y -= 0.02
            ax.text(0.05, y, "Active-site hits", fontsize=12, weight="bold")
            y -= 0.028
            for r in dfo_results:
                if r.get("active_flag"):
                    ax.text(0.06, y, f"{r['name'][:34]:<36} "
                                     f"{', '.join(r['active_hits'])}",
                            fontsize=7.5, family="monospace"); y -= 0.020

        if any(r.get("ncs_unique") for r in dfo_results):
            y -= 0.02
            ax.text(0.05, y, "NCS-unique events  (seen in only one molecule)",
                    fontsize=12, weight="bold"); y -= 0.028
            for r in dfo_results:
                nu = r.get("ncs_unique", 0)
                if nu:
                    ax.text(0.06, y, f"{r['name'][:34]:<36} {nu} peak(s) in a single copy",
                            fontsize=7.5, family="monospace"); y -= 0.020
        pdf.savefig(fig); plt.close(fig)
    print(f"  recap: {recap_path}  (+ {csv_path.name})")


# ═════════════════════════════════════════════════════════════════════════
# § G. Per-map analysis orchestration (dFo maps AND SVD component maps)
# ═════════════════════════════════════════════════════════════════════════
def analyze_map(map_mtz: Path, st, model_path: Path, sigma,
                radius: float, merge_radius: float, n_peaks: int,
                reports_dir: Path, make_images: bool,
                active_site: list[int] | None = None,
                display_sigma: float = DISPLAY_SIGMA,
                display_near: float = DISPLAY_NEAR) -> dict:
    """Full peak analysis for one FoFo map. Returns a recap dict."""
    name = map_mtz.stem
    grid, rms, voxel_vol = load_grid_from_mtz(map_mtz)

    # Per-map significance cutoff: 'auto' -> noise-based, else the given number.
    cutoff = resolve_sigma_cutoff(sigma, grid, rms, st)

    peaks = find_peaks(grid, rms, cutoff, voxel_vol)
    peaks = merge_peaks(peaks, merge_radius, st)      # folds symmetry copies
    peaks = peaks[:MAX_TABLE_PEAKS]                    # safety cap only
    assign_nearest_residue(peaks, st, radius)         # symmetry-correct + volume-scaled
    annotate_ncs(peaks, st, max(merge_radius, 2.0))   # group same event across NCS copies

    tag = "auto" if (isinstance(sigma, str) and sigma.lower() == "auto") else "fixed"
    print(f"  peaks: {name}  cutoff={cutoff:.2f} sigma ({tag})  "
          f"{len(peaks)} peaks, imaging top {min(n_peaks, len(peaks))}")

    # Extra heuristics
    snr_info = compute_snr(grid, rms, st, peaks)
    hits     = active_site_hits(peaks, st, active_site, radius, cutoff)

    map_report = reports_dir / name
    map_report.mkdir(parents=True, exist_ok=True)
    write_peak_table(peaks, map_report / "peaks.csv",
                     map_report / "peaks.txt", name)     # ALL peaks above cutoff

    # Figures: only the strongest n_peaks (the table above keeps everything).
    image_peaks = peaks[:n_peaks]
    images: list[Path] = []
    if make_images and image_peaks:
        ccp4 = map_report / f"{name}.ccp4"
        write_ccp4(grid, ccp4)
        images = render_peak_images(model_path, ccp4, image_peaks, rms,
                                    display_sigma, map_report / "img",
                                    display_near=display_near)

    build_map_pdf(name, peaks, images, reports_dir / f"{name}_peaks.pdf",
                  display_sigma=display_sigma)

    flag = "  ***ACTIVE-SITE***" if hits else ""
    print(f"    S/N={snr_info['snr']:.1f} (signal={snr_info['signal']:.1f}sigma, "
          f"noise={snr_info['noise_floor']:.2f}sigma){flag}")

    best = peaks[0] if peaks else None
    ncs_unique = sum(1 for p in peaks if p.ncs_unique)
    return {
        "name": name, "n_peaks": len(peaks), "peaks": peaks,
        "best_sigma": best.sigma if best else 0.0,
        "best_res": best.short_label if best else "-",
        "snr": snr_info["snr"], "noise_floor": snr_info["noise_floor"],
        "active_flag": bool(hits),
        "active_hits": [f"{p.short_label}@{d}A" for p, d in hits],
        "cutoff": cutoff, "ncs_unique": ncs_unique,
    }


# ═════════════════════════════════════════════════════════════════════════
# § H. Feature 4 — SLURM cluster submission (optional mode only)
# ═════════════════════════════════════════════════════════════════════════
def submit_cluster(work_dir: Path, mtz_files: list[Path], ref: Path,
                   reslim: float, args) -> int:
    """
    Submit the series to SLURM, reusing the TR-FRX_autoPROC_cbf.sh pattern:
    one job per source MTZ (diffmap + peak analysis), then a final SVD+recap
    job chained with --dependency=afterok. Non-interactive: reslim is fixed.
    """
    jobs_dir = work_dir / "pipeline_jobs"
    jobs_dir.mkdir(exist_ok=True)
    script = Path(__file__).resolve()
    py     = sys.executable

    common = (f"--dir {work_dir} --ref {ref} --sigma {args.sigma} "
              f"--radius {args.radius} --n-peaks {args.n_peaks} "
              f"--high-res {reslim}")
    if args.skip_images:
        common += " --skip-images"
    if args.active_site:
        common += " --active-site " + " ".join(str(r) for r in args.active_site)

    part = SLURM_PARTITION
    dep_ids = []
    for mtz in mtz_files:
        if mtz == ref:
            continue
        js = jobs_dir / f"peaks_{mtz.stem}.sh"
        js.write_text(
            "#!/bin/bash\nset -euo pipefail\n"
            "module load phenix ccp4 pymol 2>/dev/null || true\n"
            "unset PYTHONHOME PYTHONPATH DISPLAY || true\n"  # no X11 in batch jobs
            "export MPLBACKEND=Agg\n"
            f"cd {work_dir}\n"
            f"{py} {script} {common} --only {mtz.name}\n"
        )
        if args.dry_run:
            print(f"  [dry-run] sbatch {js}")
            continue
        out = subprocess.run(["sbatch", "-p", part, "-n", "1", "-c", "8",
                              "--mem", "16000", str(js)],
                             text=True, stdout=subprocess.PIPE).stdout
        jid = out.strip().split()[-1] if out.strip() else None
        if jid:
            dep_ids.append(jid)
            print(f"  submitted {js.name} -> job {jid}")

    # Final SVD + recap job, depends on all diffmap jobs
    final = jobs_dir / "svd_recap.sh"
    final.write_text(
        "#!/bin/bash\nset -euo pipefail\n"
        "module load phenix ccp4 pymol 2>/dev/null || true\n"
        f"cd {work_dir}\n"
        f"{py} {script} {common} --svd-only\n"
    )
    if args.dry_run:
        print(f"  [dry-run] sbatch --dependency=afterok:<all> {final}")
        return 0
    dep = f"--dependency=afterok:{':'.join(dep_ids)}" if dep_ids else ""
    cmd = ["sbatch", "-p", part, "-n", "1", "-c", "8", "--mem", "24000"]
    if dep:
        cmd.append(dep)
    cmd.append(str(final))
    subprocess.run(cmd, check=False)
    print(f"  submitted {final.name} (afterok chain). Jobs dir: {jobs_dir}")
    return 0


# ═════════════════════════════════════════════════════════════════════════
# § I. main
# ═════════════════════════════════════════════════════════════════════════
def _legacy_pipeline_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="TR-FRX pipeline: Fo-Fo diff maps + SVD + peak analysis.",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)

    # Preserved from the source scripts
    parser.add_argument("--ref",       type=Path, metavar="MTZ",
                        help="Reference MTZ (default: smallest trailing index).")
    parser.add_argument("--model",     type=Path, metavar="PDB",
                        help="PDB model (default: auto-detect).")
    parser.add_argument("--dir",       type=Path, default=Path("."), metavar="DIR",
                        help="Working directory (default: current).")
    parser.add_argument("--high-res",  type=float, default=None, metavar="A",
                        help="Fixed high-res cutoff for ALL maps "
                             "(default: interactive / worst-common).")
    parser.add_argument("--time-step", type=float, default=None, metavar="MS",
                        help="Time interval in ms between datasets (SVD X-axis).")
    parser.add_argument("--skip-svd",  action="store_true", help="Diffmaps only.")
    parser.add_argument("--dry-run",   action="store_true",
                        help="Print the plan without running anything.")

    # New — kept small on purpose
    parser.add_argument("--sigma",   default=SIGMA_DEFAULT, metavar="N|auto",
                        help="Peak cutoff: 'auto' (noise-based, per map — default) "
                             "or a fixed number of sigma, e.g. --sigma 4.")
    parser.add_argument("--radius",  type=float, default=RADIUS_DEFAULT, metavar="A",
                        help=f"Base nearest-residue search radius (default {RADIUS_DEFAULT}); "
                             "scales up with peak volume.")
    parser.add_argument("--n-peaks", type=int, default=N_PEAKS_DEFAULT, metavar="N",
                        help=f"How many peaks to IMAGE per map (default {N_PEAKS_DEFAULT}); "
                             "the peaks.csv table keeps all peaks above the cutoff.")
    parser.add_argument("--display-sigma", type=float, default=DISPLAY_SIGMA, metavar="N",
                        help=f"Contour level (|sigma|) for the mesh in the figures "
                             f"(default {DISPLAY_SIGMA}); independent of peak detection.")
    parser.add_argument("--display-near", type=float, default=DISPLAY_NEAR, metavar="A",
                        help=f"Show residues within this many A of the peak in each "
                             f"figure (default {DISPLAY_NEAR}).")
    parser.add_argument("--cluster", action="store_true",
                        help="Submit to SLURM instead of running locally.")
    parser.add_argument("--skip-images", action="store_true",
                        help="Skip PyMOL peak images (tables + PDFs still made).")
    parser.add_argument("--active-site", type=int, nargs="+", default=None,
                        metavar="RESNUM",
                        help="Residue numbers of interest (e.g. --active-site "
                             "82 88 148 151 175). Flags maps with a >=sigma peak "
                             "within --radius of them. Off unless specified.")

    # Advanced / internal (used by cluster jobs; kept out of the way)
    parser.add_argument("--merge-radius", type=float, default=None,
                        help=argparse.SUPPRESS)
    parser.add_argument("--only", metavar="MTZ", default=None,
                        help=argparse.SUPPRESS)          # process one source MTZ
    parser.add_argument("--svd-only", action="store_true",
                        help=argparse.SUPPRESS)          # SVD + SVD-map peaks + recap
    return parser


def run_pipeline(args) -> int:
    work_dir = args.dir.resolve()
    out_dfo  = work_dir / "output_dfo"
    out_svd  = work_dir / "output_svd"
    # Where phenix diffmap logs go (gathered logs/ when the orchestrator sets it).
    logs_dir = Path(args.logs_dir).resolve() if getattr(args, "logs_dir", None) else out_dfo
    make_images = not args.skip_images

    try:
        mtz_files = collect_source_mtz(work_dir)
        model     = args.model.resolve() if args.model else find_pdb(work_dir)
    except RuntimeError as e:
        print(f"ERROR: {e}")
        return 1

    ref = args.ref.resolve() if args.ref else mtz_files[0]
    if not ref.is_file():
        print(f"ERROR: Reference MTZ not found: {ref}")
        return 1

    print("Checking dependencies...")
    if not check_dependencies(need_images=make_images):
        print("\nERROR: missing required dependencies listed above. Aborting.")
        return 1
    if make_images and not check_pymol():
        print("NOTE: PyMOL not found — continuing with --skip-images behaviour.")
        make_images = False
    print("Dependencies OK.\n")

    print(f"Working directory : {work_dir}")
    print(f"Reference MTZ     : {ref.name}")
    print(f"Model PDB         : {model.name}")
    print(f"MTZ files ({len(mtz_files)}):")
    for mtz in mtz_files:
        print(f"  {mtz.name}{' <- reference' if mtz == ref else ''}")
    if args.dry_run:
        print("DRY RUN -- nothing will be written or executed.")
    print()

    try:
        labels = detect_fsig_labels(ref)
    except RuntimeError as e:
        print(f"ERROR: {e}")
        return 1
    print(f"Labels detected: {labels}\n")

    # ── SVD-only mode (cluster final job) ────────────────────────────────
    if args.svd_only:
        rc = _run_svd_and_peaks(out_dfo, out_svd, model, args, make_images)
        # This is the last job of a --cluster run: build the HTML report now
        # that all per-timepoint jobs have finished.
        if not getattr(args, "no_html", False):
            try:
                print(f"\nHTML report -> {build_report_for_analysis(work_dir)}")
            except Exception as e:
                print(f"  HTML report skipped ({e})")
        return rc

    # ── Resolution (Feature 1) ───────────────────────────────────────────
    if args.high_res is not None:
        reslim = args.high_res
        print(f"Resolution (user-specified): {reslim} A")
    else:
        res_values = resolution_table(mtz_files)
        if not res_values:
            print("ERROR: could not detect resolution for any MTZ. Use --high-res.")
            return 1
        # Prompt on any interactive terminal — including when SUBMITTING to the
        # cluster (the chosen resolution is then baked into every job via
        # --high-res). The per-map jobs themselves pass --high-res, so they never
        # reach this branch and never block on input.
        interactive = sys.stdin.isatty() and args.only is None
        reslim, mtz_files = prompt_resolution_choice(mtz_files, res_values, interactive)
        if ref not in mtz_files and args.only is None:
            mtz_files = [ref] + mtz_files       # keep reference available
    print(f"  -> All difference maps use {reslim} A.\n")

    merge_radius = args.merge_radius if args.merge_radius else max(MERGE_RADIUS_MIN, 0.0)

    # ── Cluster mode (optional) ──────────────────────────────────────────
    if args.cluster:
        return submit_cluster(work_dir, mtz_files, ref, reslim, args)

    if not args.dry_run:
        out_dfo.mkdir(exist_ok=True)

    st = None if args.dry_run else load_structure(model)
    reports_dfo = out_dfo / "reports"
    dfo_results: list[dict] = []

    # ── Diffmaps + per-map peak analysis ─────────────────────────────────
    targets = [args.only] if args.only else None
    ok = skipped = failed = 0
    for mtz in mtz_files:
        if mtz == ref:
            print(f"Skipping reference: {mtz.name}")
            continue
        if targets and mtz.name not in targets:
            continue

        existing = list(out_dfo.glob(f"dFo_{mtz.stem}_*.mtz"))
        print(f"\n=== {mtz.name}  minus  {ref.name} ===")
        try:
            if existing:
                print(f"  diffmap exists, reusing: {existing[0].name}")
                skipped += 1
            else:
                compute_diffmap(mtz, ref, model, labels, reslim, args.dry_run,
                                out_dfo, log_dir=logs_dir)
                ok += 1
            if not args.dry_run:
                produced = out_dfo / f"dFo_{mtz.stem}_1.mtz"
                if produced.exists():
                    result = analyze_map(produced, st, model, args.sigma, args.radius,
                                         merge_radius, args.n_peaks, reports_dfo,
                                         make_images, active_site=args.active_site,
                                         display_sigma=args.display_sigma,
                                         display_near=getattr(args, "display_near", DISPLAY_NEAR))
                    # Attach map-quality info for the recap checklist
                    try:
                        iso = cell_isomorphism(mtz, ref)
                        result["iso"] = iso
                        result["iso_status"] = isomorphism_status(iso)
                        result["resolution"] = detect_resolution_limit(mtz)
                        result["scaling"] = parse_fobs_agreement(
                            logs_dir / f"{mtz.stem}-{ref.stem}_diffmap.log")
                    except Exception as e:
                        print(f"  (quality check skipped: {e})")
                    dfo_results.append(result)
        # E2: one bad timepoint must not abort the series (or block the report).
        except Exception as e:
            print(f"  ERROR ({mtz.name}): {e}")
            failed += 1
        print("-----------------------------------")

    print(f"\n  Computed: {ok}  Reused: {skipped}  Failed: {failed}\n")

    if args.dry_run:
        print("[dry-run] would run SVD, SVD-map peak analysis, and build recap.")
        return 0

    # ── SVD + SVD-map peaks + recap ──────────────────────────────────────
    if args.only:
        return 0 if failed == 0 else 2   # cluster per-map job stops here
    if args.skip_svd:
        if dfo_results:
            build_series_recap(dfo_results, None,
                               out_dfo / "series_recap.pdf", args.radius)
        return 0 if failed == 0 else 2

    rc = _run_svd_and_peaks(out_dfo, out_svd, model, args, make_images,
                            st=st, dfo_results=dfo_results)
    return rc if failed == 0 else 2


def _run_svd_and_peaks(out_dfo: Path, out_svd: Path, model: Path, args,
                       make_images: bool, st=None, dfo_results=None) -> int:
    """Run SVD, analyse the SVD component maps, then build the series recap."""
    print("Running SVD...")
    try:
        sv_values = run_svd(out_dfo, out_svd, args.time_step)
    except Exception as e:
        print(f"SVD ERROR: {e}")
        return 2

    if st is None:
        st = load_structure(model)
    merge_radius = args.merge_radius if args.merge_radius else MERGE_RADIUS_MIN

    # NOTE: peak finding / imaging is intentionally NOT run on the SVD component
    # maps (lSV_*). Peaks are analysed only on the Fo-Fo difference maps.

    # Recap is built from the per-timepoint dFo results; on the svd-only cluster
    # path those aren't in memory, so rebuild them from the existing dFo maps.
    if dfo_results is None:
        dfo_results = _recap_from_existing(out_dfo, st, model, args, merge_radius)
    if dfo_results:
        build_series_recap(dfo_results, sv_values,
                           out_svd / "series_recap.pdf", args.radius)
    return 0


def _recap_from_existing(out_dfo: Path, st, model: Path, args,
                         merge_radius: float) -> list[dict]:
    """Light re-scan of existing dFo maps for the recap (svd-only path)."""
    results = []
    reports_dfo = out_dfo / "reports"
    for dfo in sorted(out_dfo.glob("dFo_*_1.mtz")):
        try:
            results.append(
                analyze_map(dfo, st, model, args.sigma, args.radius,
                            merge_radius, args.n_peaks, reports_dfo, False,
                            active_site=args.active_site,
                            display_sigma=args.display_sigma))
        except Exception as e:
            print(f"  recap: skipping {dfo.name} ({e})")
    return results


# ═════════════════════════════════════════════════════════════════════════
# § J. Stages 1-3 — copy, clean, dimple
#      Merged from copy_files.py and fetch_clean_mtz.py; conventions kept
#      verbatim. These run locally, before the § I analysis pipeline.
# ═════════════════════════════════════════════════════════════════════════
SUBFOLDER_PREFIX = "autoPROC_"

# Files copied verbatim from every autoPROC_* chunk (copy_files.py default set).
DEFAULT_COPY_FILES = [
    "staraniso_alldata-unique.mtz",
    "summary.html",
    "summary.tar.gz",
    "truncate-unique.mtz",
    "XDS_ASCII.HKL",
    "truncate.log",
]

# The two structure-factor files the user can choose between; the choice drives
# BOTH which file dimple refines against and which file is cleaned for the maps.
TARGET_FILES = {
    "staraniso": "staraniso_alldata-unique.mtz",
    "truncate":  "truncate-unique.mtz",
}

RANGE_DIR_RE = re.compile(r"^(\d+)[\s_-]+(\d+)$")   # 1-300, 1_300, 1 - 300
CAD_EXE      = "cad"
MTZDUMP_EXE  = "mtzdmp"
DIMPLE_EXE   = "dimple"
SLURM_PARTITION = "nice"   # SLURM partition used by --cluster (dimple srun + diffmap sbatch)


@dataclass
class CopyResult:
    subfolder: str
    copied: list = field(default_factory=list)
    missing: list = field(default_factory=list)
    skipped: list = field(default_factory=list)


def _process_subfolder(subfolder: Path, dest_root: Path, files_to_copy: list,
                       dry_run: bool, resume: bool = False) -> CopyResult:
    """Copy the selected files from one autoPROC_* chunk (from copy_files.py).

    autoPROC_1_300 -> <dest_root>/1_300/autoPROC_1_300/<files>. With *resume*,
    a destination file that already exists is left untouched.
    """
    result = CopyResult(subfolder=subfolder.name)
    name = subfolder.name
    range_name = name[len(SUBFOLDER_PREFIX):] if name.startswith(SUBFOLDER_PREFIX) else name
    dest_dir = dest_root / range_name / subfolder.name
    if not dry_run:
        dest_dir.mkdir(parents=True, exist_ok=True)
    for filename in files_to_copy:
        src = subfolder / filename
        if src.is_file():
            dst = dest_dir / filename
            if resume and dst.is_file():
                result.skipped.append(filename)
            elif not dry_run:
                shutil.copy2(src, dst)
                result.copied.append(filename)
            else:
                result.copied.append(filename)
        else:
            result.missing.append(filename)
    return result


def copy_chunks(source: Path, dest_root: Path, files_to_copy: list,
                dry_run: bool, workers: int, resume: bool = False) -> Path:
    """Copy every autoPROC_* chunk under *source* into *dest_root* and return it."""
    print("\n── Stage 1: copy chunks ──────────────────────────")
    print(f"  Source      : {source}")
    print(f"  Destination : {dest_root}")
    if dry_run:
        print("  (dry-run — nothing written)")
    else:
        dest_root.mkdir(parents=True, exist_ok=True)

    subfolders = sorted(p for p in source.iterdir()
                        if p.is_dir() and p.name.startswith(SUBFOLDER_PREFIX))
    if not subfolders:
        raise RuntimeError(f"No '{SUBFOLDER_PREFIX}*' subfolders found in {source}")

    results, failures = [], []
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(_process_subfolder, sf, dest_root, files_to_copy,
                             dry_run, resume): sf
                   for sf in subfolders}
        for fut in as_completed(futures):
            sf = futures[fut]
            try:
                results.append(fut.result())
            except Exception as exc:
                failures.append(f"{sf.name}: {exc}")
    copied = sum(len(r.copied) for r in results)
    missing = sum(len(r.missing) for r in results)
    skipped = sum(len(r.skipped) for r in results)
    extra = f"   already present: {skipped}" if resume else ""
    print(f"  Chunks: {len(results)}   files copied: {copied}   missing: {missing}{extra}")
    for f in failures:
        print(f"  FAILED {f}")
    return dest_root


def _natural_int(s: str, default: int = 10 ** 12) -> int:
    m = re.search(r"\d+", s)
    return int(m.group(0)) if m else default


def iter_range_subfolders(root: Path) -> list:
    """Range subfolders (1-300, 1_300) sorted by start index (fetch_clean_mtz.py)."""
    ranged = []
    for child in root.iterdir():
        if child.is_dir():
            m = RANGE_DIR_RE.match(child.name)
            if m:
                ranged.append((int(m.group(1)), child))
    ranged.sort(key=lambda t: t[0])
    return [p for _, p in ranged]


def iter_numeric_subfolders(root: Path) -> list:
    """Numeric subfolders (1, 2, ...) sorted numerically (fetch_clean_mtz.py)."""
    nums = [c for c in root.iterdir() if c.is_dir() and c.name.isdigit()]
    nums.sort(key=lambda x: int(x.name))
    return nums


def ordered_timepoints(root: Path) -> list:
    """Timepoint folders in acquisition order — range folders first, else numeric."""
    return iter_range_subfolders(root) or iter_numeric_subfolders(root)


def find_target_mtz(folder: Path, target_name: str):
    """Locate *target_name* in the folder's autoPROC_* chunk with the lowest index
    (fetch_clean_mtz.py.find_mtz_in_folder, parametrised by the file name)."""
    autoproc = [d for d in folder.iterdir()
                if d.is_dir() and d.name.startswith(SUBFOLDER_PREFIX)]
    candidates = []
    for ap in autoproc:
        mtz = ap / target_name
        if mtz.is_file():
            candidates.append((_natural_int(ap.name), mtz))
    if not candidates:
        return None
    candidates.sort(key=lambda t: t[0])
    return candidates[0][1]


def _run_checked(cmd: list, *, input_text=None) -> subprocess.CompletedProcess:
    """subprocess.run with helpful CCP4 error messages (fetch_clean_mtz.py._run)."""
    try:
        return subprocess.run(cmd, input=input_text, text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    except FileNotFoundError as e:
        raise RuntimeError(f"Required program not found: {cmd[0]}. "
                           f"Is CCP4 set up and {cmd[0]} in your PATH?") from e
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}\n"
                           f"Exit code: {e.returncode}\n--- stdout ---\n{e.stdout}\n"
                           f"--- stderr ---\n{e.stderr}") from e


def _mtz_labels(mtz_path: Path) -> list:
    """Column labels parsed from mtzdmp output (fetch_clean_mtz.py)."""
    cp = _run_checked([MTZDUMP_EXE, str(mtz_path)])
    raw = (cp.stdout or "") + "\n" + (cp.stderr or "")
    labels, in_block = [], False
    for line in raw.splitlines():
        s = line.strip()
        if s.startswith("* Column Labels"):
            in_block = True
            if ":" in s:
                after = s.split(":", 1)[1].strip()
                if after:
                    labels.extend(after.split())
            continue
        if in_block:
            if s.startswith("* ") and not s.startswith("* Column Labels"):
                break
            if not s or (s.startswith("<") and s.endswith(">")):
                continue
            labels.extend(s.split())
    seen, uniq = set(), []
    for lab in labels:
        if lab not in seen:
            uniq.append(lab)
            seen.add(lab)
    return uniq


def clean_mtz(src_mtz: Path, out_mtz: Path, cached_labels=None) -> list:
    """CAD-clean an MTZ to exactly HKL + F SIGF FreeR_flag (fetch_clean_mtz.py).

    Returns the labels used, so the caller can cache them across a series that
    shares one column layout; falls back to per-file mtzdmp if they don't match.
    """
    required = {"F", "SIGF", "FreeR_flag"}
    labels = cached_labels
    if labels is not None and not required.issubset(set(labels)):
        labels = None
    if labels is None:
        labels = _mtz_labels(src_mtz)
    missing = [r for r in sorted(required) if r not in set(labels)]
    if missing:
        raise RuntimeError(
            f"Cannot clean {src_mtz.name}: required columns not found: "
            f"{', '.join(missing)}. Available: {', '.join(labels) or '(none)'}")
    cad_in = ("LABIN FILE 1 E1=F E2=SIGF E3=FreeR_flag\n"
              "LABOUT FILE 1 E1=F E2=SIGF E3=FreeR_flag\n"
              "END\n")
    _run_checked([CAD_EXE, "HKLIN1", str(src_mtz), "HKLOUT", str(out_mtz)],
                 input_text=cad_in)
    return labels


def clean_all_timepoints(copied_root: Path, dfo_dir: Path, base_name: str,
                         target_name: str, dry_run: bool, resume: bool = False) -> int:
    """Clean the chosen MTZ from every timepoint into dfo_dir as
    <base_name>_<timepoint>.mtz. Returns the number of files written."""
    print(f"\n── Stage 2: clean MTZ into {dfo_dir.name}/ ─────────────")
    subdirs = ordered_timepoints(copied_root)
    if not subdirs:
        raise RuntimeError(f"No range (1-300) or numeric (1,2,...) timepoint folders "
                           f"under {copied_root}")
    if not dry_run:
        dfo_dir.mkdir(parents=True, exist_ok=True)
    base = base_name
    written = missing = failed = skipped = 0
    cached = None
    for sd in subdirs:
        label = sd.name.strip().replace(" ", "")
        mtz = find_target_mtz(sd, target_name)
        if mtz is None:
            print(f"  WARNING: no {target_name} for timepoint {label}")
            missing += 1
            continue
        out_path = dfo_dir / f"{base}_{label}.mtz"
        if dry_run:
            print(f"  [dry-run] clean {mtz} -> {out_path.name}")
            continue
        if resume and out_path.is_file():
            print(f"  skip {label} (already cleaned)")
            skipped += 1
            written += 1
            continue
        try:
            with tempfile.TemporaryDirectory(prefix="mtz_clean_") as td:
                tmp = Path(td) / out_path.name
                used = clean_mtz(mtz, tmp, cached_labels=cached)
                cached = cached or used
                shutil.copy2(tmp, out_path)
            print(f"  cleaned {label} -> {out_path.name}")
            written += 1
        except Exception as e:
            print(f"  ERROR cleaning {label}: {e}")
            failed += 1
    extra = f"   already present: {skipped}" if resume else ""
    print(f"  written: {written}   missing: {missing}   failed: {failed}{extra}")
    if written == 0:
        raise RuntimeError("No MTZ files were cleaned — cannot continue.")
    return written


def parse_dimple_r(dimple_log: Path):
    """(R, R_free) from dimple's final refinement, read from dimple.log.

    dimple.log is a configparser-style file; the last refinement section carries
    'overall_r' and 'free_r'. Taking the last occurrence gives the final values.
    """
    r = rfree = None
    if dimple_log.is_file():
        for line in dimple_log.read_text(errors="ignore").splitlines():
            s = line.strip()
            if s.startswith("overall_r:"):
                try:
                    r = float(s.split(":", 1)[1])
                except ValueError:
                    pass
            elif s.startswith("free_r:"):
                try:
                    rfree = float(s.split(":", 1)[1])
                except ValueError:
                    pass
    return r, rfree


def run_dimple(raw_mtz: Path, model: Path, dimple_dir: Path, dry_run: bool,
               resume: bool = False, max_r: float = 0.4, launcher=None) -> Path:
    """Run CCP4 dimple on the reference timepoint; return the refined final.pdb.

    *launcher* is an optional command prefix (e.g. an ``srun`` invocation) so the
    run can execute on a cluster compute node while its output still streams live,
    exactly like a local run. After refinement the final R is read from dimple.log
    and, if it exceeds *max_r*, a clear warning is printed (warn-only).
    """
    launcher = list(launcher or [])
    final_pdb = dimple_dir / "final.pdb"
    print("\n── Stage 3: dimple (reference timepoint) ─────────")
    print(f"  MTZ   : {raw_mtz}")
    print(f"  model : {model}")
    print(f"  out   : {dimple_dir}")
    if dry_run:
        print("  [dry-run] dimple not executed"
              + (f" (would run via: {' '.join(launcher)})" if launcher else ""))
        return final_pdb
    if resume and final_pdb.is_file():
        print(f"  skip — {final_pdb.name} already exists")
        return final_pdb
    dimple_dir.parent.mkdir(parents=True, exist_ok=True)
    # Stream dimple's output live (prefixed) so the user can see it working;
    # PYTHONUNBUFFERED nudges the CCP4 python tools to flush progress promptly.
    where = "on the cluster (srun)" if launcher else "locally"
    print(f"  running dimple {where} (live output below; a few minutes)...")
    print("  " + "-" * 48)
    import os
    import shlex
    env = {**os.environ, "PYTHONUNBUFFERED": "1"}
    dimple_cmd = [DIMPLE_EXE, str(raw_mtz), str(model), str(dimple_dir)]
    # CCP4's scratch dir ($CCP4_SCR) must exist or refmac warns ("No such
    # directory: $CCP4_SCR, refmac shall not work!"). Locally we can just create
    # it; under srun it must be created ON the compute node, so we wrap the call.
    ccp4_scr = env.get("CCP4_SCR")
    if launcher:
        pre = 'mkdir -p "$CCP4_SCR" 2>/dev/null || true; ' if ccp4_scr else ""
        cmd = launcher + ["bash", "-c", pre + "exec " + shlex.join(dimple_cmd)]
    else:
        if ccp4_scr:
            try:
                os.makedirs(ccp4_scr, exist_ok=True)
            except OSError:
                pass
        cmd = dimple_cmd
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True,
                                bufsize=1, env=env)
    except FileNotFoundError as e:
        raise RuntimeError(f"dimple not found in PATH ({DIMPLE_EXE}).") from e
    for line in proc.stdout:
        print(f"  | {line.rstrip()}")
    rc = proc.wait()
    print("  " + "-" * 48)
    if rc != 0:
        raise RuntimeError(f"dimple failed (exit {rc}); see the output above "
                           f"and {dimple_dir / 'dimple.log'}.")
    if not final_pdb.is_file():
        raise RuntimeError(f"dimple finished but {final_pdb} was not produced.")
    print(f"  -> {final_pdb}")

    # B1: report the refinement R and warn (only) if the reference fits poorly.
    r, rfree = parse_dimple_r(dimple_dir / "dimple.log")
    if r is not None:
        rfree_txt = f"{rfree:.3f}" if rfree is not None else "n/a"
        print(f"  dimple R = {r:.3f}   R_free = {rfree_txt}")
        if r > max_r:
            print(f"  ⚠ WARNING: dimple R ({r:.3f}) exceeds {max_r:.2f} — the reference "
                  f"model may fit the data poorly, which weakens every difference map. "
                  f"Continuing anyway.")
    else:
        print("  (could not read R from dimple.log)")
    return final_pdb


def prompt_name(cli_name, default: str) -> str:
    """Dataset name for this run — used for the copied folder AND the cleaned-MTZ
    prefix (<name>_<timepoint>.mtz). CLI flag if given, else asked interactively."""
    if cli_name:
        return cli_name.strip()
    if not sys.stdin.isatty():
        print(f"Non-interactive input: using dataset name '{default}'.")
        return default
    while True:
        ans = input(f"Dataset name for this run [Enter = {default}]: ").strip()
        name = ans or default
        # Keep it filesystem-safe (it becomes a folder name and a filename prefix).
        if name and "/" not in name and name not in (".", ".."):
            return name
        print("  Please enter a valid name (no '/').")


def prompt_file_choice(cli_choice) -> str:
    """truncate vs staraniso — the CLI flag if given, else ask (default truncate)."""
    if cli_choice:
        return cli_choice
    if not sys.stdin.isatty():
        print("Non-interactive input: defaulting to truncate.")
        return "truncate"
    while True:
        ans = input("Which structure factors? [t]runcate / [s]taraniso "
                    "(default t): ").strip().lower()
        if ans in ("", "t", "truncate"):
            return "truncate"
        if ans in ("s", "staraniso"):
            return "staraniso"
        print("  Please answer 't' or 's'.")


def _dry_run_plan(chunks: Path, base_name: str, analysis: Path, dimple_dir: Path,
                  target_name: str, user_model: Path) -> None:
    """Print the Stage 2-4 plan from the SOURCE chunks (nothing is copied yet)."""
    tps = source_timepoints(chunks)

    print(f"\n── Stage 2: clean MTZ into {analysis}/ ─────────────")
    for r, d in tps:
        print(f"  [dry-run] clean {d / target_name} -> {base_name}_{r}.mtz")

    print("\n── Stage 3: dimple (reference timepoint) ─────────")
    if tps:
        r0, d0 = tps[0]
        print(f"  [dry-run] dimple {d0 / target_name}  +  {user_model}")
        print(f"  [dry-run]   -> {dimple_dir}")
        print(f"  [dry-run] copy final.pdb -> {analysis / 'final.pdb'}")

    print(f"\n── Stage 4: TR-FRX analysis pipeline in {analysis} ────")
    print("  [dry-run] would run Fo-Fo diff maps, SVD and peak analysis here.")


# ═════════════════════════════════════════════════════════════════════════
# § J-bis. Reproducibility, validation, reference choice and the HTML report
# ═════════════════════════════════════════════════════════════════════════
STAGE_ORDER = ["copy", "clean", "dimple", "analyze"]


class _Tee:
    """Duplicate a text stream to a log file (for a full run transcript)."""

    def __init__(self, stream, log):
        self._stream, self._log = stream, log

    def write(self, s):
        self._stream.write(s)
        if not self._log.closed:
            self._log.write(s)
        return len(s)

    def flush(self):
        self._stream.flush()
        if not self._log.closed:
            self._log.flush()

    def isatty(self):
        return getattr(self._stream, "isatty", lambda: False)()


def setup_tee(log_path: Path):
    """Send everything printed from here on to *log_path* as well as the console."""
    log = open(log_path, "a", encoding="utf-8")
    sys.stdout = _Tee(sys.stdout, log)
    sys.stderr = _Tee(sys.stderr, log)
    print(f"  logging to  -> {log_path}")
    return log


def close_tee(log) -> None:
    """Restore the real stdout/stderr, then close the log file (avoids a flush
    on a closed file at interpreter shutdown)."""
    if log is None:
        return
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__
    try:
        log.close()
    except Exception:
        pass


def write_run_config(path: Path, cfg: dict) -> None:
    """Record exactly how this run was launched (reproducibility)."""
    lines = ["TR-FRX full pipeline — run configuration", "=" * 44]
    for k, v in cfg.items():
        lines.append(f"{k:<18}: {v}")
    path.write_text("\n".join(lines) + "\n")
    print(f"  run config  -> {path}")


def validate_model_pdb(model: Path) -> None:
    """Warn (don't fail) if the model lacks a CRYST1 record."""
    try:
        text = model.read_text(errors="ignore")
    except Exception:
        return
    if not any(ln.startswith("CRYST1") for ln in text.splitlines()):
        print(f"  WARNING: {model.name} has no CRYST1 record — symmetry-dependent "
              f"steps (PyMOL symexp, phasing) may misbehave.")


def validate_reference_mtz(raw_mtz: Path) -> None:
    """Fail fast if the chosen SF file lacks the columns CAD cleaning needs."""
    required = {"F", "SIGF", "FreeR_flag"}
    labels = _mtz_labels(raw_mtz)
    missing = required - set(labels)
    if missing:
        raise RuntimeError(
            f"{raw_mtz.name} is missing column(s) {', '.join(sorted(missing))} "
            f"required for cleaning. Found: {', '.join(labels) or '(none)'}. "
            f"If this file uses FP/SIGFP, pick the other --file or relabel it.")
    print(f"  labels OK   : {' '.join(labels)}  ({raw_mtz.name})")


def source_timepoints(chunks: Path) -> list:
    """[(label, chunk_dir), ...] from the SOURCE autoPROC_* chunks, in order."""
    dirs = sorted(p for p in chunks.iterdir()
                  if p.is_dir() and p.name.startswith(SUBFOLDER_PREFIX))
    tps = []
    for d in dirs:
        r = d.name[len(SUBFOLDER_PREFIX):]
        tps.append((_natural_int(r), r, d))
    tps.sort(key=lambda t: t[0])
    return [(r, d) for _, r, d in tps]


def resolve_reference_timepoint(timepoints: list, ref_label):
    """Pick the reference timepoint folder (default = first / lowest index)."""
    if not ref_label:
        return timepoints[0]
    want = ref_label.strip().replace(" ", "")
    for tp in timepoints:
        if tp.name.strip().replace(" ", "") == want:
            return tp
    avail = ", ".join(tp.name for tp in timepoints)
    raise RuntimeError(f"--ref-timepoint '{ref_label}' not found. Available: {avail}")


def next_run_dir(base: Path, resume: bool) -> Path:
    """Pick the run_NN analysis folder: reuse the latest on --resume, else the next.

    Copy/clean/dimple live at the dataset level (shared); each analysis rerun is a
    fresh run_NN so results are kept side-by-side (delete old ones yourself).
    """
    existing = sorted(d for d in base.glob("run_*")
                      if d.is_dir() and d.name[4:].isdigit()) if base.exists() else []
    if resume and existing:
        return existing[-1]
    nums = [int(d.name[4:]) for d in existing]
    return base / f"run_{(max(nums) + 1) if nums else 1:02d}"


# ── HTML report (self-contained: images embedded as base64, opens over SSH) ──
def _img_data_uri(png_path: Path) -> str:
    import base64
    return "data:image/png;base64," + base64.b64encode(png_path.read_bytes()).decode("ascii")


def _read_csv_rows(path: Path):
    import csv
    if not path or not path.is_file():
        return [], []
    with path.open(newline="") as f:
        rows = list(csv.reader(f))
    return (rows[0], rows[1:]) if rows else ([], [])


def _html_table(headers, rows, max_rows=None) -> str:
    import html as H
    if not headers:
        return "<p class='muted'>(no data)</p>"
    shown = rows if max_rows is None else rows[:max_rows]
    th = "".join(f"<th>{H.escape(str(h))}</th>" for h in headers)
    trs = "".join("<tr>" + "".join(f"<td>{H.escape(str(c))}</td>" for c in r) + "</tr>"
                  for r in shown)
    more = ""
    if max_rows is not None and len(rows) > max_rows:
        more = f"<p class='muted'>… {len(rows) - max_rows} more row(s) — see CSV/PDF.</p>"
    return f"<table><thead><tr>{th}</tr></thead><tbody>{trs}</tbody></table>{more}"


_HTML_CSS = """
* { box-sizing: border-box; }
body { font-family: -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif;
       margin: 0; padding: 0 0 4rem; color: #1a1a1a; background: #f6f7f9; }
header { background: #10243e; color: #fff; padding: 1.2rem 1.6rem; }
header h1 { margin: 0 0 .2rem; font-size: 1.4rem; }
header .meta { font-size: .85rem; opacity: .85; }
main { max-width: 1100px; margin: 0 auto; padding: 0 1.2rem; }
section { background: #fff; border: 1px solid #e2e6ea; border-radius: 8px;
          margin: 1.2rem 0; padding: 1rem 1.2rem; }
h2 { font-size: 1.1rem; margin: .2rem 0 .8rem; border-bottom: 2px solid #10243e;
     padding-bottom: .3rem; }
.tablewrap { overflow-x: auto; }
table { border-collapse: collapse; width: 100%; font-size: .82rem; }
th, td { border: 1px solid #e2e6ea; padding: 3px 7px; text-align: right; white-space: nowrap; }
th { background: #eef1f4; }
td:first-child, th:first-child { text-align: left; }
.muted { color: #6b7280; font-size: .8rem; }
.imgrid { display: grid; grid-template-columns: repeat(auto-fill, minmax(300px, 1fr));
          gap: .8rem; margin-top: .8rem; }
.imgrid figure { margin: 0; border: 1px solid #e2e6ea; border-radius: 6px; padding: .4rem; }
.imgrid img { width: 100%; height: auto; border-radius: 4px; }
.imgrid figcaption { font-size: .75rem; color: #374151; margin-top: .3rem; }
details { margin: .5rem 0; }
summary { cursor: pointer; font-weight: 600; padding: .3rem 0; }
.plot img { max-width: 100%; height: auto; }
.tocwrap { font-size: .9rem; margin: 1rem 0 0; }
.tocwrap a { color: #10243e; text-decoration: none; margin-right: 1rem; }
"""


def build_html_report(dfo: Path, out_html: Path, name: str, meta: dict) -> Path:
    """One self-contained HTML page: timepoint summary, SVD, persistence, and a
    collapsible per-map section with the peak tables and images embedded inline.
    Reads the artifacts the pipeline already writes (stdlib only)."""
    import html as H
    from datetime import datetime

    dfo_dir = dfo / "output_dfo"
    svd_dir = dfo / "output_svd"
    reports = dfo_dir / "reports"
    recap_dir = svd_dir if (svd_dir / "series_recap.csv").is_file() else dfo_dir

    def _wrap(t):
        return f"<div class='tablewrap'>{t}</div>"

    parts = []

    rc_head, rc_rows = _read_csv_rows(recap_dir / "series_recap.csv")
    parts.append("<section id='summary'><h2>Timepoint summary</h2>"
                 + _wrap(_html_table(rc_head, rc_rows)) + "</section>")

    mq_head, mq_rows = _read_csv_rows(recap_dir / "map_quality.csv")
    if mq_head:
        parts.append("<section id='quality'><h2>Map quality</h2>"
                     + _wrap(_html_table(mq_head, mq_rows)) + "</section>")

    sv_head, sv_rows = _read_csv_rows(svd_dir / "singular_values.csv")
    plot_png = svd_dir / "rSV_plot.png"
    plot_html = (f"<div class='plot'><img src='{_img_data_uri(plot_png)}' "
                 f"alt='rSV plot'></div>") if plot_png.is_file() else ""
    if sv_head or plot_html:
        parts.append(f"<section id='svd'><h2>SVD</h2>{plot_html}"
                     + (_wrap(_html_table(sv_head, sv_rows)) if sv_head else "")
                     + "</section>")

    ps_head, ps_rows = _read_csv_rows(recap_dir / "series_recap_persistence.csv")
    if ps_head:
        parts.append("<section id='persist'><h2>Persistent sites</h2>"
                     + _wrap(_html_table(ps_head, ps_rows)) + "</section>")

    if reports.is_dir():
        map_blocks = []
        for md in sorted(d for d in reports.iterdir() if d.is_dir()):
            ph, pr = _read_csv_rows(md / "peaks.csv")
            figs, n_imgs = "", 0
            for sub in ("img", "img2d"):     # 3D cages and/or 2D sections
                d = md / sub
                if not d.is_dir():
                    continue
                for p in sorted(d.glob("peak_*.png")):
                    figs += (f"<figure><img src='{_img_data_uri(p)}' "
                             f"alt='{H.escape(sub + '/' + p.stem)}'>"
                             f"<figcaption>{H.escape(sub + '/' + p.stem)}</figcaption></figure>")
                    n_imgs += 1
            grid = f"<div class='imgrid'>{figs}</div>" if figs else ""
            map_blocks.append(
                f"<details><summary>{H.escape(md.name)} "
                f"({len(pr)} peaks, {n_imgs} images)</summary>"
                f"{_wrap(_html_table(ph, pr, max_rows=40))}{grid}</details>")
        if map_blocks:
            parts.append("<section id='maps'><h2>Per-map peaks &amp; images</h2>"
                         + "".join(map_blocks) + "</section>")

    toc = ("<div class='tocwrap'><a href='#summary'>Summary</a>"
           "<a href='#svd'>SVD</a><a href='#persist'>Persistence</a>"
           "<a href='#maps'>Maps</a></div>")
    metabits = " · ".join(f"{k}: {v}" for k, v in meta.items())
    doc = (f"<!doctype html><html><head><meta charset='utf-8'>"
           f"<meta name='viewport' content='width=device-width, initial-scale=1'>"
           f"<title>{H.escape(name)} — TR-FRX report</title><style>{_HTML_CSS}</style></head>"
           f"<body><header><h1>{H.escape(name)} — TR-FRX report</h1>"
           f"<div class='meta'>{H.escape(metabits)} · generated "
           f"{datetime.now().strftime('%Y-%m-%d %H:%M')}</div></header>"
           f"<main>{toc}{''.join(parts)}</main></body></html>")
    out_html.write_text(doc, encoding="utf-8")
    return out_html


def build_report_for_analysis(analysis: Path, meta_extra: dict = None) -> Path:
    """Build report.html for a finished analysis dir, deriving the run folder and
    dataset name from the path <output>/<name>/run_NN/analysis. Used by the
    cluster final job so --cluster runs still produce a report."""
    analysis = Path(analysis).resolve()
    run_dir  = analysis.parent
    name     = analysis.parent.parent.name or run_dir.name
    meta     = {"run": run_dir.name}
    if meta_extra:
        meta.update(meta_extra)
    model = analysis / "final.pdb"
    if model.is_file():
        meta.setdefault("model", model.name)
    return build_html_report(analysis, run_dir / "report.html", name, meta)


# ═════════════════════════════════════════════════════════════════════════
# § K. main — the merged end-to-end orchestration
# ═════════════════════════════════════════════════════════════════════════
def main() -> int:
    argv = sys.argv[1:]

    # Advanced / cluster re-entry: the original single-stage pipeline CLI is
    # preserved verbatim and selected whenever a pipeline-only flag is present.
    # (submit_cluster re-invokes this script with --dir/--only/--svd-only.)
    if "--dir" in argv or "--only" in argv or "--svd-only" in argv:
        return run_pipeline(_legacy_pipeline_parser().parse_args(argv))

    parser = argparse.ArgumentParser(
        description="TR-FRX end-to-end pipeline: copy autoPROC chunks -> clean "
                    "MTZ -> dimple (first timepoint) -> Fo-Fo maps + SVD + peaks.",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    parser.add_argument("model",  type=Path, nargs="?", help="Reference PDB model.")
    parser.add_argument("chunks", type=Path, nargs="?",
                        help="Folder of autoPROC_* chunks (from the autoPROC step).")
    parser.add_argument("output", type=Path, nargs="?", help="Output folder.")
    parser.add_argument("--name", default=None, metavar="NAME",
                        help="Dataset name for this run (output subfolder + cleaned "
                             "MTZ prefix). Asked interactively if omitted.")
    parser.add_argument("--file", choices=["staraniso", "truncate"], default=None,
                        help="Structure-factor file for dimple + cleaning "
                             "(asked interactively if omitted).")
    parser.add_argument("--workers", type=int, default=4, metavar="N",
                        help="Parallel copy threads (default 4).")
    parser.add_argument("--resume", action="store_true",
                        help="Skip work already done (copied files, cleaned MTZs, "
                             "dimple final.pdb).")
    parser.add_argument("--skip-dimple", action="store_true",
                        help="Skip dimple; use the provided model directly for the maps.")
    parser.add_argument("--max-r", type=float, default=0.4, metavar="R",
                        help="Warn if dimple's final R exceeds this (default 0.4; "
                             "warn-only, never stops).")
    parser.add_argument("--ref-timepoint", default=None, metavar="LABEL",
                        help="Timepoint used as reference / dimpled (e.g. 1_300). "
                             "Default: the lowest-index timepoint.")
    parser.add_argument("--start-stage", choices=STAGE_ORDER, default="copy",
                        help="First stage to run (default copy).")
    parser.add_argument("--stop-after", choices=STAGE_ORDER, default="analyze",
                        help="Last stage to run (default analyze).")
    parser.add_argument("--no-html", action="store_true",
                        help="Do not build the HTML summary report.")
    # Forwarded to the analysis stage (same meaning as the pipeline CLI)
    parser.add_argument("--sigma", default=SIGMA_DEFAULT, metavar="N|auto",
                        help="Peak cutoff: 'auto' (default) or a number of sigma.")
    parser.add_argument("--radius", type=float, default=RADIUS_DEFAULT, metavar="A",
                        help=f"Base nearest-residue search radius (default {RADIUS_DEFAULT}).")
    parser.add_argument("--n-peaks", type=int, default=N_PEAKS_DEFAULT, metavar="N",
                        help=f"Peaks imaged per map (default {N_PEAKS_DEFAULT}).")
    parser.add_argument("--display-sigma", type=float, default=DISPLAY_SIGMA, metavar="N",
                        help=f"Mesh contour level in figures (default {DISPLAY_SIGMA}).")
    parser.add_argument("--display-near", type=float, default=DISPLAY_NEAR, metavar="A",
                        help=f"Show residues within this many A of the peak "
                             f"(default {DISPLAY_NEAR}).")
    parser.add_argument("--time-step", type=float, default=None, metavar="MS",
                        help="Time interval in ms between datasets (SVD X-axis).")
    parser.add_argument("--high-res", type=float, default=None, metavar="A",
                        help="Fixed high-res cutoff for ALL maps (else interactive).")
    parser.add_argument("--skip-images", action="store_true",
                        help="Skip PyMOL peak images (tables + PDFs still made).")
    parser.add_argument("--skip-svd", action="store_true", help="Diffmaps only.")
    parser.add_argument("--active-site", type=int, nargs="+", default=None,
                        metavar="RESNUM", help="Residue numbers of interest.")
    parser.add_argument("--cluster", action="store_true",
                        help="Submit the analysis stage to SLURM.")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print the plan for every stage without running it.")
    args = parser.parse_args()

    if not (args.model and args.chunks and args.output):
        parser.error("MODEL, CHUNKS and OUTPUT are required "
                     "(or pass --dir for pipeline-only mode).")

    user_model = args.model.expanduser().resolve()
    chunks     = args.chunks.expanduser().resolve()
    output     = args.output.expanduser().resolve()

    if not user_model.is_file():
        print(f"ERROR: model PDB not found: {user_model}")
        return 1
    if not chunks.is_dir():
        print(f"ERROR: chunks folder not found: {chunks}")
        return 1

    # Stage 1-3 tools (analysis tools are checked later by the pipeline itself).
    for exe in (CAD_EXE, MTZDUMP_EXE, DIMPLE_EXE):
        if not shutil.which(exe):
            print(f"ERROR: '{exe}' not in PATH — set up CCP4 (module load ccp4).")
            return 1

    name = prompt_name(args.name, chunks.name)
    choice = prompt_file_choice(args.file)
    target_name = TARGET_FILES[choice]

    start = STAGE_ORDER.index(args.start_stage) + 1
    stop  = STAGE_ORDER.index(args.stop_after) + 1
    if start > stop:
        print(f"ERROR: --start-stage ({args.start_stage}) is after "
              f"--stop-after ({args.stop_after}).")
        return 1

    print(f"\nDataset name     : {name}")
    print(f"Structure factors: {choice}  ({target_name})")
    print(f"Stages           : {args.start_stage} .. {args.stop_after}"
          + ("   [resume]" if args.resume else ""))

    logging.basicConfig(level=logging.INFO, format="%(message)s")

    # ---- Output layout (numbered reruns; copy/clean/dimple shared) ----
    #   <output>/<name>/autoproc_copy/    (Stage 1, shared across reruns)
    #   <output>/<name>/run_NN/analysis/  (cleaned MTZ, final.pdb, output_dfo/, output_svd/)
    #   <output>/<name>/run_NN/dimple/, pipeline_<stamp>.log, report.html, run_config.txt
    base        = output / name
    copied_root = base / "autoproc_copy"
    run_dir     = next_run_dir(base, args.resume)
    analysis    = run_dir / "analysis"
    dimple_dir  = run_dir / "dimple"
    if not args.dry_run:
        for d in (base, run_dir, analysis):
            d.mkdir(parents=True, exist_ok=True)

    # Reproducibility: tee the console to a log in the run folder + record config.
    logf = None
    if not args.dry_run:
        from datetime import datetime
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        logf = setup_tee(run_dir / f"pipeline_{stamp}.log")
        print(f"Run folder       : {run_dir}")
        write_run_config(run_dir / "run_config.txt", {
            "timestamp": stamp, "name": name, "run": run_dir.name,
            "structure_factors": choice, "target_file": target_name,
            "model": user_model, "chunks": chunks, "output": output,
            "stages": f"{args.start_stage}..{args.stop_after}",
            "resume": args.resume, "skip_dimple": args.skip_dimple,
            "ref_timepoint": args.ref_timepoint or "(first)",
            "sigma": args.sigma, "radius": args.radius, "n_peaks": args.n_peaks,
            "high_res": args.high_res, "time_step": args.time_step,
            "cad": shutil.which(CAD_EXE), "dimple": shutil.which(DIMPLE_EXE),
            "phenix": shutil.which("phenix.fobs_minus_fobs_map"),
        })

    ref_mtz = None
    model_for_pipeline = user_model
    try:
        # ---- Fail-fast validation (before any heavy work) ----
        validate_model_pdb(user_model)
        src_tps = source_timepoints(chunks)
        if not src_tps and start <= 2:
            print(f"ERROR: no '{SUBFOLDER_PREFIX}*' chunks in {chunks}.")
            return 1
        if src_tps and not args.dry_run:
            want = args.ref_timepoint.strip().replace(" ", "") if args.ref_timepoint else None
            if want:
                ref_src = next((d for r, d in src_tps if r.replace(" ", "") == want), None)
                if ref_src is None:
                    avail = ", ".join(r for r, _ in src_tps)
                    print(f"ERROR: --ref-timepoint '{args.ref_timepoint}' not found. "
                          f"Available: {avail}")
                    return 1
            else:
                ref_src = src_tps[0][1]
            ref_raw = ref_src / target_name
            if ref_raw.is_file():
                validate_reference_mtz(ref_raw)
            else:
                print(f"  WARNING: reference {target_name} not found in {ref_src.name}; "
                      f"skipping label validation.")

        # ---- Stage 1: copy (shared; skip files already present) ----
        if start <= 1 <= stop:
            copied_root = copy_chunks(chunks, copied_root, DEFAULT_COPY_FILES,
                                      args.dry_run, args.workers,
                                      resume=args.resume or copied_root.exists())
        else:
            print(f"\n(skip Stage 1 copy — using existing {copied_root})")

        if args.dry_run:
            _dry_run_plan(chunks, name, analysis, dimple_dir, target_name, user_model)
            return 0

        timepoints = ordered_timepoints(copied_root) if copied_root.is_dir() else []

        # ---- Stage 2: clean into run_NN/analysis ----
        if start <= 2 <= stop:
            clean_all_timepoints(copied_root, analysis, name, target_name,
                                 args.dry_run, args.resume)

        # Reference timepoint (drives dimple + the analysis --ref)
        ref_tp = resolve_reference_timepoint(timepoints, args.ref_timepoint) if timepoints else None
        if ref_tp is not None:
            cand = analysis / f"{name}_{ref_tp.name.strip().replace(' ', '')}.mtz"
            ref_mtz = cand if cand.is_file() else None

        # ---- Stage 3: dimple (reference timepoint) -> run_NN/dimple ----
        if start <= 3 <= stop and not args.skip_dimple:
            if ref_tp is None:
                print("ERROR: no timepoint folders found for dimple.")
                return 1
            raw_mtz = find_target_mtz(ref_tp, target_name)
            if raw_mtz is None:
                print(f"ERROR: no {target_name} in reference timepoint {ref_tp.name}.")
                return 1
            # D1: with --cluster, run dimple on a compute node via srun (blocking +
            # streamed, so it looks just like a local run), not on the login node.
            dimple_launcher = None
            if args.cluster:
                if shutil.which("srun"):
                    dimple_launcher = ["srun", "-p", SLURM_PARTITION, "-n", "1",
                                       "-c", "4", "--mem", "8000", "-J", "dimple"]
                else:
                    print("  NOTE: --cluster set but 'srun' not found; "
                          "running dimple locally instead.")
            final_pdb = run_dimple(raw_mtz, user_model, dimple_dir,
                                   args.dry_run, args.resume, max_r=args.max_r,
                                   launcher=dimple_launcher)
            analysis.mkdir(parents=True, exist_ok=True)
            shutil.copy2(final_pdb, analysis / "final.pdb")
            print(f"  copied final.pdb -> {analysis / 'final.pdb'}")
        elif args.skip_dimple and start <= 3 <= stop:
            print("\n(skip Stage 3 dimple — using the provided model directly)")

        # Model the analysis will use
        if args.skip_dimple:
            model_for_pipeline = user_model
        elif (analysis / "final.pdb").is_file():
            model_for_pipeline = analysis / "final.pdb"
        else:
            print(f"  WARNING: {analysis / 'final.pdb'} not found; using provided model.")
            model_for_pipeline = user_model

        if not (start <= 4 <= stop):
            print("\n(stopping before Stage 4 analysis, as requested)")
            if logf:
                close_tee(logf)
            return 0
    except RuntimeError as e:
        print(f"ERROR: {e}")
        if logf:
            close_tee(logf)
        return 1

    # ---- Stage 4: the Fo-Fo / SVD / peak pipeline, run inside run_NN/analysis ----
    print(f"\n── Stage 4: TR-FRX analysis pipeline in {analysis} ────")
    pargs = _legacy_pipeline_parser().parse_args([])   # pipeline defaults
    pargs.dir           = analysis
    pargs.model         = model_for_pipeline
    pargs.ref           = ref_mtz
    pargs.sigma         = args.sigma
    pargs.radius        = args.radius
    pargs.n_peaks       = args.n_peaks
    pargs.display_sigma = args.display_sigma
    pargs.display_near  = args.display_near
    pargs.time_step     = args.time_step
    pargs.high_res      = args.high_res
    pargs.skip_images   = args.skip_images
    pargs.skip_svd      = args.skip_svd
    pargs.active_site   = args.active_site
    pargs.cluster       = args.cluster
    rc = run_pipeline(pargs)

    # HTML summary report (self-contained; embeds the peak images for SSH viewing).
    # E2: build it whenever we ran the analysis, even if some timepoints failed.
    if not args.cluster and not args.no_html:
        try:
            html_path = build_html_report(analysis, run_dir / "report.html", name, {
                "run": run_dir.name,
                "structure factors": choice,
                "model": model_for_pipeline.name,
                "reference": (ref_mtz.name if ref_mtz else "(auto)"),
            })
            print(f"\nHTML report -> {html_path}")
        except Exception as e:
            print(f"  HTML report skipped ({e})")
    elif args.cluster and not args.no_html:
        print(f"\n(cluster) SLURM jobs submitted. Once they finish, the final job "
              f"writes:\n    {run_dir / 'report.html'}")

    if logf:
        close_tee(logf)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
