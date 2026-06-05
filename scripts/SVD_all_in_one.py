#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_series_diffmaps.py — Compute Fo-Fo difference maps for one series, then run SVD.

For each MTZ in the working directory (except the reference and already-computed
dFo_* maps), runs phenix.fobs_minus_fobs_map against the reference MTZ.
SVD analysis is then run directly — no external script needed.

Supported MTZ trailing-index patterns:
    _N.mtz          e.g. CaMDH_073_137_2.mtz
    _start-end.mtz  e.g. CaMDH_073_275_1-250.mtz
    _start_end.mtz  e.g. CaMDH_073_275_1_300.mtz

All output goes into ./output/ (created automatically).

Usage:
    ./run_series_diffmaps.py                         # auto-detect ref and model
    ./run_series_diffmaps.py --ref sample_1_300.mtz  # explicit reference
    ./run_series_diffmaps.py --model refined.pdb     # explicit model
    ./run_series_diffmaps.py --skip-svd              # diffmaps only
    ./run_series_diffmaps.py --dry-run               # print plan, do nothing
"""


from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

# ─────────────────────────────────────────────
# Startup dependency check
# ─────────────────────────────────────────────
def check_dependencies() -> bool:
    ok = True

    # External programs
    programs = {
        "mtzdmp":                     "module load ccp4",
        "cad":                        "module load ccp4",
        "phenix.fobs_minus_fobs_map": "module load phenix",
    }
    for prog, fix in programs.items():
        if not shutil.which(prog):
            print(f"  MISSING program : {prog}  ->  {fix}")
            ok = False

    # Required Python packages
    packages = {
        "numpy":  "pip install numpy",
        "pandas": "pip install pandas",
        "gemmi":  "pip install gemmi   (or: conda install -c conda-forge gemmi)",
    }
    for pkg, fix in packages.items():
        try:
            __import__(pkg)
        except ImportError:
            print(f"  MISSING package  : {pkg}  ->  {fix}")
            ok = False

    # Optional packages (warn, do not block)
    optional = {
        "dask":    "pip install dask[array]  (SVD uses numpy fallback if absent)",
        "scipy":   "pip install scipy        (needed for SV0 monoexp fit)",
        "seaborn": "pip install seaborn      (needed for SVD plots)",
    }
    for pkg, note in optional.items():
        try:
            __import__(pkg)
        except ImportError:
            print(f"  optional missing : {pkg}  ->  {note}")

    return ok


# ─────────────────────────────────────────────
# MTZ index patterns (tried in order)
# ─────────────────────────────────────────────
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


# ─────────────────────────────────────────────
# CCP4 helpers
# ─────────────────────────────────────────────
def _run_silent(cmd: list) -> str:
    try:
        cp = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
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
    """
    Extract the high-resolution limit.

    Tries gemmi first (fast, no subprocess, always available when SVD is used),
    then falls back to mtzdmp for environments where gemmi is not installed.
    """
    # Primary: gemmi reads directly from the binary header — no subprocess needed
    try:
        import gemmi
        res = gemmi.read_mtz_file(str(mtz)).resolution_high()
        if res and res > 0:
            return round(res, 4)
    except Exception:
        pass

    # Fallback: mtzdmp (replicates the original bash logic)
    txt   = _run_silent(["mtzdmp", str(mtz)])
    lines = txt.splitlines()

    for i, line in enumerate(lines):
        if "Resolution Range" in line or "Resolution range" in line:
            for bline in lines[i:i+4]:
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


def collect_source_mtz(directory: Path) -> list:
    """Return source MTZ files sorted by index. Excludes dFo_* files."""
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


# ─────────────────────────────────────────────
# Difference map
# ─────────────────────────────────────────────
def compute_diffmap(target, ref, model, labels, reslim, dry_run, out_dir):
    stem     = target.stem
    out_dir.mkdir(parents=True, exist_ok=True)
    eff_path = out_dir / f"dFo-Fo_{stem}.eff"
    log_path = out_dir / f"{stem}-{ref.stem}_diffmap.log"
    prefix   = out_dir / f"dFo_{stem}"

    eff_content = (
        f"f_obs_1_file_name = {target}\n"
        f"f_obs_1_label = {labels}\n"
        f"f_obs_2_file_name = {ref}\n"
        f"f_obs_2_label = {labels}\n"
        f"high_resolution = {reslim}\n"
        f"low_resolution = 10.0\n"
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


# ─────────────────────────────────────────────
# SVD
# ─────────────────────────────────────────────
def run_svd(dfo_dir: Path, out_dir: Path, time_step_ms: float | None = None) -> None:
    """
    Scan work_dir for dFo_*.mtz, build a map matrix, run SVD, write output.

    X-axis (time in ms) from dataset number n, depending on series kind:
      *_550_* : n=1 -> 0 ms,  n=2 -> 550 ms,  n=3 -> 1100 ms, ...
      *_275_* : n=1 -> 0 ms,  n=2 -> 275 ms,  n=3 -> 550 ms,  ...
      *_137_* : n=1,2 -> 0 ms; n>=3 -> floor((n-2)*137.5) ms
    """
    import numpy as np
    import pandas as pd
    import gemmi

    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
        HAVE_PLOT = True
        HAVE_SNS  = True
    except Exception as _e:
        print(f"  SVD: plotting unavailable ({_e}). Is the svd_pipeline env active?")
        HAVE_PLOT = False
        HAVE_SNS  = False

    try:
        from scipy.optimize import curve_fit
        HAVE_SCIPY = True
    except Exception:
        HAVE_SCIPY = False

    PLOT_SV_LIMIT = 7
    DO_FIT_SV0    = True
    FIT_MODEL     = "decay"     # "decay" or "rise"

    # ── Series helpers ───────────────────────
    def dataset_number(name):
        ints = re.findall(r"\d+", name)
        return int(ints[-2]) if len(ints) >= 2 else None

    def to_x(n, path, step):
        """Return X axis value: time in ms if step given, else the dFo filename."""
        if step is not None:
            return (n - 1) * step
        return path.stem   # dFo filename without extension

    # ── MTZ <-> matrix ───────────────────────
    _grid_size = None   # set from the first map, reused for all others

    def mtz_to_matrix(path):
        nonlocal _grid_size
        mtz  = gemmi.read_mtz_file(str(path))
        if _grid_size is None:
            # Determine grid from first map, then lock it for all maps
            grid      = mtz.transform_f_phi_to_map("FoFo", "PHFc")
            arr       = np.array(grid, copy=False)
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

    # ── Collect dFo_*.mtz ────────────────────
    entries = []
    for p in sorted(dfo_dir.glob("dFo_*.mtz")):
        n = dataset_number(p.stem)
        if n is None:
            print(f"  SVD: skipping {p.name} (cannot parse dataset number)")
            continue
        entries.append((n, p))

    if not entries:
        print("  SVD: no dFo_*.mtz files found — skipping.")
        return

    step  = time_step_ms
    x_label = f"Time (ms), step={step} ms" if step is not None else "Dataset"
    print(f"  SVD: X axis = {x_label}, {len(entries)} maps")

    df = pd.DataFrame({
        "n":      [e[0] for e in entries],
        "x_num":  [float(e[0] - 1) if step is None else float((e[0]-1)*step)
                   for e in entries],                          # always numeric for numpy
        "x_label":[e[1].stem if step is None
                   else f"{(e[0]-1)*step} ms"
                   for e in entries],                         # display label
        "path":   [e[1] for e in entries],
    }).sort_values("n").reset_index(drop=True)

    out_dir.mkdir(exist_ok=True)

    # ── Build matrix A ───────────────────────
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

    # ── SVD ─────────────────────────────────
    try:
        import dask.array as da
        U, S, VT = da.linalg.svd(da.array(A))
        U  = np.array(U); S = np.array(S); VT = np.array(VT)
        print("  SVD: using dask backend")
    except ImportError:
        print("  SVD: dask not found, using numpy.linalg.svd")
        U, S, VT = np.linalg.svd(A, full_matrices=False)

    sigmat = np.zeros((n_files, n_files), dtype=np.float32)
    for i in range(len(S)):
        sigmat[i, i] = S[i]
    US = np.array(np.matrix(U) @ np.matrix(sigmat))

    # ── Write maps ───────────────────────────
    template = df.path[0]
    for rank in range(n_files):
        matrix_to_mtz(U[:, rank].reshape(shape),  template, out_dir / f"lSV_{rank}.mtz")
        matrix_to_mtz(US[:, rank].reshape(shape), template, out_dir / f"SVscaled_lSV_{rank}.mtz")

    if n_files >= 2:
        Ar1 = np.array(US @ VT.T[1, :])
        matrix_to_mtz(Ar1.reshape(shape), df.path[1], out_dir / f"reformed_{df.path[1].stem}.mtz")

    # ── CSV ──────────────────────────────────
    t_index   = np.array(df.x_num, dtype=float)
    x_labels  = list(df.x_label)
    scaled_tf = (sigmat @ VT)[:n_files, :]
    index_label = "time_ms" if step is not None else "dataset"
    pd.DataFrame({f"rSV{i}": scaled_tf[i] for i in range(n_files)},
                 index=x_labels).to_csv(out_dir / "rSV.csv", index_label=index_label)

    # Print singular values (weights)
    print(f"  SVD: singular values (relative weight):")
    s_total = S.sum()
    for i, s in enumerate(S):
        bar = "#" * int(30 * s / s_total)
        print(f"    SV{i}: {s:.4g}  ({100*s/s_total:.1f}%)  {bar}")

    print(f"  SVD: wrote lSV/SVscaled maps and rSV.csv -> {out_dir}")

    # ── Plot ─────────────────────────────────
    if not HAVE_PLOT:
        print("  SVD: matplotlib not available — skipping plots.")
        return

    def monoexp(t, A, tau, C):
        return A * np.exp(-t / tau) + C if FIT_MODEL == "decay" \
               else A * (1 - np.exp(-t / tau)) + C

    svlim   = min(PLOT_SV_LIMIT, n_files)
    if HAVE_SNS:
        palette = [sns.color_palette("bright", n_colors=svlim)[i] for i in range(svlim)]
    else:
        palette = [plt.cm.tab10(i / 10) for i in range(svlim)]

    order    = np.argsort(t_index)
    t        = t_index[order]
    t_labels = [x_labels[i] for i in order]

    fig, ax = plt.subplots(figsize=(10, 6))

    for i in range(svlim):
        yi = np.array(scaled_tf[i], dtype=float)[order]
        ax.plot(range(len(t)), yi, "o--", color=palette[i],
                markersize=8, linewidth=2, alpha=0.7, label=f"SV{i}")

        # SV0 monoexp fit
        if i == 0 and DO_FIT_SV0 and HAVE_SCIPY and step is not None and len(t) >= 4:
            try:
                C0   = float(np.mean(yi[-max(1, len(yi) // 10):]))
                A0   = float(yi[0] - C0) if FIT_MODEL == "decay" else float(yi[-1] - C0)
                tau0 = float(max((t.max() - t.min()) / 3.0, 1.0))
                popt, _ = curve_fit(monoexp, t, yi, p0=[A0, tau0, C0],
                                    bounds=([-np.inf, 1e-6, -np.inf], [np.inf, np.inf, np.inf]),
                                    maxfev=20000)
                tfit = np.linspace(float(t.min()), float(t.max()), 300)
                # Map tfit back to plot x coordinates
                xfit = np.interp(tfit, t, range(len(t)))
                ax.plot(xfit, monoexp(tfit, *popt), "-", color=palette[0],
                        linewidth=2.5, alpha=0.95,
                        label=f"SV0 fit (tau={popt[1]:.2f} ms)")
                with open(out_dir / "SV0_monoexp_fit.txt", "w") as f:
                    f.write(f"model {FIT_MODEL}\nA {popt[0]}\ntau_ms {popt[1]}\nC {popt[2]}\n")
                print(f"  SVD: SV0 fit -> A={popt[0]:.4g}, tau={popt[1]:.4g} ms, C={popt[2]:.4g}")
            except Exception as e:
                print(f"  SVD: SV0 fit failed ({e})")

    ax.set_xticks(range(len(t_labels)))
    ax.set_xticklabels(t_labels, rotation=45, ha="right", fontsize=8)
    ax.set_xlabel(x_label)
    ax.set_ylabel("Magnitude")
    ax.set_title("Right Singular Vectors")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    plot_path = out_dir / "rSV_plot.png"
    fig.savefig(plot_path, dpi=150)
    print(f"  SVD: plot saved -> {plot_path}")
    plt.show()


# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────
def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compute Fo-Fo difference maps then run SVD — all in one script.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--ref",      type=Path, metavar="MTZ",
                        help="Reference MTZ (default: smallest trailing index).")
    parser.add_argument("--model",    type=Path, metavar="PDB",
                        help="PDB model (default: auto-detect).")
    parser.add_argument("--dir",      type=Path, default=Path("."), metavar="DIR",
                        help="Working directory (default: current).")
    parser.add_argument("--high-res", type=float, default=None, metavar="A",
                        help="Fixed high-resolution cutoff in Angstroms for ALL maps. "
                             "Default: use the reference MTZ resolution (recommended for SVD).")
    parser.add_argument("--time-step", type=float, default=None, metavar="MS",
                        help="Time interval in ms between consecutive datasets (e.g. 550, 275, 137). "
                             "Used only for SVD plot X-axis labelling. "
                             "Default: auto-detect from filename.")
    parser.add_argument("--skip-svd", action="store_true",
                        help="Skip SVD after computing maps.")
    parser.add_argument("--dry-run",  action="store_true",
                        help="Print plan without running anything.")
    args = parser.parse_args()

    work_dir = args.dir.resolve()
    out_dir  = work_dir / "output_dfo"

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
    if not check_dependencies():
        print()
        print("ERROR: missing dependencies listed above. Aborting.")
        return 1
    print("All required dependencies found.\n")

    print(f"Working directory : {work_dir}")
    print(f"Reference MTZ     : {ref.name}")
    print(f"Model PDB         : {model.name}")
    print(f"Output directory  : {out_dir}")
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

    # Resolution: use --high-res if given, otherwise use the WORST resolution
    # across all MTZs (most conservative common denominator — required for SVD
    # grid consistency and to avoid extrapolating beyond any dataset's actual data).
    if args.high_res is not None:
        reslim = args.high_res
        print(f"Resolution (user-specified)  : {reslim} A")
    else:
        res_values = {}
        for mtz in mtz_files:
            try:
                res_values[mtz.name] = detect_resolution_limit(mtz)
            except RuntimeError as e:
                print(f"WARNING: could not detect resolution for {mtz.name}: {e}")
        if not res_values:
            print("ERROR: could not detect resolution for any MTZ. Use --high-res.")
            return 1
        reslim = max(res_values.values())  # max Angstrom value = worst resolution
        worst_mtz = max(res_values, key=res_values.get)
        print(f"Resolution limits per MTZ:")
        for name, r in sorted(res_values.items(), key=lambda x: x[1]):
            print(f"  {r:.3f} A  {name}")
        print(f"Resolution (worst, used for all) : {reslim} A  <- {worst_mtz}")
    print(f"  -> All difference maps will use this resolution.\n")

    if not args.dry_run:
        out_dir.mkdir(exist_ok=True)

    total_ok = total_skipped = total_failed = 0

    for mtz in mtz_files:
        if mtz == ref:
            print(f"Skipping reference: {mtz.name}")
            continue

        existing = list(out_dir.glob(f"dFo_{mtz.stem}_*.mtz"))
        if existing:
            print(f"Already computed, skipping: {mtz.name}")
            total_skipped += 1
            continue

        print(f"\n=== {mtz.name}  minus  {ref.name} ===")
        try:
            compute_diffmap(mtz, ref, model, labels, reslim, args.dry_run, out_dir)
            total_ok += 1
        except RuntimeError as e:
            print(f"  ERROR: {e}")
            total_failed += 1
        print("-----------------------------------")

    print()
    print("-" * 48)
    print(f"  Computed : {total_ok}")
    print(f"  Skipped  : {total_skipped}")
    print(f"  Failed   : {total_failed}")
    print()

    if not args.skip_svd and not args.dry_run:
        print("Running SVD...")
        try:
            run_svd(out_dir, work_dir / "output_svd", args.time_step)
        except Exception as e:
            print(f"SVD ERROR: {e}")
            return 2
    elif args.dry_run:
        print("[dry-run] would run SVD")

    return 0 if total_failed == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())