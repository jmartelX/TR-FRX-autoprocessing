#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
diffmaps.py — Fo-Fo difference map calculation (batch or single-pair).

BATCH MODE (default):
    Runs phenix.fobs_minus_fobs_map for every MTZ in a directory against
    the one with the smallest trailing index (the reference).

    Supported MTZ trailing-index patterns:
        _N.mtz          e.g. CaMDH_073_137_2.mtz
        _start-end.mtz  e.g. CaMDH_073_275_1-250.mtz
        _start_end.mtz  e.g. CaMDH_073_275_1_300.mtz

    PDB selection:
        - one  .pdb found  → use it
        - many .pdb found  → use model.pdb

    Usage:
        ./diffmaps.py                             # run in current directory
        ./diffmaps.py --dir /path/to/mtz          # specify directory
        ./diffmaps.py --high-res 1.8 --low-res 8.0
        ./diffmaps.py --dry-run                   # print plan without running

SIMPLE MODE (--simple):
    Compute a single diff map between two explicit MTZ files.

    Usage:
        ./diffmaps.py --simple model.pdb ref.mtz final.mtz
        ./diffmaps.py --simple model.pdb ref.mtz final.mtz /output/dir
        ./diffmaps.py --simple model.pdb ref.mtz final.mtz --high-res 1.8
        ./diffmaps.py --simple model.pdb ref.mtz final.mtz --dry-run
"""

from __future__ import annotations

import re
import sys
import subprocess
import argparse
from pathlib import Path

# ─────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────
HIGH_RES_DEFAULT = 2.0
LOW_RES_DEFAULT  = 10.0

# Trailing-index patterns, tried in order. Group 1 is always the sort key.
_INDEX_PATTERNS = [
    re.compile(r"_(\d+)-(\d+)\.mtz$"),   # _start-end.mtz
    re.compile(r"_(\d+)_(\d+)\.mtz$"),   # _start_end.mtz
    re.compile(r"_(\d+)\.mtz$"),          # _N.mtz
]


# ─────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────
def mtz_sort_index(mtz: Path) -> int | None:
    """Return the leading numeric sort key for an MTZ file, or None if unrecognised."""
    name = mtz.name
    for pat in _INDEX_PATTERNS:
        m = pat.search(name)
        if m:
            return int(m.group(1))
    return None


def _run_silent(cmd: list[str], input_text: str = "") -> str:
    """Run a command, return combined stdout+stderr. Never raises."""
    try:
        cp = subprocess.run(
            cmd, input=input_text, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        return cp.stdout + cp.stderr
    except FileNotFoundError:
        return ""


def detect_fsig_labels(mtz: Path) -> str:
    """
    Inspect an MTZ and return 'FP,SIGFP' or 'F,SIGF'.
    Tries mtzdmp first, then phenix.mtz.dump.
    Raises RuntimeError if neither pair is found.
    """
    import shutil

    if shutil.which("mtzdmp"):
        txt = _run_silent(["mtzdmp", str(mtz)])
    else:
        txt = _run_silent(["phenix.mtz.dump", str(mtz)])

    words = set(txt.split())
    if "FP" in words and "SIGFP" in words:
        return "FP,SIGFP"
    if "F" in words and "SIGF" in words:
        return "F,SIGF"

    raise RuntimeError(
        f"Could not find FP/SIGFP or F/SIGF in: {mtz}\n"
        f"Available tokens (first 80): {' '.join(list(words)[:80])}"
    )


def find_pdb(directory: Path) -> Path:
    """Return the PDB file to use, applying the one-or-model.pdb rule."""
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


def collect_mtz_files(directory: Path) -> list[Path]:
    """Return MTZ files sorted by trailing index, smallest first."""
    all_mtz = list(directory.glob("*.mtz"))
    if not all_mtz:
        raise RuntimeError("No .mtz files found in directory.")

    indexed = []
    unrecognised = []
    for mtz in all_mtz:
        idx = mtz_sort_index(mtz)
        if idx is not None:
            indexed.append((idx, mtz))
        else:
            unrecognised.append(mtz.name)

    if not indexed:
        raise RuntimeError(
            "No MTZ files with a recognised trailing index (_N, _start-end, _start_end) found.\n"
            f"Files present: {', '.join(m.name for m in all_mtz)}"
        )

    if unrecognised:
        print(f"Note: ignoring {len(unrecognised)} unrecognised MTZ file(s): {', '.join(unrecognised)}")

    indexed.sort(key=lambda t: t[0])
    return [mtz for _, mtz in indexed]


def write_eff(
    eff_path: Path,
    mtz_other: Path,
    labels_other: str,
    mtz_ref: Path,
    labels_ref: str,
    pdb: Path,
    high_res: float,
    low_res: float,
) -> None:
    eff_path.write_text(
        f"f_obs_1_file_name = {mtz_other}\n"
        f"f_obs_1_label = {labels_other}\n"
        f"\n"
        f"f_obs_2_file_name = {mtz_ref}\n"
        f"f_obs_2_label = {labels_ref}\n"
        f"\n"
        f"high_resolution = {high_res}\n"
        f"low_resolution = {low_res}\n"
        f"\n"
        f"sigma_cutoff = 3.0\n"
        f"phase_source = {pdb}\n"
        f"ignore_non_isomorphous_unit_cells = True\n"
        f"\n"
        f"advanced {{\n"
        f"  multiscale = True\n"
        f"}}\n"
    )


# ─────────────────────────────────────────────
# Single-pair runner (shared by both modes)
# ─────────────────────────────────────────────
def run_single_pair(
    pdb: Path,
    ref_mtz: Path,
    other_mtz: Path,
    out_dir: Path,
    high_res: float,
    low_res: float,
    dry_run: bool,
    known_ref_labels: str | None = None,
    known_other_labels: str | None = None,
) -> int:
    """
    Compute one dFo map (other_mtz minus ref_mtz).
    Returns 0 on success, 1 on failure.

    *known_ref_labels* and *known_other_labels* can be pre-supplied by the
    caller to skip redundant mtzdmp/phenix.mtz.dump calls when processing
    many files with identical column layouts.
    """
    try:
        ref_labels   = known_ref_labels   or detect_fsig_labels(ref_mtz)
        other_labels = known_other_labels or detect_fsig_labels(other_mtz)
    except RuntimeError as e:
        print(f"ERROR detecting labels: {e}")
        return 1

    print(f"Reference MTZ  : {ref_mtz.name}  [{ref_labels}]")
    print(f"Other MTZ      : {other_mtz.name}  [{other_labels}]")
    print(f"PDB            : {pdb.name}")
    print(f"Resolution     : {low_res} – {high_res} Å")
    print(f"Output dir     : {out_dir}")
    if dry_run:
        print("DRY RUN — nothing will be written or executed.")
        return 0

    out_dir.mkdir(parents=True, exist_ok=True)

    prefix   = out_dir / f"dFo_{other_mtz.stem}-{ref_mtz.stem}"
    eff_file = prefix.with_suffix(".eff")
    log_file = prefix.with_suffix(".log")

    write_eff(
        eff_path=eff_file,
        mtz_other=other_mtz,
        labels_other=other_labels,
        mtz_ref=ref_mtz,
        labels_ref=ref_labels,
        pdb=pdb,
        high_res=high_res,
        low_res=low_res,
    )

    try:
        with log_file.open("w") as log:
            subprocess.run(
                ["phenix.fobs_minus_fobs_map", str(eff_file),
                 f"file_name_prefix={prefix}", "job_id=1"],
                stdout=log, stderr=log, check=True,
            )
        print(f"  MTZ : {prefix}_1.mtz")
        print(f"  MAP : {prefix}_1.map")
        print(f"  LOG : {log_file.name}")
        return 0
    except subprocess.CalledProcessError:
        print(f"  ERROR: phenix failed — see {log_file.name}")
        return 1
    except FileNotFoundError:
        print("  ERROR: phenix.fobs_minus_fobs_map not found in PATH.")
        return 1


# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────
def main() -> int:
    parser = argparse.ArgumentParser(
        description="Fo-Fo difference maps via phenix.fobs_minus_fobs_map.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # ── Simple mode ──────────────────────────
    parser.add_argument(
        "--simple", nargs="+", metavar="ARG",
        help=(
            "Single-pair mode. Arguments: model.pdb ref.mtz final.mtz [outdir]\n"
            "  model.pdb  — phase source\n"
            "  ref.mtz    — reference dataset (f_obs_2)\n"
            "  final.mtz  — dataset to compare (f_obs_1)\n"
            "  outdir     — output directory (default: ./output)"
        ),
    )

    # ── Batch mode ───────────────────────────
    parser.add_argument(
        "--dir", type=Path, default=Path("."),
        metavar="DIR",
        help="[batch] Directory containing MTZ (and PDB) files. (default: .)",
    )

    # ── Shared ───────────────────────────────
    parser.add_argument("--high-res", type=float, default=HIGH_RES_DEFAULT, metavar="Å")
    parser.add_argument("--low-res",  type=float, default=LOW_RES_DEFAULT,  metavar="Å")
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print the plan without running phenix.",
    )

    args = parser.parse_args()

    # ── Dispatch ────────────────────────────
    if args.simple is not None:
        return _run_simple(args)
    else:
        return _run_batch(args)


def _run_simple(args) -> int:
    """Handle --simple model.pdb ref.mtz final.mtz [outdir]."""
    positional = args.simple
    if len(positional) < 3 or len(positional) > 4:
        print(
            "ERROR: --simple requires 3 or 4 arguments:\n"
            "  model.pdb  ref.mtz  final.mtz  [outdir]"
        )
        return 1

    pdb_path   = Path(positional[0]).resolve()
    ref_mtz    = Path(positional[1]).resolve()
    other_mtz  = Path(positional[2]).resolve()
    out_dir    = Path(positional[3]).resolve() if len(positional) == 4 else Path("output").resolve()

    for label, path in [("PDB", pdb_path), ("ref MTZ", ref_mtz), ("final MTZ", other_mtz)]:
        if not path.exists():
            print(f"ERROR: {label} not found: {path}")
            return 1

    return run_single_pair(
        pdb=pdb_path,
        ref_mtz=ref_mtz,
        other_mtz=other_mtz,
        out_dir=out_dir,
        high_res=args.high_res,
        low_res=args.low_res,
        dry_run=args.dry_run,
    )


def _run_batch(args) -> int:
    """Handle batch mode (original behaviour)."""
    work_dir = args.dir.resolve()

    try:
        mtz_files = collect_mtz_files(work_dir)
        pdb       = find_pdb(work_dir)
    except RuntimeError as e:
        print(f"ERROR: {e}")
        return 1

    ref_mtz = mtz_files[0]

    print(f"MTZ files detected ({len(mtz_files)}):")
    for mtz in mtz_files:
        marker = " ← reference" if mtz == ref_mtz else ""
        print(f"  {mtz.name}{marker}")
    print()

    try:
        ref_labels = detect_fsig_labels(ref_mtz)
    except RuntimeError as e:
        print(f"ERROR: {e}")
        return 1

    out_dir = work_dir / "output"
    if not args.dry_run:
        out_dir.mkdir(exist_ok=True)

    print(f"Reference MTZ    : {ref_mtz.name}")
    print(f"Reference labels : {ref_labels}")
    print(f"PDB              : {pdb.name}")
    print(f"Resolution       : {args.low_res} – {args.high_res} Å")
    print(f"Output directory : {out_dir}")
    if args.dry_run:
        print("DRY RUN — nothing will be written or executed.")
    print()

    total_ok = 0
    total_failed = 0

    # Cache other-file labels after the first detection so we don't re-run
    # mtzdmp/phenix.mtz.dump for every file. All files from the same autoPROC
    # batch share identical column layouts. Falls back to per-file detection
    # automatically (known_other_labels=None) until the first file is processed.
    _cached_other_labels: str | None = None

    for mtz in mtz_files:
        if mtz == ref_mtz:
            print(f"Skipping reference: {mtz.name}")
            continue

        print(f"=== {mtz.name}  minus  {ref_mtz.name} ===")

        rc = run_single_pair(
            pdb=pdb,
            ref_mtz=ref_mtz,
            other_mtz=mtz,
            out_dir=out_dir,
            high_res=args.high_res,
            low_res=args.low_res,
            dry_run=args.dry_run,
            known_ref_labels=ref_labels,
            known_other_labels=_cached_other_labels,
        )

        # After the first successful run, cache the other-file labels.
        if rc == 0 and _cached_other_labels is None:
            try:
                _cached_other_labels = detect_fsig_labels(mtz)
                print(f"  (labels cached from {mtz.name}: {_cached_other_labels})")
            except RuntimeError:
                pass  # non-fatal — will re-detect next time
        if rc == 0:
            total_ok += 1
        else:
            total_failed += 1
            if rc == 1:
                return 1  # phenix not in PATH — abort early

        print("-----------------------------------")

    print()
    print("All done!")
    print(f"  Completed : {total_ok}")
    print(f"  Failed    : {total_failed}")
    print(f"  Output    : {out_dir}")
    return 0 if total_failed == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())