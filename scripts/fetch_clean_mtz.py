#!/usr/bin/env python3
"""
fetch_clean_mtz.py

Usage:
  python3 fetch_clean_mtz.py /path/to/CaMDH_* /path/to/destination/all_mtz
  python3 fetch_clean_mtz.py ../CaMDH_047 /tmp/all_mtz

What it does:
- Finds staraniso_alldata-unique.mtz inside autoPROC_* folders
- Writes a cleaned MTZ per range/numeric folder:
    <CaMDH_basename>_<folderName>.mtz

Cleaning (CCP4 CAD):
- Keeps only: H K L F SIGF FreeR_flag
- Output labels are renamed to exactly: F SIGF FreeR_flag
- HKL are always carried by CAD automatically.

Prereqs:
- CCP4 programs in PATH: mtzdmp, cad
"""

from __future__ import annotations

import re
import sys
import shutil
import tempfile
import subprocess
from pathlib import Path
from typing import Optional


TARGET_MTZ_NAME = "staraniso_alldata-unique.mtz"
RANGE_DIR_RE = re.compile(r"^(\d+)[\s_-]+(\d+)$")  # accepts 1-300, 1_300, 1 - 300

CAD_EXE = "cad"
MTZDUMP_EXE = "mtzdmp"


def is_glob_pattern(s: str) -> bool:
    return any(ch in s for ch in ["*", "?", "["])


def natural_int(s: str, default: int = 10**12) -> int:
    m = re.search(r"\d+", s)
    return int(m.group(0)) if m else default


def expand_input(pattern_or_dir: str) -> list[Path]:
    p = Path(pattern_or_dir).expanduser()

    if is_glob_pattern(pattern_or_dir):
        import glob
        matches = [Path(m).expanduser() for m in glob.glob(pattern_or_dir)]
        folders = [m for m in matches if m.is_dir()]
        return sorted(folders, key=lambda x: x.name)

    if p.is_dir():
        return [p.resolve()]

    return []


def iter_range_subfolders(camh_dir: Path) -> list[Path]:
    ranged = []
    for child in camh_dir.iterdir():
        if not child.is_dir():
            continue
        m = RANGE_DIR_RE.match(child.name)
        if m:
            start = int(m.group(1))
            ranged.append((start, child))
    ranged.sort(key=lambda t: t[0])
    return [p for _, p in ranged]


def iter_numeric_subfolders(camh_dir: Path) -> list[Path]:
    numeric_dirs = []
    for child in camh_dir.iterdir():
        if child.is_dir() and child.name.isdigit():
            numeric_dirs.append(child)
    numeric_dirs.sort(key=lambda x: int(x.name))
    return numeric_dirs


def find_mtz_in_folder(folder: Path) -> Optional[Path]:
    autoproc_dirs = [d for d in folder.iterdir() if d.is_dir() and d.name.startswith("autoPROC_")]

    candidates: list[tuple[int, Path]] = []
    for ap in autoproc_dirs:
        mtz = ap / TARGET_MTZ_NAME
        if mtz.is_file():
            candidates.append((natural_int(ap.name), mtz))

    if not candidates:
        return None

    candidates.sort(key=lambda t: t[0])
    return candidates[0][1]


def _run(cmd: list[str], *, input_text: str | None = None) -> subprocess.CompletedProcess:
    try:
        return subprocess.run(
            cmd,
            input=input_text,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
        )
    except FileNotFoundError as e:
        raise RuntimeError(
            f"Required program not found: {cmd[0]}\n"
            f"Make sure CCP4 is set up and {cmd[0]} is in your PATH."
        ) from e
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"Command failed: {' '.join(cmd)}\n"
            f"Exit code: {e.returncode}\n"
            f"--- stdout ---\n{e.stdout}\n"
            f"--- stderr ---\n{e.stderr}\n"
        ) from e


def mtz_column_labels_from_mtzdmp(mtz_path: Path) -> list[str]:
    """
    Parse the '* Column Labels :' block from mtzdmp output.

    Your mtzdmp output is HTML-wrapped, but the text block contains:
      * Column Labels :
        H K L ... FreeR_flag

    This parser:
    - reads stdout+stderr
    - finds the '* Column Labels' line
    - consumes subsequent lines until the next '* ' section header
    """
    cp = _run([MTZDUMP_EXE, str(mtz_path)])
    raw = (cp.stdout or "") + "\n" + (cp.stderr or "")
    lines = raw.splitlines()

    labels: list[str] = []
    in_block = False

    for line in lines:
        s = line.strip()

        # Start of block
        if s.startswith("* Column Labels"):
            in_block = True
            # Might have labels on same line after ':'
            if ":" in s:
                after = s.split(":", 1)[1].strip()
                if after:
                    labels.extend(after.split())
            continue

        if in_block:
            # Stop at next section header
            if s.startswith("* ") and not s.startswith("* Column Labels"):
                break
            if not s:
                continue
            # Skip HTML-ish tags if any sneak in
            if s.startswith("<") and s.endswith(">"):
                continue
            # Accumulate tokens
            labels.extend(s.split())

    # De-dup preserve order
    seen = set()
    uniq = []
    for lab in labels:
        if lab not in seen:
            uniq.append(lab)
            seen.add(lab)

    return uniq


def clean_mtz_with_cad(
    src_mtz: Path,
    out_mtz: Path,
    cached_labels: list[str] | None = None,
) -> list[str]:
    """Clean an MTZ with CAD. Returns the labels used (for caching by the caller).

    If *cached_labels* are provided they are tried first; if the required
    columns are missing the function automatically falls back to running
    mtzdmp on this specific file.
    """
    required = {"F", "SIGF", "FreeR_flag"}

    labels = cached_labels
    if labels is not None and not required.issubset(set(labels)):
        # Cached labels don't have what we need — read this file directly
        print(f"  Note: cached labels incomplete for {src_mtz.name}, re-reading with mtzdmp.")
        labels = None

    if labels is None:
        labels = mtz_column_labels_from_mtzdmp(src_mtz)

    missing = [r for r in sorted(required) if r not in set(labels)]
    if missing:
        raise RuntimeError(
            f"Cannot clean {src_mtz} because required columns were not found: {', '.join(missing)}\n"
            f"Available labels: {', '.join(labels) if labels else '(none parsed)'}"
        )

    cad_in = (
        "LABIN FILE 1 E1=F E2=SIGF E3=FreeR_flag\n"
        "LABOUT FILE 1 E1=F E2=SIGF E3=FreeR_flag\n"
        "END\n"
    )

    _run([CAD_EXE, "HKLIN1", str(src_mtz), "HKLOUT", str(out_mtz)], input_text=cad_in)
    return labels


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print(f"Usage: {Path(argv[0]).name} /path/to/CaMDH_*  /path/to/destination/all_mtz")
        return 1

    input_arg = argv[1]
    dest_arg = argv[2]

    camh_folders = expand_input(input_arg)
    if not camh_folders:
        print(f"ERROR: No input folders found for: {input_arg}")
        return 1

    dest_dir = Path(dest_arg).expanduser().resolve()
    dest_dir.mkdir(parents=True, exist_ok=True)

    print("Destination folder:")
    print(f"  {dest_dir}")
    print()

    total_written = 0
    total_missing = 0
    total_failed = 0

    # Cache column labels after the first successful mtzdmp call so we don't
    # re-run mtzdmp for every file (they all come from the same autoPROC output
    # and share identical column layouts).  Falls back to per-file detection
    # automatically if the cached labels don't match a particular file.
    _cached_labels: list[str] | None = None

    for camh in camh_folders:
        camh = camh.resolve()
        camh_name = camh.name

        subdirs = iter_range_subfolders(camh)
        if not subdirs:
            subdirs = iter_numeric_subfolders(camh)

        if not subdirs:
            print(f"WARNING: No range folders (1-250,...) or numeric folders (1,2,...) found in: {camh}")
            continue

        for sd in subdirs:
            label = sd.name.strip().replace(" ", "")
            mtz_path = find_mtz_in_folder(sd)

            if mtz_path is None:
                print(f"WARNING: Missing {TARGET_MTZ_NAME} in {sd} (no autoPROC_* contains it)")
                total_missing += 1
                continue

            out_name = f"{camh_name}_{label}.mtz"
            out_path = dest_dir / out_name

            try:
                with tempfile.TemporaryDirectory(prefix="mtz_clean_") as td:
                    tmp_out = Path(td) / out_name
                    used_labels = clean_mtz_with_cad(mtz_path, tmp_out, cached_labels=_cached_labels)
                    if _cached_labels is None:
                        _cached_labels = used_labels
                        print(f"  Labels cached from first file: {' '.join(used_labels)}")
                    shutil.copy2(tmp_out, out_path)

                print(f"Cleaned: {mtz_path} -> {out_path}  (HKL + F SIGF FreeR_flag)")
                total_written += 1
            except Exception as e:
                print(f"ERROR: Failed to clean {mtz_path}")
                print(f"  Reason: {e}")
                total_failed += 1

    print()
    print("Done!")
    print(f"Cleaned+written: {total_written}")
    print(f"Missing source MTZ: {total_missing}")
    print(f"Failed cleaning: {total_failed}")
    print(f"Output: {dest_dir}")
    return 0 if total_failed == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))