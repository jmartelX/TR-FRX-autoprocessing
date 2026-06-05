#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
copy_filesV3.py — Copy selected files from autoPROC_* subfolders.

Each autoPROC_XXX folder is copied into a parent folder named XXX inside
the destination. For example:
    autoPROC_1_300/   → <dest>/<source>_copy/1_300/
    autoPROC_301_600/ → <dest>/<source>_copy/301_600/

Usage:
    ./copy_filesV3.py <source> <destination> [options]

Examples:
    ./copy_filesV3.py /data/run1 /backup
    ./copy_filesV3.py /data/run1 /backup --dry-run
    ./copy_filesV3.py /data/run1 /backup --workers 8
    ./copy_filesV3.py /data/run1 /backup --log-file copy.log --verbose
"""

from __future__ import annotations

import argparse
import logging
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path

# ─────────────────────────────────────────────
# Hardcoded subfolder prefix
# ─────────────────────────────────────────────
SUBFOLDER_PREFIX = "autoPROC_"

# ─────────────────────────────────────────────
# Default list of files to search for
# ─────────────────────────────────────────────
DEFAULT_FILES = [
    "staraniso_alldata-unique.mtz",
    "summary.html",
    "summary.tar.gz",
    "truncate-unique.mtz",
    "XDS_ASCII.HKL",
    "truncate.log",
]


# ─────────────────────────────────────────────
# Result tracking
# ─────────────────────────────────────────────
@dataclass
class CopyResult:
    subfolder: str
    copied: list[str] = field(default_factory=list)
    missing: list[str] = field(default_factory=list)


# ─────────────────────────────────────────────
# Core logic
# ─────────────────────────────────────────────
def process_subfolder(
    subfolder: Path,
    dest_root: Path,
    files_to_copy: list[str],
    dry_run: bool,
) -> CopyResult:
    """Copy selected files from one subfolder into the destination.

    The autoPROC_ prefix is stripped from the subfolder name to form the
    destination directory, e.g. autoPROC_1_300 -> dest_root/1_300/.
    """
    result = CopyResult(subfolder=subfolder.name)
    range_name = subfolder.name.removeprefix(SUBFOLDER_PREFIX)
    dest_dir = dest_root / range_name / subfolder.name

    if not dry_run:
        dest_dir.mkdir(parents=True, exist_ok=True)

    for filename in files_to_copy:
        src = subfolder / filename
        if src.is_file():
            if not dry_run:
                shutil.copy2(src, dest_dir / filename)
            result.copied.append(filename)
            logging.debug("  [%s] %s → copied", subfolder.name, filename)
        else:
            result.missing.append(filename)
            logging.debug("  [%s] %s → MISSING", subfolder.name, filename)

    return result


def copy_structure(
    source: Path,
    output: Path,
    files_to_copy: list[str],
    dry_run: bool,
    workers: int,
) -> None:
    # -- Resolve destination --
    # Use "_copy" suffix only when outputting into the same directory as the source,
    # to avoid a name collision. Otherwise just reuse the source folder name.
    same_parent = output.resolve() == source.parent.resolve()
    folder_name = f"{source.name}_copy" if same_parent else source.name
    dest_root = output / folder_name
    if dry_run:
        logging.info("Dry-run mode -- nothing will be written.")
    logging.info("Source      : %s", source)
    logging.info("Destination : %s", dest_root)
    logging.info("Prefix      : %s*", SUBFOLDER_PREFIX)
    logging.info("Workers     : %d", workers)

    if not dry_run:
        dest_root.mkdir(parents=True, exist_ok=True)

    # -- Collect matching subfolders --
    subfolders = sorted(
        p for p in source.iterdir()
        if p.is_dir() and p.name.startswith(SUBFOLDER_PREFIX)
    )

    if not subfolders:
        logging.warning("No subfolders found matching prefix '%s*'.", SUBFOLDER_PREFIX)
        return

    logging.info("Found %d matching subfolder(s).\n", len(subfolders))

    # ── Process (possibly in parallel) ───────
    results: list[CopyResult] = []

    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(process_subfolder, sf, dest_root, files_to_copy, dry_run): sf
            for sf in subfolders
        }
        for future in as_completed(futures):
            try:
                results.append(future.result())
            except Exception as exc:
                logging.error("Unexpected error processing %s: %s", futures[future], exc)

    # Sort results by subfolder name for deterministic output
    results.sort(key=lambda r: r.subfolder)

    # ── Summary ──────────────────────────────
    total_copied = sum(len(r.copied) for r in results)
    total_missing = sum(len(r.missing) for r in results)

    print()
    if dry_run:
        print("── DRY RUN SUMMARY ─────────────────────────────")
    else:
        print(f"✅ Copy complete → {dest_root}")
        print("── SUMMARY ─────────────────────────────────────")

    print(f"   Subfolders processed : {len(results)}")
    print(f"   Files copied         : {total_copied}")
    print(f"   Files missing        : {total_missing}")

    if total_missing:
        print("\n❌ Missing files:")
        for r in results:
            for f in r.missing:
                print(f"   {r.subfolder} → {f}")
    else:
        print("\n🎉 All files were found and copied successfully.")


# ─────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Copy selected files from autoPROC_* subfolders. "
            "Each autoPROC_XXX folder is placed under a parent named XXX."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("source", type=Path, help="Source directory.")
    parser.add_argument("destination", type=Path, help="Output parent directory.")
    parser.add_argument(
        "--files",
        nargs="+",
        metavar="FILE",
        default=DEFAULT_FILES,
        help="Override the list of filenames to copy.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Simulate the copy without writing anything.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        metavar="N",
        help="Number of parallel copy threads. (default: 4)",
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable DEBUG-level logging.",
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        metavar="PATH",
        help="Also write log output to this file.",
    )
    return parser


def setup_logging(verbose: bool, log_file: Path | None) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    handlers: list[logging.Handler] = [logging.StreamHandler(sys.stdout)]
    if log_file:
        handlers.append(logging.FileHandler(log_file, encoding="utf-8"))
    logging.basicConfig(
        level=level,
        format="%(levelname)-8s %(message)s",
        handlers=handlers,
    )


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    setup_logging(args.verbose, args.log_file)

    # Validate paths
    if not args.source.is_dir():
        logging.error("Source directory not found: %s", args.source)
        sys.exit(1)
    if not args.destination.is_dir():
        logging.error("Destination directory not found: %s", args.destination)
        sys.exit(1)

    copy_structure(
        source=args.source.resolve(),
        output=args.destination.resolve(),
        files_to_copy=args.files,
        dry_run=args.dry_run,
        workers=args.workers,
    )


if __name__ == "__main__":
    main()
