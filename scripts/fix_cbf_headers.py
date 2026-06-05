#!/usr/bin/env python3
"""
Fix broken CBF headers for PILATUS3 6M data from ESRF BM07.

Two cases are handled:
  Case A: file has _Misc.* block AND _array_data.header_contents
          → just remove the _Misc.* block (the metadata is already there)
  Case B: file has _Misc.* block but NO _array_data.header_contents
          → remove the _Misc.* block AND inject a reconstructed header
          (requires experiment geometry parameters below)

Usage:
    # Preview without modifying files
    python fix_cbf_headers.py --dry-run --pattern "dataset_1_####.cbf" --first 1 --last 3000

    # Apply fix
    python fix_cbf_headers.py --pattern "dataset_1_####.cbf" --first 1 --last 3000
"""

import os
import re
import sys
import argparse


# ─────────────────────────────────────────────────────────────────────────────
# EXPERIMENT PARAMETERS  (fill in from a working CBF header)
# ─────────────────────────────────────────────────────────────────────────────
DETECTOR        = "PILATUS3 6M, S/N 60-0128, ESRF BM07"
PIXEL_SIZE_M    = 172e-6          # metres
SENSOR_M        = 0.001000        # metres  (1 mm silicon)
OVERLOAD        = 1048500         # Count_cutoff
NX              = 2463
NY              = 2527
POLARIZATION    = 0.99

# These vary per dataset — update before running:
WAVELENGTH_A    = 0.979470        # Angstroms
DISTANCE_M      = 0.389991        # metres
BEAM_X_PIX      = 1223.31         # pixels
BEAM_Y_PIX      = 1289.47         # pixels
ANGLE_INCREMENT = 0.1000          # degrees per image  ← set for your dataset
START_ANGLE_IMG1 = 0.0            # start angle of image #1 (degrees)
FIRST_IMAGE_NUM  = 1              # number of the first image in the series
EXPOSURE_TIME   = 0.050000        # seconds
EXPOSURE_PERIOD = 0.050950        # seconds
THRESHOLD_EV    = 6329
# ─────────────────────────────────────────────────────────────────────────────


HEADER_TEMPLATE = """\
# Detector: {detector}
# {date}

# Silicon sensor, thickness {sensor_m:.6f} m

# Pixel_size {px_m} m x {px_m} m
# N_oscillations 1
# Oscillation_axis omega
# Chi 0.0000 deg.
# Phi -9999.0000 deg.
# Kappa -9999.0000 deg.
# Alpha 0.0000 deg.
# Polarization {pol}
# Detector_2theta 0.0000 deg.
# Angle_increment {angle_inc:.4f} deg.
# Beam_xy ({beam_x:.2f}, {beam_y:.2f}) pixels
# Detector_Voffset 0.0000 m
# Detector_distance {dist_m:.6f} m
# Wavelength {wl:.6f} A
# Flat_field: (nil)
# Excluded_pixels:  badpix_mask.tif
# Threshold_setting {threshold} eV
# Count_cutoff {overload}
# Tau = 0 s
# Exposure_period {exp_period:.6f} s
# Exposure_time {exp_time:.6f} s
# Start_angle {start_angle:.4f} deg.
"""


def image_number_from_path(path):
    """Extract the trailing image number from a filename like dataset_1_0042.cbf."""
    stem = os.path.splitext(os.path.basename(path))[0]
    m = re.search(r"(\d+)$", stem)
    if m:
        return int(m.group(1))
    return None


def start_angle_for_image(img_num):
    offset = img_num - FIRST_IMAGE_NUM
    return START_ANGLE_IMG1 + offset * ANGLE_INCREMENT


def build_header(img_num):
    import datetime
    date_str = datetime.datetime.now().strftime("%Y/%b/%d %H:%M:%S")
    px_str = f"{PIXEL_SIZE_M}".rstrip("0")  # e.g. "0.000172"
    return HEADER_TEMPLATE.format(
        detector=DETECTOR,
        date=date_str,
        sensor_m=SENSOR_M,
        px_m=px_str,
        pol=POLARIZATION,
        angle_inc=ANGLE_INCREMENT,
        beam_x=BEAM_X_PIX,
        beam_y=BEAM_Y_PIX,
        dist_m=DISTANCE_M,
        wl=WAVELENGTH_A,
        threshold=THRESHOLD_EV,
        overload=OVERLOAD,
        exp_period=EXPOSURE_PERIOD,
        exp_time=EXPOSURE_TIME,
        start_angle=start_angle_for_image(img_num),
    )


def fix_cbf(filepath, dry_run=False):
    with open(filepath, "rb") as f:
        data = f.read()

    # Nothing to do if no _Misc block
    if b"_Misc." not in data[:4096]:
        return "clean"

    has_header_contents = b"_array_data.header_contents" in data[:8192]

    # ── Locate the CIF data block name (first line, e.g. "data_image_0") ──
    first_newline = data.find(b"\n")
    if first_newline == -1:
        return "skip"
    block_line = data[: first_newline + 1]

    if has_header_contents:
        # Case A: strip _Misc block, keep existing header_contents
        marker = b"_array_data.header_contents"
        pos = data.find(marker)
        new_data = block_line + b"\n" + data[pos:]
    else:
        # Case B: strip _Misc block, inject reconstructed header
        marker = b"_array_data.data"
        pos = data.find(marker)
        if pos == -1:
            return "skip"

        img_num = image_number_from_path(filepath)
        if img_num is None:
            print(f"  WARNING: cannot extract image number from {filepath}, using image 1")
            img_num = FIRST_IMAGE_NUM

        header_txt = build_header(img_num)
        header_block = (
            b"_array_data.header_contents\n"
            b";\n"
            + header_txt.encode("ascii")
            + b";\n\n"
        )
        new_data = block_line + b"\n" + header_block + data[pos:]

    if new_data == data:
        return "clean"

    case = "A" if has_header_contents else "B"
    removed = len(data) - len(new_data)

    if dry_run:
        print(f"  DRY-RUN [{case}] {filepath}  ({removed:+d} bytes)")
        return "would_fix"

    with open(filepath, "wb") as f:
        f.write(new_data)
    return case


def expand_pattern(pattern, first, last):
    digits = pattern.count("#")
    fmt = pattern.replace("#" * digits, f"{{:0{digits}d}}")
    return [fmt.format(i) for i in range(first, last + 1)]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("files", nargs="*", help="CBF files to fix")
    parser.add_argument("--pattern", help="CBF template with #### (e.g. data_1_####.cbf)")
    parser.add_argument("--first", type=int, default=1)
    parser.add_argument("--last",  type=int, default=1)
    parser.add_argument("--dry-run", action="store_true",
                        help="Show what would be done without modifying files")
    args = parser.parse_args()

    if args.pattern:
        files = expand_pattern(args.pattern, args.first, args.last)
    elif args.files:
        files = args.files
    else:
        parser.error("Provide file arguments or --pattern/--first/--last")

    counts = {"A": 0, "B": 0, "clean": 0, "skip": 0, "missing": 0, "would_fix": 0}

    for path in files:
        if not os.path.exists(path):
            counts["missing"] += 1
            continue
        result = fix_cbf(path, dry_run=args.dry_run)
        counts[result] = counts.get(result, 0) + 1
        if result in ("A", "B"):
            print(f"  FIXED [case {result}] {path}")

    print()
    if args.dry_run:
        print(f"Would fix: {counts['would_fix']} file(s)")
    else:
        print(f"Fixed case A (had header_contents): {counts['A']}")
        print(f"Fixed case B (missing header_contents): {counts['B']}")
    print(f"Already clean: {counts['clean']}  |  Skipped: {counts['skip']}  |  Not found: {counts['missing']}")


if __name__ == "__main__":
    main()
