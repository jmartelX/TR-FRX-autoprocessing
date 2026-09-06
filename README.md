# TR-FRX autoprocessing

A reproducible command-line pipeline for processing **Time-Resolved Functional
Rotation Crystallography (TR-FRX)** data — from raw diffraction images all the
way to **Fo–Fo difference maps**, their **singular-value decomposition (SVD)**
and an automated **peak analysis**.

TR-FRX captures real-time structural snapshots of ligand binding and enzymatic
catalysis in **single protein crystals at room temperature**, using standard
rotation-based X-ray data collection on conventional beamline infrastructure —
no specialized sample-delivery hardware required. A reaction is triggered by
dispensing nanoliter-scale droplets of ligand or substrate onto a crystal
mounted on a conventional holder, and diffraction is recorded continuously
during crystal rotation. The continuous data collection is then **split into a
series of complete subdatasets, each corresponding to one time point**, spanning
timescales from ms to minutes. (Temporal resolution = the time needed to
collect one complete subdataset, so higher-symmetry space groups allow faster
time points.)

This repository covers the **downstream data-processing** half of TR-FRX: it
turns that series of per-time-point subdatasets into **Fo–Fo / isomorphous
difference maps** — Fo(t) – Fo(reference) — the **SVD components** that summarise
the structural changes across the time series, and a per-timepoint **peak
analysis** with a self-contained HTML report. It wraps and orchestrates
well-established crystallography software (autoPROC, CCP4, PHENIX, PyMOL) with a
small set of Python/shell tools so an entire time-resolved dataset can be
processed with **two commands**.

> Throughout this documentation, **"chunk" = "subdataset" = one time window** —
> a complete sub-dataset extracted from the continuous rotation collection.

See the associated paper for the method itself ([Citing this
work](#citing-this-work)).

---

## Pipeline at a glance

Two scripts take you from raw images to difference maps, SVD and peaks:

| Step | Script | What it does |
|:---:|---|---|
| 0 | `scripts/TR-FRX_autoPROC.sh` | Submits **chunked autoPROC jobs** to SLURM (auto-detects **CBF** or **HDF5**) → `autoPROC_*` chunk folders + consolidated reports |
| 1+ | `scripts/trfrx_full_pipeline.py` | **One command**: copy chunks → clean MTZs → `dimple` the reference → **Fo–Fo difference maps** → **SVD** → **peak analysis** → **HTML report** |

Each script is self-documenting — run it with `-h`/`--help`, or read the
docstring at the top of the file.

---

## Requirements

**External software** (on your `PATH` / loaded as modules):

| Tool | Used by | Provides |
|---|---|---|
| [autoPROC](https://www.globalphasing.com/autoproc/) | step 0 | `process` |
| SLURM | step 0, `--cluster` | `sbatch`, `srun` |
| [CCP4](https://www.ccp4.ac.uk/) | clean + dimple + reindex | `cad`, `mtzdmp`, `dimple`, `pointless` |
| [PHENIX](https://phenix-online.org/) | difference maps | `phenix.fobs_minus_fobs_map` |
| [PyMOL](https://pymol.org/) | peak figures | `pymol` |

**Python ≥ 3.10**, with: `numpy`, `pandas`, `scipy`, `gemmi`, `matplotlib`,
`seaborn`, `dask`, `h5py` (Eiger header reading). `setup_env.sh` installs all
of them.

### Environment setup

A helper loads every external tool as a module and activates the Python
environment. **Source** it once per shell (don't execute it):

```bash
source scripts/setup_env.sh
```

> [!IMPORTANT]
> **Use all scripts in this env.** `setup_env.sh` is needed for the autoPROC step (Step 0), and the pipeline (Step 1+).

Make the scripts executable once after cloning:

```bash
chmod +x scripts/*.py scripts/*.sh
```

---

## Full walkthrough

Throughout, replace `SAMPLE` / paths / templates with your own dataset.

### Step 0 — Submit chunked autoPROC jobs

> Run **from the directory containing your images** (a CBF template or an Eiger
> HDF5 master file — the mode is detected automatically).

Edit the variables at the top of `scripts/TR-FRX_autoPROC.sh` (the exact names
are documented in the script header):

| Variable | Description | Example |
|---|---|---|
| `RUN_ID` | Dataset run number | `006` |
| `IMAGE_TEMPLATE` | Image template (`####` = frame no.) or HDF5 master | `SAMPLE_1_####.cbf` |
| `FIRST_IMG` / `LAST_IMG` | Total frame range | `1` / `3000` |
| `REF_FIRST_IMG` / `REF_LAST_IMG` | Reference chunk range | `1` / `300` |
| `CHUNK_SIZE` | Frames per subsequent chunk | `300` |
| `SYMM` | Space group (edit for your crystal) | `I41` |
| `CELL` | Unit-cell parameters (edit for your crystal) | `114 114 118 90 90 90` |
| `SLURM_PARTITION` | SLURM partition | `nice` |
| `SLURM_CPUS` / `SLURM_MEM` | Resources per job | `24` / `24000` |
| `FRAME_TIME_MS` / `OSC_PER_IMAGE` | Frame period (ms) / oscillation (°); blank = read from the image header | `` / `` |
| `T0_IMAGE` / `T0_SECONDS` | Time origin; blank `T0_IMAGE` = t = 0 at the end of the reference chunk | `` / `0` |

Then:

```bash
cd /data/images/SAMPLE
./path/to/TR-FRX-autoprocessing/scripts/TR-FRX_autoPROC.sh
```

> [!TIP]
> **Separate input / output folders?** Instead of `cd`-ing into the image
> folder, pass them as arguments:
>
> ```bash
> ./TR-FRX_autoPROC.sh <input_dir> <output_dir>
> ```

It submits one SLURM job for the **reference** chunk (e.g. 1–300), then one job
per subsequent chunk (301–600, 601–900, …), each depending on the reference
job. Output lands in `autoproc_chunks/`, and a final job regroups the
statistics into `autoproc_chunks/reports/`:

```
autoproc_chunks/
    autoPROC_1_300/       ← reference
    autoPROC_301_600/
    autoPROC_601_900/
    ...
    reports/              ← consolidated STARANISO / truncate statistics
    reports_qc/           ← QC provenance (qc_fixed.tsv, reindex_standalone.tsv)
```

`reports/` also holds `time_windows.csv`, the real per-dataset acquisition time
that Step 1+ uses for the SVD time axis.

#### Automatic recovery and QC

Two things happen without you asking:

- **Failed `-ref` scaling.** If a chunk fails to scale against the reference
  (typically non-isomorphous drift on an otherwise good window) it is retried
  **standalone**, then reindexed onto the reference with CCP4 `pointless` so the
  Fo−Fo subtraction stays valid. Such windows are listed in
  `reports_qc/reindex_standalone.tsv` and hatched blue in the report PDFs.
- **Outlier `Rmerge`.** Once every chunk has finished, a QC job re-runs any
  window whose overall `Rmerge` is an outlier, excluding dead leading frames and
  keeping the trim with the lowest `Rmerge` at acceptable completeness. The
  superseded run is kept as `autoPROC_<first>_<last>.qc_orig_backup`, the change
  is recorded in `reports_qc/qc_fixed.tsv`, and the window is shaded amber in the
  report PDFs. Tune with the `QC_*` variables; `QC_ENABLE=0` turns it off.

The QC pass verifies its own effect: it compares the number of images used in
the merge before and after each trim, and rejects a trim that did not actually
reduce it (which would mean `EXCLUDE_DATA_RANGE` had no effect in your autoPROC
build). If the count cannot be parsed the trim is still accepted, but
`reports_qc/qc_refix.log` records the exclusion as unverified.

> [!NOTE]
> `.qc_orig_backup` folders are deliberately **skipped** by Step 1+ — they hold
> superseded pre-QC data. The pipeline prints what it skipped.

#### Regenerating the reports only

```bash
./TR-FRX_autoPROC.sh --stats /data/processed/SAMPLE
```

Add `--max-image N` to keep only chunks ending at or before image `N` (a chunk
straddling `N` is excluded, never truncated); those reports get an `_upto<N>`
suffix so the full-series report is never overwritten.

> [!CAUTION]
> **Wait for all jobs to finish before continuing.**

### Step 1+ — The full pipeline (one command)

**USE IN A TERMINAL WITH THE ENV ACTIVATED**

Everything downstream is a single script. Point it at your **model**, the
**autoproc_chunks** folder from Step 0, and an **output** folder:

```bash
python path/to/scripts/trfrx_full_pipeline.py  model.pdb  /data/images/SAMPLE/autoproc_chunks  ./
```

It asks for a **dataset name** and whether to use the **staraniso** or
**truncate** MTZ, then runs, in order:

1. copy the `autoPROC_*` chunks into the output folder,
2. CAD-clean the chosen MTZ (staraniso **or** truncate) for every timepoint,
3. `dimple` the reference timepoint against your model → `final.pdb`
   (its final R is reported; a warning is printed if R > `--max-r`, default 0.4),
4. `phenix.fobs_minus_fobs_map` for every timepoint against the reference,
5. **SVD** of the difference-map series + **peak finding** + PyMOL figures + a
   self-contained **`report.html`**.

**Output layout** — numbered reruns; copy/clean/dimple are shared across reruns:

```
/data/processed/SAMPLE/
    autoproc_copy/                  # copied chunks (created once)
        reports/                    # statistics, damage PDFs, time_windows.csv
        reports_qc/                 # QC provenance, carried over with the chunks
    run_01/
        analysis/                   # cleaned MTZ, final.pdb, output_dfo/, output_svd/
        dimple/                     # dimple run + its logs
        pipeline_<stamp>.log        # full console transcript of this run
        report.html                 # ← open this
        run_config.txt
    run_02/  ...                    # each rerun keeps its own folder
```

**Handy options** (run `scripts/trfrx_full_pipeline.py -h` for the full list):

```bash
--cluster                  # dimple via srun + diffmaps via sbatch (SLURM)
--name D1                  # skip the dataset-name prompt
--file staraniso|truncate  # skip the structure-factor-file prompt
--high-res 2.5             # force ONE uniform resolution (A) for every map + SVD
--dry-run                  # print the plan for every stage, do nothing
--resume                   # reuse the latest run_NN, skip finished work
--start-stage analyze      # re-run only the analysis (retune) into a new run_NN
--ref-timepoint 1_300      # choose which timepoint is the reference / dimpled
--display-sigma 3.0        # peak-figure contour level (sigma)
--display-near 5.0         # residues shown around each peak (A)

```

> [!IMPORTANT]
> **On a cluster:** `--cluster` runs `dimple` on a compute node via `srun`
> (streamed live) and submits the difference maps as parallel `sbatch` jobs; the
> final job writes `report.html` once every timepoint has finished.

### Choosing the map resolution (auto vs uniform)

Before computing the difference maps, the pipeline measures each timepoint's
**honest resolution** — the resolution where the data is still spherically
~90 % complete — instead of trusting the header high-resolution limit (which,
for STARANISO files, is the optimistic *tip* of the anisotropic ellipsoid). It
then asks for a **resolution mode**:

| Mode | How to pick | Every map computed at | Best for |
|---|---|---|---|
| **uniform** | press **Enter** (worst non-outlier), type a number, or `--high-res X` | one shared cutoff | **comparing** peaks / volumes across timepoints (SVD, kinetics) |
| **auto (per-map)** | type **`A`** at the prompt | its own honest resolution | **visualising** the sharpest map at a single timepoint |

> [!IMPORTANT]
> Peak height, width and **integrated volume scale with resolution**, so maps at
> different resolutions are **not comparable across timepoints**. For any
> cross-timepoint or kinetic analysis (peak-volume integration, the SVD) use
> **uniform**. Reserve **auto** for eyeballing individual high-resolution
> timepoints. `Enter` and headless/`--cluster` non-interactive runs default to
> uniform for this reason.

**SVD is always unbiased.** Whatever mode you pick, the SVD truncates every map
to one common resolution and automatically **drops resolution outliers** — a
lone, much-lower-resolution timepoint (e.g. from radiation-damage decay) that
would blur every component — from the decomposition; that timepoint's individual
map is still produced. The `report.html` map-quality table lists each
timepoint's nominal vs honest resolution, the resolution its map was computed at,
whether it entered the SVD, and a short quality verdict
(outlier / anisotropic / weak dataset agreement).

---

## Quick reference

```bash
# ── One shell, environment activated ─────────────────────
source scripts/setup_env.sh

# 0. edit variables, then submit autoPROC from the image folder
cd /data/images/SAMPLE
bash /path/to/scripts/TR-FRX_autoPROC.sh          # → wait for SLURM

# 1+. after all SLURM jobs finish, run the full pipeline
python scripts/trfrx_full_pipeline.py  model.pdb  /data/images/SAMPLE/autoproc_chunks  /data/processed
#     → /data/processed/SAMPLE/run_01/report.html
```

---

## Citing this work

If this pipeline is useful in your research, please cite the associated paper.
A [`CITATION.cff`](CITATION.cff) is included, so GitHub also shows a **"Cite
this repository"** button on the repo page.

> 📄 **Associated paper:**
>
> Martel, J. M. J., Caramello, N., Coquille, S., Mathieu, E., Petit, L.,
> Jacquet, P., Appolaire, A., Leonarski, F., Olieric, V., Wang, M., Madern, D.,
> Royant, A. & Engilberge, S. (2026). *Time-resolved functional rotation
> crystallography reveals protein dynamics and catalysis.* **bioRxiv**.
> doi: [10.64898/2026.04.24.718481](https://doi.org/10.64898/2026.04.24.718481)

---

## License

Released under the [MIT License](LICENSE) — free to use, modify and
redistribute with attribution.

## Acknowledgements

Parts of this codebase were developed with the assistance of
[Claude Code](https://claude.com/claude-code) (Anthropic).

## Contact

Julien Martel · [@jmartelX](https://github.com/jmartelX)
