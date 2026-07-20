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
| [CCP4](https://www.ccp4.ac.uk/) | clean + dimple | `cad`, `mtzdmp`, `dimple` |
| [PHENIX](https://phenix-online.org/) | difference maps | `phenix.fobs_minus_fobs_map` |
| [PyMOL](https://pymol.org/) | peak figures | `pymol` |

**Python ≥ 3.10**, with: `numpy`, `pandas`, `scipy`, `gemmi`, `matplotlib`,
`seaborn`, `dask`.

### Environment setup

A helper loads every external tool as a module and activates the Python
environment. **Source** it once per shell (don't execute it):

```bash
source scripts/setup_env.sh
```

> [!IMPORTANT]
> **Use two terminals.** `setup_env.sh` is only needed for the pipeline
> (Step 1+). Run the autoPROC step (Step 0) in its **own, plain** terminal, and
> the pipeline in a **separate** terminal where you have sourced `setup_env.sh`
> — the pipeline environment (PHENIX/CCP4/PyMOL modules + virtualenv) can
> otherwise interfere with autoPROC's own module setup.

Make the scripts executable once after cloning:

```bash
chmod +x scripts/*.py scripts/*.sh
```

---

## Full walkthrough

Throughout, replace `CaMDH_001` / paths / templates with your own dataset.

### Step 0 — Submit chunked autoPROC jobs

> Run **from the directory containing your images** (a CBF template or an Eiger
> HDF5 master file — the mode is detected automatically).

Edit the variables at the top of `scripts/TR-FRX_autoPROC.sh` (the exact names
are documented in the script header):

| Variable | Description | Example |
|---|---|---|
| `RUN_ID` | Dataset run number | `001` |
| `IMAGE_TEMPLATE` | Image template (`####` = frame no.) or HDF5 master | `CaMDH_001_1_####.cbf` |
| `FIRST_IMG` / `LAST_IMG` | Total frame range | `1` / `3000` |
| `REF_FIRST_IMG` / `REF_LAST_IMG` | Reference chunk range | `1` / `300` |
| `CHUNK_SIZE` | Frames per subsequent chunk | `300` |
| `SYMM` | Space group | `P3121` |
| `CELL` | Unit-cell parameters | `107 107 103 90 90 120` |
| `SLURM_PARTITION` | SLURM partition | `nice` |
| `SLURM_CPUS` / `SLURM_MEM` | Resources per job | `24` / `24000` |

Then:

```bash
cd /data/images/CaMDH_001
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
```

> [!CAUTION]
> **Wait for all jobs to finish before continuing.**

### Step 1+ — The full pipeline (one command)

**USE IN A TERMINAL WITH THE ENV ACTIVATED**

Everything downstream is a single script. Point it at your **model**, the
**autoproc_chunks** folder from Step 0, and an **output** folder:

```bash
python path/to/scripts/trfrx_full_pipeline.py  model.pdb  /data/images/PfuGRHPR_006/autoproc_chunks  ./
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
/data/processed/PfuGRHPR_006/
    autoproc_copy/                  # copied chunks (created once)
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

---

## Quick reference

```bash
# ── Terminal 1 · plain shell ─────────────────────────────
# 0. edit variables, then submit autoPROC from the image folder
cd /data/images/PfuGRHPR_006
bash /path/to/scripts/TR-FRX_autoPROC.sh          # → wait for SLURM

# ── Terminal 2 · pipeline environment ────────────────────
source scripts/setup_env.sh
python scripts/trfrx_full_pipeline.py  model.pdb  /data/images/PfuGRHPR_006/autoproc_chunks  /data/processed
#     → /data/processed/PfuGRHPR_006/run_01/report.html
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

## Contact

Julien Martel · [@jmartelX](https://github.com/jmartelX)
