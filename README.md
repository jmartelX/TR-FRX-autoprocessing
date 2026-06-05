# TR-FRX autoprocessing

A reproducible command-line pipeline for processing **Time-Resolved Functional
Rotation Crystallography (TR-FRX)** data — from raw CBF diffraction images all
the way to **Fo–Fo difference maps** and their **singular-value decomposition
(SVD)**.

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
difference maps** — Fo(t) – Fo(reference) — and the **SVD components** that
summarise the structural changes across the time series. It wraps and
orchestrates well-established crystallography software (autoPROC, CCP4, PHENIX)
with a small set of Python/shell tools so an entire time-resolved dataset can be
processed with a handful of commands.

> Throughout this documentation, **"chunk" = "subdataset" = one time point** —
> a complete sub-dataset extracted from the continuous rotation collection.

See the associated paper for the method itself ([Citing this
work](#citing-this-work)).

---

## Pipeline at a glance

| Step | Script | What it does |
|:---:|---|---|
| 0 | `scripts/TR-FRX_autoPROC_cbf.sh` | Submits **chunked autoPROC jobs** to SLURM → produces `autoPROC_*` folders |
| 0b | `scripts/fix_cbf_headers.py` | *(optional)* Repairs broken/missing CBF headers before processing |
| 1 | `scripts/copy_files.py` | Collects the useful files from each `autoPROC_*` folder into a clean tree |
| 2 | `scripts/fetch_clean_mtz.py` | Finds and **cleans** the MTZs (keeps only `H K L F SIGF FreeR_flag` via CCP4 `cad`) |
| 3 | `scripts/diffmaps.py` | Computes **Fo–Fo difference maps** for every dataset against a reference (PHENIX) |
| 4 | `scripts/SVD_all_in_one.py` | Runs **SVD** on the difference-map series to extract time-resolved components |

Each script is self-documenting — run it with `-h`/`--help`, or read the
docstring at the top of the file.

---

## Requirements

**External software** (must be available on your `PATH` / loaded as modules):

| Tool | Used by | Provides |
|---|---|---|
| [autoPROC](https://www.globalphasing.com/autoproc/) | step 0 | `process` |
| SLURM | step 0 | `sbatch` (cluster job submission) |
| [CCP4](https://www.ccp4.ac.uk/) | step 2 | `cad`, `mtzdmp` |
| [PHENIX](https://phenix-online.org/) | steps 3–4 | `phenix.fobs_minus_fobs_map`, `phenix.python` |

On an HPC system these are typically loaded with environment modules, e.g.:

```bash
module load autoPROC ccp4 phenix
```

**Python ≥ 3.10**, with: `numpy`, `pandas`, `scipy`, `gemmi`, `matplotlib`,
`seaborn`, `wxmplot`, `dask`.

### Quick environment setup

A helper is provided to build a virtual environment with all Python
dependencies (it reuses `phenix.python` if available, so `gemmi`/`numpy` come
for free):

```bash
bash setup_svd_env.sh
# then, in any new session:
source ~/.venv/svd_pipeline/bin/activate
```

Make the scripts executable once after cloning:

```bash
chmod +x scripts/*.py scripts/*.sh
```

---

## Full walkthrough

Throughout, replace `CaMDH_073` / paths / `PfuGRHPR_006` with your own dataset.

### Step 0 — Submit chunked autoPROC jobs

> Run **from the directory containing your CBF images.**

Edit the variables at the top of `scripts/TR-FRX_autoPROC_cbf.sh`:

| Variable | Description | Example |
|---|---|---|
| `RUN_ID` | Dataset run number | `006` |
| `CBF_TEMPLATE` | CBF filename template (`####` = frame number) | `PfuGRHPR_006_1_####.cbf` |
| `FIRST_IMG` / `LAST_IMG` | Total frame range | `1` / `3000` |
| `REF_FIRST_IMG` / `REF_LAST_IMG` | Reference chunk range | `1` / `300` |
| `CHUNK_SIZE` | Frames per subsequent chunk | `300` |
| `SYMM` | Space group | `I41` |
| `CELL` | Unit-cell parameters | `114 114 118 90 90 90` |
| `SLURM_PARTITION` | SLURM partition | `nice` |
| `SLURM_CPUS` / `SLURM_MEM` | Resources per job | `24` / `24000` |

Then:

```bash
cd /data/images/PfuGRHPR_006
bash /path/to/TR-FRX-autoprocessing/scripts/TR-FRX_autoPROC_cbf.sh
```

It submits one SLURM job for the **reference** chunk (e.g. 1–300), then one
job per subsequent chunk (301–600, 601–900, …), each depending on the
reference job. Output lands in `autoproc_chunks/`:

```
autoproc_chunks/
    autoPROC_1_300/       ← reference
    autoPROC_301_600/
    autoPROC_601_900/
    ...
```

Wait for all jobs to finish before continuing.

> **CBF header problems?** If autoPROC complains about missing detector
> metadata, repair the headers first with
> `scripts/fix_cbf_headers.py --help` (preview with `--dry-run`).

### Step 1 — Collect the useful files

```bash
./scripts/copy_files.py <source_dir> <destination_dir>
# e.g.
./scripts/copy_files.py /data/images/PfuGRHPR_006/autoproc_chunks /data/processed
```

Scans for `autoPROC_*` folders and copies the key outputs
(`staraniso_alldata-unique.mtz`, `truncate-unique.mtz`, `summary.html`,
`report.pdf`, `XDS_ASCII.HKL`, …) into a structured destination, one folder per
chunk.

### Step 2 — Clean & collect the MTZs

Create a working folder and run the cleaner from inside it:

```bash
mkdir /data/processed/PfuGRHPR_006/dFo
cd    /data/processed/PfuGRHPR_006/dFo
../../../TR-FRX-autoprocessing/scripts/fetch_clean_mtz.py /data/processed/PfuGRHPR_006 .
```

For each chunk it finds `staraniso_alldata-unique.mtz`, strips it down to
`H K L F SIGF FreeR_flag` with CCP4 `cad`, and writes a named MTZ:

```
dFo/
    PfuGRHPR_006_1_300.mtz       ← will serve as reference (lowest index)
    PfuGRHPR_006_301_600.mtz
    PfuGRHPR_006_601_900.mtz
```

### Step 3 — Add a model, compute Fo–Fo difference maps

Drop your refined model into the `dFo` folder (named `model.pdb` if there are
several PDBs), then:

```bash
cd /data/processed/PfuGRHPR_006/dFo
/path/to/TR-FRX-autoprocessing/scripts/diffmaps.py
```

Runs `phenix.fobs_minus_fobs_map` for every MTZ against the
lowest-index reference. Results go to `dFo/output/`:

```
dFo/output/
    dFo_PfuGRHPR_006_301_600-PfuGRHPR_006_1_300_1.mtz
    dFo_PfuGRHPR_006_301_600-PfuGRHPR_006_1_300_1.map
    ...
```

Useful options:

```bash
--dir /path/to/dFo   # run from elsewhere
--high-res 1.8       # high-resolution cutoff (Å)
--low-res  8.0       # low-resolution  cutoff (Å)
--dry-run            # print the plan without running PHENIX
```

### Step 4 — SVD of the difference-map series

```bash
/path/to/TR-FRX-autoprocessing/scripts/SVD_all_in_one.py --help
```

Performs a singular-value decomposition across the Fo–Fo difference-map series
to separate time-resolved structural signal from noise. See the script's
`--help` for inputs and options.

---

## Quick reference

```bash
# 0. edit variables, submit from the CBF folder
cd /data/images/PfuGRHPR_006
bash TR-FRX_autoPROC_cbf.sh                 # → wait for SLURM

# 1. collect
./copy_files.py /data/images/PfuGRHPR_006/autoproc_chunks /data/processed

# 2. clean MTZs
mkdir /data/processed/PfuGRHPR_006/dFo && cd $_
./fetch_clean_mtz.py /data/processed/PfuGRHPR_006 .

# 3. add model.pdb, then difference maps
./diffmaps.py

# 4. SVD
./SVD_all_in_one.py
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
