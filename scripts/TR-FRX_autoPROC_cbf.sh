#!/bin/bash

set -euo pipefail

RUN_ID=006

# --- À lancer depuis le dossier contenant les images CBF ---
IMAGES_DIR="$(pwd)"
CBF_TEMPLATE="PfuGRHPR_006_1_####.cbf"

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

module load autoPROC

# --- Tous les traitements seront créés dans un sous-dossier du dossier courant ---
PROCESS_DIR="${IMAGES_DIR}/autoproc_chunks"
mkdir -p "$PROCESS_DIR"

FIRST_OUTDIR=""
FIRST_JOBID=""

# --- Traitement du jeu de référence : images 1 à 110 ---
REF_JOB_BASENAME="${REF_FIRST_IMG}-${REF_LAST_IMG}_autoproc"
REF_JOB_SCRIPT="${PROCESS_DIR}/${REF_JOB_BASENAME}.sh"
REF_OUTDIR="${PROCESS_DIR}/autoPROC_${REF_FIRST_IMG}_${REF_LAST_IMG}"

cat > "$REF_JOB_SCRIPT" <<EOF
#!/bin/bash
set -euo pipefail
module load autoPROC
cd "$PROCESS_DIR"
process -Id "$RUN_ID,$IMAGES_DIR,$CBF_TEMPLATE,$REF_FIRST_IMG,$REF_LAST_IMG" symm="$SYMM" cell="$CELL" AutoProcScale_RunStaraniso=yes -d "$REF_OUTDIR"
EOF

sbatch_output=$(sbatch -p "$SLURM_PARTITION" -n 1 -c "$SLURM_CPUS" --mem="$SLURM_MEM" "$REF_JOB_SCRIPT")
FIRST_JOBID=$(echo "$sbatch_output" | awk '{print $4}')
FIRST_OUTDIR="$REF_OUTDIR"

# --- Traitement des jeux suivants par blocs de 220 images ---
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
process -Id "$RUN_ID,$IMAGES_DIR,$CBF_TEMPLATE,$i,$j" -ref "$FIRST_OUTDIR/truncate-unique.mtz" AutoProcScale_RunStaraniso=yes -d "$OUTDIR"
EOF

    sbatch -p "$SLURM_PARTITION" -n 1 -c "$SLURM_CPUS" --mem="$SLURM_MEM" --dependency=afterok:"$FIRST_JOBID" "$JOB_SCRIPT"

    i=$((j + 1))
done

echo "Traitements soumis. Résultats dans : $PROCESS_DIR"
