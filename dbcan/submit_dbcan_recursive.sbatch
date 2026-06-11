#!/bin/bash -l
#SBATCH -J dbcan_submitter
#SBATCH -p commons
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#SBATCH -o /projects/alvarez/scripts/logs/dbcan_submitter.%j.out
#SBATCH -e /projects/alvarez/scripts/logs/dbcan_submitter.%j.err

set -Eeuo pipefail

# Required variables must be provided via --export
LIST=${LIST:?LIST missing}
WORKER=${WORKER:?WORKER missing}
START=${START:?START missing}
N=${N:?N missing}
CHUNK=${CHUNK:-200}      # Number of tasks per batch (adjust based on cluster/QOS; 200 is generally safe)
MAXPAR=${MAXPAR:-30}     # Maximum concurrent tasks within each batch (%30)
PARTITION=${PARTITION:-commons}

END=$(( START + CHUNK - 1 ))
if [[ $END -ge $N ]]; then
  END=$(( N - 1 ))
fi

RANGE="${START}-${END}%${MAXPAR}"

echo "[INFO] Submit array: ${RANGE}"

jid=$(sbatch -p "$PARTITION" --array="$RANGE" "$WORKER" | awk '{print $4}')

echo "[INFO] Array job ID=${jid}"

NEXT=$(( END + 1 ))

if [[ $NEXT -lt $N ]]; then
  echo "[INFO] Schedule next submitter with START=${NEXT}"

  sbatch -p "$PARTITION" \
    --dependency=afterok:${jid} \
    --export=ALL,LIST="$LIST",WORKER="$WORKER",START="$NEXT",N="$N",CHUNK="$CHUNK",MAXPAR="$MAXPAR",PARTITION="$PARTITION" \
    $0
else
  echo "[OK] All batches have been scheduled."
fi
