#!/bin/bash -l
#SBATCH -J dbcan_submitter
#SBATCH -p commons
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#SBATCH -o /projects/alvarez/scripts/logs/dbcan_submitter.%j.out
#SBATCH -e /projects/alvarez/scripts/logs/dbcan_submitter.%j.err

set -Eeuo pipefail

# 必要变量通过 --export 传入
LIST=${LIST:?LIST missing}
WORKER=${WORKER:?WORKER missing}
START=${START:?START missing}
N=${N:?N missing}
CHUNK=${CHUNK:-200}      # 每批任务数（按你们QOS调；200较稳）
MAXPAR=${MAXPAR:-30}     # 每批内部并发（%30）
PARTITION=${PARTITION:-commons}

END=$(( START + CHUNK - 1 ))
if [[ $END -ge $N ]]; then END=$(( N - 1 )); fi
RANGE="${START}-${END}%${MAXPAR}"

echo "[INFO] submit array: ${RANGE}"
jid=$(sbatch -p "$PARTITION" --array="$RANGE" "$WORKER" | awk '{print $4}')
echo "[INFO] array jobid=${jid}"

NEXT=$(( END + 1 ))
if [[ $NEXT -lt $N ]]; then
  echo "[INFO] schedule next submitter for NEXT=${NEXT}"
  sbatch -p "$PARTITION" \
    --dependency=afterok:${jid} \
    --export=ALL,LIST="$LIST",WORKER="$WORKER",START="$NEXT",N="$N",CHUNK="$CHUNK",MAXPAR="$MAXPAR",PARTITION="$PARTITION" \
    $0
else
  echo "[OK] all batches scheduled."
fi
