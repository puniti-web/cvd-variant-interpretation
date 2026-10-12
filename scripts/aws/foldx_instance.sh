#!/bin/bash
# EC2 user-data for the Phase 5 FoldX run (Ubuntu 24.04, x86_64).
#
# Expects in S3 (private bucket, same account as the instance role):
#   s3://$BUCKET/foldx5_Linux_0.zip           FoldX binary (licensed, never public)
#   s3://$BUCKET/$PREFIX/work.tar.gz          output of `foldx_ddg.py prepare` + foldx_ddg.py
# Writes back:
#   s3://$BUCKET/$PREFIX/results/             one CSV per batch, synced every 5 min
#   s3://$BUCKET/$PREFIX/repaired/            repaired chunks, so a restart skips RepairPDB
#   s3://$BUCKET/$PREFIX/logs/shard-$SHARD.log
#   s3://$BUCKET/$PREFIX/done/shard-$SHARD    marker written when the shard finishes
#
# A restarted instance resumes from what is already in S3. The instance
# shuts itself down at the end; launch it with shutdown behaviour "terminate".
set -euo pipefail

BUCKET="${BUCKET:-hermes-biodb-artifacts-448016109107}"
PREFIX="${PREFIX:-cvd-foldx}"
SHARD="${SHARD:-0}"
SHARD_WEIGHTS="${SHARD_WEIGHTS:-1}"
AUTO_SHUTDOWN="${AUTO_SHUTDOWN:-1}"

ROOT=/opt/foldx
LOG=$ROOT/shard-$SHARD.log
mkdir -p $ROOT && cd $ROOT
exec > >(tee -a $LOG) 2>&1
echo "start $(date -u) shard=$SHARD weights=$SHARD_WEIGHTS"

if ! command -v aws >/dev/null; then
  apt-get update -q && apt-get install -yq unzip curl
  curl -sSL https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip -o awscliv2.zip
  unzip -q awscliv2.zip && ./aws/install && rm -rf aws awscliv2.zip
fi
command -v unzip >/dev/null || { apt-get update -q && apt-get install -yq unzip; }

S3="s3://$BUCKET/$PREFIX"
aws s3 cp --only-show-errors "s3://$BUCKET/foldx5_Linux_0.zip" .
unzip -oq foldx5_Linux_0.zip -d bin && chmod +x bin/foldx_*
FOLDX=$(ls $ROOT/bin/foldx_* | head -1)
aws s3 cp --only-show-errors "$S3/work.tar.gz" . && tar xzf work.tar.gz
aws s3 sync --only-show-errors "$S3/results/" work/results/
aws s3 sync --only-show-errors "$S3/repaired/" work/repaired/

sync_up() {
  aws s3 sync --only-show-errors work/results/ "$S3/results/" --exclude "*.part"
  aws s3 sync --only-show-errors work/repaired/ "$S3/repaired/"
  aws s3 cp --only-show-errors $LOG "$S3/logs/shard-$SHARD.log"
}
( while true; do sleep 300; sync_up || true; done ) &
SYNC_PID=$!

python3 foldx_ddg.py run --foldx "$FOLDX" --work-dir work --workers "$(nproc)" \
  --shard-weights "$SHARD_WEIGHTS" --shard "$SHARD" || echo "run exited with $?"

kill $SYNC_PID || true
sync_up
echo "finished $(date -u)" | aws s3 cp - "$S3/done/shard-$SHARD"
echo "finished $(date -u)"
sync_up
[ "$AUTO_SHUTDOWN" = "1" ] && shutdown -h now
