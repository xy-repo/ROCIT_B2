#!/usr/bin/env bash
#SBATCH --job-name=squeezemeta_db
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=3-00:00:00
#SBATCH --output=/group/sbms004/yxia/GUT/logs/squeezemeta_db_%j.out
#SBATCH --error=/group/sbms004/yxia/GUT/logs/squeezemeta_db_%j.err

set -eo pipefail

mkdir -p /group/sbms004/yxia/GUT/logs

source ~/miniconda3/etc/profile.d/conda.sh
conda activate SqueezeMeta

DB_DIR="/group/sbms004/yxia/GUT/db/SqueezeMeta_DB"

mkdir -p "${DB_DIR}"

echo "=========================================="
echo "[+] Starting SqueezeMeta database download"
echo "DB_DIR: ${DB_DIR}"
echo "Date: $(date)"
echo "Disk space:"
df -h /group/sbms004/yxia/GUT/db
echo "=========================================="

download_databases.pl "${DB_DIR}"

echo "=========================================="
echo "[+] Database download finished"
echo "Date: $(date)"
echo "Database directory:"
ls -lh "${DB_DIR}"
echo "Disk space:"
df -h /group/sbms004/yxia/GUT/db
echo "=========================================="

test_install.pl