#!/usr/bin/env bash
#SBATCH --job-name=gtdbtk-r232
#SBATCH --output=/group/sbms004/yxia/GUT/logs/gtdbtk_r232_%j.out
#SBATCH --error=/group/sbms004/yxia/GUT/logs/gtdbtk_r232_%j.err
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

set -eo pipefail

DB_DIR="/group/sbms004/yxia/GUT/db/release232"
LOG_DIR="/group/sbms004/yxia/GUT/logs"
TARBALL="${DB_DIR}/gtdbtk_data.tar.gz"

URL="https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz"

mkdir -p "${DB_DIR}"
mkdir -p "${LOG_DIR}"

source ~/miniconda3/etc/profile.d/conda.sh

set +u
conda activate SqueezeMeta18
set -u

cd "${DB_DIR}"

echo "Downloading GTDB-Tk database..."
wget -c -O "${TARBALL}" "${URL}"

echo "Extracting..."
tar -xvzf "${TARBALL}" -C "${DB_DIR}"

export GTDBTK_DATA_PATH="${DB_DIR}"

echo "Checking database files..."
find "${GTDBTK_DATA_PATH}" -name "markers.bin" | head
find "${GTDBTK_DATA_PATH}" -name "*.sketch" | head

echo "Running GTDB-Tk check_install..."
gtdbtk check_install

echo "Done."
echo "Use this path in SqueezeMeta:"
echo "GTDBTK_DB=\"${GTDBTK_DATA_PATH}\""