#!/bin/bash
#SBATCH --job-name=agrf_download
#SBATCH --output=agrf_download.out
#SBATCH --error=agrf_download.err
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

# Load conda environment if needed
# source ~/miniconda3/etc/profile.d/conda.sh
# conda activate your_env

# Credentials
USER="MarkNicol"
PASS="hoThip-5mokxu-diqmap"
HOST="agrf-data.agrf.org.au"
REMOTE_DIR="files/AGRF_NXGSQCAGRF25110222-1_23K7VWLT3"
LOCAL_DIR="$PWD/AGRF_download"

mkdir -p "$LOCAL_DIR"

# Run lftp
lftp -u "$USER","$PASS" sftp://$HOST <<EOF
set sftp:auto-confirm yes
set net:max-retries 2
set net:timeout 20
set xfer:clobber yes

mirror --verbose \
       --continue \
       --parallel=4 \
       $REMOTE_DIR $LOCAL_DIR

bye
EOF
