#!/usr/bin/env python3

import os
import sys
from math import ceil

# -------- CONFIGURATION --------
output_root = "/group/sbms004/yxia/GUT/sample_list"
input_dir   = "/group/sbms004/yxia/GUT/cleaned_fastq"
out_prefix  = "sample_list"
samples_per_group = 1
max_groups = 300
# --------------------------------

os.makedirs(output_root, exist_ok=True)

if not os.path.isdir(input_dir):
    print(f"ERROR: Directory does not exist: {input_dir}")
    sys.exit(1)

sample_map = {}

for fname in sorted(os.listdir(input_dir)):
    if not fname.endswith(".fastq.gz"):
        continue

    if fname.endswith("_R1_clean.fastq.gz"):
        sample_id = fname.replace("_R1_clean.fastq.gz", "")
        tag = "R1"
    elif fname.endswith("_R2_clean.fastq.gz"):
        sample_id = fname.replace("_R2_clean.fastq.gz", "")
        tag = "R2"
    else:
        continue

    sample_map.setdefault(sample_id, {})[tag] = fname

valid_samples = sorted(
    sample_id for sample_id, reads in sample_map.items()
    if "R1" in reads and "R2" in reads
)

if not valid_samples:
    print(f"ERROR: No valid paired clean FASTQ files found in {input_dir}")
    print("Expected file names like:")
    print("  SAMPLE_R1_clean.fastq.gz")
    print("  SAMPLE_R2_clean.fastq.gz")
    sys.exit(1)

print(f"Found {len(valid_samples)} valid paired samples.")

# 1) Write full sample list first
all_samples_file = os.path.join(output_root, "all_samples.tsv")

with open(all_samples_file, "w") as out:
    for sample_id in valid_samples:
        fnames = sample_map[sample_id]
        out.write(f"{sample_id}\t{fnames['R1']}\tpair1\n")
        out.write(f"{sample_id}\t{fnames['R2']}\tpair2\n")

print(f"Written full sample list: {all_samples_file}")
print(f"Total lines written: {len(valid_samples) * 2}")

# 2) Split into grouped sample lists
n_groups = ceil(len(valid_samples) / samples_per_group)

if n_groups > max_groups:
    print(f"WARNING: {len(valid_samples)} samples require {n_groups} groups.")
    print(f"Only writing first {max_groups} groups because max_groups={max_groups}.")
    n_groups = max_groups

for grp in range(n_groups):
    chunk = valid_samples[
        grp * samples_per_group : (grp + 1) * samples_per_group
    ]

    out_file = os.path.join(output_root, f"{out_prefix}_{grp + 1}.tsv")

    with open(out_file, "w") as out:
        for sample_id in chunk:
            fnames = sample_map[sample_id]
            out.write(f"{sample_id}\t{fnames['R1']}\tpair1\n")
            out.write(f"{sample_id}\t{fnames['R2']}\tpair2\n")

    print(
        f"Written group {grp + 1}: "
        f"{len(chunk)} samples, {len(chunk) * 2} lines -> {out_file}"
    )

print("Done.")