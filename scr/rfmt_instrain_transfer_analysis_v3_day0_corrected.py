#!/usr/bin/env python3
"""
RFMT inStrain transfer analysis

Purpose
-------
Use per-RFMT inStrain genomeWide_compare.tsv tables plus per-sample bintables
and a cleaned RFMT metadata sample list to classify donor-to-recipient strain
transfer at MAG level and summarize successfully transmitted species.

Expected inputs
---------------
1) Clean sample metadata table, e.g.
   RFMT_analysis_sample_long.tsv
   Required columns: RFMT_ID, sample_name_expected, role

2) inStrain final root, e.g.
   /group/sbms004/yxia/GUT/inStrain_final
   Expected compare file example:
   /group/sbms004/yxia/GUT/inStrain_final/RFMT015/compare_out_with_stb/output/compare_out_with_stb_genomeWide_compare.tsv

3) MAG root, e.g.
   /group/sbms004/yxia/GUT/MAG_per_sample
   Expected bintable example:
   /group/sbms004/yxia/GUT/MAG_per_sample/sample_list_13/RFMT013_5430/RFMT013_5430/results/18.RFMT013_5430.bintable

Main outputs
------------
- mag_transfer_calls.tsv
- species_transmission_summary.tsv
- rfmt_transfer_summary.tsv
- bintable_mag_annotation.tsv
- missing_inputs.tsv

Notes
-----
The genome name in genomeWide_compare.tsv is expected to look like:
  RFMT015_5423_RFMT015_5423.metabat2.10
The bintable Bin ID is expected to look like:
  RFMT015_5423.metabat2.10
This script creates the join key:
  genome = sample_name + "_" + Bin ID
"""

from __future__ import annotations

import argparse
import glob
import os
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd


# Timeline interpretation used by this script
# recipient_d0 is PRE-Donor1 baseline, not a post-transfer sample.
# recipient_pre_donor2 is after Donor1 and immediately before Donor2 when available.
DONOR1_DIRECT_POST_ROLES = ["recipient_pre_donor2"]
DONOR1_PERSISTENCE_POST_DONOR2_ROLES = ["recipient_d2", "recipient_d5", "recipient_d7"]
DONOR2_POST_ROLES = ["recipient_d2", "recipient_d5", "recipient_d7"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Classify Donor1/Donor2 RFMT strain transfer with recipient_d0 as pre-Donor1 baseline using inStrain genomeWide_compare.tsv and bintables."
    )
    p.add_argument(
        "--metadata",
        required=True,
        help="Clean RFMT sample_long TSV, e.g. RFMT_analysis_sample_long.tsv",
    )
    p.add_argument(
        "--instrain-root",
        default="/group/sbms004/yxia/GUT/inStrain_final",
        help="Root containing RFMTxxx/compare_out_with_stb/output/*genomeWide_compare.tsv",
    )
    p.add_argument(
        "--mag-root",
        default="/group/sbms004/yxia/GUT/MAG_per_sample",
        help="Root containing sample_list_*/SAMPLE/SAMPLE/results/*.bintable",
    )
    p.add_argument(
        "--outdir",
        required=True,
        help="Output directory",
    )
    p.add_argument(
        "--min-compared-bases",
        type=int,
        default=50000,
        help="High-confidence minimum compared_bases_count",
    )
    p.add_argument(
        "--min-percent-compared",
        type=float,
        default=5.0,
        help="High-confidence minimum percent_compared. If values are 0-1 in your table, the script auto-converts to percent.",
    )
    p.add_argument(
        "--min-conani",
        type=float,
        default=99.99,
        help="Minimum conANI for donor-like call",
    )
    p.add_argument(
        "--min-popani",
        type=float,
        default=99.90,
        help="Minimum popANI for donor-like call when popANI is present",
    )
    p.add_argument(
        "--snp-rate-margin",
        type=float,
        default=1e-5,
        help="Donor SNP rate must be this much lower than baseline SNP rate to call donor closer when both exist",
    )
    p.add_argument(
        "--include-limited-groups",
        action="store_true",
        help="Also analyze RFMT groups without donor samples where possible. Default keeps only groups with >=2 samples and at least one donor.",
    )
    return p.parse_args()


def read_tsv_auto(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, low_memory=False)


def to_num(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def normalize_sample_name(x: object) -> str:
    """Convert paths or profile folder names to RFMTxxx_#### sample IDs when possible."""
    if pd.isna(x):
        return ""
    x = str(x).strip()
    x = os.path.basename(x.rstrip("/"))
    x = re.sub(r"_IS$", "", x)
    m = re.search(r"(RFMT\d{3}_\d+)", x)
    return m.group(1) if m else x


def make_full_sample_id(rfmt: str, sample: object) -> str:
    """
    Convert sample IDs to the folder naming style used by your MAG output.

    Examples
    --------
    RFMT015 + 5423          -> RFMT015_5423
    RFMT015 + RFMT015_5423  -> RFMT015_5423
    """
    if pd.isna(sample):
        return ""
    sample = str(sample).strip()
    if sample.endswith(".0"):
        sample = sample[:-2]
    if sample.startswith(rfmt + "_"):
        return sample
    if sample.startswith("RFMT"):
        return sample
    return f"{rfmt}_{sample}"


def find_compare_file(instrain_root: str, rfmt: str) -> Optional[str]:
    patterns = [
        f"{instrain_root}/{rfmt}/compare_out_with_stb/output/*genomeWide_compare.tsv",
        f"{instrain_root}/{rfmt}/compare_out_with_stb/output/*genome*compare*.tsv",
        f"{instrain_root}/{rfmt}/compare_out_with_stb/**/*genomeWide_compare.tsv",
        f"{instrain_root}/{rfmt}/compare_out_with_stb/**/*genome*compare*.tsv",
    ]
    hits: List[str] = []
    for pat in patterns:
        hits.extend(glob.glob(pat, recursive=True))
    hits = sorted(set(hits))
    if not hits:
        return None
    # Prefer explicit genomeWide file.
    for h in hits:
        if "genomewide" in os.path.basename(h).lower():
            return h
    return hits[0]


def find_bintable(mag_root: str, rfmt: str, sample_name: object) -> Optional[str]:
    """
    Find a sample bintable while ignoring the sample_list_* number.

    Example target
    --------------
    /group/sbms004/yxia/GUT/MAG_per_sample/sample_list_13/RFMT013_5430/RFMT013_5430/results/18.RFMT013_5430.bintable

    This function searches sample_list_* and accepts either full sample names
    such as RFMT013_5430 or sample numbers such as 5430.
    """
    full_sample = make_full_sample_id(rfmt, sample_name)
    patterns = [
        f"{mag_root}/sample_list_*/{full_sample}/{full_sample}/results/*.bintable",
        f"{mag_root}/sample_list_*/{full_sample}/{full_sample}/results/*{full_sample}*.bintable",
        f"{mag_root}/**/{full_sample}/{full_sample}/results/*.bintable",
    ]
    hits: List[str] = []
    for pat in patterns:
        hits.extend(glob.glob(pat, recursive=True))
    hits = sorted(set(hits))
    if not hits:
        return None
    # Prefer non-empty file.
    for h in hits:
        try:
            if os.path.getsize(h) > 0:
                return h
        except OSError:
            continue
    return hits[0]


def gtdb_species(tax: object) -> str:
    """Extract species from GTDB-Tk taxonomy string if possible."""
    if pd.isna(tax):
        return "Unclassified"
    text = str(tax)
    if text.strip() == "":
        return "Unclassified"
    parts = [p.strip() for p in text.split(";")]
    sp = [p for p in parts if p.startswith("s__")]
    if sp:
        val = sp[-1].replace("s__", "").strip()
        return val if val else "Unclassified_species"
    # fallback: last non-empty token
    nonempty = [p for p in parts if p]
    return nonempty[-1] if nonempty else "Unclassified"


def gtdb_genus(tax: object) -> str:
    if pd.isna(tax):
        return "Unclassified"
    parts = [p.strip() for p in str(tax).split(";")]
    gen = [p for p in parts if p.startswith("g__")]
    if gen:
        val = gen[-1].replace("g__", "").strip()
        return val if val else "Unclassified_genus"
    return "Unclassified"


def read_bintable_auto(path: str | Path) -> pd.DataFrame:
    """
    Read a bintable and auto-detect whether the real header is on line 1 or line 2.

    Your bintables may need:
        pd.read_csv(path, sep="\t", header=1, dtype=str)
    because the first line can be metadata and the second line contains:
        Bin ID, Method, Tax, Tax GTDB-Tk, Length, ...
    """
    path = str(path)
    attempts = []
    for header in [0, 1]:
        try:
            df = pd.read_csv(path, sep="\t", header=header, dtype=str, low_memory=False)
            attempts.append((header, list(df.columns)))
            if "Bin ID" in df.columns:
                return df
        except Exception as e:
            attempts.append((header, f"READ_ERROR: {e}"))
    raise ValueError(f"No 'Bin ID' column in bintable: {path}; attempts={attempts}")


def read_one_bintable(path: str, sample_name: str, rfmt: str) -> pd.DataFrame:
    sample_name = make_full_sample_id(rfmt, sample_name)
    df = read_bintable_auto(path)
    if "Bin ID" not in df.columns:
        raise ValueError(f"No 'Bin ID' column in bintable: {path}")

    df["RFMT_ID"] = rfmt
    df["sample_name"] = sample_name
    df["bintable_path"] = path
    df["Bin ID"] = df["Bin ID"].astype(str)
    df["genome"] = sample_name + "_" + df["Bin ID"]

    # Add stable taxonomy fields.
    tax_col = "Tax GTDB-Tk" if "Tax GTDB-Tk" in df.columns else ("Tax" if "Tax" in df.columns else None)
    if tax_col:
        df["taxonomy_source"] = tax_col
        df["taxonomy"] = df[tax_col]
    else:
        df["taxonomy_source"] = "missing"
        df["taxonomy"] = np.nan
    df["genus"] = df["taxonomy"].map(gtdb_genus)
    df["species"] = df["taxonomy"].map(gtdb_species)

    # Normalize useful numeric columns if present.
    for col in ["Length", "GC perc", "Num contigs", "Disparity", "Completeness", "Contamination", "Strain het"]:
        if col in df.columns:
            df[col] = to_num(df[col])

    cov_col = f"Coverage {sample_name}"
    tpm_col = f"TPM {sample_name}"
    df["donor_sample_coverage_from_bintable"] = to_num(df[cov_col]) if cov_col in df.columns else np.nan
    df["donor_sample_TPM_from_bintable"] = to_num(df[tpm_col]) if tpm_col in df.columns else np.nan
    return df


def load_bintables(metadata: pd.DataFrame, mag_root: str) -> Tuple[pd.DataFrame, List[dict]]:
    rows = []
    missing = []
    for _, r in metadata[["RFMT_ID", "sample_name_expected"]].drop_duplicates().iterrows():
        rfmt = r["RFMT_ID"]
        sample = make_full_sample_id(rfmt, r["sample_name_expected"])
        bt = find_bintable(mag_root, rfmt, sample)
        if bt is None:
            missing.append({"RFMT_ID": rfmt, "sample_name": sample, "missing_type": "bintable", "path_checked": mag_root})
            continue
        try:
            rows.append(read_one_bintable(bt, sample, rfmt))
        except Exception as e:
            missing.append({"RFMT_ID": rfmt, "sample_name": sample, "missing_type": "bintable_read_error", "path_checked": bt, "error": str(e)})
    if not rows:
        return pd.DataFrame(), missing
    return pd.concat(rows, ignore_index=True, sort=False), missing


def standardize_compare(df: pd.DataFrame, rfmt: str, path: str) -> pd.DataFrame:
    df = df.copy()
    df["RFMT_ID"] = rfmt
    df["compare_path"] = path

    # Column aliases across inStrain versions.
    colmap = {}
    lower = {c.lower(): c for c in df.columns}

    def pick(names: Iterable[str]) -> Optional[str]:
        for n in names:
            if n in df.columns:
                return n
            if n.lower() in lower:
                return lower[n.lower()]
        return None

    genome_col = pick(["genome", "Genome", "scaffold", "scaffold_name"])
    name1_col = pick(["name1", "sample1", "Sample1", "sample_1", "profile1"])
    name2_col = pick(["name2", "sample2", "Sample2", "sample_2", "profile2"])

    if genome_col is None or name1_col is None or name2_col is None:
        raise ValueError(
            f"Could not identify genome/name1/name2 columns in {path}. Columns: {list(df.columns)}"
        )

    colmap[genome_col] = "genome"
    colmap[name1_col] = "name1"
    colmap[name2_col] = "name2"
    df = df.rename(columns=colmap)

    for c in ["name1", "name2"]:
        df[c] = df[c].map(normalize_sample_name)

    for c in [
        "conANI",
        "popANI",
        "consensus_SNPs",
        "population_SNPs",
        "compared_bases_count",
        "percent_compared",
        "coverage_overlap",
    ]:
        if c in df.columns:
            df[c] = to_num(df[c])
        else:
            df[c] = np.nan

    # If ANI is 0-1, convert to percent.
    for c in ["conANI", "popANI"]:
        if c in df.columns and df[c].dropna().between(0, 1.5).all() and len(df[c].dropna()) > 0:
            df[c] = df[c] * 100.0

    # If percent_compared is 0-1, convert to percent.
    if df["percent_compared"].dropna().between(0, 1.5).all() and len(df["percent_compared"].dropna()) > 0:
        df["percent_compared"] = df["percent_compared"] * 100.0

    df["pair_key"] = df.apply(
        lambda r: tuple(sorted([str(r["name1"]), str(r["name2"])])), axis=1
    )
    df["consensus_snp_rate"] = df["consensus_SNPs"] / df["compared_bases_count"]
    df["population_snp_rate"] = df["population_SNPs"] / df["compared_bases_count"]
    return df


def load_compares(metadata: pd.DataFrame, instrain_root: str) -> Tuple[pd.DataFrame, List[dict]]:
    rows = []
    missing = []
    for rfmt in sorted(metadata["RFMT_ID"].dropna().unique()):
        path = find_compare_file(instrain_root, rfmt)
        if path is None:
            missing.append({"RFMT_ID": rfmt, "missing_type": "genomeWide_compare", "path_checked": f"{instrain_root}/{rfmt}/compare_out_with_stb"})
            continue
        try:
            raw = read_tsv_auto(path)
            rows.append(standardize_compare(raw, rfmt, path))
        except Exception as e:
            missing.append({"RFMT_ID": rfmt, "missing_type": "compare_read_error", "path_checked": path, "error": str(e)})
    if not rows:
        return pd.DataFrame(), missing
    return pd.concat(rows, ignore_index=True, sort=False), missing


def get_samples_by_role(meta_rfmt: pd.DataFrame, role: str) -> List[str]:
    return sorted(meta_rfmt.loc[meta_rfmt["role"] == role, "sample_name_expected"].dropna().astype(str).unique())


def pair_subset(compare_rfmt: pd.DataFrame, genome: str, s1: str, s2: str) -> pd.DataFrame:
    key = tuple(sorted([s1, s2]))
    sub = compare_rfmt[(compare_rfmt["genome"] == genome) & (compare_rfmt["pair_key"] == key)].copy()
    return sub


def choose_best_pair(compare_rfmt: pd.DataFrame, genome: str, s1: str, s2: str) -> Optional[pd.Series]:
    sub = pair_subset(compare_rfmt, genome, s1, s2)
    if sub.empty:
        return None
    # If duplicate rows exist, choose the row with largest compared_bases_count.
    sub = sub.sort_values("compared_bases_count", ascending=False, na_position="last")
    return sub.iloc[0]


def passes_similarity(row: Optional[pd.Series], args: argparse.Namespace) -> bool:
    if row is None:
        return False
    vals = row
    if pd.isna(vals.get("compared_bases_count", np.nan)) or vals["compared_bases_count"] < args.min_compared_bases:
        return False
    if pd.isna(vals.get("percent_compared", np.nan)) or vals["percent_compared"] < args.min_percent_compared:
        return False
    if pd.isna(vals.get("conANI", np.nan)) or vals["conANI"] < args.min_conani:
        return False
    # If popANI is present, require it. If entirely missing for a row, do not fail.
    if not pd.isna(vals.get("popANI", np.nan)) and vals["popANI"] < args.min_popani:
        return False
    return True


def dist(row: Optional[pd.Series]) -> float:
    if row is None:
        return np.inf
    x = row.get("consensus_snp_rate", np.nan)
    if pd.isna(x):
        # fallback from conANI if SNP fields are missing
        ca = row.get("conANI", np.nan)
        if pd.isna(ca):
            return np.inf
        return max(0.0, 1.0 - float(ca) / 100.0)
    return float(x)


def row_metrics(prefix: str, row: Optional[pd.Series]) -> Dict[str, object]:
    fields = [
        "conANI",
        "popANI",
        "consensus_SNPs",
        "population_SNPs",
        "compared_bases_count",
        "percent_compared",
        "coverage_overlap",
        "consensus_snp_rate",
        "population_snp_rate",
    ]
    out = {}
    for f in fields:
        out[f"{prefix}_{f}"] = np.nan if row is None else row.get(f, np.nan)
    return out


def classify_call(
    donor_row: Optional[pd.Series],
    baseline_post_row: Optional[pd.Series],
    donor_baseline_row: Optional[pd.Series],
    donor_label: str,
    baseline_label: str,
    args: argparse.Namespace,
) -> Tuple[str, str]:
    donor_high = passes_similarity(donor_row, args)
    base_high = passes_similarity(baseline_post_row, args)
    donor_base_high = passes_similarity(donor_baseline_row, args)

    if donor_row is None:
        return "no_donor_post_comparison", "no donor-post in genomeWide_compare"

    if not donor_high:
        if base_high:
            return "recipient_or_previous_state_persistent", "baseline-post passes; donor-post does not"
        return "not_high_confidence", "donor-post fails high-confidence thresholds"

    donor_d = dist(donor_row)
    base_d = dist(baseline_post_row)

    if baseline_post_row is None:
        return f"{donor_label}_engrafted_baseline_missing", "donor-post passes; no baseline-post comparison"

    if donor_base_high and base_high:
        if donor_d + args.snp_rate_margin < base_d:
            return f"{donor_label}_engrafted_but_already_related_baseline", "donor-post closer than baseline-post; donor-baseline also high"
        return "already_shared_or_baseline_persistent", "donor-baseline and baseline-post both pass"

    if base_high and not (donor_d + args.snp_rate_margin < base_d):
        return "recipient_or_previous_state_persistent", "baseline-post passes and is not farther than donor-post"

    if donor_d + args.snp_rate_margin < base_d or not base_high:
        return f"{donor_label}_engrafted", "donor-post passes and is closer/stronger than baseline-post"

    return "ambiguous", "donor-post passes but source assignment is ambiguous"


def build_transfer_calls(metadata: pd.DataFrame, compare: pd.DataFrame, annotations: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    ann_cols = ["genome", "sample_name", "Bin ID", "taxonomy", "genus", "species", "Completeness", "Contamination", "Length"]
    ann = annotations[[c for c in ann_cols if c in annotations.columns]].drop_duplicates("genome") if not annotations.empty else pd.DataFrame(columns=["genome"])
    calls = []

    for rfmt in sorted(metadata["RFMT_ID"].dropna().unique()):
        meta_r = metadata[metadata["RFMT_ID"] == rfmt].copy()
        samples = sorted(meta_r["sample_name_expected"].dropna().astype(str).unique())
        if len(samples) < 2:
            continue
        if not args.include_limited_groups and not any(meta_r["role"].isin(["donor1", "donor2"])):
            continue
        cmp_r = compare[compare["RFMT_ID"] == rfmt].copy()
        if cmp_r.empty:
            continue

        genomes = sorted(cmp_r["genome"].dropna().astype(str).unique())
        donor1s = get_samples_by_role(meta_r, "donor1")
        donor2s = get_samples_by_role(meta_r, "donor2")
        recipient_d0s = get_samples_by_role(meta_r, "recipient_d0")
        pre_donor2s = get_samples_by_role(meta_r, "recipient_pre_donor2")

        donor_configs = []

        # Donor1 analysis:
        # recipient_d0 is the pre-Donor1 baseline.
        # recipient_pre_donor2, when present, is the cleanest post-Donor1/pre-Donor2 sample.
        # day2/day5/day7 are after Donor2, so Donor1-like calls there are interpreted as persistence after Donor2.
        for d in donor1s:
            baseline_samples = recipient_d0s
            if pre_donor2s:
                donor_configs.append((
                    "donor1",
                    d,
                    baseline_samples,
                    DONOR1_DIRECT_POST_ROLES,
                    "donor1_direct_engraftment_pre_donor2",
                ))
            donor_configs.append((
                "donor1",
                d,
                baseline_samples,
                DONOR1_PERSISTENCE_POST_DONOR2_ROLES,
                "donor1_persistence_after_donor2",
            ))

        # Donor2 analysis:
        # recipient_pre_donor2 is the preferred immediate baseline.
        # If missing, recipient_d0 is used only as a weak baseline and is labelled accordingly.
        for d in donor2s:
            baseline_samples = pre_donor2s if pre_donor2s else recipient_d0s
            donor_configs.append((
                "donor2",
                d,
                baseline_samples,
                DONOR2_POST_ROLES,
                "donor2_direct_engraftment",
            ))

        for donor_label, donor_sample, baseline_samples, post_roles, analysis_context in donor_configs:
            baseline_sample = baseline_samples[0] if baseline_samples else None
            if donor_label == "donor1":
                baseline_role = "recipient_d0_pre_donor1" if baseline_sample else "missing"
            else:
                baseline_role = "recipient_pre_donor2" if pre_donor2s else ("recipient_d0_pre_donor1_weak_for_donor2" if baseline_sample else "missing")

            for post_role in post_roles:
                post_samples = get_samples_by_role(meta_r, post_role)
                for post_sample in post_samples:
                    if post_sample == donor_sample:
                        continue
                    for genome in genomes:
                        donor_post = choose_best_pair(cmp_r, genome, donor_sample, post_sample)
                        baseline_post = choose_best_pair(cmp_r, genome, baseline_sample, post_sample) if baseline_sample else None
                        donor_baseline = choose_best_pair(cmp_r, genome, donor_sample, baseline_sample) if baseline_sample else None

                        status, reason = classify_call(donor_post, baseline_post, donor_baseline, donor_label, baseline_role, args)
                        call = {
                            "RFMT_ID": rfmt,
                            "donor_label": donor_label,
                            "donor_sample": donor_sample,
                            "baseline_sample": baseline_sample,
                            "baseline_role": baseline_role,
                            "post_sample": post_sample,
                            "post_role": post_role,
                            "analysis_context": analysis_context,
                            "genome": genome,
                            "classification": status,
                            "classification_reason": reason,
                        }
                        call.update(row_metrics("donor_post", donor_post))
                        call.update(row_metrics("baseline_post", baseline_post))
                        call.update(row_metrics("donor_baseline", donor_baseline))
                        calls.append(call)

    out = pd.DataFrame(calls)
    if out.empty:
        return out
    out = out.merge(ann, on="genome", how="left")
    out["is_engrafted"] = out["classification"].str.contains("engrafted", na=False)
    out["is_high_confidence_engrafted"] = out["classification"].isin(["donor1_engrafted", "donor2_engrafted"])
    out["is_direct_donor1_engraftment"] = (
        (out["donor_label"] == "donor1")
        & (out["post_role"] == "recipient_pre_donor2")
        & out["is_high_confidence_engrafted"]
    )
    out["is_donor1_persistent_after_donor2"] = (
        (out["donor_label"] == "donor1")
        & (out["post_role"].isin(["recipient_d2", "recipient_d5", "recipient_d7"]))
        & out["is_engrafted"]
    )
    out["is_direct_donor2_engraftment"] = (
        (out["donor_label"] == "donor2")
        & (out["post_role"].isin(["recipient_d2", "recipient_d5", "recipient_d7"]))
        & out["is_high_confidence_engrafted"]
    )
    out["is_persistent_to_d7"] = out["is_engrafted"] & (out["post_role"] == "recipient_d7")
    return out


def summarize_species(calls: pd.DataFrame) -> pd.DataFrame:
    if calls.empty:
        return pd.DataFrame()
    df = calls.copy()
    df["species"] = df["species"].fillna("Unclassified")
    # One RFMT + donor + species unit, to prevent inflation by duplicate MAGs.
    unit_cols = ["RFMT_ID", "donor_label", "species"]
    unit = df.groupby(unit_cols, dropna=False).agg(
        n_MAG_rows=("genome", "count"),
        n_unique_MAGs=("genome", "nunique"),
        any_engrafted=("is_engrafted", "max"),
        any_high_conf_engrafted=("is_high_confidence_engrafted", "max"),
        any_persistent_to_d7=("is_persistent_to_d7", "max"),
        best_conANI=("donor_post_conANI", "max"),
        best_popANI=("donor_post_popANI", "max"),
        max_compared_bases=("donor_post_compared_bases_count", "max"),
        best_percent_compared=("donor_post_percent_compared", "max"),
    ).reset_index()

    summary = unit.groupby(["species"], dropna=False).agg(
        n_RFMT_donor_units=("RFMT_ID", "count"),
        n_RFMTs=("RFMT_ID", "nunique"),
        n_units_engrafted=("any_engrafted", "sum"),
        n_units_high_conf_engrafted=("any_high_conf_engrafted", "sum"),
        n_units_persistent_to_d7=("any_persistent_to_d7", "sum"),
        best_conANI=("best_conANI", "max"),
        best_popANI=("best_popANI", "max"),
        max_compared_bases=("max_compared_bases", "max"),
    ).reset_index()

    d1 = unit[unit["donor_label"] == "donor1"].groupby("species")["any_engrafted"].sum().rename("n_donor1_engrafted_units")
    d2 = unit[unit["donor_label"] == "donor2"].groupby("species")["any_engrafted"].sum().rename("n_donor2_engrafted_units")
    summary = summary.merge(d1, on="species", how="left").merge(d2, on="species", how="left")
    summary[["n_donor1_engrafted_units", "n_donor2_engrafted_units"]] = summary[["n_donor1_engrafted_units", "n_donor2_engrafted_units"]].fillna(0).astype(int)
    summary["transmission_rate_among_tested_units"] = summary["n_units_engrafted"] / summary["n_RFMT_donor_units"].replace(0, np.nan)
    summary = summary.sort_values(["n_units_high_conf_engrafted", "n_units_engrafted", "n_units_persistent_to_d7"], ascending=False)
    return summary


def summarize_rfmt(calls: pd.DataFrame) -> pd.DataFrame:
    if calls.empty:
        return pd.DataFrame()
    df = calls.copy()
    out = df.groupby(["RFMT_ID", "donor_label"], dropna=False).agg(
        n_rows=("genome", "count"),
        n_MAGs_tested=("genome", "nunique"),
        n_engrafted_rows=("is_engrafted", "sum"),
        n_high_conf_engrafted_rows=("is_high_confidence_engrafted", "sum"),
        n_persistent_to_d7_rows=("is_persistent_to_d7", "sum"),
        n_species_tested=("species", "nunique"),
    ).reset_index()

    # Count unique MAGs and species that are engrafted.
    eng = df[df["is_engrafted"]]
    mag_eng = eng.groupby(["RFMT_ID", "donor_label"])["genome"].nunique().rename("n_MAGs_engrafted")
    sp_eng = eng.groupby(["RFMT_ID", "donor_label"])["species"].nunique().rename("n_species_engrafted")
    d7 = df[df["is_persistent_to_d7"]].groupby(["RFMT_ID", "donor_label"])["genome"].nunique().rename("n_MAGs_persistent_to_d7")
    out = out.merge(mag_eng, on=["RFMT_ID", "donor_label"], how="left").merge(sp_eng, on=["RFMT_ID", "donor_label"], how="left").merge(d7, on=["RFMT_ID", "donor_label"], how="left")
    for c in ["n_MAGs_engrafted", "n_species_engrafted", "n_MAGs_persistent_to_d7"]:
        out[c] = out[c].fillna(0).astype(int)
    out["MAG_engraftment_rate"] = out["n_MAGs_engrafted"] / out["n_MAGs_tested"].replace(0, np.nan)
    out["species_engraftment_rate"] = out["n_species_engrafted"] / out["n_species_tested"].replace(0, np.nan)
    return out.sort_values(["RFMT_ID", "donor_label"])


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    metadata = read_tsv_auto(args.metadata)
    required = {"RFMT_ID", "sample_name_expected", "role"}
    missing_cols = required - set(metadata.columns)
    if missing_cols:
        raise ValueError(f"Metadata missing columns: {sorted(missing_cols)}")

    # Remove rows without sample, normalize sample names, ignore one-sample groups by group size.
    metadata = metadata.dropna(subset=["RFMT_ID", "sample_name_expected", "role"]).copy()
    metadata["RFMT_ID"] = metadata["RFMT_ID"].astype(str)
    metadata["sample_name_expected"] = metadata.apply(
        lambda r: make_full_sample_id(str(r["RFMT_ID"]), r["sample_name_expected"]), axis=1
    )
    metadata["role"] = metadata["role"].astype(str)
    group_n = metadata.groupby("RFMT_ID")["sample_name_expected"].nunique()
    keep_rfmts = group_n[group_n >= 2].index
    metadata = metadata[metadata["RFMT_ID"].isin(keep_rfmts)].copy()

    print(f"Loaded metadata rows after >=2-sample filtering: {len(metadata)}")
    print(f"RFMT groups retained: {metadata['RFMT_ID'].nunique()}")

    annotations, missing_bintables = load_bintables(metadata, args.mag_root)
    print(f"Loaded bintable annotation rows: {len(annotations)}")

    compare, missing_compares = load_compares(metadata, args.instrain_root)
    print(f"Loaded genomeWide compare rows: {len(compare)}")

    missing = pd.DataFrame(missing_bintables + missing_compares)
    missing.to_csv(outdir / "missing_inputs.tsv", sep="\t", index=False)

    if not annotations.empty:
        annotations.to_csv(outdir / "bintable_mag_annotation.tsv", sep="\t", index=False)

    if compare.empty:
        raise SystemExit("No compare tables loaded. Check --instrain-root and compare output paths.")

    calls = build_transfer_calls(metadata, compare, annotations, args)
    calls.to_csv(outdir / "mag_transfer_calls.tsv", sep="\t", index=False)
    print(f"Wrote MAG transfer calls: {len(calls)}")

    species = summarize_species(calls)
    species.to_csv(outdir / "species_transmission_summary.tsv", sep="\t", index=False)
    print(f"Wrote species summary rows: {len(species)}")

    rfmt = summarize_rfmt(calls)
    rfmt.to_csv(outdir / "rfmt_transfer_summary.tsv", sep="\t", index=False)
    print(f"Wrote RFMT summary rows: {len(rfmt)}")

    # A compact high-confidence table for quick inspection.
    if not calls.empty:
        hc = calls[calls["is_engrafted"]].copy()
        hc = hc.sort_values(["RFMT_ID", "donor_label", "post_role", "species", "donor_post_compared_bases_count"], ascending=[True, True, True, True, False])
        hc.to_csv(outdir / "engrafted_MAGs_only.tsv", sep="\t", index=False)
        print(f"Wrote engrafted-only table rows: {len(hc)}")

    print("Done.")


if __name__ == "__main__":
    main()
