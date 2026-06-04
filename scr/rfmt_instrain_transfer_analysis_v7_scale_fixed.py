#!/usr/bin/env python3
"""
RFMT strain-transmission analysis with:
1) MAG quality filtering
2) correct donor-available engraftment denominator
3) antibiotic metadata merge
4) per-sample MAG counts
5) TPM matrix for top MAG visualization

Main denominator:
    MAG_engraftment_rate = n_engrafted donor MAGs / n_quality_filtered donor MAGs available in donor bintable

Required inputs:
    --metadata RFMT_analysis_sample_long.tsv
    --instrain-root /group/sbms004/yxia/GUT/inStrain_final
    --mag-root /group/sbms004/yxia/GUT/MAG_per_sample
    --antibiotics RFMT_antibiotic_indicators_used.tsv
"""

import argparse
import re
from pathlib import Path
from typing import Optional, Dict, List, Tuple

import numpy as np
import pandas as pd


# -----------------------------
# General helpers
# -----------------------------

def read_tsv(path: Path, dtype=str) -> pd.DataFrame:
    if not path.exists():
        print(f"[MISS] {path}")
        return pd.DataFrame()
    print(f"[READ] {path}")
    return pd.read_csv(path, sep="\t", dtype=dtype)


def write_tsv(df: pd.DataFrame, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    print(f"[WRITE] {path} rows={len(df)}")


def clean_numeric(series):
    return pd.to_numeric(series, errors="coerce")


def guess_col(df: pd.DataFrame, candidates: List[str], required: bool = False) -> Optional[str]:
    for c in candidates:
        if c in df.columns:
            return c
    if required:
        raise ValueError(f"Cannot find any of {candidates}. Available columns: {list(df.columns)}")
    return None


def norm_text(x) -> str:
    if pd.isna(x):
        return ""
    return str(x).strip()


def make_full_sample_id(rfmt: str, sample: str) -> str:
    sample = norm_text(sample)
    if sample.endswith(".0"):
        sample = sample[:-2]
    if sample.startswith(rfmt + "_"):
        return sample
    if sample.startswith("RFMT"):
        return sample
    return f"{rfmt}_{sample}"


def parse_rfmt_from_sample(sample: str) -> str:
    m = re.match(r"(RFMT\d{3})", str(sample))
    return m.group(1) if m else ""


def extract_sample_id_from_instrain_name(x: str) -> str:
    """
    Robustly recover sample IDs like RFMT015_5423 from inStrain name1/name2 values.

    Handles examples such as:
      RFMT015_5423
      RFMT015_5423_IS
      /path/to/RFMT015_5423_IS
      RFMT015_5423.sorted.bam
      RFMT015_5423.profile
    """
    s = str(x).strip()
    m = re.search(r"(RFMT\d{3}_\d+)", s)
    if m:
        return m.group(1)
    name = Path(s).name
    name = name.replace("_IS", "")
    name = re.sub(r"\.sorted\.bam$|\.bam$|\.profile$|\.IS$", "", name)
    return name


def genome_aliases(genome: str) -> set:
    """
    Robust genome matching between donor bintable names and inStrain genome names.

    Possible forms:
      donor bintable-derived: RFMT015_5423_RFMT015_5423.metabat2.10
      inStrain genome:        RFMT015_5423_RFMT015_5423.metabat2.10
      or sometimes:           RFMT015_5423.metabat2.10

    Return likely aliases so exact naming differences do not kill all comparisons.
    """
    g = str(genome).strip()
    aliases = {g}

    m = re.match(r"^(RFMT\d{3}_\d+)_(RFMT\d{3}_\d+\..+)$", g)
    if m:
        sample, binid = m.group(1), m.group(2)
        aliases.add(binid)
        aliases.add(f"{sample}_{binid}")

    m2 = re.match(r"^(RFMT\d{3}_\d+\..+)$", g)
    if m2:
        binid = m2.group(1)
        sample = binid.split(".")[0]
        aliases.add(f"{sample}_{binid}")

    return aliases


def short_species(x):
    if pd.isna(x) or str(x).strip() == "":
        return "Unknown"
    x = str(x).strip()
    if ";" in x:
        parts = [p.strip() for p in x.split(";")]
        sp = [p for p in parts if p.startswith("s__")]
        if sp:
            return sp[-1].replace("s__", "").strip() or "Unknown"
        gen = [p for p in parts if p.startswith("g__")]
        if gen:
            return gen[-1].replace("g__", "").strip() + " sp."
    if x.startswith("s__"):
        return x.replace("s__", "", 1).strip() or "Unknown"
    return x


# -----------------------------
# Metadata roles
# -----------------------------

ROLE_ORDER = {
    "donor1": 0,
    "recipient_day0": 1,
    "recipient_pre_donor2": 2,
    "donor2": 3,
    "recipient_d2": 4,
    "recipient_d5": 5,
    "recipient_d7": 6,
}

POST_ROLES_DONOR1_DIRECT = ["recipient_pre_donor2"]
POST_ROLES_DONOR1_LATE = ["recipient_d2", "recipient_d5", "recipient_d7"]
POST_ROLES_DONOR2 = ["recipient_d2", "recipient_d5", "recipient_d7"]


def canonical_role(x: str) -> str:
    """
    Convert metadata role labels into canonical roles.

    Important:
    Check recipient_pre_donor2 BEFORE donor2.
    Otherwise strings like 'recipient_pre_donor2' are incorrectly classified as donor2.
    """
    s = norm_text(x).lower()
    s2 = s.replace("-", "_").replace(" ", "_")

    # Recipient baselines / timepoints first.
    if "pre_donor2" in s2 or "predonor2" in s2 or "pre_donor_2" in s2 or s2 == "recipient_pre_donor2":
        return "recipient_pre_donor2"

    if "day0" in s2 or s2.endswith("_d0") or s2 == "d0" or "recipient_d0" in s2 or s2 == "recipient_day0":
        return "recipient_day0"

    if "day2" in s2 or s2.endswith("_d2") or s2 == "d2" or "recipient_d2" in s2:
        return "recipient_d2"

    if "day5" in s2 or s2.endswith("_d5") or s2 == "d5" or "recipient_d5" in s2:
        return "recipient_d5"

    if "day7" in s2 or s2.endswith("_d7") or s2 == "d7" or "recipient_d7" in s2:
        return "recipient_d7"

    # Donors after recipient timepoints, so pre_donor2 is not swallowed by donor2.
    if s2 in {"donor1", "donor_1", "dnr1", "dnr_1", "d1"} or "donor1" in s2 or "donor_1" in s2 or "dnr#1" in s:
        return "donor1"

    if s2 in {"donor2", "donor_2", "dnr2", "dnr_2"} or "donor2" in s2 or "donor_2" in s2 or "dnr#2" in s:
        return "donor2"

    return s2

def load_metadata(metadata_path: Path) -> pd.DataFrame:
    meta = pd.read_csv(metadata_path, sep="\t", dtype=str)

    rfmt_col = guess_col(meta, ["RFMT_ID", "rfmt", "RFMT", "Group", "group"], required=True)
    sample_col = guess_col(meta, ["sample_id", "SampleID", "sample", "Sample", "sample_number", "Sample Number"], required=True)
    role_col = guess_col(meta, ["role", "sample_role", "Role", "timepoint", "Timepoint"], required=True)

    out = meta.copy()
    out["RFMT_ID"] = out[rfmt_col].astype(str).str.strip()
    out["sample_raw"] = out[sample_col].astype(str).str.strip()
    out["sample_name"] = [make_full_sample_id(r, s) for r, s in zip(out["RFMT_ID"], out["sample_raw"])]
    out["role"] = out[role_col].map(canonical_role)
    out["role_order"] = out["role"].map(ROLE_ORDER).fillna(999).astype(int)

    out = out[out["RFMT_ID"].str.match(r"RFMT\d{3}", na=False)].copy()
    out = out[out["sample_name"].str.match(r"RFMT\d{3}_", na=False)].copy()

    # Drop duplicate sample-role rows.
    out = out.drop_duplicates(subset=["RFMT_ID", "sample_name", "role"])
    return out[["RFMT_ID", "sample_name", "role", "role_order"]].sort_values(["RFMT_ID", "role_order", "sample_name"])


# -----------------------------
# File discovery and table parsing
# -----------------------------

def find_compare_file(instrain_root: Path, rfmt: str) -> Optional[Path]:
    compare_dir = instrain_root / rfmt / "compare_out_with_stb" / "output"
    if not compare_dir.exists():
        return None
    candidates = sorted(compare_dir.glob("*genomeWide_compare.tsv"))
    if not candidates:
        candidates = sorted(compare_dir.glob("*genome*compare*.tsv"))
    return candidates[0] if candidates else None


def find_bintable(mag_root: Path, sample_name: str) -> Optional[Path]:
    candidates = sorted(mag_root.glob(f"sample_list_*/{sample_name}/{sample_name}/results/*.bintable"))
    return candidates[0] if candidates else None


def read_bintable(path: Path) -> pd.DataFrame:
    # Some SqueezeMeta bintables have one descriptive line before the header.
    tries = []
    for header in [0, 1]:
        try:
            df = pd.read_csv(path, sep="\t", header=header, dtype=str)
            tries.append((header, list(df.columns)))
            if "Bin ID" in df.columns:
                return df
        except Exception:
            continue
    raise ValueError(f"Could not read bintable with 'Bin ID': {path}; tried {tries}")


def normalize_bin_to_genome(sample_name: str, bin_id: str) -> str:
    return f"{sample_name}_{str(bin_id).strip()}"


def parse_bintables(meta: pd.DataFrame,
                    mag_root: Path,
                    min_completeness: float,
                    max_contamination: float) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Returns:
      mag_annot: one row per MAG from all bintables
      donor_available: quality-passed donor MAG denominator
      sample_mag_counts: per-sample MAG count and TPM summary
      tpm_long: long TPM matrix rows from bintables
    """
    all_mag_rows = []
    donor_avail_rows = []
    sample_count_rows = []
    tpm_rows = []
    missing = []

    for _, row in meta.iterrows():
        rfmt = row["RFMT_ID"]
        sample = row["sample_name"]
        role = row["role"]

        bintable = find_bintable(mag_root, sample)
        if bintable is None:
            missing.append({"RFMT_ID": rfmt, "sample_name": sample, "role": role, "missing": "bintable"})
            continue

        try:
            bt = read_bintable(bintable)
        except Exception as e:
            missing.append({"RFMT_ID": rfmt, "sample_name": sample, "role": role, "missing": f"bad_bintable:{e}"})
            continue

        # Required/optional columns
        comp_col = "Completeness" if "Completeness" in bt.columns else None
        cont_col = "Contamination" if "Contamination" in bt.columns else None
        length_col = "Length" if "Length" in bt.columns else None
        tax_col = "Tax GTDB-Tk" if "Tax GTDB-Tk" in bt.columns else ("Tax" if "Tax" in bt.columns else None)

        bt["Completeness_num"] = clean_numeric(bt[comp_col]) if comp_col else np.nan
        bt["Contamination_num"] = clean_numeric(bt[cont_col]) if cont_col else np.nan
        bt["Length_num"] = clean_numeric(bt[length_col]) if length_col else np.nan

        bt["quality_pass"] = (
            (bt["Completeness_num"] >= min_completeness) &
            (bt["Contamination_num"] <= max_contamination)
        )

        tpm_cols = [c for c in bt.columns if c.startswith("TPM ")]
        cov_cols = [c for c in bt.columns if c.startswith("Coverage ")]

        n_total = len(bt)
        n_quality = int(bt["quality_pass"].sum())
        sum_tpm_quality = 0.0

        for _, b in bt.iterrows():
            bin_id = norm_text(b["Bin ID"])
            genome = normalize_bin_to_genome(sample, bin_id)
            tax = norm_text(b[tax_col]) if tax_col else ""
            species = short_species(tax)

            tpm_self = np.nan
            tpm_self_col = f"TPM {sample}"
            if tpm_self_col in bt.columns:
                tpm_self = pd.to_numeric(b.get(tpm_self_col), errors="coerce")
                if bool(b["quality_pass"]) and pd.notna(tpm_self):
                    sum_tpm_quality += float(tpm_self)

            annot_row = {
                "RFMT_ID": rfmt,
                "sample_name": sample,
                "role": role,
                "genome": genome,
                "Bin ID": bin_id,
                "taxonomy": tax,
                "species": species,
                "Completeness": b["Completeness_num"],
                "Contamination": b["Contamination_num"],
                "Length": b["Length_num"],
                "quality_pass": bool(b["quality_pass"]),
                "bintable": str(bintable),
                "TPM_self": tpm_self,
            }
            all_mag_rows.append(annot_row)

            if role in {"donor1", "donor2"} and bool(b["quality_pass"]):
                donor_avail_rows.append({
                    "RFMT_ID": rfmt,
                    "donor_label": role,
                    "donor_sample": sample,
                    **annot_row
                })

            # Long TPM values. If only one TPM column exists, it still works.
            for tpm_col in tpm_cols:
                tpm_sample = tpm_col.replace("TPM ", "", 1).strip()
                tpm_val = pd.to_numeric(b.get(tpm_col), errors="coerce")
                tpm_rows.append({
                    "RFMT_ID": parse_rfmt_from_sample(tpm_sample) or rfmt,
                    "source_bintable_sample": sample,
                    "sample_name": tpm_sample,
                    "source_role": role,
                    "genome": genome,
                    "species": species,
                    "TPM": tpm_val,
                    "quality_pass": bool(b["quality_pass"]),
                })

        sample_count_rows.append({
            "RFMT_ID": rfmt,
            "sample_name": sample,
            "role": role,
            "role_order": ROLE_ORDER.get(role, 999),
            "n_MAGs_total": n_total,
            "n_MAGs_quality": n_quality,
            "sum_TPM_quality_self": sum_tpm_quality,
            "bintable": str(bintable),
        })

    mag_annot = pd.DataFrame(all_mag_rows)
    donor_available = pd.DataFrame(donor_avail_rows)
    sample_mag_counts = pd.DataFrame(sample_count_rows)
    tpm_long = pd.DataFrame(tpm_rows)
    missing_df = pd.DataFrame(missing)
    return mag_annot, donor_available, sample_mag_counts, tpm_long, missing_df


# -----------------------------
# inStrain comparison parsing
# -----------------------------

def load_compare_tables(meta: pd.DataFrame, instrain_root: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    missing = []

    for rfmt in sorted(meta["RFMT_ID"].unique()):
        f = find_compare_file(instrain_root, rfmt)
        if f is None:
            missing.append({"RFMT_ID": rfmt, "missing": "genomeWide_compare"})
            continue
        try:
            df = pd.read_csv(f, sep="\t", dtype=str)
        except Exception as e:
            missing.append({"RFMT_ID": rfmt, "missing": f"bad_compare:{e}"})
            continue

        df["RFMT_ID"] = rfmt
        df["compare_file"] = str(f)
        rows.append(df)

    if rows:
        cmp = pd.concat(rows, ignore_index=True)
    else:
        cmp = pd.DataFrame()

    return cmp, pd.DataFrame(missing)


def standardize_compare(cmp: pd.DataFrame) -> pd.DataFrame:
    if cmp.empty:
        return cmp

    out = cmp.copy()
    rename = {}

    for old, new in [
        ("genome", "genome"),
        ("name1", "sample1"),
        ("name2", "sample2"),
        ("conANI", "conANI"),
        ("popANI", "popANI"),
        ("consensus_SNPs", "consensus_SNPs"),
        ("population_SNPs", "population_SNPs"),
        ("compared_bases_count", "compared_bases_count"),
        ("percent_compared", "percent_compared"),
        ("coverage_overlap", "coverage_overlap"),
    ]:
        if old in out.columns:
            rename[old] = new

    out = out.rename(columns=rename)

    required = ["RFMT_ID", "genome", "sample1", "sample2"]
    missing = [c for c in required if c not in out.columns]
    if missing:
        raise ValueError(f"Compare table missing required columns {missing}. Available: {list(out.columns)}")

    for c in ["conANI", "popANI", "consensus_SNPs", "population_SNPs", "compared_bases_count", "percent_compared", "coverage_overlap"]:
        if c in out.columns:
            out[c] = pd.to_numeric(out[c], errors="coerce")

    # Normalize inStrain profile names to real sample IDs, e.g. RFMT015_5423.
    out["sample1_raw"] = out["sample1"].astype(str)
    out["sample2_raw"] = out["sample2"].astype(str)
    out["sample1"] = out["sample1_raw"].map(extract_sample_id_from_instrain_name)
    out["sample2"] = out["sample2_raw"].map(extract_sample_id_from_instrain_name)

    return out


def get_pair_metrics(cmp: pd.DataFrame, rfmt: str, genome: str, a: str, b: str) -> Dict[str, object]:
    if cmp.empty or not a or not b:
        return {}

    aliases = genome_aliases(genome)
    sub = cmp[(cmp["RFMT_ID"] == rfmt) & (cmp["genome"].astype(str).isin(aliases))].copy()
    if sub.empty:
        return {}

    a2 = extract_sample_id_from_instrain_name(a)
    b2 = extract_sample_id_from_instrain_name(b)

    hit = sub[
        ((sub["sample1"] == a2) & (sub["sample2"] == b2)) |
        ((sub["sample1"] == b2) & (sub["sample2"] == a2))
    ]
    if hit.empty:
        return {}

    # If duplicate rows, use row with largest compared_bases_count.
    if "compared_bases_count" in hit.columns:
        hit = hit.sort_values("compared_bases_count", ascending=False)
    r = hit.iloc[0].to_dict()
    return r

def scaled_cutoff(value, cutoff, metric_name: str) -> float:
    """
    inStrain versions can report ANI / percent_compared either as:
      ANI: 0.99999 or 99.999
      percent_compared: 0.50 or 50
    This makes command-line cutoffs robust.
    """
    v = pd.to_numeric(value, errors="coerce")
    c = pd.to_numeric(cutoff, errors="coerce")
    if pd.isna(v) or pd.isna(c):
        return c

    if metric_name in {"conANI", "popANI", "percent_compared"}:
        # Output is fraction, user supplied percent.
        if v <= 1.5 and c > 1.5:
            return c / 100.0
        # Output is percent, user supplied fraction.
        if v > 1.5 and c <= 1.5:
            return c * 100.0

    return c


def pass_instrain_cutoff(metrics: Dict[str, object],
                         min_conani: float,
                         min_popani: float,
                         min_compared_bases: float,
                         min_percent_compared: float) -> bool:
    if not metrics:
        return False

    con = pd.to_numeric(metrics.get("conANI"), errors="coerce")
    pop = pd.to_numeric(metrics.get("popANI"), errors="coerce")
    bases = pd.to_numeric(metrics.get("compared_bases_count"), errors="coerce")
    pct = pd.to_numeric(metrics.get("percent_compared"), errors="coerce")

    if pd.isna(con) or pd.isna(pop) or pd.isna(bases) or pd.isna(pct):
        return False

    con_cut = scaled_cutoff(con, min_conani, "conANI")
    pop_cut = scaled_cutoff(pop, min_popani, "popANI")
    pct_cut = scaled_cutoff(pct, min_percent_compared, "percent_compared")

    return (
        con >= con_cut and
        pop >= pop_cut and
        bases >= min_compared_bases and
        pct >= pct_cut
    )


def metric_prefix(metrics: Dict[str, object], prefix: str) -> Dict[str, object]:
    fields = ["conANI", "popANI", "consensus_SNPs", "population_SNPs", "compared_bases_count", "percent_compared", "coverage_overlap"]
    out = {}
    for f in fields:
        out[f"{prefix}_{f}"] = metrics.get(f, np.nan) if metrics else np.nan

    bases = pd.to_numeric(out.get(f"{prefix}_compared_bases_count"), errors="coerce")
    cs = pd.to_numeric(out.get(f"{prefix}_consensus_SNPs"), errors="coerce")
    ps = pd.to_numeric(out.get(f"{prefix}_population_SNPs"), errors="coerce")

    out[f"{prefix}_consensus_snp_rate"] = cs / bases if pd.notna(cs) and pd.notna(bases) and bases > 0 else np.nan
    out[f"{prefix}_population_snp_rate"] = ps / bases if pd.notna(ps) and pd.notna(bases) and bases > 0 else np.nan
    return out


# -----------------------------
# Transmission calling
# -----------------------------

def build_role_map(meta: pd.DataFrame, rfmt: str) -> Dict[str, List[str]]:
    sub = meta[meta["RFMT_ID"] == rfmt]
    d = {}
    for role, ss in sub.groupby("role")["sample_name"]:
        d[role] = sorted(ss.unique())
    return d


def classify_call(donor_label: str,
                  post_role: str,
                  donor_post_pass: bool,
                  baseline_post_pass: bool,
                  donor_baseline_pass: bool,
                  has_donor_post: bool,
                  has_baseline: bool) -> Tuple[str, str, bool]:
    """
    Strict baseline-absent engraftment logic.

    Engrafted requires:
      1. donor-post passes strain cutoff
      2. donor-baseline does NOT pass cutoff
      3. baseline-post does NOT pass cutoff

    In words:
      the donor-like strain is present after transfer,
      but not already present / not already donor-like in the relevant baseline.

    Relevant baseline:
      Donor1: recipient day0, pre-Donor1 baseline
      Donor2: recipient pre-Donor2, or day0 if pre-Donor2 is unavailable

    Baseline-related calls are NOT counted as engrafted:
      donor_baseline_pass -> already_shared
      baseline_post_pass -> recipient_persistent or ambiguous
    """
    if donor_post_pass and (not donor_baseline_pass) and (not baseline_post_pass):
        if donor_label == "donor1" and post_role == "recipient_pre_donor2":
            return "donor1_engrafted_pre_donor2", "Donor1-post passes cutoff; day0 baseline does not match donor or post", True
        if donor_label == "donor1" and post_role in POST_ROLES_DONOR1_LATE:
            return "donor1_late_detected_engrafted", "Donor1-post passes cutoff after Donor2 follow-up; relevant baseline does not match", True
        if donor_label == "donor2":
            return "donor2_engrafted", "Donor2-post passes cutoff; pre-Donor2 baseline does not match donor or post", True
        return f"{donor_label}_engrafted", "donor-post passes cutoff and baseline does not explain the post-transfer strain", True

    if donor_post_pass and donor_baseline_pass:
        return f"{donor_label}_already_shared", "donor and baseline were already highly similar before this donor exposure", False

    if donor_post_pass and baseline_post_pass:
        return f"{donor_label}_ambiguous_baseline_like", "post-transfer strain is similar to both donor and baseline", False

    if baseline_post_pass and not donor_post_pass:
        return "recipient_persistent", "baseline-post passes cutoff but donor-post does not", False

    if has_donor_post:
        return "not_engrafted", "donor-post comparison available but did not pass cutoff", False

    return "not_engrafted_no_valid_post_comparison", "donor MAG available but no valid donor-post genome-wide comparison", False

def make_transfer_calls(meta: pd.DataFrame,
                        donor_available: pd.DataFrame,
                        mag_annot: pd.DataFrame,
                        cmp: pd.DataFrame,
                        min_conani: float,
                        min_popani: float,
                        min_compared_bases: float,
                        min_percent_compared: float) -> pd.DataFrame:
    if donor_available.empty:
        return pd.DataFrame()

    annotation_cols = ["genome", "Bin ID", "taxonomy", "species", "Completeness", "Contamination", "Length", "quality_pass", "TPM_self"]
    annot = mag_annot[annotation_cols].drop_duplicates("genome") if not mag_annot.empty else pd.DataFrame(columns=annotation_cols)

    calls = []
    for _, drow in donor_available.iterrows():
        rfmt = drow["RFMT_ID"]
        donor_label = drow["donor_label"]
        donor_sample = drow["donor_sample"]
        genome = drow["genome"]

        role_map = build_role_map(meta, rfmt)

        if donor_label == "donor1":
            baseline_role = "recipient_day0"
            post_roles = POST_ROLES_DONOR1_DIRECT + POST_ROLES_DONOR1_LATE
        elif donor_label == "donor2":
            baseline_role = "recipient_pre_donor2" if role_map.get("recipient_pre_donor2") else "recipient_day0"
            post_roles = POST_ROLES_DONOR2
        else:
            continue

        baseline_samples = role_map.get(baseline_role, [])
        baseline_sample = baseline_samples[0] if baseline_samples else ""

        for post_role in post_roles:
            post_samples = role_map.get(post_role, [])
            for post_sample in post_samples:
                donor_post = get_pair_metrics(cmp, rfmt, genome, donor_sample, post_sample)
                baseline_post = get_pair_metrics(cmp, rfmt, genome, baseline_sample, post_sample) if baseline_sample else {}
                donor_baseline = get_pair_metrics(cmp, rfmt, genome, donor_sample, baseline_sample) if baseline_sample else {}

                donor_post_pass = pass_instrain_cutoff(donor_post, min_conani, min_popani, min_compared_bases, min_percent_compared)
                baseline_post_pass = pass_instrain_cutoff(baseline_post, min_conani, min_popani, min_compared_bases, min_percent_compared)
                donor_baseline_pass = pass_instrain_cutoff(donor_baseline, min_conani, min_popani, min_compared_bases, min_percent_compared)

                classification, reason, is_engrafted = classify_call(
                    donor_label=donor_label,
                    post_role=post_role,
                    donor_post_pass=donor_post_pass,
                    baseline_post_pass=baseline_post_pass,
                    donor_baseline_pass=donor_baseline_pass,
                    has_donor_post=bool(donor_post),
                    has_baseline=bool(baseline_sample),
                )

                is_high_conf = bool(is_engrafted and donor_post_pass)
                is_persistent_to_d7 = bool(is_high_conf and post_role == "recipient_d7")

                row = {
                    "RFMT_ID": rfmt,
                    "donor_label": donor_label,
                    "donor_sample": donor_sample,
                    "baseline_sample": baseline_sample,
                    "baseline_role": baseline_role,
                    "post_sample": post_sample,
                    "post_role": post_role,
                    "genome": genome,
                    "classification": classification,
                    "classification_reason": reason,
                    "is_engrafted": is_engrafted,
                    "is_high_confidence_engrafted": is_high_conf,
                    "is_persistent_to_d7": is_persistent_to_d7,
                    "donor_MAG_available": True,
                }
                row.update(metric_prefix(donor_post, "donor_post"))
                row.update(metric_prefix(baseline_post, "baseline_post"))
                row.update(metric_prefix(donor_baseline, "donor_baseline"))
                calls.append(row)

    out = pd.DataFrame(calls)
    if out.empty:
        return out

    out = out.merge(annot, on="genome", how="left")
    return out


# -----------------------------
# Summaries
# -----------------------------

def summarize_rfmt(calls: pd.DataFrame, donor_available: pd.DataFrame) -> pd.DataFrame:
    rows = []

    if donor_available.empty:
        return pd.DataFrame()

    keys = donor_available[["RFMT_ID", "donor_label"]].drop_duplicates()

    for _, k in keys.iterrows():
        rfmt = k["RFMT_ID"]
        donor = k["donor_label"]
        av = donor_available[(donor_available["RFMT_ID"] == rfmt) & (donor_available["donor_label"] == donor)]
        sub = calls[(calls["RFMT_ID"] == rfmt) & (calls["donor_label"] == donor)] if not calls.empty else pd.DataFrame()

        available_genomes = set(av["genome"].dropna())
        available_species = set(av["species"].dropna())

        eng_genomes = set(sub.loc[sub["is_engrafted"].astype(bool), "genome"].dropna()) if not sub.empty else set()
        high_genomes = set(sub.loc[sub["is_high_confidence_engrafted"].astype(bool), "genome"].dropna()) if not sub.empty else set()
        d7_genomes = set(sub.loc[sub["is_persistent_to_d7"].astype(bool), "genome"].dropna()) if not sub.empty else set()

        eng_species = set(sub.loc[sub["is_engrafted"].astype(bool), "species"].dropna()) if not sub.empty and "species" in sub.columns else set()
        d7_species = set(sub.loc[sub["is_persistent_to_d7"].astype(bool), "species"].dropna()) if not sub.empty and "species" in sub.columns else set()

        tested_genomes = set(sub.loc[sub["donor_post_compared_bases_count"].notna(), "genome"].dropna()) if not sub.empty else set()

        rows.append({
            "RFMT_ID": rfmt,
            "donor_label": donor,
            "n_MAGs_available": len(available_genomes),
            "n_MAGs_tested_comparable": len(tested_genomes),
            "n_MAGs_engrafted": len(eng_genomes),
            "n_MAGs_high_conf_engrafted": len(high_genomes),
            "n_MAGs_persistent_to_d7": len(d7_genomes),
            "n_species_available": len(available_species),
            "n_species_engrafted": len(eng_species),
            "n_species_persistent_to_d7": len(d7_species),
            "MAG_engraftment_rate": len(eng_genomes) / len(available_genomes) if available_genomes else np.nan,
            "MAG_engraftment_rate_among_comparable": len(eng_genomes) / len(tested_genomes) if tested_genomes else np.nan,
            "species_engraftment_rate": len(eng_species) / len(available_species) if available_species else np.nan,
        })

    return pd.DataFrame(rows)


def summarize_species(calls: pd.DataFrame, donor_available: pd.DataFrame) -> pd.DataFrame:
    if donor_available.empty:
        return pd.DataFrame()

    # Available denominator at RFMT+donor+species level
    av_units = donor_available[["RFMT_ID", "donor_label", "species"]].drop_duplicates()
    av_units = av_units[av_units["species"].notna() & (av_units["species"] != "")].copy()

    if calls.empty:
        calls = pd.DataFrame(columns=["RFMT_ID", "donor_label", "species", "is_engrafted", "is_high_confidence_engrafted", "is_persistent_to_d7"])

    eng = calls[calls.get("is_engrafted", pd.Series(False, index=calls.index)).astype(bool)].copy() if not calls.empty else pd.DataFrame()
    high = calls[calls.get("is_high_confidence_engrafted", pd.Series(False, index=calls.index)).astype(bool)].copy() if not calls.empty else pd.DataFrame()
    d7 = calls[calls.get("is_persistent_to_d7", pd.Series(False, index=calls.index)).astype(bool)].copy() if not calls.empty else pd.DataFrame()

    rows = []
    for species, av_sp in av_units.groupby("species"):
        available_units = av_sp[["RFMT_ID", "donor_label"]].drop_duplicates()
        n_av = len(available_units)

        eng_sp = eng[eng["species"] == species] if not eng.empty else pd.DataFrame()
        high_sp = high[high["species"] == species] if not high.empty else pd.DataFrame()
        d7_sp = d7[d7["species"] == species] if not d7.empty else pd.DataFrame()

        eng_units = eng_sp[["RFMT_ID", "donor_label"]].drop_duplicates() if not eng_sp.empty else pd.DataFrame(columns=["RFMT_ID", "donor_label"])
        high_units = high_sp[["RFMT_ID", "donor_label"]].drop_duplicates() if not high_sp.empty else pd.DataFrame(columns=["RFMT_ID", "donor_label"])
        d7_units = d7_sp[["RFMT_ID", "donor_label"]].drop_duplicates() if not d7_sp.empty else pd.DataFrame(columns=["RFMT_ID", "donor_label"])

        best_con = pd.to_numeric(eng_sp.get("donor_post_conANI"), errors="coerce").max() if not eng_sp.empty else np.nan
        best_pop = pd.to_numeric(eng_sp.get("donor_post_popANI"), errors="coerce").max() if not eng_sp.empty else np.nan
        max_bases = pd.to_numeric(eng_sp.get("donor_post_compared_bases_count"), errors="coerce").max() if not eng_sp.empty else np.nan

        rows.append({
            "species": species,
            "n_RFMT_donor_units_available": n_av,
            "n_RFMTs_available": av_sp["RFMT_ID"].nunique(),
            "n_units_engrafted": len(eng_units),
            "n_units_high_conf_engrafted": len(high_units),
            "n_units_persistent_to_d7": len(d7_units),
            "n_donor1_engrafted_units": int((eng_units["donor_label"] == "donor1").sum()) if not eng_units.empty else 0,
            "n_donor2_engrafted_units": int((eng_units["donor_label"] == "donor2").sum()) if not eng_units.empty else 0,
            "transmission_rate_among_available_units": len(eng_units) / n_av if n_av else np.nan,
            "best_conANI": best_con,
            "best_popANI": best_pop,
            "max_compared_bases": max_bases,
        })

    return pd.DataFrame(rows).sort_values(["n_units_engrafted", "transmission_rate_among_available_units"], ascending=False)


def add_antibiotics_to_summary(summary: pd.DataFrame, antibiotics: pd.DataFrame) -> pd.DataFrame:
    if summary.empty or antibiotics.empty:
        return summary
    if "RFMT_ID" not in antibiotics.columns:
        return summary
    return summary.merge(antibiotics, on="RFMT_ID", how="left")



# -----------------------------
# Debug summaries
# -----------------------------

def summarize_compare_metric_ranges(cmp: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if cmp.empty:
        return pd.DataFrame()
    for col in ["conANI", "popANI", "percent_compared", "compared_bases_count"]:
        if col in cmp.columns:
            s = pd.to_numeric(cmp[col], errors="coerce").dropna()
            if len(s) == 0:
                continue
            rows.append({
                "metric": col,
                "n_non_missing": len(s),
                "min": s.min(),
                "median": s.median(),
                "max": s.max(),
                "interpretation_hint": (
                    "fraction-scale likely; cutoffs like 99.999 or 50 are auto-converted"
                    if col in {"conANI", "popANI", "percent_compared"} and s.max() <= 1.5
                    else "percent-scale likely or count metric"
                ),
            })
    return pd.DataFrame(rows)


def summarize_call_debug(calls: pd.DataFrame) -> pd.DataFrame:
    if calls.empty:
        return pd.DataFrame([{"note": "mag_transfer_calls is empty"}])

    rows = []
    for donor, sub in calls.groupby("donor_label"):
        rows.append({
            "donor_label": donor,
            "n_call_rows": len(sub),
            "n_unique_MAGs": sub["genome"].nunique() if "genome" in sub.columns else np.nan,
            "n_rows_with_donor_post_comparison": sub["donor_post_compared_bases_count"].notna().sum() if "donor_post_compared_bases_count" in sub.columns else np.nan,
            "n_rows_engrafted": sub["is_engrafted"].astype(bool).sum() if "is_engrafted" in sub.columns else np.nan,
            "n_rows_high_conf_engrafted": sub["is_high_confidence_engrafted"].astype(bool).sum() if "is_high_confidence_engrafted" in sub.columns else np.nan,
            "max_donor_post_conANI": pd.to_numeric(sub.get("donor_post_conANI"), errors="coerce").max() if "donor_post_conANI" in sub.columns else np.nan,
            "max_donor_post_popANI": pd.to_numeric(sub.get("donor_post_popANI"), errors="coerce").max() if "donor_post_popANI" in sub.columns else np.nan,
            "max_donor_post_percent_compared": pd.to_numeric(sub.get("donor_post_percent_compared"), errors="coerce").max() if "donor_post_percent_compared" in sub.columns else np.nan,
            "max_donor_post_compared_bases": pd.to_numeric(sub.get("donor_post_compared_bases_count"), errors="coerce").max() if "donor_post_compared_bases_count" in sub.columns else np.nan,
        })
    return pd.DataFrame(rows)


def summarize_matching_debug(meta: pd.DataFrame, donor_available: pd.DataFrame, cmp: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if cmp.empty or donor_available.empty:
        return pd.DataFrame([{"note": "cmp or donor_available empty"}])

    for rfmt, av in donor_available.groupby("RFMT_ID"):
        cmp_r = cmp[cmp["RFMT_ID"] == rfmt]
        meta_r = meta[meta["RFMT_ID"] == rfmt]
        donor_genomes = set(av["genome"].astype(str))
        donor_alias = set()
        for g in donor_genomes:
            donor_alias.update(genome_aliases(g))

        cmp_genomes = set(cmp_r["genome"].astype(str))
        meta_samples = set(meta_r["sample_name"].astype(str))
        cmp_samples = set(cmp_r["sample1"].astype(str)) | set(cmp_r["sample2"].astype(str))

        rows.append({
            "RFMT_ID": rfmt,
            "n_metadata_samples": len(meta_samples),
            "n_compare_samples": len(cmp_samples),
            "n_sample_overlap": len(meta_samples & cmp_samples),
            "example_metadata_samples": ",".join(sorted(list(meta_samples))[:5]),
            "example_compare_samples": ",".join(sorted(list(cmp_samples))[:5]),
            "n_donor_available_genomes": len(donor_genomes),
            "n_compare_genomes": len(cmp_genomes),
            "n_genome_overlap_with_alias": len(donor_alias & cmp_genomes),
            "example_donor_genomes": ",".join(sorted(list(donor_genomes))[:3]),
            "example_compare_genomes": ",".join(sorted(list(cmp_genomes))[:3]),
        })
    return pd.DataFrame(rows)

# -----------------------------
# Main
# -----------------------------

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--metadata", required=True)
    p.add_argument("--instrain-root", required=True)
    p.add_argument("--mag-root", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--antibiotics", default=None, help="Optional RFMT_antibiotic_indicators_used.tsv")
    p.add_argument("--min-completeness", type=float, default=50)
    p.add_argument("--max-contamination", type=float, default=10)
    p.add_argument("--min-conani", type=float, default=99.99)
    p.add_argument("--min-popani", type=float, default=99.999)
    p.add_argument("--min-compared-bases", type=float, default=50000)
    p.add_argument("--min-percent-compared", type=float, default=50)
    args = p.parse_args()

    metadata = Path(args.metadata)
    instrain_root = Path(args.instrain_root)
    mag_root = Path(args.mag_root)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    meta = load_metadata(metadata)
    write_tsv(meta, outdir / "analysis_sample_metadata.tsv")

    role_counts = (
        meta.groupby(["RFMT_ID", "role"])
        .size()
        .reset_index(name="n_samples")
        .sort_values(["RFMT_ID", "role"])
    )
    write_tsv(role_counts, outdir / "metadata_role_counts.tsv")

    antibiotics = read_tsv(Path(args.antibiotics)) if args.antibiotics else pd.DataFrame()
    if not antibiotics.empty:
        write_tsv(antibiotics, outdir / "antibiotic_indicators_used.tsv")

    mag_annot, donor_available, sample_mag_counts, tpm_long, missing_bintable = parse_bintables(
        meta, mag_root, args.min_completeness, args.max_contamination
    )

    write_tsv(mag_annot, outdir / "bintable_mag_annotation.tsv")
    write_tsv(donor_available, outdir / "donor_available_MAGs.tsv")
    write_tsv(sample_mag_counts, outdir / "sample_MAG_counts.tsv")
    write_tsv(tpm_long, outdir / "MAG_TPM_long.tsv")
    write_tsv(missing_bintable, outdir / "missing_bintables.tsv")

    filtered_out = mag_annot[~mag_annot["quality_pass"].astype(bool)].copy() if not mag_annot.empty else pd.DataFrame()
    write_tsv(filtered_out, outdir / "filtered_out_low_quality_MAGs.tsv")

    cmp_raw, missing_cmp = load_compare_tables(meta, instrain_root)
    write_tsv(missing_cmp, outdir / "missing_compare_inputs.tsv")
    cmp = standardize_compare(cmp_raw) if not cmp_raw.empty else pd.DataFrame()
    write_tsv(cmp, outdir / "all_genomeWide_compare_merged.tsv")
    write_tsv(summarize_compare_metric_ranges(cmp), outdir / "compare_metric_ranges.tsv")
    write_tsv(summarize_matching_debug(meta, donor_available, cmp), outdir / "matching_debug_summary.tsv")

    calls = make_transfer_calls(
        meta=meta,
        donor_available=donor_available,
        mag_annot=mag_annot,
        cmp=cmp,
        min_conani=args.min_conani,
        min_popani=args.min_popani,
        min_compared_bases=args.min_compared_bases,
        min_percent_compared=args.min_percent_compared,
    )
    write_tsv(calls, outdir / "mag_transfer_calls.tsv")
    write_tsv(summarize_call_debug(calls), outdir / "call_debug_summary.tsv")

    eng = calls[calls["is_engrafted"].astype(bool)].copy() if not calls.empty else pd.DataFrame()
    write_tsv(eng, outdir / "engrafted_MAGs_only.tsv")

    rfmt_summary = summarize_rfmt(calls, donor_available)
    rfmt_summary = add_antibiotics_to_summary(rfmt_summary, antibiotics)
    write_tsv(rfmt_summary, outdir / "rfmt_transfer_summary.tsv")

    species_summary = summarize_species(calls, donor_available)
    write_tsv(species_summary, outdir / "species_transmission_summary.tsv")

    params = pd.DataFrame([{
        "min_completeness": args.min_completeness,
        "max_contamination": args.max_contamination,
        "min_conANI": args.min_conani,
        "min_popANI": args.min_popani,
        "min_compared_bases": args.min_compared_bases,
        "min_percent_compared": args.min_percent_compared,
        "denominator": "quality-filtered donor MAGs from donor bintable",
        "note": "v11 fixes role parsing: recipient_pre_donor2 is classified before donor2; strict baseline-absent engraftment logic; auto-detects ANI / percent_compared scales.",
    }])
    write_tsv(params, outdir / "analysis_parameters.tsv")

    print("\nDone.")
    print(f"Output directory: {outdir}")


if __name__ == "__main__":
    main()
