#!/usr/bin/env python3
"""
Audit donor assignment for same-donor RFMT analysis.

Purpose:
  Detect cases where a post-FMT2 recipient MAG was classified as FMT1-donor-like,
  but the same post sample also matches the FMT2 donor sample. This is especially
  important because donor1 and donor2 are samples from the same donor individual.

Use case:
  RFMT031 has many FMT1-like engrafted MAGs, but Fig0 suggests FMT2 donor and day2
  are very similar. This script checks whether those FMT1-like calls are actually
  better explained by the FMT2 donor sample or should be marked as same-donor ambiguous.

Required files from analysis output:
  mag_transfer_calls.tsv
  all_genomeWide_compare_merged.tsv
  analysis_sample_metadata.tsv

Example:
  python audit_rfmt_same_donor_assignment.py \
    --indir /group/sbms004/yxia/GUT/RFMT_tran_v12_same_donor \
    --rfmt RFMT031 \
    --outdir /group/sbms004/yxia/GUT/RFMT_tran_v12_same_donor/audit_RFMT031
"""

import argparse
import re
from pathlib import Path
import pandas as pd
import numpy as np


def read_tsv(path):
    if not Path(path).exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t", dtype=str)


def to_num(x):
    return pd.to_numeric(x, errors="coerce")


def extract_sample_id(x):
    s = str(x)
    m = re.search(r"(RFMT\d{3}_\d+)", s)
    if m:
        return m.group(1)
    return Path(s).name.replace("_IS", "")


def genome_aliases(genome):
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


def scale_cutoff(value, cutoff):
    v = pd.to_numeric(value, errors="coerce")
    c = pd.to_numeric(cutoff, errors="coerce")
    if pd.isna(v) or pd.isna(c):
        return c
    if v <= 1.5 and c > 1.5:
        return c / 100.0
    if v > 1.5 and c <= 1.5:
        return c * 100.0
    return c


def get_pair(cmp_df, rfmt, genome, s1, s2):
    aliases = genome_aliases(genome)
    sub = cmp_df[(cmp_df["RFMT_ID"] == rfmt) & (cmp_df["genome"].astype(str).isin(aliases))].copy()
    if sub.empty:
        return {}

    s1 = extract_sample_id(s1)
    s2 = extract_sample_id(s2)
    hit = sub[
        ((sub["sample1"] == s1) & (sub["sample2"] == s2)) |
        ((sub["sample1"] == s2) & (sub["sample2"] == s1))
    ].copy()

    if hit.empty:
        return {}

    if "compared_bases_count" in hit.columns:
        hit["compared_bases_count_num"] = to_num(hit["compared_bases_count"])
        hit = hit.sort_values("compared_bases_count_num", ascending=False)

    return hit.iloc[0].to_dict()


def passes(row, min_conani, min_popani, min_percent_compared, min_compared_bases):
    if not row:
        return False
    con = to_num(row.get("conANI"))
    pop = to_num(row.get("popANI"))
    pct = to_num(row.get("percent_compared"))
    bases = to_num(row.get("compared_bases_count"))

    if pd.isna(con) or pd.isna(pop) or pd.isna(pct) or pd.isna(bases):
        return False

    con_cut = scale_cutoff(con, min_conani)
    pop_cut = scale_cutoff(pop, min_popani)
    pct_cut = scale_cutoff(pct, min_percent_compared)

    return con >= con_cut and pop >= pop_cut and pct >= pct_cut and bases >= min_compared_bases


def metric(row, col):
    return to_num(row.get(col)) if row else np.nan


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--indir", required=True)
    p.add_argument("--rfmt", required=True)
    p.add_argument("--outdir", default=None)
    p.add_argument("--min-conani", type=float, default=99.99)
    p.add_argument("--min-popani", type=float, default=99.999)
    p.add_argument("--min-percent-compared", type=float, default=50)
    p.add_argument("--min-compared-bases", type=float, default=50000)
    args = p.parse_args()

    indir = Path(args.indir)
    outdir = Path(args.outdir) if args.outdir else indir / f"audit_{args.rfmt}"
    outdir.mkdir(parents=True, exist_ok=True)

    calls = read_tsv(indir / "mag_transfer_calls.tsv")
    cmp_df = read_tsv(indir / "all_genomeWide_compare_merged.tsv")
    meta = read_tsv(indir / "analysis_sample_metadata.tsv")

    # Ensure normalized compare sample names.
    cmp_df["sample1"] = cmp_df["sample1"].map(extract_sample_id)
    cmp_df["sample2"] = cmp_df["sample2"].map(extract_sample_id)

    rfmt = args.rfmt
    meta_r = meta[meta["RFMT_ID"] == rfmt].copy()

    donor1 = meta_r.loc[meta_r["role"] == "donor1", "sample_name"].dropna().unique().tolist()
    donor2 = meta_r.loc[meta_r["role"] == "donor2", "sample_name"].dropna().unique().tolist()
    day0 = meta_r.loc[meta_r["role"] == "recipient_day0", "sample_name"].dropna().unique().tolist()
    pred2 = meta_r.loc[meta_r["role"] == "recipient_pre_donor2", "sample_name"].dropna().unique().tolist()

    if not donor1:
        raise ValueError(f"No donor1 sample found for {rfmt}")
    if not donor2:
        raise ValueError(f"No donor2 sample found for {rfmt}")

    donor1 = donor1[0]
    donor2 = donor2[0]
    day0 = day0[0] if day0 else ""
    pred2 = pred2[0] if pred2 else ""

    calls_r = calls[calls["RFMT_ID"] == rfmt].copy()

    post_roles = {"recipient_d2", "recipient_d5", "recipient_d7"}

    # Focus on post-FMT2 timepoints classified as donor1-like/engrafted.
    focus = calls_r[
        (calls_r["donor_label"] == "donor1") &
        (calls_r["post_role"].isin(post_roles)) &
        (calls_r["is_engrafted"].astype(str).str.lower().isin(["true", "1", "yes"]))
    ].copy()

    rows = []

    for _, r in focus.iterrows():
        genome = r["genome"]
        post_sample = r["post_sample"]
        post_role = r["post_role"]

        d1_post = get_pair(cmp_df, rfmt, genome, donor1, post_sample)
        d2_post = get_pair(cmp_df, rfmt, genome, donor2, post_sample)
        pred2_post = get_pair(cmp_df, rfmt, genome, pred2, post_sample) if pred2 else {}
        d1_d2 = get_pair(cmp_df, rfmt, genome, donor1, donor2)
        d2_pred2 = get_pair(cmp_df, rfmt, genome, donor2, pred2) if pred2 else {}
        d1_pred2 = get_pair(cmp_df, rfmt, genome, donor1, pred2) if pred2 else {}

        d1_pass = passes(d1_post, args.min_conani, args.min_popani, args.min_percent_compared, args.min_compared_bases)
        d2_pass = passes(d2_post, args.min_conani, args.min_popani, args.min_percent_compared, args.min_compared_bases)
        pred2_pass = passes(pred2_post, args.min_conani, args.min_popani, args.min_percent_compared, args.min_compared_bases)
        d1d2_pass = passes(d1_d2, args.min_conani, args.min_popani, args.min_percent_compared, args.min_compared_bases)

        d1_pop = metric(d1_post, "popANI")
        d2_pop = metric(d2_post, "popANI")
        pred2_pop = metric(pred2_post, "popANI")

        if d2_pass and not pred2_pass:
            if d1_pass and d1d2_pass:
                suggested = "same_donor_strain_detected_after_FMT2"
                reason = "post matches both FMT1 and FMT2 donor samples; FMT1/FMT2 donor samples also match"
            elif d2_pop >= d1_pop:
                suggested = "FMT2_donor_sample_associated"
                reason = "post matches FMT2 donor sample and FMT2-post popANI >= FMT1-post popANI"
            else:
                suggested = "FMT1_like_but_FMT2_also_matches"
                reason = "post matches FMT2 donor sample, but FMT1-post popANI is higher"
        elif d1_pass and not d2_pass and not pred2_pass:
            suggested = "FMT1_donor_sample_associated_late_detected"
            reason = "post matches FMT1 donor sample but not FMT2 donor sample or pre-FMT2 baseline"
        elif pred2_pass:
            suggested = "baseline_like_not_new_FMT2"
            reason = "post matches pre-FMT2 baseline"
        else:
            suggested = "uncertain_or_low_support"
            reason = "does not clearly pass FMT2/post or baseline logic"

        rows.append({
            "RFMT_ID": rfmt,
            "genome": genome,
            "species": r.get("species", ""),
            "post_sample": post_sample,
            "post_role": post_role,
            "original_classification": r.get("classification", ""),
            "original_is_engrafted": r.get("is_engrafted", ""),
            "D1_post_pass": d1_pass,
            "D2_post_pass": d2_pass,
            "preD2_post_pass": pred2_pass,
            "D1_D2_pass": d1d2_pass,
            "D1_post_conANI": metric(d1_post, "conANI"),
            "D1_post_popANI": d1_pop,
            "D1_post_percent_compared": metric(d1_post, "percent_compared"),
            "D1_post_compared_bases": metric(d1_post, "compared_bases_count"),
            "D2_post_conANI": metric(d2_post, "conANI"),
            "D2_post_popANI": d2_pop,
            "D2_post_percent_compared": metric(d2_post, "percent_compared"),
            "D2_post_compared_bases": metric(d2_post, "compared_bases_count"),
            "preD2_post_conANI": metric(pred2_post, "conANI"),
            "preD2_post_popANI": pred2_pop,
            "preD2_post_percent_compared": metric(pred2_post, "percent_compared"),
            "preD2_post_compared_bases": metric(pred2_post, "compared_bases_count"),
            "D1_D2_popANI": metric(d1_d2, "popANI"),
            "D2_preD2_popANI": metric(d2_pred2, "popANI"),
            "D1_preD2_popANI": metric(d1_pred2, "popANI"),
            "suggested_assignment": suggested,
            "audit_reason": reason,
        })

    audit = pd.DataFrame(rows)

    out = outdir / f"{rfmt}_postFMT2_donor_assignment_audit.csv"
    audit.to_csv(out, index=False)
    print(f"[WRITE] {out} rows={len(audit)}")

    if not audit.empty:
        summary = (
            audit.groupby("suggested_assignment")
            .size()
            .reset_index(name="n_MAG_calls")
            .sort_values("n_MAG_calls", ascending=False)
        )
    else:
        summary = pd.DataFrame(columns=["suggested_assignment", "n_MAG_calls"])

    out2 = outdir / f"{rfmt}_postFMT2_donor_assignment_audit_summary.csv"
    summary.to_csv(out2, index=False)
    print(f"[WRITE] {out2}")

    print("\nSummary:")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
