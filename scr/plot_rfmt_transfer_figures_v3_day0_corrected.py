#!/usr/bin/env python3
"""
Plot RFMT transfer results with the correct timeline:
  recipient_d0 = pre-Donor1 baseline
  recipient_pre_donor2 = post-Donor1 / immediately pre-Donor2
  recipient_d2/d5/d7 = post-Donor2 follow-up

Input directory must contain:
  mag_transfer_calls.tsv
  species_transmission_summary.tsv
  rfmt_transfer_summary.tsv
"""

from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        print(f"[MISS] {path}")
        return pd.DataFrame()
    print(f"[READ] {path}")
    return pd.read_csv(path, sep="\t", dtype=str, low_memory=False)


def to_numeric(df, cols):
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def to_bool_series(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.lower().isin(["true", "1", "yes", "y"])


def savefig(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.savefig(path.with_suffix(".pdf"))
    print(f"[WRITE] {path}")
    plt.close()


def short_species(x):
    if pd.isna(x) or str(x).strip() == "":
        return "Unknown"
    x = str(x).strip()
    if ";" in x:
        parts = [p.strip() for p in x.split(";")]
        sp = [p for p in parts if p.startswith("s__")]
        if sp:
            val = sp[-1].replace("s__", "").strip()
            return val or "Unknown"
    if x.startswith("s__"):
        return x.replace("s__", "", 1).strip() or "Unknown"
    return x


def role_to_timeline(role):
    """Return a label consistent with the corrected biology."""
    if pd.isna(role):
        return "unknown"
    s = str(role).lower()
    if "pre_donor2" in s or "predonor2" in s or "pre-donor2" in s:
        return "preDonor2"
    if "d7" in s or "day7" in s:
        return "day7"
    if "d5" in s or "day5" in s:
        return "day5"
    if "d2" in s or "day2" in s:
        return "day2"
    if "d0" in s or "day0" in s:
        return "day0_preDonor1"
    return str(role)


def classification_simple(x):
    if pd.isna(x):
        return "Uncertain"
    raw = str(x).strip()
    s = raw.lower()
    if s == "":
        return "Uncertain"
    if "donor2" in s and ("engraft" in s or "like" in s or "transmit" in s):
        return "Donor2-like"
    if "donor1" in s and ("engraft" in s or "like" in s or "transmit" in s):
        return "Donor1-like"
    if "replace" in s:
        return "Replacement"
    if "recipient" in s and ("persist" in s or "resident" in s or "like" in s or "previous" in s):
        return "Recipient-like"
    if "shared" in s:
        return "Already shared"
    if "not" in s or "no_donor" in s:
        return "Not transmitted"
    if "uncertain" in s or "low" in s or "ambiguous" in s:
        return "Uncertain"
    return raw


def ensure_required(df, cols, table_name):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        print(f"[SKIP] {table_name}: missing columns {missing}")
        print("Available columns:", list(df.columns))
        return False
    return True


def plot_rfmt_engraftment_counts(rfmt_summary, outdir):
    required = ["RFMT_ID", "donor_label", "n_MAGs_engrafted", "MAG_engraftment_rate"]
    if rfmt_summary.empty or not ensure_required(rfmt_summary, required, "rfmt_transfer_summary.tsv"):
        return
    df = to_numeric(rfmt_summary.copy(), ["n_MAGs_engrafted", "MAG_engraftment_rate"])
    count_df = df.pivot_table(index="RFMT_ID", columns="donor_label", values="n_MAGs_engrafted", aggfunc="sum", fill_value=0).sort_index()
    if count_df.empty:
        return
    plt.figure(figsize=(max(8, 0.45 * len(count_df)), 5))
    x = np.arange(len(count_df)); bottom = np.zeros(len(count_df))
    for donor in count_df.columns:
        vals = count_df[donor].astype(float).values
        plt.bar(x, vals, bottom=bottom, label=str(donor))
        bottom += vals
    plt.xticks(x, count_df.index, rotation=90)
    plt.xlabel("RFMT group")
    plt.ylabel("Number of engrafted MAG calls")
    plt.title("Donor-derived MAG engraftment across RFMT groups")
    plt.legend(title="Donor", frameon=False)
    savefig(outdir / "Fig1_RFMT_engrafted_MAG_counts_by_donor.png")

    rate_df = df.pivot_table(index="RFMT_ID", columns="donor_label", values="MAG_engraftment_rate", aggfunc="mean", fill_value=0).sort_index()
    plt.figure(figsize=(max(8, 0.45 * len(rate_df)), 5))
    x = np.arange(len(rate_df)); width = 0.8 / max(1, len(rate_df.columns))
    for i, donor in enumerate(rate_df.columns):
        plt.bar(x + i*width, rate_df[donor].astype(float).values, width=width, label=str(donor))
    plt.xticks(x + width*(len(rate_df.columns)-1)/2, rate_df.index, rotation=90)
    plt.xlabel("RFMT group")
    plt.ylabel("MAG engraftment rate")
    plt.title("MAG engraftment rate across RFMT groups")
    plt.legend(title="Donor", frameon=False)
    savefig(outdir / "Fig1b_RFMT_MAG_engraftment_rate_by_donor.png")


def plot_species_transmission_frequency(species_summary, outdir, top_n=25):
    required = ["species", "n_donor1_engrafted_units", "n_donor2_engrafted_units", "n_units_engrafted", "transmission_rate_among_tested_units"]
    if species_summary.empty or not ensure_required(species_summary, required, "species_transmission_summary.tsv"):
        return
    df = species_summary.copy()
    df["species_short"] = df["species"].map(short_species)
    df = to_numeric(df, ["n_donor1_engrafted_units", "n_donor2_engrafted_units", "n_units_engrafted", "transmission_rate_among_tested_units", "n_units_persistent_to_d7"])
    df = df[df["n_units_engrafted"].fillna(0) > 0].sort_values("n_units_engrafted", ascending=False).head(top_n).iloc[::-1]
    if df.empty:
        print("[SKIP] no engrafted species")
        return
    y = np.arange(len(df))
    d1 = df["n_donor1_engrafted_units"].fillna(0).astype(float).values
    d2 = df["n_donor2_engrafted_units"].fillna(0).astype(float).values
    plt.figure(figsize=(9, max(5, 0.32 * len(df))))
    plt.barh(y, d1, label="Donor1")
    plt.barh(y, d2, left=d1, label="Donor2")
    plt.yticks(y, df["species_short"])
    plt.xlabel("Number of RFMT donor units with successful transmission")
    plt.ylabel("Species")
    plt.title(f"Top {top_n} transmitted species")
    plt.legend(frameon=False)
    savefig(outdir / "Fig2_transmitted_species_frequency.png")

    df_rate = df.sort_values("transmission_rate_among_tested_units", ascending=True)
    plt.figure(figsize=(9, max(5, 0.32 * len(df_rate))))
    plt.barh(np.arange(len(df_rate)), df_rate["transmission_rate_among_tested_units"].fillna(0).astype(float).values)
    plt.yticks(np.arange(len(df_rate)), df_rate["species_short"])
    plt.xlabel("Transmission rate among tested donor units")
    plt.ylabel("Species")
    plt.title(f"Transmission rate of top {top_n} transmitted species")
    savefig(outdir / "Fig2b_species_transmission_rate.png")


def plot_corrected_timeline(mag_calls, outdir):
    required = ["donor_label", "post_role", "is_engrafted", "is_high_confidence_engrafted"]
    if mag_calls.empty or not ensure_required(mag_calls, required, "mag_transfer_calls.tsv"):
        return
    df = mag_calls.copy()
    df["time"] = df["post_role"].map(role_to_timeline)
    df["engrafted"] = to_bool_series(df["is_engrafted"])
    df["high_conf"] = to_bool_series(df["is_high_confidence_engrafted"])

    # Correct timeline: day0 is baseline, not plotted as post-transfer. Donor1 can be evaluated at preDonor2 and later persistence.
    order_by_donor = {
        "donor1": ["preDonor2", "day2", "day5", "day7"],
        "donor2": ["day2", "day5", "day7"],
    }
    for flag, title, fname in [
        ("engrafted", "Engrafted MAG calls over corrected timeline", "Fig3_corrected_timeline_engrafted_MAGs.png"),
        ("high_conf", "High-confidence engrafted MAG calls over corrected timeline", "Fig3b_corrected_timeline_high_conf_MAGs.png"),
    ]:
        sub = df[df[flag]].copy()
        if sub.empty:
            print(f"[SKIP] no rows for {flag}")
            continue
        plt.figure(figsize=(7, 4))
        for donor, order in order_by_donor.items():
            dsub = sub[sub["donor_label"] == donor]
            if dsub.empty:
                continue
            counts = dsub.groupby("time").size().reindex(order).fillna(0)
            plt.plot(order, counts.values, marker="o", label=donor)
        plt.xlabel("Recipient sample relative to donor exposure")
        plt.ylabel("Number of MAG calls")
        plt.title(title)
        plt.legend(title="Donor", frameon=False)
        savefig(outdir / fname)


def plot_classification_heatmap(mag_calls, outdir, top_n=60):
    required = ["RFMT_ID", "donor_label", "species", "classification"]
    if mag_calls.empty or not ensure_required(mag_calls, required, "mag_transfer_calls.tsv"):
        return
    df = mag_calls.copy()
    df["species_short"] = df["species"].map(short_species)
    df["classification_simple"] = df["classification"].map(classification_simple)
    df["col_label"] = df["RFMT_ID"].astype(str) + "_" + df["donor_label"].astype(str)
    priority = {"Donor2-like":7, "Donor1-like":6, "Replacement":5, "Recipient-like":4, "Already shared":3, "Uncertain":2, "Not transmitted":1, "Absent":0}
    event_classes = ["Donor1-like", "Donor2-like", "Replacement", "Recipient-like", "Already shared"]
    score = df[df["classification_simple"].isin(event_classes)].groupby("species_short").size().sort_values(ascending=False).head(top_n)
    if score.empty:
        print("[SKIP] no classified species events for heatmap")
        return
    df = df[df["species_short"].isin(score.index)].copy()
    df["priority"] = df["classification_simple"].map(priority).fillna(2)
    idx = df.sort_values("priority").groupby(["species_short", "col_label"]).tail(1)
    mat = idx.pivot(index="species_short", columns="col_label", values="classification_simple").fillna("Absent")
    num = mat.replace(priority).apply(pd.to_numeric, errors="coerce").fillna(priority["Uncertain"]).astype(float)
    plt.figure(figsize=(max(10, 0.35 * num.shape[1]), max(6, 0.25 * num.shape[0])))
    plt.imshow(num.to_numpy(), aspect="auto", interpolation="nearest")
    plt.xticks(np.arange(num.shape[1]), num.columns, rotation=90)
    plt.yticks(np.arange(num.shape[0]), num.index)
    plt.xlabel("RFMT donor unit")
    plt.ylabel("Species")
    plt.title("Species-level strain-source classification across RFMT donor units")
    cbar = plt.colorbar()
    cbar.set_ticks(list(priority.values()))
    cbar.set_ticklabels(list(priority.keys()))
    savefig(outdir / "Fig4_species_classification_heatmap.png")


def plot_donor_vs_baseline_similarity(mag_calls, outdir):
    required = ["donor_label", "donor_post_conANI", "baseline_post_conANI", "is_engrafted"]
    if mag_calls.empty or not ensure_required(mag_calls, required, "mag_transfer_calls.tsv"):
        return
    df = to_numeric(mag_calls.copy(), ["donor_post_conANI", "baseline_post_conANI"])
    df = df.dropna(subset=["donor_post_conANI", "baseline_post_conANI"])
    if df.empty:
        print("[SKIP] no rows with both donor_post_conANI and baseline_post_conANI")
        return
    df["engrafted"] = to_bool_series(df["is_engrafted"])
    plt.figure(figsize=(5.6, 5.6))
    for lab, sub in [("Not engrafted", df[~df["engrafted"]]), ("Engrafted", df[df["engrafted"]])]:
        if not sub.empty:
            plt.scatter(sub["baseline_post_conANI"], sub["donor_post_conANI"], alpha=0.55, label=lab)
    lo = min(df["baseline_post_conANI"].min(), df["donor_post_conANI"].min())
    hi = max(df["baseline_post_conANI"].max(), df["donor_post_conANI"].max())
    plt.plot([lo, hi], [lo, hi], linestyle="--")
    plt.xlabel("Baseline vs post conANI")
    plt.ylabel("Donor vs post conANI")
    plt.title("Donor-post similarity versus correct baseline-post similarity")
    plt.legend(frameon=False)
    savefig(outdir / "Fig5_donor_post_vs_baseline_post_conANI.png")


def plot_species_persistence(species_summary, outdir, top_n=25):
    required = ["species", "n_units_engrafted", "n_units_persistent_to_d7"]
    if species_summary.empty or not ensure_required(species_summary, required, "species_transmission_summary.tsv"):
        return
    df = species_summary.copy()
    df["species_short"] = df["species"].map(short_species)
    df = to_numeric(df, ["n_units_engrafted", "n_units_persistent_to_d7"])
    df = df[df["n_units_engrafted"].fillna(0) > 0]
    if df.empty:
        return
    df = df.sort_values(["n_units_persistent_to_d7", "n_units_engrafted"], ascending=False).head(top_n).iloc[::-1]
    y = np.arange(len(df))
    plt.figure(figsize=(9, max(5, 0.32 * len(df))))
    plt.barh(y, df["n_units_engrafted"].fillna(0).astype(float), label="Engrafted")
    plt.barh(y, df["n_units_persistent_to_d7"].fillna(0).astype(float), label="Persistent to day7")
    plt.yticks(y, df["species_short"])
    plt.xlabel("Number of RFMT donor units")
    plt.ylabel("Species")
    plt.title("Persistence of transmitted species to day7")
    plt.legend(frameon=False)
    savefig(outdir / "Fig6_species_persistence_to_day7.png")


def write_templates(outdir):
    txt = """Corrected interpretation note\n\nRecipient day0 is treated as the pre-Donor1 baseline and is not plotted as a post-transfer timepoint. Donor1 direct engraftment is best evaluated at recipient pre-Donor2 when that sample is available. Donor1-like calls at day2/day5/day7 indicate persistence after Donor2 exposure. Donor2 engraftment is evaluated at day2/day5/day7 using recipient pre-Donor2 as the preferred baseline; if pre-Donor2 is missing, day0 is only a weak baseline.\n\nFigure 1: MAG-level donor engraftment varied across RFMT groups. Counts are stratified by donor source.\nFigure 2: Successful transmission was concentrated in a subset of species, summarized after MAG-level strain evidence.\nFigure 3: Corrected timeline plot. Donor1 is evaluated at preDonor2 and later persistence timepoints; Donor2 is evaluated at day2/day5/day7.\nFigure 4: Species-level source classifications across RFMT donor units.\nFigure 5: Donor-post similarity compared against the correct baseline-post similarity.\nFigure 6: A subset of transmitted species persisted to day7.\n"""
    path = outdir / "figure_result_sentence_templates_corrected_timeline.txt"
    path.write_text(txt)
    print(f"[WRITE] {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", default="/group/sbms004/yxia/GUT/RFMT_tran_v2")
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--top-species", type=int, default=25)
    ap.add_argument("--top-heatmap-rows", type=int, default=60)
    args = ap.parse_args()
    indir = Path(args.indir)
    outdir = Path(args.outdir) if args.outdir else indir / "figures_corrected_timeline"
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"[INDIR] {indir}")
    print(f"[OUTDIR] {outdir}")
    mag = read_tsv(indir / "mag_transfer_calls.tsv")
    sp = read_tsv(indir / "species_transmission_summary.tsv")
    rfmt = read_tsv(indir / "rfmt_transfer_summary.tsv")
    plot_rfmt_engraftment_counts(rfmt, outdir)
    plot_species_transmission_frequency(sp, outdir, args.top_species)
    plot_corrected_timeline(mag, outdir)
    plot_classification_heatmap(mag, outdir, args.top_heatmap_rows)
    plot_donor_vs_baseline_similarity(mag, outdir)
    plot_species_persistence(sp, outdir, args.top_species)
    write_templates(outdir)
    print("Done. Figures are in:", outdir)


if __name__ == "__main__":
    main()
