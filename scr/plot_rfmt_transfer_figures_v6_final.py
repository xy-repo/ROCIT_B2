#!/usr/bin/env python3
"""
Plot RFMT transfer results with:
- antibiotic annotation strips reusable across RFMT-level plots
- sample MAG count heatmap
- corrected donor-available engraftment rate
- species plots with n_units_engrafted >= 2
- species-level and MAG-level classification heatmaps
- top50 MAG TPM heatmap
"""

import argparse
from pathlib import Path
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


# -------------------------
# Helpers
# -------------------------

def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        print(f"[MISS] {path}")
        return pd.DataFrame()
    print(f"[READ] {path}")
    return pd.read_csv(path, sep="\t", dtype=str)


def to_num(df, cols):
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def to_bool(s):
    return s.astype(str).str.strip().str.lower().isin(["true", "1", "yes", "y"])


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


def savefig(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.savefig(path.with_suffix(".pdf"))
    print(f"[WRITE] {path}")
    plt.close()


def role_label(role):
    mapping = {
        "donor1": "Donor1",
        "recipient_day0": "Day0\npre-D1",
        "recipient_pre_donor2": "Pre-D2",
        "donor2": "Donor2",
        "recipient_d2": "Day2",
        "recipient_d5": "Day5",
        "recipient_d7": "Day7",
    }
    return mapping.get(str(role), str(role))


def classification_simple(x):
    if pd.isna(x):
        return "Uncertain"
    s = str(x).strip().lower()
    if "donor2" in s and "engraft" in s:
        return "Donor2-like"
    if "donor1" in s and "engraft" in s:
        return "Donor1-like"
    if "late_detected" in s:
        return "Donor1-like"
    if "recipient" in s and "persistent" in s:
        return "Recipient-like"
    if "already_shared" in s or "shared" in s:
        return "Already shared"
    if "ambiguous" in s:
        return "Uncertain"
    if "not_engrafted" in s or "not" in s:
        return "Not transmitted"
    if "replace" in s:
        return "Replacement"
    if "uncertain" in s:
        return "Uncertain"
    return str(x)


# -------------------------
# Antibiotics
# -------------------------

DEFAULT_ABX_STAGES = [
    "FMT1_48h_prior",
    "FMT1_at_FMT",
    "FMT1_during_engraftment",
    "FMT1_D7",
    "FMT2_at_FMT",
    "FMT2_during_engraftment",
    "FMT2_D7",
]

ABX_LABELS = {
    "FMT1_48h_prior": "F1 48h",
    "FMT1_at_FMT": "F1 at FMT",
    "FMT1_during_engraftment": "F1 engraft",
    "FMT1_D7": "F1 D7",
    "FMT2_at_FMT": "F2 at FMT",
    "FMT2_during_engraftment": "F2 engraft",
    "FMT2_D7": "F2 D7",
}

ABX_COLORS = {
    "FMT1_48h_prior": "#8dd3c7",
    "FMT1_at_FMT": "#ffffb3",
    "FMT1_during_engraftment": "#bebada",
    "FMT1_D7": "#fb8072",
    "FMT2_at_FMT": "#80b1d3",
    "FMT2_during_engraftment": "#fdb462",
    "FMT2_D7": "#b3de69",
}


def load_antibiotics(indir: Path, explicit_path=None):
    path = Path(explicit_path) if explicit_path else indir / "antibiotic_indicators_used.tsv"
    abx = read_tsv(path)
    if abx.empty:
        return abx
    if "RFMT_ID" not in abx.columns:
        print("[WARN] antibiotic table lacks RFMT_ID")
        return pd.DataFrame()

    for st in DEFAULT_ABX_STAGES:
        if st in abx.columns:
            abx[st] = pd.to_numeric(abx[st], errors="coerce").fillna(0).astype(int)
    return abx


def draw_abx_strip(ax, rfmts, abx):
    """
    Draw antibiotic exposure strip. Colored square = yes for that stage.
    White = no/unknown.
    """
    if abx.empty:
        ax.axis("off")
        return

    stages = [s for s in DEFAULT_ABX_STAGES if s in abx.columns]
    if not stages:
        ax.axis("off")
        return

    abx_idx = abx.set_index("RFMT_ID")
    mat = np.zeros((len(stages), len(rfmts)))
    for j, r in enumerate(rfmts):
        if r in abx_idx.index:
            for i, st in enumerate(stages):
                mat[i, j] = int(abx_idx.loc[r, st])

    ax.set_xlim(-0.5, len(rfmts) - 0.5)
    ax.set_ylim(len(stages) - 0.5, -0.5)

    for i, st in enumerate(stages):
        for j, r in enumerate(rfmts):
            color = ABX_COLORS[st] if mat[i, j] == 1 else "white"
            ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, facecolor=color, edgecolor="0.75", linewidth=0.4))

    ax.set_yticks(range(len(stages)))
    ax.set_yticklabels([ABX_LABELS.get(s, s) for s in stages], fontsize=7)
    ax.set_xticks(range(len(rfmts)))
    ax.set_xticklabels([])
    ax.tick_params(length=0)
    ax.set_ylabel("Antibiotics", fontsize=8)
    for spine in ax.spines.values():
        spine.set_visible(False)


def abx_legend_handles():
    handles = [Patch(facecolor=ABX_COLORS[s], edgecolor="0.5", label=ABX_LABELS.get(s, s)) for s in DEFAULT_ABX_STAGES]
    handles.append(Patch(facecolor="white", edgecolor="0.5", label="No/unknown"))
    return handles


# -------------------------
# Plots
# -------------------------

def plot_rfmt_counts_and_rates(rfmt_summary, abx, outdir):
    if rfmt_summary.empty:
        return
    required = {"RFMT_ID", "donor_label", "n_MAGs_engrafted", "MAG_engraftment_rate"}
    if not required.issubset(rfmt_summary.columns):
        print("[SKIP] rfmt summary missing columns", required - set(rfmt_summary.columns))
        return

    df = rfmt_summary.copy()
    df = to_num(df, ["n_MAGs_engrafted", "MAG_engraftment_rate", "n_MAGs_available"])
    rfmts = sorted(df["RFMT_ID"].dropna().unique())

    # Counts
    count_df = df.pivot_table(index="RFMT_ID", columns="donor_label", values="n_MAGs_engrafted", aggfunc="sum", fill_value=0).reindex(rfmts)
    fig = plt.figure(figsize=(max(9, 0.5 * len(rfmts)), 6.2))
    gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[4.5, 1.3], hspace=0.05)
    ax = fig.add_subplot(gs[0, 0])
    ax_abx = fig.add_subplot(gs[1, 0], sharex=ax)

    x = np.arange(len(rfmts))
    bottom = np.zeros(len(rfmts))
    for donor in count_df.columns:
        vals = count_df[donor].astype(float).values
        ax.bar(x, vals, bottom=bottom, label=str(donor))
        bottom += vals

    ax.set_xticks(x)
    ax.set_xticklabels(rfmts, rotation=90)
    ax.set_ylabel("Engrafted MAGs")
    ax.set_title("Donor-derived MAG engraftment counts")
    ax.legend(frameon=False, title="Donor")
    draw_abx_strip(ax_abx, rfmts, abx)
    fig.legend(handles=abx_legend_handles(), loc="upper right", bbox_to_anchor=(1.02, 0.98), fontsize=7, frameon=False)
    savefig(outdir / "Fig1_RFMT_engrafted_MAG_counts_with_antibiotics.png")

    # Rates
    rate_df = df.pivot_table(index="RFMT_ID", columns="donor_label", values="MAG_engraftment_rate", aggfunc="mean", fill_value=0).reindex(rfmts)
    fig = plt.figure(figsize=(max(9, 0.5 * len(rfmts)), 6.2))
    gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[4.5, 1.3], hspace=0.05)
    ax = fig.add_subplot(gs[0, 0])
    ax_abx = fig.add_subplot(gs[1, 0], sharex=ax)

    width = 0.8 / max(1, len(rate_df.columns))
    for i, donor in enumerate(rate_df.columns):
        ax.bar(x + i * width, rate_df[donor].astype(float).values, width=width, label=str(donor))
    offset = width * (len(rate_df.columns) - 1) / 2
    ax.set_xticks(x + offset)
    ax.set_xticklabels(rfmts, rotation=90)
    ax.set_ylabel("Engraftment rate\n(engrafted / donor-available)")
    ax.set_title("Corrected MAG engraftment rate")
    ax.legend(frameon=False, title="Donor")
    draw_abx_strip(ax_abx, rfmts, abx)
    fig.legend(handles=abx_legend_handles(), loc="upper right", bbox_to_anchor=(1.02, 0.98), fontsize=7, frameon=False)
    savefig(outdir / "Fig1b_corrected_MAG_engraftment_rate_with_antibiotics.png")


def plot_sample_mag_count_heatmap(sample_counts, abx, outdir):
    if sample_counts.empty:
        return
    required = {"RFMT_ID", "role", "n_MAGs_quality"}
    if not required.issubset(sample_counts.columns):
        print("[SKIP] sample_MAG_counts missing columns", required - set(sample_counts.columns))
        return

    df = sample_counts.copy()
    df = to_num(df, ["n_MAGs_quality", "n_MAGs_total"])
    roles = ["donor1", "recipient_day0", "recipient_pre_donor2", "donor2", "recipient_d2", "recipient_d5", "recipient_d7"]
    rfmts = sorted(df["RFMT_ID"].dropna().unique())

    mat = pd.DataFrame(np.nan, index=roles, columns=rfmts)
    for _, row in df.iterrows():
        if row["role"] in mat.index and row["RFMT_ID"] in mat.columns:
            mat.loc[row["role"], row["RFMT_ID"]] = row["n_MAGs_quality"]

    fig = plt.figure(figsize=(max(9, 0.5 * len(rfmts)), 6.5))
    gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[4.8, 1.2], hspace=0.05)
    ax = fig.add_subplot(gs[0, 0])
    ax_abx = fig.add_subplot(gs[1, 0], sharex=ax)

    im = ax.imshow(mat.astype(float).values, aspect="auto", interpolation="nearest")
    ax.set_xticks(np.arange(len(rfmts)))
    ax.set_xticklabels(rfmts, rotation=90)
    ax.set_yticks(np.arange(len(roles)))
    ax.set_yticklabels([role_label(r) for r in roles])
    ax.set_title("Number of quality-filtered MAGs per sample")
    ax.set_ylabel("Sample role / timepoint")
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("MAG count\n(completeness ≥50%, contamination ≤10%)")

    draw_abx_strip(ax_abx, rfmts, abx)
    fig.legend(handles=abx_legend_handles(), loc="upper right", bbox_to_anchor=(1.02, 0.98), fontsize=7, frameon=False)
    savefig(outdir / "Fig0_sample_MAG_count_heatmap_with_antibiotics.png")


def plot_species_frequency(species_summary, outdir, min_events=2, top_n=50):
    if species_summary.empty:
        return
    required = {"species", "n_units_engrafted", "n_donor1_engrafted_units", "n_donor2_engrafted_units", "n_RFMT_donor_units_available", "transmission_rate_among_available_units"}
    if not required.issubset(species_summary.columns):
        print("[SKIP] species summary missing columns", required - set(species_summary.columns))
        return

    df = species_summary.copy()
    df = to_num(df, ["n_units_engrafted", "n_donor1_engrafted_units", "n_donor2_engrafted_units", "n_RFMT_donor_units_available", "transmission_rate_among_available_units"])
    df = df[df["n_units_engrafted"] >= min_events].copy()
    if df.empty:
        print("[SKIP] no species with n_units_engrafted >=", min_events)
        return

    df["species_label"] = df["species"].map(short_species)
    df["frac_label"] = df["species_label"] + " (" + df["n_units_engrafted"].astype(int).astype(str) + "/" + df["n_RFMT_donor_units_available"].astype(int).astype(str) + ")"
    df = df.sort_values("n_units_engrafted", ascending=False).head(top_n).iloc[::-1]

    y = np.arange(len(df))
    d1 = df["n_donor1_engrafted_units"].fillna(0).astype(float).values
    d2 = df["n_donor2_engrafted_units"].fillna(0).astype(float).values

    plt.figure(figsize=(9, max(5, 0.32 * len(df))))
    plt.barh(y, d1, label="Donor1")
    plt.barh(y, d2, left=d1, label="Donor2")
    plt.yticks(y, df["frac_label"])
    plt.xlabel("Number of RFMT-donor combinations with successful transmission")
    plt.ylabel("Species (engrafted / donor-available)")
    plt.title(f"Transmitted species with ≥{min_events} successful RFMT-donor combinations")
    plt.legend(frameon=False)
    savefig(outdir / "Fig2_transmitted_species_frequency_min2.png")

    df_rate = df.sort_values("transmission_rate_among_available_units", ascending=True)
    plt.figure(figsize=(9, max(5, 0.32 * len(df_rate))))
    plt.barh(np.arange(len(df_rate)), df_rate["transmission_rate_among_available_units"].astype(float).values)
    plt.yticks(np.arange(len(df_rate)), df_rate["frac_label"])
    plt.xlabel("Transmission rate among donor-available RFMT-donor combinations")
    plt.ylabel("Species (engrafted / donor-available)")
    plt.title(f"Transmission rate for species with ≥{min_events} successful transmissions")
    savefig(outdir / "Fig2b_species_transmission_rate_min2_corrected_denominator.png")


def plot_corrected_timeline(mag_calls, outdir):
    if mag_calls.empty:
        return
    required = {"donor_label", "post_role", "is_engrafted", "is_high_confidence_engrafted"}
    if not required.issubset(mag_calls.columns):
        print("[SKIP] mag calls missing columns", required - set(mag_calls.columns))
        return

    df = mag_calls.copy()
    df["is_engrafted_bool"] = to_bool(df["is_engrafted"])
    df["is_high_bool"] = to_bool(df["is_high_confidence_engrafted"])

    order_roles = ["recipient_pre_donor2", "recipient_d2", "recipient_d5", "recipient_d7"]
    role_labels = [role_label(r) for r in order_roles]

    for flag, title, fname in [
        ("is_engrafted_bool", "Engrafted MAG calls across corrected timeline", "Fig3_corrected_timeline_engrafted_MAGs.png"),
        ("is_high_bool", "High-confidence engrafted MAG calls across corrected timeline", "Fig3b_corrected_timeline_high_conf_MAGs.png"),
    ]:
        sub = df[df[flag]].copy()
        if sub.empty:
            continue
        count_df = sub.groupby(["donor_label", "post_role"]).size().reset_index(name="n")
        plt.figure(figsize=(6.5, 4))
        for donor, ddf in count_df.groupby("donor_label"):
            vals = ddf.set_index("post_role").reindex(order_roles)["n"].fillna(0).values
            plt.plot(role_labels, vals, marker="o", label=donor)
        plt.xlabel("Recipient timepoint")
        plt.ylabel("MAG calls")
        plt.title(title)
        plt.legend(frameon=False)
        savefig(outdir / fname)


def plot_classification_heatmap(mag_calls, outdir, level="species", top_n=60):
    if mag_calls.empty:
        return
    required = {"RFMT_ID", "donor_label", "classification"}
    if not required.issubset(mag_calls.columns):
        print("[SKIP] mag calls missing columns", required - set(mag_calls.columns))
        return

    df = mag_calls.copy()
    df["classification_simple"] = df["classification"].map(classification_simple)
    df["col_label"] = df["RFMT_ID"].astype(str) + "_" + df["donor_label"].astype(str)

    if level == "species":
        if "species" not in df.columns:
            return
        df["row_label"] = df["species"].map(short_species)
        fname = "Fig4_species_classification_heatmap.png"
        title = "Species-level classification heatmap"
    else:
        if "genome" not in df.columns:
            return
        # MAG label includes species short + compact genome ending
        df["row_label"] = df["genome"].astype(str)
        if "species" in df.columns:
            df["row_label"] = df["species"].map(short_species) + " | " + df["genome"].astype(str)
        fname = "Fig4b_MAG_classification_heatmap.png"
        title = "MAG-level classification heatmap"

    priority = {
        "Donor2-like": 7,
        "Donor1-like": 6,
        "Replacement": 5,
        "Recipient-like": 4,
        "Already shared": 3,
        "Uncertain": 2,
        "Not transmitted": 1,
        "Absent": 0,
    }

    event_classes = ["Donor1-like", "Donor2-like", "Replacement", "Recipient-like", "Already shared"]
    score = df[df["classification_simple"].isin(event_classes)].groupby("row_label").size().sort_values(ascending=False).head(top_n)
    if score.empty:
        return

    df = df[df["row_label"].isin(score.index)].copy()
    df["priority"] = df["classification_simple"].map(priority).fillna(2)
    idx = df.sort_values("priority").groupby(["row_label", "col_label"]).tail(1)

    mat = idx.pivot(index="row_label", columns="col_label", values="classification_simple").fillna("Absent")

    class_to_num = {
        "Absent": 0,
        "Not transmitted": 1,
        "Uncertain": 2,
        "Already shared": 3,
        "Recipient-like": 4,
        "Replacement": 5,
        "Donor1-like": 6,
        "Donor2-like": 7,
    }

    num = mat.replace(class_to_num)
    num = num.apply(pd.to_numeric, errors="coerce").fillna(2).astype(float)

    plt.figure(figsize=(max(10, 0.35 * num.shape[1]), max(6, 0.25 * num.shape[0])))
    plt.imshow(num.to_numpy(), aspect="auto", interpolation="nearest")
    plt.xticks(np.arange(num.shape[1]), num.columns, rotation=90)
    plt.yticks(np.arange(num.shape[0]), num.index)
    plt.title(title)
    plt.xlabel("RFMT-donor combination")
    plt.ylabel(level.capitalize())
    cbar = plt.colorbar()
    cbar.set_ticks(list(class_to_num.values()))
    cbar.set_ticklabels(list(class_to_num.keys()))
    savefig(outdir / fname)


def plot_donor_vs_baseline_conani(mag_calls, outdir):
    if mag_calls.empty:
        return
    required = {"donor_post_conANI", "baseline_post_conANI", "is_engrafted"}
    if not required.issubset(mag_calls.columns):
        return

    df = mag_calls.copy()
    df = to_num(df, ["donor_post_conANI", "baseline_post_conANI"])
    df = df.dropna(subset=["donor_post_conANI", "baseline_post_conANI"])
    if df.empty:
        return

    eng = to_bool(df["is_engrafted"])

    plt.figure(figsize=(5.5, 5.5))
    plt.scatter(df.loc[~eng, "baseline_post_conANI"], df.loc[~eng, "donor_post_conANI"], alpha=0.25, label="Not engrafted")
    plt.scatter(df.loc[eng, "baseline_post_conANI"], df.loc[eng, "donor_post_conANI"], alpha=0.75, label="Engrafted")
    mn = min(df["baseline_post_conANI"].min(), df["donor_post_conANI"].min())
    mx = max(df["baseline_post_conANI"].max(), df["donor_post_conANI"].max())
    plt.plot([mn, mx], [mn, mx], linestyle="--")
    plt.xlabel("Baseline vs post conANI")
    plt.ylabel("Donor vs post conANI")
    plt.title("Donor-post versus baseline-post strain similarity")
    plt.legend(frameon=False)
    savefig(outdir / "Fig5_donor_post_vs_baseline_post_conANI.png")


def plot_species_persistence(species_summary, outdir, min_events=2, top_n=50):
    if species_summary.empty:
        return
    required = {"species", "n_units_engrafted", "n_units_persistent_to_d7"}
    if not required.issubset(species_summary.columns):
        return
    df = species_summary.copy()
    df = to_num(df, ["n_units_engrafted", "n_units_persistent_to_d7"])
    df = df[df["n_units_engrafted"] >= min_events].copy()
    if df.empty:
        return
    df["species_label"] = df["species"].map(short_species)
    df = df.sort_values(["n_units_persistent_to_d7", "n_units_engrafted"], ascending=False).head(top_n).iloc[::-1]

    y = np.arange(len(df))
    plt.figure(figsize=(9, max(5, 0.32 * len(df))))
    plt.barh(y, df["n_units_engrafted"].astype(float).values, label="Engrafted")
    plt.barh(y, df["n_units_persistent_to_d7"].fillna(0).astype(float).values, label="Persistent to D7")
    plt.yticks(y, df["species_label"])
    plt.xlabel("Number of RFMT-donor combinations")
    plt.ylabel("Species")
    plt.title(f"Persistence to D7 among species with ≥{min_events} transmissions")
    plt.legend(frameon=False)
    savefig(outdir / "Fig6_species_persistence_to_day7_min2.png")


def plot_top50_mag_tpm(tpm_long, sample_meta, outdir, top_n=50):
    if tpm_long.empty:
        print("[SKIP] TPM long table empty")
        return
    required = {"sample_name", "genome", "TPM"}
    if not required.issubset(tpm_long.columns):
        print("[SKIP] TPM table missing", required - set(tpm_long.columns))
        return

    df = tpm_long.copy()
    df = to_num(df, ["TPM"])
    if "quality_pass" in df.columns:
        df = df[df["quality_pass"].astype(str).str.lower().isin(["true", "1", "yes"])]
    df = df.dropna(subset=["TPM"])
    if df.empty:
        return

    # Determine top MAGs by max TPM across samples.
    top = df.groupby("genome")["TPM"].max().sort_values(ascending=False).head(top_n).index
    df = df[df["genome"].isin(top)].copy()

    sample_order = sample_meta.sort_values(["RFMT_ID", "role_order", "sample_name"])["sample_name"].drop_duplicates().tolist()
    sample_order = [s for s in sample_order if s in set(df["sample_name"])]

    if not sample_order:
        sample_order = sorted(df["sample_name"].unique())

    mat = df.pivot_table(index="genome", columns="sample_name", values="TPM", aggfunc="max", fill_value=0)
    mat = mat.reindex(index=top, columns=sample_order, fill_value=0)

    vals = np.log10(mat.astype(float).values + 1)

    plt.figure(figsize=(max(10, 0.28 * len(sample_order)), max(8, 0.22 * len(mat))))
    im = plt.imshow(vals, aspect="auto", interpolation="nearest")
    plt.xticks(np.arange(len(sample_order)), sample_order, rotation=90, fontsize=6)
    plt.yticks(np.arange(len(mat.index)), mat.index, fontsize=6)
    plt.xlabel("Sample")
    plt.ylabel("MAG")
    plt.title(f"Top {top_n} MAG TPM distribution across samples")
    cbar = plt.colorbar(im)
    cbar.set_label("log10(TPM + 1)")
    savefig(outdir / "Fig7_top50_MAG_TPM_heatmap.png")


def write_notes(outdir):
    txt = """Figure notes

Antibiotic annotation:
Colored square indicates antibiotic exposure recorded for the RFMT at the indicated stage.
White indicates no/unknown. The same annotation scheme is used for RFMT-level plots.

Corrected MAG engraftment rate:
MAG_engraftment_rate = n_MAGs_engrafted / n_MAGs_available.
n_MAGs_available is defined from quality-filtered donor MAGs in donor bintables
(completeness >= 50%, contamination <= 10% by default).

Figure 4 species-level heatmap:
Each cell is one species x RFMT-donor combination.
If multiple MAGs from one species occur in the same cell, the strongest MAG-level classification is shown.

Figure 4b MAG-level heatmap:
Each row is a specific MAG, not collapsed by species.

Figure 2:
Only species with n_units_engrafted >= 2 are shown.
Labels include engrafted / donor-available denominators.
"""
    path = outdir / "figure_notes_v6.txt"
    path.write_text(txt)
    print(f"[WRITE] {path}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--indir", required=True)
    p.add_argument("--antibiotics", default=None, help="Optional antibiotic indicator TSV. Default: <indir>/antibiotic_indicators_used.tsv")
    p.add_argument("--outdir", default=None)
    p.add_argument("--top-species", type=int, default=50)
    p.add_argument("--top-heatmap-rows", type=int, default=60)
    p.add_argument("--top-mag-tpm", type=int, default=50)
    p.add_argument("--min-species-events", type=int, default=2)
    args = p.parse_args()

    indir = Path(args.indir)
    outdir = Path(args.outdir) if args.outdir else indir / "figures_v6_final"
    outdir.mkdir(parents=True, exist_ok=True)

    mag_calls = read_tsv(indir / "mag_transfer_calls.tsv")
    species_summary = read_tsv(indir / "species_transmission_summary.tsv")
    rfmt_summary = read_tsv(indir / "rfmt_transfer_summary.tsv")
    sample_counts = read_tsv(indir / "sample_MAG_counts.tsv")
    tpm_long = read_tsv(indir / "MAG_TPM_long.tsv")
    sample_meta = read_tsv(indir / "analysis_sample_metadata.tsv")
    abx = load_antibiotics(indir, args.antibiotics)

    plot_sample_mag_count_heatmap(sample_counts, abx, outdir)
    plot_rfmt_counts_and_rates(rfmt_summary, abx, outdir)
    plot_species_frequency(species_summary, outdir, min_events=args.min_species_events, top_n=args.top_species)
    plot_corrected_timeline(mag_calls, outdir)
    plot_classification_heatmap(mag_calls, outdir, level="species", top_n=args.top_heatmap_rows)
    plot_classification_heatmap(mag_calls, outdir, level="MAG", top_n=args.top_heatmap_rows)
    plot_donor_vs_baseline_conani(mag_calls, outdir)
    plot_species_persistence(species_summary, outdir, min_events=args.min_species_events, top_n=args.top_species)
    plot_top50_mag_tpm(tpm_long, sample_meta, outdir, top_n=args.top_mag_tpm)
    write_notes(outdir)

    print("\nDone. Figures written to:")
    print(outdir)


if __name__ == "__main__":
    main()
