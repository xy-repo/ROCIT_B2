#!/usr/bin/env python3
"""
RFMT transfer plotting script v7

Fixes compared with previous version:
1. Antibiotic indicators are read from *_indicator columns, e.g. FMT1_at_FMT_indicator.
2. RFMT column labels are shown under antibiotic annotation strips.
3. Antibiotic annotation strips are added to RFMT-level plots and sample MAG-count heatmap.
4. A standalone antibiotic heatmap and legend are generated.
5. Species figure keeps species with n_units_engrafted >= 2.
6. Includes species-level and MAG-level classification heatmaps.
7. Includes top50 MAG TPM heatmap.

Input directory should contain:
    mag_transfer_calls.tsv
    species_transmission_summary.tsv
    rfmt_transfer_summary.tsv
    sample_MAG_counts.tsv
    MAG_TPM_long.tsv
    analysis_sample_metadata.tsv

Antibiotic table:
    RFMT_antibiotic_indicators_used.tsv
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


# -------------------------
# Basic helpers
# -------------------------

def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        print(f"[MISS] {path}")
        return pd.DataFrame()
    print(f"[READ] {path}")
    return pd.read_csv(path, sep="\t", dtype=str)


def to_num(df: pd.DataFrame, cols):
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def to_bool(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.lower().isin(["true", "1", "yes", "y"])


def savefig(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
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
            return sp[-1].replace("s__", "").strip() or "Unknown"
        gen = [p for p in parts if p.startswith("g__")]
        if gen:
            return gen[-1].replace("g__", "").strip() + " sp."

    if x.startswith("s__"):
        return x.replace("s__", "", 1).strip() or "Unknown"

    return x


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
# Antibiotic annotation
# -------------------------

ABX_STAGE_INFO = [
    ("FMT1_48h_prior", "FMT1_48h_prior_indicator", "F1 48h"),
    ("FMT1_at_FMT", "FMT1_at_FMT_indicator", "F1 at FMT"),
    ("FMT1_during_engraftment", "FMT1_during_engraftment_indicator", "F1 engraft"),
    ("FMT1_D7", "FMT1_D7_indicator", "F1 D7"),
    ("FMT2_at_FMT", "FMT2_at_FMT_indicator", "F2 at FMT"),
    ("FMT2_during_engraftment", "FMT2_during_engraftment_indicator", "F2 engraft"),
    ("FMT2_D7", "FMT2_D7_indicator", "F2 D7"),
]

ABX_COLORS = {
    "FMT1_48h_prior": "#66c2a5",
    "FMT1_at_FMT": "#fc8d62",
    "FMT1_during_engraftment": "#8da0cb",
    "FMT1_D7": "#e78ac3",
    "FMT2_at_FMT": "#a6d854",
    "FMT2_during_engraftment": "#ffd92f",
    "FMT2_D7": "#e5c494",
}


def load_antibiotics(indir: Path, explicit_path=None) -> pd.DataFrame:
    path = Path(explicit_path) if explicit_path else indir / "antibiotic_indicators_used.tsv"
    abx = read_tsv(path)

    if abx.empty:
        return abx

    if "RFMT_ID" not in abx.columns:
        print("[WARN] antibiotic table lacks RFMT_ID")
        return pd.DataFrame()

    # Create clean 0/1 stage columns using *_indicator when present.
    for stage, ind_col, _label in ABX_STAGE_INFO:
        if ind_col in abx.columns:
            abx[stage + "__plot"] = pd.to_numeric(abx[ind_col], errors="coerce").fillna(0).astype(int)
        elif stage in abx.columns:
            # Fallback if only text stage exists.
            abx[stage + "__plot"] = abx[stage].astype(str).str.strip().str.lower().isin(["yes", "y", "1", "true"]).astype(int)
        else:
            abx[stage + "__plot"] = 0

    return abx


def antibiotic_matrix(abx: pd.DataFrame, rfmts):
    stages = []
    labels = []

    for stage, _ind_col, label in ABX_STAGE_INFO:
        col = stage + "__plot"
        if col in abx.columns:
            stages.append((stage, col))
            labels.append(label)

    if not stages or abx.empty:
        return np.zeros((0, len(rfmts))), [], []

    abx_idx = abx.set_index("RFMT_ID")
    mat = np.zeros((len(stages), len(rfmts)), dtype=int)

    for j, r in enumerate(rfmts):
        if r in abx_idx.index:
            for i, (_stage, col) in enumerate(stages):
                val = abx_idx.loc[r, col]
                if isinstance(val, pd.Series):
                    val = val.iloc[0]
                mat[i, j] = int(pd.to_numeric(val, errors="coerce")) if pd.notna(pd.to_numeric(val, errors="coerce")) else 0

    stage_names = [x[0] for x in stages]
    return mat, stage_names, labels


def draw_abx_strip(ax, rfmts, abx, show_xlabels=True):
    """
    Draw antibiotic exposure strip below RFMT-level plots.
    Colored square = antibiotic yes for that stage.
    White square = no/unknown.
    """
    mat, stages, labels = antibiotic_matrix(abx, rfmts)

    if mat.shape[0] == 0:
        ax.axis("off")
        return

    nstage, nrfmt = mat.shape

    ax.set_xlim(-0.5, nrfmt - 0.5)
    ax.set_ylim(nstage - 0.5, -0.5)

    for i, stage in enumerate(stages):
        for j in range(nrfmt):
            color = ABX_COLORS.get(stage, "black") if mat[i, j] == 1 else "white"
            ax.add_patch(
                plt.Rectangle(
                    (j - 0.5, i - 0.5),
                    1,
                    1,
                    facecolor=color,
                    edgecolor="0.70",
                    linewidth=0.45,
                )
            )

    ax.set_yticks(range(nstage))
    ax.set_yticklabels(labels, fontsize=8)

    ax.set_xticks(range(nrfmt))
    if show_xlabels:
        ax.set_xticklabels(rfmts, rotation=90, fontsize=8)
    else:
        ax.set_xticklabels([])

    ax.tick_params(length=0)
    ax.set_ylabel("Antibiotics", fontsize=9)

    for spine in ax.spines.values():
        spine.set_visible(False)


def antibiotic_legend_handles():
    handles = [
        Patch(facecolor=ABX_COLORS[stage], edgecolor="0.5", label=label)
        for stage, _ind, label in ABX_STAGE_INFO
    ]
    handles.append(Patch(facecolor="white", edgecolor="0.5", label="No / unknown"))
    return handles


def plot_antibiotic_heatmap(abx: pd.DataFrame, rfmts, outdir: Path):
    if abx.empty:
        print("[SKIP] antibiotic table empty")
        return

    mat, stages, labels = antibiotic_matrix(abx, rfmts)
    if mat.shape[0] == 0:
        print("[SKIP] no antibiotic indicator columns found")
        return

    fig, ax = plt.subplots(figsize=(max(9, 0.45 * len(rfmts)), 2.8))

    ax.set_xlim(-0.5, len(rfmts) - 0.5)
    ax.set_ylim(len(stages) - 0.5, -0.5)

    for i, stage in enumerate(stages):
        for j in range(len(rfmts)):
            color = ABX_COLORS.get(stage, "black") if mat[i, j] == 1 else "white"
            ax.add_patch(
                plt.Rectangle(
                    (j - 0.5, i - 0.5),
                    1,
                    1,
                    facecolor=color,
                    edgecolor="0.70",
                    linewidth=0.45,
                )
            )

    ax.set_xticks(range(len(rfmts)))
    ax.set_xticklabels(rfmts, rotation=90, fontsize=8)
    ax.set_yticks(range(len(stages)))
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_title("Antibiotic exposure indicators by RFMT")
    ax.tick_params(length=0)

    ax.legend(
        handles=antibiotic_legend_handles(),
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        frameon=False,
        fontsize=8,
    )

    savefig(outdir / "FigABX_antibiotic_exposure_heatmap.png")


# -------------------------
# RFMT-level plots with antibiotic strips
# -------------------------

def plot_rfmt_counts_and_rates(rfmt_summary: pd.DataFrame, abx: pd.DataFrame, outdir: Path):
    if rfmt_summary.empty:
        print("[SKIP] rfmt_transfer_summary.tsv empty")
        return

    required = {"RFMT_ID", "donor_label", "n_MAGs_engrafted", "MAG_engraftment_rate"}
    if not required.issubset(rfmt_summary.columns):
        print("[SKIP] rfmt summary missing columns:", required - set(rfmt_summary.columns))
        return

    df = rfmt_summary.copy()
    df = to_num(df, ["n_MAGs_engrafted", "MAG_engraftment_rate", "n_MAGs_available"])
    rfmts = sorted(df["RFMT_ID"].dropna().unique())

    plot_antibiotic_heatmap(abx, rfmts, outdir)

    # Counts
    count_df = (
        df.pivot_table(
            index="RFMT_ID",
            columns="donor_label",
            values="n_MAGs_engrafted",
            aggfunc="sum",
            fill_value=0,
        )
        .reindex(rfmts)
    )

    fig = plt.figure(figsize=(max(9, 0.5 * len(rfmts)), 7.2))
    gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[4.5, 2.2], hspace=0.08)

    ax = fig.add_subplot(gs[0, 0])
    ax_abx = fig.add_subplot(gs[1, 0], sharex=ax)

    x = np.arange(len(rfmts))
    bottom = np.zeros(len(rfmts))

    for donor in count_df.columns:
        vals = count_df[donor].astype(float).values
        ax.bar(x, vals, bottom=bottom, label=str(donor))
        bottom += vals

    ax.set_xticks(x)
    ax.set_xticklabels([])
    ax.set_ylabel("Engrafted MAGs")
    ax.set_title("Donor-derived MAG engraftment counts")
    ax.legend(frameon=False, title="Donor", loc="upper left")

    draw_abx_strip(ax_abx, rfmts, abx, show_xlabels=True)

    fig.legend(
        handles=antibiotic_legend_handles(),
        loc="upper right",
        bbox_to_anchor=(1.02, 0.98),
        fontsize=7,
        frameon=False,
    )

    savefig(outdir / "Fig1_RFMT_engrafted_MAG_counts_with_antibiotics.png")

    # Rates
    rate_df = (
        df.pivot_table(
            index="RFMT_ID",
            columns="donor_label",
            values="MAG_engraftment_rate",
            aggfunc="mean",
            fill_value=0,
        )
        .reindex(rfmts)
    )

    fig = plt.figure(figsize=(max(9, 0.5 * len(rfmts)), 7.2))
    gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[4.5, 2.2], hspace=0.08)

    ax = fig.add_subplot(gs[0, 0])
    ax_abx = fig.add_subplot(gs[1, 0], sharex=ax)

    width = 0.8 / max(1, len(rate_df.columns))

    for i, donor in enumerate(rate_df.columns):
        ax.bar(x + i * width, rate_df[donor].astype(float).values, width=width, label=str(donor))

    offset = width * (len(rate_df.columns) - 1) / 2
    ax.set_xticks(x + offset)
    ax.set_xticklabels([])
    ax.set_ylabel("Engraftment rate\nengrafted / donor-available")
    ax.set_title("Corrected MAG engraftment rate")
    ax.legend(frameon=False, title="Donor", loc="upper left")

    draw_abx_strip(ax_abx, rfmts, abx, show_xlabels=True)

    fig.legend(
        handles=antibiotic_legend_handles(),
        loc="upper right",
        bbox_to_anchor=(1.02, 0.98),
        fontsize=7,
        frameon=False,
    )

    savefig(outdir / "Fig1b_corrected_MAG_engraftment_rate_with_antibiotics.png")


def plot_sample_mag_count_heatmap(sample_counts: pd.DataFrame, abx: pd.DataFrame, outdir: Path):
    if sample_counts.empty:
        print("[SKIP] sample_MAG_counts.tsv empty")
        return

    required = {"RFMT_ID", "role", "n_MAGs_quality"}
    if not required.issubset(sample_counts.columns):
        print("[SKIP] sample_MAG_counts missing columns:", required - set(sample_counts.columns))
        return

    df = sample_counts.copy()
    df = to_num(df, ["n_MAGs_quality", "n_MAGs_total"])

    roles = [
        "donor1",
        "recipient_day0",
        "recipient_pre_donor2",
        "donor2",
        "recipient_d2",
        "recipient_d5",
        "recipient_d7",
    ]

    rfmts = sorted(df["RFMT_ID"].dropna().unique())

    mat = pd.DataFrame(np.nan, index=roles, columns=rfmts)

    for _, row in df.iterrows():
        r = row["RFMT_ID"]
        role = row["role"]
        if role in mat.index and r in mat.columns:
            mat.loc[role, r] = row["n_MAGs_quality"]

    fig = plt.figure(figsize=(max(9, 0.5 * len(rfmts)), 7.5))
    gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[5.2, 2.2], hspace=0.08)

    ax = fig.add_subplot(gs[0, 0])
    ax_abx = fig.add_subplot(gs[1, 0], sharex=ax)

    im = ax.imshow(mat.astype(float).values, aspect="auto", interpolation="nearest")

    ax.set_xticks(np.arange(len(rfmts)))
    ax.set_xticklabels([])
    ax.set_yticks(np.arange(len(roles)))
    ax.set_yticklabels([role_label(r) for r in roles])
    ax.set_ylabel("Sample role / timepoint")
    ax.set_title("Number of quality-filtered MAGs per sample")

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("MAG count\nCompleteness ≥50%, Contamination ≤10%")

    draw_abx_strip(ax_abx, rfmts, abx, show_xlabels=True)

    fig.legend(
        handles=antibiotic_legend_handles(),
        loc="upper right",
        bbox_to_anchor=(1.02, 0.98),
        fontsize=7,
        frameon=False,
    )

    savefig(outdir / "Fig0_sample_MAG_count_heatmap_with_antibiotics.png")


# -------------------------
# Other result plots
# -------------------------

def plot_species_frequency(species_summary: pd.DataFrame, outdir: Path, min_events=2, top_n=50):
    if species_summary.empty:
        return

    required = {
        "species",
        "n_units_engrafted",
        "n_donor1_engrafted_units",
        "n_donor2_engrafted_units",
        "n_RFMT_donor_units_available",
        "transmission_rate_among_available_units",
    }

    if not required.issubset(species_summary.columns):
        print("[SKIP] species summary missing columns:", required - set(species_summary.columns))
        return

    df = species_summary.copy()
    df = to_num(
        df,
        [
            "n_units_engrafted",
            "n_donor1_engrafted_units",
            "n_donor2_engrafted_units",
            "n_RFMT_donor_units_available",
            "transmission_rate_among_available_units",
        ],
    )

    df = df[df["n_units_engrafted"] >= min_events].copy()

    if df.empty:
        print(f"[SKIP] no species with n_units_engrafted >= {min_events}")
        return

    df["species_label"] = df["species"].map(short_species)
    df["frac_label"] = (
        df["species_label"]
        + " ("
        + df["n_units_engrafted"].astype(int).astype(str)
        + "/"
        + df["n_RFMT_donor_units_available"].astype(int).astype(str)
        + ")"
    )

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


def plot_corrected_timeline(mag_calls: pd.DataFrame, outdir: Path):
    if mag_calls.empty:
        return

    required = {"donor_label", "post_role", "is_engrafted", "is_high_confidence_engrafted"}
    if not required.issubset(mag_calls.columns):
        print("[SKIP] mag calls missing columns:", required - set(mag_calls.columns))
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
            print(f"[SKIP] no rows for {title}")
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


def plot_classification_heatmap(mag_calls: pd.DataFrame, outdir: Path, level="species", top_n=60):
    if mag_calls.empty:
        return

    required = {"RFMT_ID", "donor_label", "classification"}
    if not required.issubset(mag_calls.columns):
        print("[SKIP] mag calls missing columns:", required - set(mag_calls.columns))
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
        ylabel = "Species"
    else:
        if "genome" not in df.columns:
            return
        df["row_label"] = df["genome"].astype(str)
        if "species" in df.columns:
            df["row_label"] = df["species"].map(short_species) + " | " + df["genome"].astype(str)
        fname = "Fig4b_MAG_classification_heatmap.png"
        title = "MAG-level classification heatmap"
        ylabel = "MAG"

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

    score = (
        df[df["classification_simple"].isin(event_classes)]
        .groupby("row_label")
        .size()
        .sort_values(ascending=False)
        .head(top_n)
    )

    if score.empty:
        print(f"[SKIP] no classified events for {level} heatmap")
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

    plt.xticks(np.arange(num.shape[1]), num.columns, rotation=90, fontsize=7)
    plt.yticks(np.arange(num.shape[0]), num.index, fontsize=7)
    plt.title(title)
    plt.xlabel("RFMT-donor combination")
    plt.ylabel(ylabel)

    cbar = plt.colorbar()
    cbar.set_ticks(list(class_to_num.values()))
    cbar.set_ticklabels(list(class_to_num.keys()))

    savefig(outdir / fname)


def plot_donor_vs_baseline_conani(mag_calls: pd.DataFrame, outdir: Path):
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


def plot_species_persistence(species_summary: pd.DataFrame, outdir: Path, min_events=2, top_n=50):
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


def plot_top50_mag_tpm(tpm_long: pd.DataFrame, sample_meta: pd.DataFrame, outdir: Path, top_n=50):
    if tpm_long.empty:
        print("[SKIP] TPM long table empty")
        return

    required = {"sample_name", "genome", "TPM"}

    if not required.issubset(tpm_long.columns):
        print("[SKIP] TPM table missing:", required - set(tpm_long.columns))
        return

    df = tpm_long.copy()
    df = to_num(df, ["TPM"])

    if "quality_pass" in df.columns:
        df = df[df["quality_pass"].astype(str).str.lower().isin(["true", "1", "yes"])]

    df = df.dropna(subset=["TPM"])

    if df.empty:
        return

    top = df.groupby("genome")["TPM"].max().sort_values(ascending=False).head(top_n).index
    df = df[df["genome"].isin(top)].copy()

    if not sample_meta.empty and {"RFMT_ID", "role_order", "sample_name"}.issubset(sample_meta.columns):
        sample_meta2 = sample_meta.copy()
        sample_meta2["role_order"] = pd.to_numeric(sample_meta2["role_order"], errors="coerce")
        sample_order = sample_meta2.sort_values(["RFMT_ID", "role_order", "sample_name"])["sample_name"].drop_duplicates().tolist()
        sample_order = [s for s in sample_order if s in set(df["sample_name"])]
    else:
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


def write_notes(outdir: Path):
    txt = """Figure notes v7

Antibiotic annotation:
The plot reads *_indicator columns from RFMT_antibiotic_indicators_used.tsv.
Colored square = antibiotic exposure recorded as yes.
White square = no or unknown.
RFMT column labels are shown below the antibiotic strip.

Corrected MAG engraftment rate:
MAG_engraftment_rate = n_MAGs_engrafted / n_MAGs_available.
n_MAGs_available is defined from quality-filtered donor MAGs in donor bintables.

Sample MAG count heatmap:
Rows are sample roles/timepoints.
Columns are RFMT groups.
Cells show number of quality-filtered MAGs.

Figure 4:
Species-level collapsed heatmap. If multiple MAGs from one species occur in one RFMT-donor combination, the strongest MAG-level classification is shown.

Figure 4b:
MAG-level heatmap. Each row is a specific MAG.

Figure 2:
Only species with n_units_engrafted >= 2 are shown.
Labels show engrafted / donor-available denominator.
"""
    path = outdir / "figure_notes_v7.txt"
    path.write_text(txt)
    print(f"[WRITE] {path}")


# -------------------------
# Main
# -------------------------

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--indir", required=True)
    parser.add_argument("--antibiotics", default=None, help="Optional antibiotic indicator TSV. Default: <indir>/antibiotic_indicators_used.tsv")
    parser.add_argument("--outdir", default=None)
    parser.add_argument("--top-species", type=int, default=50)
    parser.add_argument("--top-heatmap-rows", type=int, default=60)
    parser.add_argument("--top-mag-tpm", type=int, default=50)
    parser.add_argument("--min-species-events", type=int, default=2)

    args = parser.parse_args()

    indir = Path(args.indir)
    outdir = Path(args.outdir) if args.outdir else indir / "figures_v7_antibiotics_fixed"
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"[INDIR]  {indir}")
    print(f"[OUTDIR] {outdir}")

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
