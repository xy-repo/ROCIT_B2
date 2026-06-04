#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        print(f"[MISS] {path}")
        return pd.DataFrame()

    print(f"[READ] {path}")
    return pd.read_csv(path, sep="\t", dtype=str)


def to_numeric(df: pd.DataFrame, cols):
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def to_bool_series(s: pd.Series) -> pd.Series:
    return (
        s.astype(str)
        .str.strip()
        .str.lower()
        .isin(["true", "1", "yes", "y"])
    )


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
        species_parts = [p for p in parts if p.startswith("s__")]
        if species_parts:
            sp = species_parts[-1].replace("s__", "").strip()
            return sp if sp else "Unknown"

    if x.startswith("s__"):
        return x.replace("s__", "", 1).strip() or "Unknown"

    return x


def role_to_day(role):
    if pd.isna(role):
        return "unknown"

    s = str(role).lower()

    if "d7" in s or "day7" in s:
        return "day7"
    if "d5" in s or "day5" in s:
        return "day5"
    if "d2" in s or "day2" in s:
        return "day2"
    if "pre" in s:
        return "preDonor2"
    if "d0" in s or "day0" in s:
        return "day0"

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
    if "recipient" in s and ("persist" in s or "like" in s or "resident" in s):
        return "Recipient-like"
    if "shared" in s:
        return "Already shared"
    if "replace" in s:
        return "Replacement"
    if "not" in s:
        return "Not transmitted"
    if "uncertain" in s or "low" in s:
        return "Uncertain"

    return raw


def ensure_required(df: pd.DataFrame, cols, table_name: str) -> bool:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        print(f"[SKIP] {table_name}: missing columns {missing}")
        print(f"Available columns: {list(df.columns)}")
        return False
    return True


def plot_rfmt_engraftment_counts(rfmt_summary: pd.DataFrame, outdir: Path):
    if rfmt_summary.empty:
        print("[SKIP] rfmt_transfer_summary.tsv is empty")
        return

    required = [
        "RFMT_ID",
        "donor_label",
        "n_MAGs_engrafted",
        "n_species_engrafted",
        "MAG_engraftment_rate",
        "species_engraftment_rate",
    ]

    if not ensure_required(rfmt_summary, required, "rfmt_transfer_summary.tsv"):
        return

    df = rfmt_summary.copy()
    df = to_numeric(
        df,
        [
            "n_MAGs_tested",
            "n_MAGs_engrafted",
            "n_species_engrafted",
            "MAG_engraftment_rate",
            "species_engraftment_rate",
        ],
    )

    count_df = (
        df.pivot_table(
            index="RFMT_ID",
            columns="donor_label",
            values="n_MAGs_engrafted",
            aggfunc="sum",
            fill_value=0,
        )
        .sort_index()
    )

    if count_df.empty:
        print("[SKIP] No RFMT engraftment counts to plot")
        return

    plt.figure(figsize=(max(8, 0.45 * len(count_df)), 5))

    x = np.arange(len(count_df))
    bottom = np.zeros(len(count_df))

    for donor in count_df.columns:
        vals = count_df[donor].astype(float).values
        plt.bar(x, vals, bottom=bottom, label=str(donor))
        bottom += vals

    plt.xticks(x, count_df.index, rotation=90)
    plt.xlabel("RFMT group")
    plt.ylabel("Number of engrafted MAGs")
    plt.title("Donor-derived MAG engraftment across RFMT groups")
    plt.legend(title="Donor", frameon=False)

    savefig(outdir / "Fig1_RFMT_engrafted_MAG_counts_by_donor.png")

    rate_df = (
        df.pivot_table(
            index="RFMT_ID",
            columns="donor_label",
            values="MAG_engraftment_rate",
            aggfunc="mean",
            fill_value=0,
        )
        .sort_index()
    )

    if not rate_df.empty:
        plt.figure(figsize=(max(8, 0.45 * len(rate_df)), 5))

        x = np.arange(len(rate_df))
        width = 0.8 / max(1, len(rate_df.columns))

        for i, donor in enumerate(rate_df.columns):
            vals = rate_df[donor].astype(float).values
            plt.bar(x + i * width, vals, width=width, label=str(donor))

        offset = width * (len(rate_df.columns) - 1) / 2
        plt.xticks(x + offset, rate_df.index, rotation=90)
        plt.xlabel("RFMT group")
        plt.ylabel("MAG engraftment rate")
        plt.title("MAG engraftment rate across RFMT groups")
        plt.legend(title="Donor", frameon=False)

        savefig(outdir / "Fig1b_RFMT_MAG_engraftment_rate_by_donor.png")


def plot_species_transmission_frequency(species_summary: pd.DataFrame, outdir: Path, top_n: int = 25):
    if species_summary.empty:
        print("[SKIP] species_transmission_summary.tsv is empty")
        return

    required = [
        "species",
        "n_donor1_engrafted_units",
        "n_donor2_engrafted_units",
        "n_units_engrafted",
        "n_units_high_conf_engrafted",
        "n_units_persistent_to_d7",
        "transmission_rate_among_tested_units",
    ]

    if not ensure_required(species_summary, required, "species_transmission_summary.tsv"):
        return

    df = species_summary.copy()
    df["species_short"] = df["species"].map(short_species)

    df = to_numeric(
        df,
        [
            "n_donor1_engrafted_units",
            "n_donor2_engrafted_units",
            "n_units_engrafted",
            "n_units_high_conf_engrafted",
            "n_units_persistent_to_d7",
            "transmission_rate_among_tested_units",
        ],
    )

    df["plot_total"] = df["n_units_engrafted"].fillna(0)
    df = df[df["plot_total"] > 0].copy()

    if df.empty:
        print("[SKIP] No transmitted species found")
        return

    df = df.sort_values("plot_total", ascending=False).head(top_n)
    df = df.iloc[::-1]

    y = np.arange(len(df))
    d1 = df["n_donor1_engrafted_units"].fillna(0).astype(float).values
    d2 = df["n_donor2_engrafted_units"].fillna(0).astype(float).values

    plt.figure(figsize=(9, max(5, 0.32 * len(df))))
    plt.barh(y, d1, label="Donor1")
    plt.barh(y, d2, left=d1, label="Donor2")
    plt.yticks(y, df["species_short"])
    plt.xlabel("Number of donor units with successful transmission")
    plt.ylabel("Species")
    plt.title(f"Top {top_n} transmitted species")
    plt.legend(frameon=False)

    savefig(outdir / "Fig2_transmitted_species_frequency.png")

    df_rate = df.sort_values("transmission_rate_among_tested_units", ascending=True)

    plt.figure(figsize=(9, max(5, 0.32 * len(df_rate))))
    plt.barh(
        np.arange(len(df_rate)),
        df_rate["transmission_rate_among_tested_units"].fillna(0).astype(float).values,
    )
    plt.yticks(np.arange(len(df_rate)), df_rate["species_short"])
    plt.xlabel("Transmission rate among tested donor units")
    plt.ylabel("Species")
    plt.title(f"Transmission rate of top {top_n} transmitted species")

    savefig(outdir / "Fig2b_species_transmission_rate.png")


def plot_persistence_over_time(mag_calls: pd.DataFrame, outdir: Path):
    if mag_calls.empty:
        print("[SKIP] mag_transfer_calls.tsv is empty")
        return

    required = [
        "donor_label",
        "post_role",
        "is_engrafted",
        "is_high_confidence_engrafted",
    ]

    if not ensure_required(mag_calls, required, "mag_transfer_calls.tsv"):
        return

    df = mag_calls.copy()
    df["day"] = df["post_role"].map(role_to_day)
    df["is_engrafted_bool"] = to_bool_series(df["is_engrafted"])
    df["is_high_conf_bool"] = to_bool_series(df["is_high_confidence_engrafted"])

    df = df[df["day"].isin(["day2", "day5", "day7"])].copy()

    if df.empty:
        print("[SKIP] No day2/day5/day7 rows found")
        return

    day_order = ["day2", "day5", "day7"]

    for flag_col, label, filename in [
        ("is_engrafted_bool", "Engrafted MAG calls", "Fig3_engrafted_MAGs_over_time.png"),
        ("is_high_conf_bool", "High-confidence engrafted MAG calls", "Fig3b_high_conf_engrafted_MAGs_over_time.png"),
    ]:
        sub = df[df[flag_col]].copy()

        if sub.empty:
            print(f"[SKIP] No rows for {label}")
            continue

        count_df = (
            sub.groupby(["donor_label", "day"])
            .size()
            .reset_index(name="n_calls")
        )

        plt.figure(figsize=(6, 4))

        for donor, donor_df in count_df.groupby("donor_label"):
            donor_df = donor_df.set_index("day").reindex(day_order).fillna(0)
            plt.plot(day_order, donor_df["n_calls"].values, marker="o", label=str(donor))

        plt.xlabel("Recipient post-transfer timepoint")
        plt.ylabel("Number of MAG calls")
        plt.title(label + " over time")
        plt.legend(title="Donor", frameon=False)

        savefig(outdir / filename)


def plot_classification_heatmap(mag_calls: pd.DataFrame, outdir: Path, top_n: int = 60):
    if mag_calls.empty:
        print("[SKIP] mag_transfer_calls.tsv is empty")
        return

    required = ["RFMT_ID", "donor_label", "species", "classification"]

    if not ensure_required(mag_calls, required, "mag_transfer_calls.tsv"):
        return

    df = mag_calls.copy()
    df["species_short"] = df["species"].map(short_species)
    df["classification_simple"] = df["classification"].map(classification_simple)

    df["col_label"] = df["RFMT_ID"].astype(str) + "_" + df["donor_label"].astype(str)

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
        .groupby("species_short")
        .size()
        .sort_values(ascending=False)
        .head(top_n)
    )

    if score.empty:
        print("[SKIP] No classified species events for heatmap")
        return

    df = df[df["species_short"].isin(score.index)].copy()
    df["priority"] = df["classification_simple"].map(priority).fillna(2)

    idx = df.sort_values("priority").groupby(["species_short", "col_label"]).tail(1)

    mat = (
        idx.pivot(index="species_short", columns="col_label", values="classification_simple")
        .fillna("Absent")
    )

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
    num = num.apply(pd.to_numeric, errors="coerce").fillna(class_to_num["Uncertain"])
    num = num.astype(float)

    plt.figure(figsize=(max(10, 0.35 * num.shape[1]), max(6, 0.25 * num.shape[0])))
    plt.imshow(num.to_numpy(), aspect="auto", interpolation="nearest")

    plt.xticks(np.arange(num.shape[1]), num.columns, rotation=90)
    plt.yticks(np.arange(num.shape[0]), num.index)
    plt.xlabel("RFMT donor unit")
    plt.ylabel("Species")
    plt.title("Species-level strain-source classification across RFMT donor units")

    cbar = plt.colorbar()
    cbar.set_ticks(list(class_to_num.values()))
    cbar.set_ticklabels(list(class_to_num.keys()))

    savefig(outdir / "Fig4_species_classification_heatmap.png")


def plot_donor_vs_baseline_similarity(mag_calls: pd.DataFrame, outdir: Path):
    if mag_calls.empty:
        print("[SKIP] mag_transfer_calls.tsv is empty")
        return

    required = [
        "donor_label",
        "donor_post_conANI",
        "baseline_post_conANI",
        "is_engrafted",
        "classification",
    ]

    if not ensure_required(mag_calls, required, "mag_transfer_calls.tsv"):
        return

    df = mag_calls.copy()
    df = to_numeric(df, ["donor_post_conANI", "baseline_post_conANI"])
    df = df.dropna(subset=["donor_post_conANI", "baseline_post_conANI"])

    if df.empty:
        print("[SKIP] No rows with both donor_post_conANI and baseline_post_conANI")
        return

    df["is_engrafted_bool"] = to_bool_series(df["is_engrafted"])

    plt.figure(figsize=(5.5, 5.5))

    non = df[~df["is_engrafted_bool"]]
    yes = df[df["is_engrafted_bool"]]

    if not non.empty:
        plt.scatter(
            non["baseline_post_conANI"],
            non["donor_post_conANI"],
            alpha=0.35,
            label="Not engrafted",
        )

    if not yes.empty:
        plt.scatter(
            yes["baseline_post_conANI"],
            yes["donor_post_conANI"],
            alpha=0.75,
            label="Engrafted",
        )

    min_val = min(df["baseline_post_conANI"].min(), df["donor_post_conANI"].min())
    max_val = max(df["baseline_post_conANI"].max(), df["donor_post_conANI"].max())

    plt.plot([min_val, max_val], [min_val, max_val], linestyle="--")

    plt.xlabel("Baseline vs post conANI")
    plt.ylabel("Donor vs post conANI")
    plt.title("Donor-post similarity versus baseline-post similarity")
    plt.legend(frameon=False)

    savefig(outdir / "Fig5_donor_post_vs_baseline_post_conANI.png")


def plot_species_persistence(species_summary: pd.DataFrame, outdir: Path, top_n: int = 25):
    if species_summary.empty:
        print("[SKIP] species_transmission_summary.tsv is empty")
        return

    required = [
        "species",
        "n_units_engrafted",
        "n_units_persistent_to_d7",
    ]

    if not ensure_required(species_summary, required, "species_transmission_summary.tsv"):
        return

    df = species_summary.copy()
    df["species_short"] = df["species"].map(short_species)
    df = to_numeric(df, ["n_units_engrafted", "n_units_persistent_to_d7"])

    df = df[df["n_units_engrafted"].fillna(0) > 0].copy()

    if df.empty:
        print("[SKIP] No engrafted species for persistence plot")
        return

    df["persistence_fraction"] = (
        df["n_units_persistent_to_d7"].fillna(0) / df["n_units_engrafted"].replace(0, np.nan)
    )

    df = df.sort_values(["n_units_persistent_to_d7", "n_units_engrafted"], ascending=False).head(top_n)
    df = df.iloc[::-1]

    y = np.arange(len(df))

    plt.figure(figsize=(9, max(5, 0.32 * len(df))))
    plt.barh(y, df["n_units_engrafted"].fillna(0).astype(float).values, label="Engrafted")
    plt.barh(y, df["n_units_persistent_to_d7"].fillna(0).astype(float).values, label="Persistent to day7")
    plt.yticks(y, df["species_short"])
    plt.xlabel("Number of donor units")
    plt.ylabel("Species")
    plt.title("Persistence of transmitted species to day7")
    plt.legend(frameon=False)

    savefig(outdir / "Fig6_species_persistence_to_day7.png")


def write_result_sentence_templates(outdir: Path):
    text = """Result sentence templates

Figure 1:
MAG-level donor engraftment varied across RFMT groups. Both Donor1- and Donor2-associated strain transmission events were detected, but the number and rate of engrafted MAGs differed across donor-recipient units.

Figure 2:
Successful transmission was concentrated in a subset of species. Species-level transmission was inferred only after at least one MAG assigned to that species showed donor-like strain similarity by genome-wide inStrain comparison.

Figure 3:
Longitudinal analysis of recipient post-transfer samples showed that donor-like MAG calls varied across day2, day5, and day7, distinguishing transient detection from persistence.

Figure 4:
The species-level heatmap summarizes the dominant strain-source classification across RFMT donor units, separating donor-like engraftment, recipient-like persistence, already-shared strains, and uncertain calls.

Figure 5:
Comparison of donor-post and baseline-post conANI values provides evidence for donor-associated strain replacement or engraftment when recipient post-transfer populations are closer to donor than to baseline.

Figure 6:
A subset of successfully transmitted species persisted to day7, suggesting stable colonization rather than transient detection.
"""
    path = outdir / "figure_result_sentence_templates.txt"
    path.write_text(text)
    print(f"[WRITE] {path}")


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--indir",
        default="/group/sbms004/yxia/GUT/RFMT_tran_v2",
        help="Directory containing transfer analysis output tables",
    )

    parser.add_argument(
        "--outdir",
        default=None,
        help="Output figure directory. Default: <indir>/figures",
    )

    parser.add_argument("--top-species", type=int, default=25)
    parser.add_argument("--top-heatmap-rows", type=int, default=60)

    args = parser.parse_args()

    indir = Path(args.indir)
    outdir = Path(args.outdir) if args.outdir else indir / "figures"
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"[INDIR]  {indir}")
    print(f"[OUTDIR] {outdir}")

    mag_calls = read_tsv(indir / "mag_transfer_calls.tsv")
    species_summary = read_tsv(indir / "species_transmission_summary.tsv")
    rfmt_summary = read_tsv(indir / "rfmt_transfer_summary.tsv")

    plot_rfmt_engraftment_counts(rfmt_summary, outdir)
    plot_species_transmission_frequency(species_summary, outdir, top_n=args.top_species)
    plot_persistence_over_time(mag_calls, outdir)
    plot_classification_heatmap(mag_calls, outdir, top_n=args.top_heatmap_rows)
    plot_donor_vs_baseline_similarity(mag_calls, outdir)
    plot_species_persistence(species_summary, outdir, top_n=args.top_species)

    write_result_sentence_templates(outdir)

    print("\nDone. Figures are in:")
    print(outdir)


if __name__ == "__main__":
    main()