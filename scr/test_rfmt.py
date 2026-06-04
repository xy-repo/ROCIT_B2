#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd


def find_compare_file(instrain_root: Path, rfmt: str):
    """
    Find genomeWide compare table for one RFMT.

    Expected example:
    /group/sbms004/yxia/GUT/inStrain_final/RFMT015/compare_out_with_stb/output/compare_out_with_stb_genomeWide_compare.tsv
    """
    compare_dir = instrain_root / rfmt / "compare_out_with_stb" / "output"

    if not compare_dir.exists():
        return None

    candidates = sorted(compare_dir.glob("*genomeWide_compare.tsv"))

    if not candidates:
        candidates = sorted(compare_dir.glob("*genome*compare*.tsv"))

    return candidates[0] if candidates else None


def make_full_sample_id(rfmt: str, sample: str) -> str:
    """
    Convert sample ID to full sample folder name.

    Examples:
      rfmt = RFMT015, sample = 5423        -> RFMT015_5423
      rfmt = RFMT015, sample = RFMT015_5423 -> RFMT015_5423
    """
    sample = str(sample).strip()

    # Clean Excel-style numeric strings like 5423.0
    if sample.endswith(".0"):
        sample = sample[:-2]

    if sample.startswith(rfmt + "_"):
        return sample

    if sample.startswith("RFMT"):
        return sample

    return f"{rfmt}_{sample}"


def find_bintable(mag_root: Path, rfmt: str, sample: str):
    """
    Find bintable while ignoring sample_list_* number.

    Example target:
    /group/sbms004/yxia/GUT/MAG_per_sample/sample_list_13/RFMT013_5430/RFMT013_5430/results/18.RFMT013_5430.bintable

    The sample_list_13 part is irrelevant, so this searches:
    sample_list_*/RFMT013_5430/RFMT013_5430/results/*.bintable
    """
    full_sample = make_full_sample_id(rfmt, sample)

    candidates = sorted(
        mag_root.glob(f"sample_list_*/{full_sample}/{full_sample}/results/*.bintable")
    )

    return candidates[0] if candidates else None


def normalize_bintable_to_instrain_genome(sample: str, bin_id: str) -> str:
    """
    Convert bintable Bin ID to inStrain genome name.

    Example:
      sample = RFMT015_5423
      Bin ID = RFMT015_5423.metabat2.10

    inStrain genome:
      RFMT015_5423_RFMT015_5423.metabat2.10
    """
    sample = str(sample).strip()
    bin_id = str(bin_id).strip()
    return f"{sample}_{bin_id}"


def guess_column(columns, possible_names):
    for name in possible_names:
        if name in columns:
            return name
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Test RFMT metadata, inStrain genomeWide_compare, and MAG bintable file matching."
    )

    parser.add_argument(
        "--metadata",
        required=True,
        help="Path to RFMT_analysis_sample_long.tsv",
    )

    parser.add_argument(
        "--instrain-root",
        required=True,
        help="Root folder containing RFMTxxx/compare_out_with_stb",
    )

    parser.add_argument(
        "--mag-root",
        required=True,
        help="Root folder containing sample_list_* folders",
    )

    parser.add_argument(
        "--rfmt",
        default=None,
        help="Optional: test only one RFMT, e.g. RFMT015",
    )

    parser.add_argument(
        "--out",
        default="rfmt_input_read_test_summary.tsv",
        help="Output summary TSV",
    )

    args = parser.parse_args()

    metadata = Path(args.metadata)
    instrain_root = Path(args.instrain_root)
    mag_root = Path(args.mag_root)
    out_path = Path(args.out)

    print("=== Checking input paths ===")
    print(f"metadata:      {metadata}")
    print(f"inStrain root: {instrain_root}")
    print(f"MAG root:      {mag_root}")
    print(f"output:        {out_path}")

    if not metadata.exists():
        raise FileNotFoundError(f"Metadata not found: {metadata}")

    if not instrain_root.exists():
        raise FileNotFoundError(f"inStrain root not found: {instrain_root}")

    if not mag_root.exists():
        raise FileNotFoundError(f"MAG root not found: {mag_root}")

    meta = pd.read_csv(metadata, sep="\t", dtype=str)

    print("\n=== Metadata loaded ===")
    print(f"Rows: {len(meta)}")
    print(f"Columns: {list(meta.columns)}")
    print("\nMetadata preview:")
    print(meta.head(10).to_string(index=False))

    possible_rfmt_cols = [
        "RFMT_ID",
        "rfmt_id",
        "RFMT",
        "rfmt",
        "group",
        "Group",
    ]

    possible_sample_cols = [
        "sample_id_full",
        "full_sample_id",
        "sample_full",
        "SampleID_full",
        "sample_id",
        "SampleID",
        "sample",
        "Sample",
        "sample_number",
        "Sample number",
        "Sample Number",
    ]

    rfmt_col = guess_column(meta.columns, possible_rfmt_cols)
    sample_col = guess_column(meta.columns, possible_sample_cols)

    if rfmt_col is None:
        raise ValueError(
            f"Cannot find RFMT column. Available columns: {list(meta.columns)}"
        )

    if sample_col is None:
        raise ValueError(
            f"Cannot find sample column. Available columns: {list(meta.columns)}"
        )

    print("\n=== Column detection ===")
    print(f"RFMT column:   {rfmt_col}")
    print(f"Sample column: {sample_col}")

    if args.rfmt:
        meta = meta[meta[rfmt_col] == args.rfmt].copy()

    rfmts = sorted(meta[rfmt_col].dropna().unique())

    print(f"\nRFMT groups to test: {len(rfmts)}")
    print(", ".join(rfmts[:20]) + (" ..." if len(rfmts) > 20 else ""))

    all_rows = []

    for rfmt in rfmts:
        print("\n" + "=" * 90)
        print(f"Testing {rfmt}")

        compare_file = find_compare_file(instrain_root, rfmt)

        if compare_file is None:
            print(f"[MISS] genomeWide compare file not found for {rfmt}")

            all_rows.append(
                {
                    "RFMT_ID": rfmt,
                    "status": "missing_compare",
                    "compare_file": "",
                    "n_compare_rows": 0,
                    "n_compare_genomes": 0,
                    "n_metadata_samples": 0,
                    "n_bintables_found": 0,
                    "n_bintables_missing": 0,
                    "n_bintable_genomes": 0,
                    "n_matched_genomes": 0,
                    "match_rate_vs_compare_percent": 0,
                }
            )
            continue

        print(f"[OK] compare file: {compare_file}")

        try:
            cmp_df = pd.read_csv(compare_file, sep="\t", dtype=str)
        except Exception as e:
            print(f"[FAIL] could not read compare file: {e}")

            all_rows.append(
                {
                    "RFMT_ID": rfmt,
                    "status": "bad_compare_read",
                    "compare_file": str(compare_file),
                    "n_compare_rows": 0,
                    "n_compare_genomes": 0,
                    "n_metadata_samples": 0,
                    "n_bintables_found": 0,
                    "n_bintables_missing": 0,
                    "n_bintable_genomes": 0,
                    "n_matched_genomes": 0,
                    "match_rate_vs_compare_percent": 0,
                }
            )
            continue

        print(f"[OK] compare rows: {len(cmp_df)}")
        print(f"Compare columns: {list(cmp_df.columns)[:30]}")

        genome_col = guess_column(
            cmp_df.columns,
            [
                "genome",
                "Genome",
                "genome_name",
                "Genome Name",
                "scaffold",
                "Scaffold",
            ],
        )

        if genome_col is None:
            print("[WARN] No genome/scaffold column found in compare table.")
            print(f"Columns are: {list(cmp_df.columns)}")
            compare_genomes = set()
        else:
            compare_genomes = set(cmp_df[genome_col].dropna().astype(str))
            print(f"[OK] compare genome column: {genome_col}")
            print(f"[OK] unique compare genomes: {len(compare_genomes)}")
            print("Example compare genomes:")
            for x in list(sorted(compare_genomes))[:5]:
                print(f"  {x}")

        sub_meta = meta.loc[meta[rfmt_col] == rfmt].copy()

        raw_samples = (
            sub_meta[sample_col]
            .dropna()
            .astype(str)
            .str.strip()
        )

        rfmt_samples = sorted(
            {
                make_full_sample_id(rfmt, s)
                for s in raw_samples
                if str(s).strip() not in {"", "nan", "NaN", "None"}
            }
        )

        print("\nSamples from metadata after full-ID conversion:")
        for s in rfmt_samples:
            print(f"  {s}")

        bintable_genomes = set()
        found_bintables = 0
        missing_bintables = []
        bad_bintables = []

        print("\nChecking bintables:")

        for sample in rfmt_samples:
            bintable = find_bintable(mag_root, rfmt, sample)

            if bintable is None:
                print(f"  [MISS] {sample}: no bintable found")
                missing_bintables.append(sample)
                continue

            found_bintables += 1
            print(f"  [OK] {sample}: {bintable}")

            try:
                bt = pd.read_csv(bintable, sep="\t", header=1, dtype=str)
            except Exception as e:
                print(f"  [FAIL] {sample}: could not read bintable: {e}")
                bad_bintables.append(sample)
                continue
             
            if "Bin ID" not in bt.columns:
                print(f"  [WARN] {sample}: no 'Bin ID' column.")
                print(f"         Columns: {list(bt.columns)}")
                bad_bintables.append(sample)
                continue

            for bin_id in bt["Bin ID"].dropna().astype(str):
                bintable_genomes.add(
                    normalize_bintable_to_instrain_genome(sample, bin_id)
                )

        matched = compare_genomes & bintable_genomes

        if len(compare_genomes) > 0:
            match_rate = len(matched) / len(compare_genomes) * 100
        else:
            match_rate = 0

        print("\nMatch summary:")
        print(f"  samples in metadata:       {len(rfmt_samples)}")
        print(f"  bintables found:           {found_bintables}")
        print(f"  bintables missing:         {len(missing_bintables)}")
        print(f"  bad bintables:             {len(bad_bintables)}")
        print(f"  unique bintable genomes:   {len(bintable_genomes)}")
        print(f"  unique compare genomes:    {len(compare_genomes)}")
        print(f"  matched genome names:      {len(matched)}")
        print(f"  match rate vs compare:     {match_rate:.2f}%")

        if missing_bintables:
            print("\nMissing bintable samples:")
            for s in missing_bintables[:20]:
                print(f"  {s}")

        if len(matched) == 0 and len(compare_genomes) > 0 and len(bintable_genomes) > 0:
            print("\n[WARN] No genome names matched.")
            print("This usually means the genome-name conversion rule is wrong.")
            print("Example compare genome:")
            print(f"  {next(iter(compare_genomes))}")
            print("Example bintable-derived genome:")
            print(f"  {next(iter(bintable_genomes))}")

        else:
            print("\nExample matched genomes:")
            for x in list(sorted(matched))[:5]:
                print(f"  {x}")

        all_rows.append(
            {
                "RFMT_ID": rfmt,
                "status": "ok",
                "compare_file": str(compare_file),
                "n_compare_rows": len(cmp_df),
                "n_compare_genomes": len(compare_genomes),
                "n_metadata_samples": len(rfmt_samples),
                "n_bintables_found": found_bintables,
                "n_bintables_missing": len(missing_bintables),
                "n_bad_bintables": len(bad_bintables),
                "n_bintable_genomes": len(bintable_genomes),
                "n_matched_genomes": len(matched),
                "match_rate_vs_compare_percent": round(match_rate, 3),
            }
        )

    summary = pd.DataFrame(all_rows)
    summary.to_csv(out_path, sep="\t", index=False)

    print("\n" + "=" * 90)
    print("Finished.")
    print(f"Summary written to: {out_path.resolve()}")
    print("\nSummary:")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()