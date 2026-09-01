#!/usr/bin/env python3

import argparse
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

def natural_sort_key(value):
    """
    Natural sorting so RH2 comes before RH10.
    """
    return [
        int(x) if x.isdigit() else x.lower()
        for x in re.split(r"(\d+)", str(value))
    ]


def get_dynamic_bins(values, multiplier=1.75):
    """
    Determine histogram bin number dynamically.

    Uses NumPy's automatic bin selection as a starting point,
    then increases the number of bins by the specified multiplier.
    """
    edges = np.histogram_bin_edges(
        values,
        bins="auto"
    )

    base_bins = len(edges) - 1

    return max(
        5,
        int(np.ceil(base_bins * multiplier))
    )


def load_run_info(kallisto_dir):
    """
    Find all */run_info.json files and return one row per sample.
    """
    records = []

    for json_file in kallisto_dir.glob("*/run_info.json"):

        sample = json_file.parent.name

        try:
            with open(json_file) as f:
                info = json.load(f)

            records.append(
                {
                    "sample": sample,
                    "n_processed": info["n_processed"],
                    "n_pseudoaligned": info["n_pseudoaligned"],
                    "n_unique": info["n_unique"],
                    "p_pseudoaligned": info["p_pseudoaligned"],
                    "p_unique": info["p_unique"],
                }
            )

        except (json.JSONDecodeError, KeyError) as e:
            print(
                f"WARNING: skipping {json_file}: {e}"
            )

    if not records:
        raise RuntimeError(
            f"No valid run_info.json files found under {kallisto_dir}"
        )

    df = pd.DataFrame(records)

    ordered_samples = sorted(
        df["sample"],
        key=natural_sort_key
    )

    df["sample"] = pd.Categorical(
        df["sample"],
        categories=ordered_samples,
        ordered=True
    )

    df = (
        df
        .sort_values("sample")
        .reset_index(drop=True)
    )

    return df


# ---------------------------------------------------------------------
# Report pages
# ---------------------------------------------------------------------

def add_summary_page(pdf, df):

    fig = plt.figure(
        figsize=(11, 8.5)
    )

    fig.suptitle(
        "Kallisto Quantification QC Summary",
        fontsize=20,
        fontweight="bold",
        y=0.95
    )

    total_processed = df["n_processed"].sum()
    total_pseudo = df["n_pseudoaligned"].sum()

    summary_text = (
        f"Samples processed: {len(df)}\n\n"
        f"Total fragments processed: {total_processed:,}\n"
        f"Total pseudoaligned fragments: {total_pseudo:,}\n\n"
        f"Pseudoalignment rate\n"
        f"  Median: {df['p_pseudoaligned'].median():.2f}%\n"
        f"  Mean:   {df['p_pseudoaligned'].mean():.2f}%\n"
        f"  Range:  {df['p_pseudoaligned'].min():.2f}% - "
        f"{df['p_pseudoaligned'].max():.2f}%\n\n"
        f"Unique pseudoalignment rate\n"
        f"  Median: {df['p_unique'].median():.2f}%\n"
        f"  Mean:   {df['p_unique'].mean():.2f}%\n"
        f"  Range:  {df['p_unique'].min():.2f}% - "
        f"{df['p_unique'].max():.2f}%\n\n"
        f"Fragments processed per sample\n"
        f"  Median: {df['n_processed'].median() / 1e6:.2f} million\n"
        f"  Range:  {df['n_processed'].min() / 1e6:.2f} - "
        f"{df['n_processed'].max() / 1e6:.2f} million"
    )

    fig.text(
        0.08,
        0.82,
        summary_text,
        fontsize=13,
        va="top",
        family="monospace"
    )

    metrics = [
        "n_processed",
        "n_pseudoaligned",
        "n_unique",
        "p_pseudoaligned",
        "p_unique",
    ]

    summary = (
        df[metrics]
        .describe()
        .loc[
            ["mean", "50%", "min", "max"]
        ]
        .T
    )

    summary.columns = [
        "Mean",
        "Median",
        "Minimum",
        "Maximum"
    ]

    formatted = []

    for metric, row in summary.iterrows():

        if metric.startswith("n_"):
            values = [
                f"{x:,.0f}"
                for x in row
            ]

        else:
            values = [
                f"{x:.2f}"
                for x in row
            ]

        formatted.append(
            [metric] + values
        )

    ax = fig.add_axes(
        [0.48, 0.30, 0.47, 0.48]
    )

    ax.axis("off")

    table = ax.table(
        cellText=formatted,
        colLabels=[
            "Metric",
            "Mean",
            "Median",
            "Minimum",
            "Maximum"
        ],
        loc="center",
        cellLoc="center"
    )

    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.6)

    pdf.savefig(
        fig,
        bbox_inches="tight"
    )

    plt.close(fig)


def add_fragment_counts_plot(pdf, df):

    plot_df = (
        df
        .sort_values("n_processed")
        .reset_index(drop=True)
    )

    x = np.arange(
        len(plot_df)
    )

    fig, ax = plt.subplots(
        figsize=(11, 7)
    )

    ax.plot(
        x,
        plot_df["n_processed"] / 1e6,
        marker="o",
        label="Processed"
    )

    ax.plot(
        x,
        plot_df["n_pseudoaligned"] / 1e6,
        marker="o",
        label="Pseudoaligned"
    )

    ax.plot(
        x,
        plot_df["n_unique"] / 1e6,
        marker="o",
        label="Unique"
    )

    ax.set_xticks(x)

    ax.set_xticklabels(
        plot_df["sample"],
        rotation=90
    )

    ax.set_ylabel(
        "Fragments (millions)"
    )

    ax.set_xlabel(
        "Sample"
    )

    ax.set_title(
        "Sequencing depth and Kallisto assignment"
    )

    ax.legend(
        frameon=False
    )

    ax.grid(
        axis="y",
        alpha=0.25
    )

    plt.tight_layout()

    pdf.savefig(fig)
    plt.close(fig)


def add_distribution_plots(pdf, df):

    # ---------------------------------------------------------
    # Pseudoalignment distribution
    # ---------------------------------------------------------

    pseudo_bins = get_dynamic_bins(
        df["p_pseudoaligned"]
    )

    fig, ax = plt.subplots(
        figsize=(10, 6)
    )

    ax.hist(
        df["p_pseudoaligned"],
        bins=pseudo_bins,
        edgecolor="black"
    )

    ax.axvline(
        df["p_pseudoaligned"].median(),
        linestyle="--",
        label=(
            f"Median = "
            f"{df['p_pseudoaligned'].median():.1f}%"
        )
    )

    ax.set_xlabel(
        "Pseudoalignment (%)"
    )

    ax.set_ylabel(
        "Number of samples"
    )

    ax.set_title(
        "Distribution of pseudoalignment rates"
    )

    ax.legend(
        frameon=False
    )

    plt.tight_layout()

    pdf.savefig(fig)
    plt.close(fig)


    # ---------------------------------------------------------
    # Unique pseudoalignment distribution
    # ---------------------------------------------------------

    unique_bins = get_dynamic_bins(
        df["p_unique"]
    )

    fig, ax = plt.subplots(
        figsize=(10, 6)
    )

    ax.hist(
        df["p_unique"],
        bins=unique_bins,
        edgecolor="black"
    )

    ax.axvline(
        df["p_unique"].median(),
        linestyle="--",
        label=(
            f"Median = "
            f"{df['p_unique'].median():.1f}%"
        )
    )

    ax.set_xlabel(
        "Unique pseudoalignment (%)"
    )

    ax.set_ylabel(
        "Number of samples"
    )

    ax.set_title(
        "Distribution of unique pseudoalignment rates"
    )

    ax.legend(
        frameon=False
    )

    plt.tight_layout()

    pdf.savefig(fig)
    plt.close(fig)


def add_scatter_plot(pdf, df):

    fig, ax = plt.subplots(
        figsize=(8, 7)
    )

    ax.scatter(
        df["p_pseudoaligned"],
        df["p_unique"],
        s=50
    )

    ax.set_xlabel(
        "Pseudoaligned (%)"
    )

    ax.set_ylabel(
        "Unique (%)"
    )

    ax.set_title(
        "Pseudoalignment rate vs. unique assignment"
    )

    ax.grid(
        alpha=0.25
    )

    # Label samples at the extremes of pseudoalignment rate
    pseudo_low = df.nsmallest(
        min(3, len(df)),
        "p_pseudoaligned"
    )

    pseudo_high = df.nlargest(
        min(3, len(df)),
        "p_pseudoaligned"
    )

    label_df = (
        pd.concat(
            [
                pseudo_low,
                pseudo_high
            ]
        )
        .drop_duplicates()
    )

    for _, row in label_df.iterrows():

        ax.annotate(
            str(row["sample"]),
            (
                row["p_pseudoaligned"],
                row["p_unique"]
            ),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8
        )

    plt.tight_layout()

    pdf.savefig(fig)
    plt.close(fig)


def add_table_pages(
    pdf,
    df,
    rows_per_page=25
):

    table_df = df.copy()

    table_df["n_processed"] = (
        table_df["n_processed"]
        .map(
            lambda x: f"{x:,}"
        )
    )

    table_df["n_pseudoaligned"] = (
        table_df["n_pseudoaligned"]
        .map(
            lambda x: f"{x:,}"
        )
    )

    table_df["n_unique"] = (
        table_df["n_unique"]
        .map(
            lambda x: f"{x:,}"
        )
    )

    table_df["p_pseudoaligned"] = (
        table_df["p_pseudoaligned"]
        .map(
            lambda x: f"{x:.1f}"
        )
    )

    table_df["p_unique"] = (
        table_df["p_unique"]
        .map(
            lambda x: f"{x:.1f}"
        )
    )

    table_df.columns = [
        "Sample",
        "Processed",
        "Pseudoaligned",
        "Unique",
        "Pseudoaligned (%)",
        "Unique (%)"
    ]

    n_pages = math.ceil(
        len(table_df) /
        rows_per_page
    )

    for page in range(n_pages):

        start = (
            page *
            rows_per_page
        )

        end = min(
            start + rows_per_page,
            len(table_df)
        )

        subset = table_df.iloc[
            start:end
        ]

        fig, ax = plt.subplots(
            figsize=(11, 8.5)
        )

        ax.axis("off")

        ax.set_title(
            (
                "Per-sample Kallisto statistics "
                f"({start + 1}-{end} of {len(table_df)})"
            ),
            fontsize=15,
            pad=20
        )

        table = ax.table(
            cellText=subset.values,
            colLabels=subset.columns,
            loc="center",
            cellLoc="center"
        )

        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(1, 1.35)

        for col in range(
            len(subset.columns)
        ):
            table[
                (0, col)
            ].set_text_props(
                weight="bold"
            )

        plt.tight_layout()

        pdf.savefig(fig)
        plt.close(fig)


# ---------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------

def make_report(
    kallisto_dir,
    output_pdf,
    output_table
):

    df = load_run_info(
        kallisto_dir
    )

    # Save machine-readable summary table
    df.to_csv(
        output_table,
        sep="\t",
        index=False
    )

    with PdfPages(
        output_pdf
    ) as pdf:

        add_summary_page(
            pdf,
            df
        )

        add_fragment_counts_plot(
            pdf,
            df
        )

        add_distribution_plots(
            pdf,
            df
        )

        add_scatter_plot(
            pdf,
            df
        )

        add_table_pages(
            pdf,
            df
        )

        metadata = pdf.infodict()

        metadata["Title"] = (
            "Kallisto Quantification QC Report"
        )

        metadata["Subject"] = (
            "Summary of Kallisto run_info.json statistics"
        )

    print(
        f"Processed {len(df)} samples"
    )

    print(
        f"PDF report: {output_pdf}"
    )

    print(
        f"Summary table: {output_table}"
    )


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Generate a cohort-level PDF QC report "
            "from Kallisto run_info.json files."
        )
    )

    parser.add_argument(
        "kallisto_dir",
        type=Path,
        help=(
            "Directory containing one Kallisto "
            "subdirectory per sample"
        )
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "kallisto_qc_report.pdf"
        ),
        help="Output PDF"
    )

    parser.add_argument(
        "--table",
        type=Path,
        default=Path(
            "kallisto_qc_summary.tsv"
        ),
        help="Output TSV summary"
    )

    args = parser.parse_args()

    make_report(
        kallisto_dir=args.kallisto_dir,
        output_pdf=args.output,
        output_table=args.table
    )


if __name__ == "__main__":
    main()
