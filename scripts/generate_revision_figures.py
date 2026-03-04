#!/usr/bin/env python3
"""
Generate revision figures/tables for PB-S-1341 from benchmark CSV outputs.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Generate revision figures/tables from benchmark CSVs.")
    p.add_argument("--gap-csv", type=Path, required=True, help="Combined gap_sweep CSV.")
    p.add_argument("--site-csv", type=Path, required=True, help="Combined site_summary CSV.")
    p.add_argument("--out-dir", type=Path, required=True, help="Output directory under revision_docs.")
    return p


def ensure_out_dirs(base: Path) -> dict[str, Path]:
    dirs = {
        "figures": base / "figures",
        "tables": base / "tables",
        "figure_data": base / "figure_data",
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)
    return dirs


def save_figure2(gap_df: pd.DataFrame, out_fig: Path, out_data: Path) -> None:
    fig_data = gap_df[
        [
            "dataset",
            "gap_size",
            "runtime_sec",
            "median_sn_ratio",
            "off_target_q95",
            "off_target_n",
            "detected_rate",
        ]
    ].copy()
    fig_data.sort_values(["dataset", "gap_size"], inplace=True)
    fig_data.to_csv(out_data, index=False)

    colors = {"Rf1A": "#1f77b4", "CRP1": "#ff7f0e"}
    datasets = list(fig_data["dataset"].drop_duplicates())

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), dpi=200)

    ax0 = axes[0]
    for ds in datasets:
        sub = fig_data[fig_data["dataset"] == ds].sort_values("gap_size")
        ax0.plot(
            sub["gap_size"],
            sub["median_sn_ratio"],
            marker="o",
            linewidth=2.2,
            label=ds,
            color=colors.get(ds, None),
        )
    ax0.set_xlabel("Gap Size")
    ax0.set_ylabel("Median S/N")
    ax0.set_title("A. Gap Size vs S/N")
    ax0.grid(alpha=0.25, linestyle="--")
    ax0.legend(frameon=False)

    ax1 = axes[1]
    for ds in datasets:
        sub = fig_data[fig_data["dataset"] == ds].sort_values("gap_size")
        ax1.plot(
            sub["gap_size"],
            sub["runtime_sec"],
            marker="s",
            linewidth=2.2,
            label=ds,
            color=colors.get(ds, None),
        )
    ax1.set_xlabel("Gap Size")
    ax1.set_ylabel("Runtime (sec)")
    ax1.set_title("B. Gap Size vs Runtime")
    ax1.grid(alpha=0.25, linestyle="--")

    fig.suptitle("Figure 2. Gap-size simulation: S/N and runtime trade-off", fontsize=13, y=1.02)
    fig.tight_layout()
    fig.savefig(out_fig, bbox_inches="tight")
    plt.close(fig)


def save_figure3(site_df: pd.DataFrame, out_fig: Path, out_data: Path) -> None:
    fig_data = site_df[
        [
            "dataset",
            "gap_size",
            "off_target_ge_validated",
            "rank",
            "rank_percentile",
            "sn_ratio",
            "validated_score",
        ]
    ].copy()
    fig_data.sort_values(["dataset", "gap_size"], inplace=True)
    fig_data.to_csv(out_data, index=False)

    colors = {"Rf1A": "#1f77b4", "CRP1": "#ff7f0e"}
    datasets = list(fig_data["dataset"].drop_duplicates())

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), dpi=200)

    ax0 = axes[0]
    for ds in datasets:
        sub = fig_data[fig_data["dataset"] == ds].sort_values("gap_size")
        ax0.plot(
            sub["gap_size"],
            sub["off_target_ge_validated"],
            marker="o",
            linewidth=2.2,
            label=ds,
            color=colors.get(ds, None),
        )
    ax0.set_xlabel("Gap Size")
    ax0.set_ylabel("Off-target windows >= validated score")
    ax0.set_yscale("log")
    ax0.set_title("A. Potential off-target burden")
    ax0.grid(alpha=0.25, linestyle="--")
    ax0.legend(frameon=False)

    ax1 = axes[1]
    for ds in datasets:
        sub = fig_data[fig_data["dataset"] == ds].sort_values("gap_size")
        ax1.plot(
            sub["gap_size"],
            sub["rank"],
            marker="s",
            linewidth=2.2,
            label=ds,
            color=colors.get(ds, None),
        )
    ax1.set_xlabel("Gap Size")
    ax1.set_ylabel("Validated-site rank (lower is better)")
    ax1.invert_yaxis()
    ax1.set_title("B. Validated-site rank")
    ax1.grid(alpha=0.25, linestyle="--")

    fig.suptitle("Figure 3. Off-target analysis and validated-site ranking", fontsize=13, y=1.02)
    fig.tight_layout()
    fig.savefig(out_fig, bbox_inches="tight")
    plt.close(fig)


def _format_value(col: str, val: object) -> str:
    if pd.isna(val):
        return ""
    if col in {"runtime_sec"}:
        return f"{float(val):.6f}"
    if col in {"validated_score", "sn_ratio", "rank_percentile", "off_target_q95"}:
        return f"{float(val):.4f}"
    if col in {"detected_rate"}:
        return f"{float(val):.3f}"
    if isinstance(val, float):
        return f"{val:.4f}"
    return str(val)


def write_markdown_table(df: pd.DataFrame, path: Path) -> None:
    cols = list(df.columns)
    header = "| " + " | ".join(cols) + " |"
    sep = "| " + " | ".join(["---"] * len(cols)) + " |"
    lines = [header, sep]
    for _, row in df.iterrows():
        vals = [_format_value(c, row[c]) for c in cols]
        lines.append("| " + " | ".join(vals) + " |")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def save_table1(gap_df: pd.DataFrame, site_df: pd.DataFrame, out_csv: Path, out_md: Path, out_data: Path) -> None:
    merged = pd.merge(
        site_df,
        gap_df[
            [
                "dataset",
                "gap_size",
                "runtime_sec",
                "off_target_q95",
                "detected_rate",
                "total_windows",
                "off_target_n",
            ]
        ],
        on=["dataset", "gap_size"],
        how="left",
    )
    table = merged[
        [
            "dataset",
            "gap_size",
            "validated_score",
            "sn_ratio",
            "rank",
            "rank_percentile",
            "off_target_ge_validated",
            "off_target_q95",
            "off_target_n",
            "runtime_sec",
            "detected_rate",
            "hit_found",
        ]
    ].copy()
    table.sort_values(["dataset", "gap_size"], inplace=True)
    table.to_csv(out_csv, index=False)
    table.to_csv(out_data, index=False)
    write_markdown_table(table, out_md)


def main() -> int:
    args = build_parser().parse_args()
    out_dirs = ensure_out_dirs(args.out_dir)

    gap_df = pd.read_csv(args.gap_csv)
    site_df = pd.read_csv(args.site_csv)

    save_figure2(
        gap_df=gap_df,
        out_fig=out_dirs["figures"] / "Figure2_gap_sweep_tradeoff.png",
        out_data=out_dirs["figure_data"] / "Figure2_gap_sweep_tradeoff_data.csv",
    )
    save_figure3(
        site_df=site_df,
        out_fig=out_dirs["figures"] / "Figure3_offtarget_rank.png",
        out_data=out_dirs["figure_data"] / "Figure3_offtarget_rank_data.csv",
    )
    save_table1(
        gap_df=gap_df,
        site_df=site_df,
        out_csv=out_dirs["tables"] / "Table1_benchmark_summary.csv",
        out_md=out_dirs["tables"] / "Table1_benchmark_summary.md",
        out_data=out_dirs["figure_data"] / "Table1_benchmark_summary_data.csv",
    )

    print(f"[DONE] Figure 2 -> {out_dirs['figures'] / 'Figure2_gap_sweep_tradeoff.png'}")
    print(f"[DONE] Figure 3 -> {out_dirs['figures'] / 'Figure3_offtarget_rank.png'}")
    print(f"[DONE] Table 1 CSV -> {out_dirs['tables'] / 'Table1_benchmark_summary.csv'}")
    print(f"[DONE] Table 1 MD  -> {out_dirs['tables'] / 'Table1_benchmark_summary.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
