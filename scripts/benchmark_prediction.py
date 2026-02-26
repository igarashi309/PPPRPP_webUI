#!/usr/bin/env python3
"""
Benchmark pipeline for reproducible PPPRPP prediction evaluation.

Outputs:
  - data/benchmark/benchmark_summary.csv
  - data/benchmark/gap_sweep.csv
"""

from __future__ import annotations

import argparse
import json
import pickle
import re
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd
from Bio import SeqIO


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from ppr_analysis import core, io as ppr_io, motif_generation  # noqa: E402


@dataclass
class KnownSite:
    site_id: str
    start: int  # 0-based inclusive
    end: int    # 0-based exclusive
    strand: Optional[str]  # "Forward" / "Reverse" / None
    chrom: str = ""


def normalize_method(value: str) -> str:
    key = value.strip().lower()
    mapping = {
        "yan": "yan",
        "p-type": "yan",
        "p_type": "yan",
        "kobayashi": "kobayashi",
        "pls-type": "kobayashi",
        "pls_type": "kobayashi",
    }
    if key not in mapping:
        raise ValueError(f"Unknown method: {value}")
    return mapping[key]


def normalize_strand(value: object) -> Optional[str]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    token = str(value).strip()
    if not token:
        return None
    low = token.lower()
    if low in {".", "?", "na", "none", "both"}:
        return None
    if low in {"+", "forward", "fwd", "plus", "1"}:
        return "Forward"
    if low in {"-", "reverse", "rev", "minus", "-1"}:
        return "Reverse"
    return None


def parse_gap_values(max_gap: int, gap_values: Optional[str]) -> list[int]:
    if gap_values:
        parsed = []
        for token in gap_values.split(","):
            token = token.strip()
            if not token:
                continue
            parsed.append(int(token))
        if not parsed:
            raise ValueError("--gap-values is empty.")
        values = sorted(set(parsed))
        if values[0] < 0:
            raise ValueError("Gap values must be >= 0.")
        return values
    if max_gap < 0:
        raise ValueError("--max-gap must be >= 0.")
    return list(range(0, max_gap + 1))


def read_first_fasta_record(path: Path) -> tuple[str, str]:
    record = next(SeqIO.parse(str(path), "fasta"))
    seq_id = str(record.id)
    sequence = re.sub(r"\s+", "", str(record.seq)).upper()
    if not sequence:
        raise ValueError(f"Empty sequence in FASTA: {path}")
    return seq_id, sequence


def load_known_sites_table(path: Path) -> pd.DataFrame:
    """Best-effort table loader supporting BED/TSV/CSV with/without header."""
    # 1) Try automatic header parsing.
    try:
        return pd.read_csv(path, sep=None, engine="python", comment="#", dtype=str)
    except Exception:
        pass
    # 2) Fallback: tab without header.
    return pd.read_csv(path, sep="\t", engine="python", comment="#", header=None, dtype=str)


def resolve_coord_mode(path: Path, mode: str) -> str:
    if mode != "auto":
        return mode
    if path.suffix.lower() == ".bed":
        return "bed"
    return "one-based-inclusive"


def _to_numeric_int(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").astype("Int64")


def load_known_sites(path: Path, coord_mode: str) -> list[KnownSite]:
    resolved = resolve_coord_mode(path, coord_mode)
    rows: list[KnownSite] = []

    if resolved == "bed":
        df = pd.read_csv(path, sep="\t", comment="#", header=None, dtype=str)
        if df.empty or df.shape[1] < 3:
            raise ValueError("BED file must have at least 3 columns (chrom, start, end).")
        for idx, row in df.iterrows():
            start = pd.to_numeric(row.iloc[1], errors="coerce")
            end = pd.to_numeric(row.iloc[2], errors="coerce")
            if pd.isna(start) or pd.isna(end):
                continue
            site_id = str(row.iloc[3]) if df.shape[1] >= 4 and pd.notna(row.iloc[3]) else f"site_{idx+1:03d}"
            strand = normalize_strand(row.iloc[5]) if df.shape[1] >= 6 else None
            start_i, end_i = int(start), int(end)
            if start_i >= end_i:
                continue
            chrom = str(row.iloc[0]) if pd.notna(row.iloc[0]) else ""
            rows.append(KnownSite(site_id=site_id, start=start_i, end=end_i, strand=strand, chrom=chrom))
        return rows

    df = load_known_sites_table(path)
    if df.empty:
        raise ValueError(f"Known-site file is empty: {path}")

    lower = {str(c).strip().lower(): c for c in df.columns}
    has_named_cols = ("start" in lower and ("end" in lower or "stop" in lower))

    if has_named_cols:
        c_start = lower["start"]
        c_end = lower.get("end", lower.get("stop"))
        c_name = lower.get("name", lower.get("id"))
        c_strand = lower.get("strand")
        c_chrom = lower.get("chrom", lower.get("chr", lower.get("seqid", "")))
        working = pd.DataFrame(
            {
                "start": _to_numeric_int(df[c_start]),
                "end": _to_numeric_int(df[c_end]),
                "name": df[c_name].astype(str) if c_name else "",
                "strand": df[c_strand] if c_strand else "",
                "chrom": df[c_chrom].astype(str) if c_chrom else "",
            }
        )
    else:
        # Fallback to positional columns.
        raw = pd.read_csv(path, sep=None, engine="python", comment="#", header=None, dtype=str)
        if raw.shape[1] < 2:
            raise ValueError("Known-site TSV/CSV must have at least start/end columns.")

        # If first row looks like header tokens, drop it and remap.
        first_tokens = [str(x).strip().lower() for x in raw.iloc[0].tolist()]
        if "start" in first_tokens and ("end" in first_tokens or "stop" in first_tokens):
            i_start = first_tokens.index("start")
            i_end = first_tokens.index("end") if "end" in first_tokens else first_tokens.index("stop")
            i_name = first_tokens.index("name") if "name" in first_tokens else None
            i_strand = first_tokens.index("strand") if "strand" in first_tokens else None
            i_chrom = first_tokens.index("chrom") if "chrom" in first_tokens else (
                first_tokens.index("chr") if "chr" in first_tokens else None
            )
            data = raw.iloc[1:].reset_index(drop=True)
        else:
            data = raw
            i_chrom = None
            i_name = 2 if raw.shape[1] >= 3 else None
            i_strand = 3 if raw.shape[1] >= 4 else None
            # Detect BED-like positional format (chrom,start,end,...)
            c0_numeric = pd.to_numeric(raw.iloc[:, 0], errors="coerce").notna().mean() > 0.8
            c1_numeric = pd.to_numeric(raw.iloc[:, 1], errors="coerce").notna().mean() > 0.8
            c2_numeric = raw.shape[1] >= 3 and (pd.to_numeric(raw.iloc[:, 2], errors="coerce").notna().mean() > 0.8)
            if (not c0_numeric) and c1_numeric and c2_numeric:
                i_chrom = 0
                i_start = 1
                i_end = 2
                i_name = 3 if raw.shape[1] >= 4 else None
                i_strand = 5 if raw.shape[1] >= 6 else (4 if raw.shape[1] >= 5 else None)
            else:
                i_start = 0
                i_end = 1

        working = pd.DataFrame(
            {
                "start": _to_numeric_int(data.iloc[:, i_start]),
                "end": _to_numeric_int(data.iloc[:, i_end]),
                "name": data.iloc[:, i_name].astype(str) if i_name is not None and i_name < data.shape[1] else "",
                "strand": data.iloc[:, i_strand] if i_strand is not None and i_strand < data.shape[1] else "",
                "chrom": data.iloc[:, i_chrom].astype(str) if i_chrom is not None and i_chrom < data.shape[1] else "",
            }
        )

    for idx, row in working.iterrows():
        if pd.isna(row["start"]) or pd.isna(row["end"]):
            continue
        start_raw = int(row["start"])
        end_raw = int(row["end"])

        if resolved == "one-based-inclusive":
            start_i = start_raw - 1
            end_i = end_raw
        elif resolved == "zero-based-half-open":
            start_i = start_raw
            end_i = end_raw
        else:
            raise ValueError(f"Unsupported coordinate mode: {resolved}")

        if start_i < 0 or start_i >= end_i:
            continue

        name = str(row["name"]).strip() if str(row["name"]).strip() else f"site_{idx+1:03d}"
        rows.append(
            KnownSite(
                site_id=name,
                start=start_i,
                end=end_i,
                strand=normalize_strand(row["strand"]),
                chrom=str(row["chrom"]).strip(),
            )
        )

    if not rows:
        raise ValueError("No valid known-site rows were parsed.")
    return rows


def discover_motifs(protein_fasta: Path, tool: str, temp_root: Path):
    out_dir = temp_root / "motif_discovery"
    out_dir.mkdir(parents=True, exist_ok=True)

    if tool == "PPRFinder":
        pkl_path = motif_generation.run_pprfinder(str(protein_fasta), str(out_dir))
    else:
        pkl_path = motif_generation.run_ps_scan_api(str(protein_fasta), str(out_dir))
    if not pkl_path:
        raise RuntimeError("Motif discovery failed: no pickle generated.")

    motif_data = ppr_io.load_ppr_data_from_pickle(pkl_path)
    with open(pkl_path, "rb") as f:
        records = pickle.load(f)
    motif_data["motif_types"] = [r.annotations.get("PPR_TYPE", "Unknown") for r in records]

    return motif_data, Path(pkl_path).stem, Path(pkl_path)


def run_prediction_once(
    motif_data: dict,
    rna_seq: str,
    ppr_name: str,
    method: str,
    use_ps_score: bool,
    gap_size: int,
    min_front: int,
    min_back: int,
) -> tuple[pd.DataFrame, float]:
    mode = "gap0" if gap_size == 0 else "varied"
    ps_key = "ps" if use_ps_score else "no_ps"
    lookup_key = f"{method}_{'gap0' if gap_size == 0 else 'gap'}_{ps_key}"
    target_analysis = {"method": method, "use_ps": use_ps_score, "gap_mode": mode}

    t0 = time.perf_counter()
    fwd = core.analyze_ppr_data(
        motif_data,
        rna_seq,
        min_front=min_front,
        min_back=min_back,
        max_gap=gap_size,
        predicted_scores=[],
        ppr_name=ppr_name,
        detailed_output=False,
        target_analysis=target_analysis,
    )
    rev_seq = core.get_reverse_complement(rna_seq)
    rev = core.analyze_ppr_data(
        motif_data,
        rev_seq,
        min_front=min_front,
        min_back=min_back,
        max_gap=gap_size,
        predicted_scores=[],
        ppr_name=ppr_name,
        detailed_output=False,
        target_analysis=target_analysis,
    )
    runtime_sec = time.perf_counter() - t0

    seq_len = len(rna_seq)
    frames = []
    for strand, source in (("Forward", fwd), ("Reverse", rev)):
        data = source.get(lookup_key, {})
        if not data or not data.get("Position_Relative"):
            continue
        df = pd.DataFrame(data)
        if df.empty:
            continue
        df["Strand"] = strand
        if strand == "Forward":
            df["Genome_Position"] = df["Position_Relative"]
        else:
            df["Genome_Position"] = seq_len - 1 - pd.to_numeric(df["Position_Relative"], errors="coerce")

        df["PPR_Binding_Score"] = pd.to_numeric(df["PPR_Binding_Score"], errors="coerce")
        df["Genome_Position"] = pd.to_numeric(df["Genome_Position"], errors="coerce")
        df["Gap_Size"] = pd.to_numeric(df["Gap_Size"], errors="coerce").fillna(0).astype(int)
        df["Segment"] = df["Segment"].astype(str)
        df["Window_Length"] = df["Segment"].str.len().astype(int)
        df["Window_End"] = df["Genome_Position"].astype(int) + 1
        if strand == "Forward":
            df["Window_Start"] = df["Genome_Position"].astype(int)
            df["Window_End"] = df["Window_Start"] + df["Window_Length"]
        else:
            df["Window_Start"] = df["Window_End"] - df["Window_Length"]
        frames.append(df)

    if not frames:
        return pd.DataFrame(), runtime_sec

    merged = pd.concat(frames, ignore_index=True)
    merged["Analysis_Key"] = lookup_key
    merged["Requested_Gap_Size"] = int(gap_size)
    return merged, runtime_sec


def interval_overlap_mask(
    starts: np.ndarray,
    ends: np.ndarray,
    site_start: int,
    site_end: int,
) -> np.ndarray:
    return (starts < site_end) & (ends > site_start)


def evaluate_predictions(
    pred_df: pd.DataFrame,
    known_sites: Iterable[KnownSite],
) -> tuple[pd.DataFrame, dict]:
    known = list(known_sites)

    if pred_df.empty:
        empty_rows = []
        for site in known:
            empty_rows.append(
                {
                    "site_id": site.site_id,
                    "site_start": site.start,
                    "site_end": site.end,
                    "site_strand": site.strand if site.strand else "Both",
                    "hit_found": False,
                    "overlap_windows": 0,
                    "validated_score": np.nan,
                    "rank": np.nan,
                    "rank_percentile": np.nan,
                    "sn_ratio": np.nan,
                    "off_target_ge_validated": np.nan,
                }
            )
        summary = {
            "total_windows": 0,
            "off_target_n": 0,
            "off_target_q95": np.nan,
            "validated_sites": len(known),
            "validated_detected": 0,
            "detected_rate": 0.0 if known else np.nan,
            "median_validated_score": np.nan,
            "median_rank_percentile": np.nan,
            "median_sn_ratio": np.nan,
        }
        return pd.DataFrame(empty_rows), summary

    work = pred_df.copy()
    work["PPR_Binding_Score"] = pd.to_numeric(work["PPR_Binding_Score"], errors="coerce").fillna(0.0)
    work["Window_Start"] = pd.to_numeric(work["Window_Start"], errors="coerce").fillna(-1).astype(int)
    work["Window_End"] = pd.to_numeric(work["Window_End"], errors="coerce").fillna(-1).astype(int)

    starts = work["Window_Start"].to_numpy()
    ends = work["Window_End"].to_numpy()
    strands = work["Strand"].to_numpy()
    scores = work["PPR_Binding_Score"].to_numpy(dtype=float)
    n = len(work)

    order = np.argsort(-scores, kind="mergesort")
    ranks = np.empty(n, dtype=int)
    ranks[order] = np.arange(1, n + 1)

    validated_mask = np.zeros(n, dtype=bool)
    for site in known:
        mask = interval_overlap_mask(starts, ends, site.start, site.end)
        if site.strand in {"Forward", "Reverse"}:
            mask &= (strands == site.strand)
        validated_mask |= mask

    off_mask = ~validated_mask
    off_scores = scores[off_mask]
    off_q95 = float(np.quantile(off_scores, 0.95)) if off_scores.size > 0 else np.nan

    rows = []
    for site in known:
        mask = interval_overlap_mask(starts, ends, site.start, site.end)
        if site.strand in {"Forward", "Reverse"}:
            mask &= (strands == site.strand)
        idxs = np.where(mask)[0]
        if idxs.size == 0:
            rows.append(
                {
                    "site_id": site.site_id,
                    "site_start": site.start,
                    "site_end": site.end,
                    "site_strand": site.strand if site.strand else "Both",
                    "hit_found": False,
                    "overlap_windows": 0,
                    "validated_score": np.nan,
                    "rank": np.nan,
                    "rank_percentile": np.nan,
                    "sn_ratio": np.nan,
                    "off_target_ge_validated": np.nan,
                }
            )
            continue

        local = idxs[np.argmax(scores[idxs])]
        validated_score = float(scores[local])
        rank = int(ranks[local])
        rank_percentile = 100.0 * (n - rank + 1) / float(n)
        sn_ratio = (validated_score / off_q95) if (not np.isnan(off_q95) and off_q95 > 0) else np.nan
        off_ge = int(np.sum(off_scores >= validated_score)) if off_scores.size > 0 else 0

        rows.append(
            {
                "site_id": site.site_id,
                "site_start": site.start,
                "site_end": site.end,
                "site_strand": site.strand if site.strand else "Both",
                "hit_found": True,
                "overlap_windows": int(idxs.size),
                "validated_score": validated_score,
                "rank": rank,
                "rank_percentile": rank_percentile,
                "sn_ratio": sn_ratio,
                "off_target_ge_validated": off_ge,
            }
        )

    per_site = pd.DataFrame(rows)
    summary = {
        "total_windows": int(n),
        "off_target_n": int(off_scores.size),
        "off_target_q95": off_q95,
        "validated_sites": int(len(known)),
        "validated_detected": int(per_site["hit_found"].sum()) if not per_site.empty else 0,
        "detected_rate": (
            float(per_site["hit_found"].mean()) if not per_site.empty else np.nan
        ),
        "median_validated_score": (
            float(per_site["validated_score"].median(skipna=True))
            if not per_site.empty
            else np.nan
        ),
        "median_rank_percentile": (
            float(per_site["rank_percentile"].median(skipna=True))
            if not per_site.empty
            else np.nan
        ),
        "median_sn_ratio": (
            float(per_site["sn_ratio"].median(skipna=True))
            if not per_site.empty
            else np.nan
        ),
    }
    return per_site, summary


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Run PPPRPP benchmark prediction pipeline.")
    p.add_argument("--protein-fasta", required=True, type=Path, help="Protein FASTA for motif discovery.")
    p.add_argument("--rna-fasta", required=True, type=Path, help="Target RNA FASTA.")
    p.add_argument("--known-sites", required=True, type=Path, help="Known-site file (BED/TSV/CSV).")
    p.add_argument("--output-dir", type=Path, default=REPO_ROOT / "data" / "benchmark", help="Output directory.")
    p.add_argument("--tool", choices=["PPRFinder", "PS_Scan"], default="PPRFinder", help="Motif discovery tool.")
    p.add_argument(
        "--method",
        choices=["yan", "kobayashi", "P-type", "PLS-type"],
        default="yan",
        help="Scoring method for target prediction.",
    )
    ps_group = p.add_mutually_exclusive_group()
    ps_group.add_argument("--use-ps-score", dest="use_ps_score", action="store_true", help="Use motif confidence scores.")
    ps_group.add_argument("--no-use-ps-score", dest="use_ps_score", action="store_false", help="Disable motif confidence scores.")
    p.set_defaults(use_ps_score=True)
    p.add_argument("--max-gap", type=int, default=3, help="Run gap sweep from 0..max-gap (inclusive).")
    p.add_argument("--gap-values", type=str, default=None, help="Comma-separated gap sizes (overrides --max-gap).")
    p.add_argument("--min-front", type=int, default=1, help="Minimum motifs before gap insertion.")
    p.add_argument("--min-back", type=int, default=1, help="Minimum motifs after gap insertion.")
    p.add_argument(
        "--known-sites-coordinate",
        choices=["auto", "bed", "zero-based-half-open", "one-based-inclusive"],
        default="auto",
        help="Coordinate mode for known-site file.",
    )
    p.add_argument("--save-predictions", action="store_true", help="Save per-gap prediction CSVs.")
    return p


def main() -> int:
    args = build_parser().parse_args()
    method = normalize_method(args.method)
    gap_values = parse_gap_values(args.max_gap, args.gap_values)

    for fp in (args.protein_fasta, args.rna_fasta, args.known_sites):
        if not fp.exists():
            raise FileNotFoundError(f"Input not found: {fp}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    print(f"[INFO] Output dir: {args.output_dir}")
    print(f"[INFO] Gap sweep: {gap_values}")

    known_sites = load_known_sites(args.known_sites, args.known_sites_coordinate)
    print(f"[INFO] Parsed known sites: {len(known_sites)}")

    rna_id, rna_seq = read_first_fasta_record(args.rna_fasta)
    print(f"[INFO] RNA loaded: id={rna_id}, length={len(rna_seq)}")

    with tempfile.TemporaryDirectory(prefix="ppprpp_benchmark_") as tmpdir:
        temp_root = Path(tmpdir)
        motif_data, ppr_name, pkl_path = discover_motifs(args.protein_fasta, args.tool, temp_root)
        print(f"[INFO] Motif discovery done: ppr={ppr_name}, motifs={len(motif_data.get('motif_labels', []))}")

        per_site_all = []
        gap_rows = []
        for gap_size in gap_values:
            print(f"[INFO] Running gap={gap_size} ...")
            pred_df, runtime_sec = run_prediction_once(
                motif_data=motif_data,
                rna_seq=rna_seq,
                ppr_name=ppr_name,
                method=method,
                use_ps_score=args.use_ps_score,
                gap_size=gap_size,
                min_front=args.min_front,
                min_back=args.min_back,
            )
            site_df, summary = evaluate_predictions(pred_df, known_sites)
            site_df["gap_size"] = int(gap_size)
            site_df["method"] = method
            site_df["use_ps_score"] = bool(args.use_ps_score)
            site_df["discovery_tool"] = args.tool
            site_df["ppr_name"] = ppr_name
            site_df["rna_id"] = rna_id
            per_site_all.append(site_df)

            gap_row = {
                "gap_size": int(gap_size),
                "runtime_sec": float(runtime_sec),
                "method": method,
                "use_ps_score": bool(args.use_ps_score),
                "discovery_tool": args.tool,
                "ppr_name": ppr_name,
                "rna_id": rna_id,
                **summary,
            }
            gap_rows.append(gap_row)

            if args.save_predictions:
                pred_path = args.output_dir / f"predictions_gap{gap_size}.csv"
                pred_df.to_csv(pred_path, index=False)
                print(f"[INFO] Saved predictions: {pred_path}")

        summary_df = pd.concat(per_site_all, ignore_index=True) if per_site_all else pd.DataFrame()
        sweep_df = pd.DataFrame(gap_rows).sort_values(by=["gap_size"], kind="mergesort")

        summary_path = args.output_dir / "benchmark_summary.csv"
        sweep_path = args.output_dir / "gap_sweep.csv"
        summary_df.sort_values(by=["gap_size", "site_id"], kind="mergesort").to_csv(summary_path, index=False)
        sweep_df.to_csv(sweep_path, index=False)

        meta = {
            "protein_fasta": str(args.protein_fasta),
            "rna_fasta": str(args.rna_fasta),
            "known_sites": str(args.known_sites),
            "tool": args.tool,
            "method": method,
            "use_ps_score": bool(args.use_ps_score),
            "gap_values": gap_values,
            "min_front": int(args.min_front),
            "min_back": int(args.min_back),
            "known_sites_coordinate": args.known_sites_coordinate,
            "ppr_pickle": str(pkl_path),
        }
        (args.output_dir / "benchmark_run_meta.json").write_text(
            json.dumps(meta, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )

        print(f"[DONE] benchmark_summary.csv -> {summary_path}")
        print(f"[DONE] gap_sweep.csv -> {sweep_path}")
        print(f"[DONE] benchmark_run_meta.json -> {args.output_dir / 'benchmark_run_meta.json'}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
