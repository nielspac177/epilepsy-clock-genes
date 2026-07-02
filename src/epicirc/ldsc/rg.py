"""Thin wrappers around LDSC (Bulik-Sullivan et al.) for h2, genetic correlation, and the
cross-trait intercept used as the shared-control ``rho`` in the heterogeneity test.

LDSC is legacy Python 2 / a pinned conda env; we shell out to it rather than reimplement. These
wrappers parse the ``.log`` output into structured results and are the only place the pipeline
touches the external binary.
"""
from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class RgResult:
    p1: str
    p2: str
    rg: float
    rg_se: float
    rg_p: float
    gcov_int: float          # cross-trait intercept ~ shared-sample confounding
    gcov_int_se: float
    h2_z_p1: float | None    # heritability z used for the power gate
    h2_z_p2: float | None


def run_rg(munged_a: str | Path, munged_b: str | Path, ld: str | Path,
           out: str | Path, ldsc: str = "ldsc.py") -> Path:
    """Run ``ldsc.py --rg`` for two munged sumstats; returns the log path to parse."""
    out = Path(out)
    subprocess.run(
        [ldsc, "--rg", f"{munged_a},{munged_b}",
         "--ref-ld-chr", str(ld), "--w-ld-chr", str(ld), "--out", str(out)],
        check=True,
    )
    return out.with_suffix(".log")


def parse_rg_log(log_path: str | Path) -> RgResult:
    """Parse the summary table LDSC writes at the end of an --rg run."""
    text = Path(log_path).read_text()
    # LDSC prints a whitespace table with a header line beginning 'p1'
    lines = text.splitlines()
    header_idx = next((i for i, ln in enumerate(lines) if ln.strip().startswith("p1")), None)
    if header_idx is None:
        raise ValueError(f"no rg summary table in {log_path}")
    cols = lines[header_idx].split()
    vals = lines[header_idx + 1].split()
    row = dict(zip(cols, vals))

    def _f(key: str) -> float:
        return float(row[key])

    return RgResult(
        p1=row["p1"], p2=row["p2"],
        rg=_f("rg"), rg_se=_f("se"), rg_p=_f("p"),
        gcov_int=_f("gcov_int"), gcov_int_se=_f("gcov_int_se"),
        h2_z_p1=_h2_z(text, which=1), h2_z_p2=_h2_z(text, which=2),
    )


def _h2_z(text: str, which: int) -> float | None:
    """Extract h2 Z-score from the per-trait 'Heritability of phenotype {which}' block."""
    m = re.search(rf"Heritability of phenotype {which}.*?Total Observed scale h2: "
                  r"([\-0-9.]+) \(([\-0-9.]+)\)", text, re.S)
    if not m:
        return None
    h2, se = float(m.group(1)), float(m.group(2))
    return h2 / se if se > 0 else None


def passes_power_gate(h2_z: float | None, threshold: float = 4.0) -> bool:
    """rg is only reported for traits whose h2 Z exceeds ``threshold`` (default 4)."""
    return h2_z is not None and h2_z >= threshold
