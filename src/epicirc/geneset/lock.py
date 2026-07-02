"""Freeze the pre-registered gene sets: write a hash + summary lock file before outcome analysis.

Run once at data lock. CI compares the live hashes to this file; any drift fails the build,
enforcing the pre-registration contract (docs/adr/0003).

    python -m epicirc.geneset.lock --config config/gene_sets --out results/geneset_lock.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from epicirc.geneset.curate import load_geneset, overlap


def build_lock(config_dir: str | Path) -> dict:
    config_dir = Path(config_dir)
    sets = {}
    for tsv in sorted(config_dir.glob("*.tsv")):
        gs = load_geneset(tsv)
        sets[gs.name] = gs

    core = sets.get("circadian_core")
    lock: dict = {"gene_sets": {}, "overlaps": {}}
    for name, gs in sets.items():
        lock["gene_sets"][name] = {
            "n_genes": len(gs),
            "sha256": gs.sha256,
            "symbols": sorted(gs.symbols),
        }
    # audit overlap of the primary set against every other set (transparency, docs/adr/0005)
    if core is not None:
        for name, gs in sets.items():
            if name == core.name:
                continue
            shared = overlap(core, gs)
            if shared:
                lock["overlaps"][f"circadian_core∩{name}"] = list(shared)
    return lock


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="config/gene_sets")
    ap.add_argument("--out", default="results/geneset_lock.json")
    args = ap.parse_args(argv)

    lock = build_lock(args.config)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(lock, indent=2, sort_keys=True) + "\n")
    n = len(lock["gene_sets"])
    print(f"Locked {n} gene sets -> {out}")
    for name, info in lock["gene_sets"].items():
        print(f"  {name}: {info['n_genes']} genes  {info['sha256'][:12]}")
    if lock["overlaps"]:
        print("Overlap audit (circadian_core vs others):")
        for k, v in lock["overlaps"].items():
            print(f"  {k}: {', '.join(v)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
