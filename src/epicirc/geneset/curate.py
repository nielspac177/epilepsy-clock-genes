"""Load, hash, validate, and control-match pre-registered gene sets.

Reproducibility contract:
* every gene set is hashed (SHA-256 over sorted symbols) -> results/geneset_lock.json
* symbols are validated against the MAGMA gene-location file when present
* negative controls are drawn size- and gene-length-matched with a fixed seed
"""
from __future__ import annotations

import csv
import hashlib
import random
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class GeneSet:
    name: str
    symbols: tuple[str, ...]

    @property
    def sha256(self) -> str:
        joined = "\n".join(sorted(s.upper() for s in self.symbols))
        return hashlib.sha256(joined.encode()).hexdigest()

    def __len__(self) -> int:
        return len(self.symbols)


def load_geneset(path: str | Path) -> GeneSet:
    """Read a gene-set TSV (comment lines start with '#', header row has 'symbol')."""
    path = Path(path)
    symbols: list[str] = []
    with path.open() as fh:
        rows = [ln for ln in fh if not ln.startswith("#")]
    reader = csv.DictReader(rows, delimiter="\t")
    if reader.fieldnames is None or "symbol" not in reader.fieldnames:
        raise ValueError(f"{path} missing 'symbol' column")
    for row in reader:
        sym = (row.get("symbol") or "").strip()
        if sym:
            symbols.append(sym.upper())
    if not symbols:
        raise ValueError(f"{path} contained no genes")
    if len(symbols) != len(set(symbols)):
        raise ValueError(f"{path} has duplicate symbols")
    return GeneSet(name=path.stem, symbols=tuple(symbols))


def validate_against_annotation(gs: GeneSet, gene_loc: str | Path) -> list[str]:
    """Return symbols in ``gs`` absent from the MAGMA .gene.loc annotation.

    An empty list means every gene maps. The pipeline treats a non-empty list as a hard error
    (a typo'd symbol silently dropped from a competitive test biases the result).
    """
    known: set[str] = set()
    with Path(gene_loc).open() as fh:
        for ln in fh:
            parts = ln.split()
            # MAGMA .gene.loc: ENTREZ CHR START STOP STRAND SYMBOL
            if len(parts) >= 6:
                known.add(parts[5].upper())
    return [s for s in gs.symbols if s not in known]


def overlap(a: GeneSet, b: GeneSet) -> tuple[str, ...]:
    return tuple(sorted(set(a.symbols) & set(b.symbols)))


def matched_negative_controls(
    gs: GeneSet,
    universe: dict[str, int],
    n_sets: int,
    seed: int,
    length_tolerance: float = 0.20,
) -> list[tuple[str, ...]]:
    """Draw ``n_sets`` random gene sets matched on size and (approx) gene length.

    ``universe`` maps symbol -> gene length (bp). Each control gene is sampled to have a length
    within ``length_tolerance`` of a corresponding target gene, giving a size- and
    length-matched null. Deterministic given ``seed``.
    """
    rng = random.Random(seed)
    pool = {s: L for s, L in universe.items() if s not in set(gs.symbols)}
    target_lengths = [universe[s] for s in gs.symbols if s in universe]
    if not target_lengths:
        raise ValueError("no target genes found in universe; check annotation join")

    controls: list[tuple[str, ...]] = []
    for _ in range(n_sets):
        chosen: list[str] = []
        used: set[str] = set()
        for L in target_lengths:
            lo, hi = L * (1 - length_tolerance), L * (1 + length_tolerance)
            candidates = [s for s, gl in pool.items() if lo <= gl <= hi and s not in used]
            if not candidates:  # widen if the band is empty
                candidates = [s for s in pool if s not in used]
            pick = rng.choice(candidates)
            chosen.append(pick)
            used.add(pick)
        controls.append(tuple(chosen))
    return controls
