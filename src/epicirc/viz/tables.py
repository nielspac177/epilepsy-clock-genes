"""JAMA/NEJM-style result tables with inline unicode CI bars, plus a matplotlib forest plot.

No such helper existed in the repo (render/ draws brain maps; viz/report.py emits HTML only).
The markdown + unicode-bar functions are pure-stdlib and unit-testable; `forest_png` needs
matplotlib at call time. Effect estimates are shown in-cell as `est [lo, hi]` per JAMA convention,
with an optional monospaced unicode bar so a column of estimates reads like a forest plot in
plain GitHub markdown (no image needed).
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class Effect:
    label: str
    est: float
    lo: float
    hi: float
    p: float | None = None
    extra: dict | None = None


def fmt_ci(est: float, lo: float, hi: float, digits: int = 2) -> str:
    return f"{est:.{digits}f} [{lo:.{digits}f}, {hi:.{digits}f}]"


def ci_bar(est: float, lo: float, hi: float, vmin: float, vmax: float,
           width: int = 21, null_value: float = 0.0) -> str:
    """A fixed-width monospaced unicode CI bar spanning [vmin, vmax].

    '│' marks the null line, '─' the CI whiskers, '●' the point estimate. Values outside the
    axis are clamped. Meant to sit in a markdown code span so columns align.
    """
    if vmax <= vmin:
        vmax = vmin + 1e-9

    def pos(x: float) -> int:
        x = min(max(x, vmin), vmax)
        return int(round((x - vmin) / (vmax - vmin) * (width - 1)))

    cells = [" "] * width
    p0 = pos(null_value)
    if 0 <= p0 < width:
        cells[p0] = "│"  # │ null line
    a, b = pos(lo), pos(hi)
    for i in range(a, b + 1):
        if cells[i] == " ":
            cells[i] = "─"  # ─
    cells[pos(est)] = "●"   # ●
    return "".join(cells)


def markdown_table(headers: list[str], rows: list[list[str]], align: list[str] | None = None) -> str:
    align = align or ["left"] * len(headers)
    sep = {"left": ":--", "right": "--:", "center": ":-:"}
    out = ["| " + " | ".join(headers) + " |",
           "| " + " | ".join(sep.get(a, ":--") for a in align) + " |"]
    for r in rows:
        out.append("| " + " | ".join(str(c) for c in r) + " |")
    return "\n".join(out)


def forest_markdown(effects: list[Effect], value_name: str = "Estimate",
                    digits: int = 2, null_value: float = 0.0, bar_width: int = 21) -> str:
    """A JAMA-style table: label | est [CI] | inline forest bar | p."""
    vals = [v for e in effects for v in (e.lo, e.hi)]
    vmin, vmax = min(vals), max(vals)
    pad = 0.08 * (vmax - vmin or 1.0)
    vmin, vmax = min(vmin, null_value) - pad, max(vmax, null_value) + pad
    headers = [value_name, f"{value_name} [95% CI]", f"Forest (null={null_value:g})", "P"]
    rows = []
    for e in effects:
        bar = "`" + ci_bar(e.est, e.lo, e.hi, vmin, vmax, bar_width, null_value) + "`"
        pstr = "" if e.p is None else (f"{e.p:.1e}" if e.p < 1e-3 else f"{e.p:.3f}")
        rows.append([e.label, fmt_ci(e.est, e.lo, e.hi, digits), bar, pstr])
    return markdown_table(headers, rows, ["left", "right", "left", "right"])


def forest_png(effects: list[Effect], out_path: str, title: str = "",
               xlabel: str = "Effect [95% CI]", null_value: float = 0.0, figsize=None) -> str:
    """Publication forest plot (matplotlib). Point + 95% CI whiskers, null reference line."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(effects)
    fig, ax = plt.subplots(figsize=figsize or (6.4, 0.45 * n + 1.2))
    ys = list(range(n, 0, -1))
    for y, e in zip(ys, effects):
        ax.plot([e.lo, e.hi], [y, y], color="#333333", lw=1.6, zorder=2)
        ax.plot([e.est], [y], "o", color="#1f4e79", ms=6, zorder=3)
    ax.axvline(null_value, color="#b03030", lw=1.0, ls="--", zorder=1)
    ax.set_yticks(ys)
    ax.set_yticklabels([e.label for e in effects])
    ax.set_ylim(0.4, n + 0.6)
    ax.set_xlabel(xlabel)
    if title:
        ax.set_title(title, fontsize=11, loc="left")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out_path
