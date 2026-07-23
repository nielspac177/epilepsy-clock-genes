"""Assemble hg19-coordinate control gene sets for the Tier-1 validation.

Each set = a symbol list intersected with reference/genes_hg19.tsv (symbol->chr,start,end). The
intersection auto-drops symbols absent from the GENCODE-v19 universe (incl. QuickGO UniProt junk).

Sets:
  positive_control  - established epilepsy genes (MUST enrich -> pipeline valid)
  ion_channel       - voltage-gated Na/K/Ca channels (neuronal-excitability positive, GGE prior)
  housekeeping      - classic constitutive genes (MUST be null -> FPR calibration)
  go_circadian      - GO:0007623 circadian rhythm + descendants (unbiased-definition robustness)
"""
from pathlib import Path

REF = Path("reference/genes_hg19.tsv")
OUT = Path("config/gene_sets")

# --- curated symbol lists --------------------------------------------------------------------
POSITIVE = ["SCN1A", "SCN1B", "SCN2A", "SCN8A", "GABRA1", "GABRG2", "GABRB3", "CACNA1H", "KCNQ2",
            "KCNQ3", "STX1B", "SLC2A1", "DEPDC5", "CHRNA4", "GABRA5", "HCN1", "GRIN2A", "LGI1",
            "CACNA1A", "KCNA1"]

# Voltage-gated channels NOT already in POSITIVE, to keep the two positive sets distinct.
ION_CHANNEL = ["SCN3A", "SCN4A", "SCN5A", "SCN9A", "SCN10A", "SCN11A", "KCNA2", "KCNB1", "KCNC1",
               "KCND2", "KCNH2", "KCNJ10", "KCNMA1", "KCNQ1", "KCNQ5", "CACNA1B", "CACNA1C",
               "CACNA1E", "CACNA1G", "CACNA1I", "CACNB4", "CACNG2"]

# Eisenberg & Levanon (2013) constitutive housekeeping genes — a null set with no epilepsy/
# circadian prior. Deliberately spread across chromosomes.
HOUSEKEEPING = ["ACTB", "GAPDH", "B2M", "HPRT1", "PGK1", "PPIA", "RPL13A", "RPLP0", "TBP", "TFRC",
                "GUSB", "YWHAZ", "SDHA", "UBC", "PPIB", "RPS18", "EEF1A1", "RPL19", "RPL27",
                "TUBB", "VPS29", "CHMP2A", "EMC7", "REEP5", "PSMB2", "PSMB4", "RAB7A", "SNRPD3",
                "VCP", "GPI", "HMBS", "POLR2A", "NACA", "TPT1", "FAU", "RPS27A", "COX6B1",
                "C1orf43", "CANX", "HSPA8"]


def load_ref():
    m = {}
    with REF.open() as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            try:
                m[f[0]] = (f[1], int(f[2]), int(f[3]))
            except ValueError:
                continue
    return m


def go_symbols():
    raw = Path("reference/quickgo_circadian_raw.tsv")
    syms = set()
    with raw.open() as fh:
        next(fh)
        for ln in fh:
            f = ln.split("\t")
            if len(f) > 2:
                syms.add(f[2].strip())
    return sorted(syms)


def write_set(name, symbols, ref):
    present, missing = [], []
    for s in symbols:
        if s in ref:
            present.append(s)
        else:
            missing.append(s)
    path = OUT / f"{name}_hg19.tsv"
    with path.open("w") as fh:
        fh.write("symbol\tchr\tstart_hg19\tend_hg19\n")
        for s in present:
            c, a, b = ref[s]
            fh.write(f"{s}\t{c}\t{a}\t{b}\n")
    print(f"{name:16s} {len(present):3d} genes  (dropped {len(missing)}: {missing})")


if __name__ == "__main__":
    ref = load_ref()
    write_set("positive_control", POSITIVE, ref)
    write_set("ion_channel", ION_CHANNEL, ref)
    write_set("housekeeping", HOUSEKEEPING, ref)
    write_set("go_circadian", go_symbols(), ref)
