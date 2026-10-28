"""Shared code for the kmerseek Botryllus-vs-human immune and HLA notebooks (232, 233).

- load_immune_genes: human genes under GO:0002376 "immune system process" (QuickGO, cached)
- with_gene_columns: human gene symbol and protein length parsed from a GENCODE header
- HlaParts: which part of an HLA protein (UniProt features) a match covers
- KNOWN_IMMUNE_FAMILY_PATTERN: gene-symbol pattern for immune families a biologist names first
- PART_COLOUR, PART_MARKER: one colour and marker per HLA part, the same in every figure
- same_residue_vs_chance, same_residue_vs_chance_in_class: how many aligned positions carry the
  same amino acid, against the number two unrelated sequences would share
- match_line: a '|' under each position where the two sequences carry the same amino acid
"""

from __future__ import annotations

import importlib.util
import json
import urllib.parse
import urllib.request
from pathlib import Path

import polars as pl

import pubfig as pf

DATA_DIR = Path(__file__).with_name("data")
GO_IMMUNE_TSV = DATA_DIR / "go0002376_immune_system_process_human_quickgo.tsv"
UNIPROT_HLA_DIR = DATA_DIR / "uniprot_hla"
QUICKGO_URL = (
    "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?"
    "goId=GO:0002376&goUsage=descendants&goUsageRelationships=is_a,part_of,occurs_in"
    "&taxonId=9606&taxonUsage=exact&geneProductType=protein&selectedFields=symbol,goId"
    "&downloadLimit=2000000"
)

CLASS_I = {"HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G"}
NOT_MAPPED = "not mapped (GENCODE isoform differs from UniProt; residues not found)"

PART_COLOUR = {
    "peptide-binding groove": pf.OKABE_ITO["vermillion"],
    "Ig-like domain": pf.OKABE_ITO["blue"],
    "transmembrane helix": pf.OKABE_ITO["orange"],
    "signal peptide": pf.OKABE_ITO["reddish_purple"],
    "cytoplasmic tail": pf.OKABE_ITO["bluish_green"],
    "other": pf.OKABE_ITO["sky_blue"],
    NOT_MAPPED: pf.GREY,
}
PART_MARKER = {part: ("x" if part == NOT_MAPPED else "o") for part in PART_COLOUR}

# HLA, antigen processing, Toll-like and NOD-like receptors, complement, cytokines and their
# receptors, interferons, TNF family, NF-kB, IRF, TRAF, JAK/STAT, lymphocyte kinases, CD
# antigens, KIR/LILR, Fc receptors, immunoglobulin and T-cell receptor genes, RAG and AIRE.
KNOWN_IMMUNE_FAMILY_PATTERN = (
    r"^(HLA-|B2M$|TAP[12]$|TAPBP|PSMB(8|9|10)$|TLR\d|NLR|NOD[12]$|C[1-9][A-Z]?\d?$|CF[BDHI]$|MASP|"
    r"IL\d|IL\d+R|CCL\d|CXCL\d|CCR\d|CXCR\d|IFN[ABGLKWE]|TNF$|TNFSF|TNFRSF|LTA$|LTB$|NFKB|REL[AB]?$|IRF\d|TRAF\d|MYD88|IRAK\d|"
    r"JAK\d|STAT\d|LCK$|ZAP70|SYK$|BTK$|CD\d|KIR\d|LILR|FCGR|FCER|IG[HKL][VCDJ]|TR[ABGD][VCJ]|RAG[12]$|AIRE$|CARD11$)"
)


def load_immune_genes() -> list[str]:
    """Human gene symbols annotated to GO:0002376 or a term below it, downloaded once."""
    if not GO_IMMUNE_TSV.exists():
        req = urllib.request.Request(QUICKGO_URL, headers={"Accept": "text/tsv"})
        GO_IMMUNE_TSV.write_bytes(urllib.request.urlopen(req, timeout=300).read())
    go_immune = pl.read_csv(GO_IMMUNE_TSV, separator="\t")
    assert go_immune.columns == ["SYMBOL", "GO TERM"], go_immune.columns  # not an outage page
    return sorted(set(go_immune["SYMBOL"]))


def with_gene_columns(frame: pl.LazyFrame | pl.DataFrame) -> pl.LazyFrame | pl.DataFrame:
    """Add human_gene and human_protein_length (fields 7 and 8 of the GENCODE header) and
    strip the trailing space kmerseek leaves on some query names."""
    fields = pl.col("target_name").str.split("|")
    return frame.with_columns(
        pl.col("query_name").str.strip_chars(),
        fields.list.get(6).alias("human_gene"),
        fields.list.get(7).cast(pl.Int64).alias("human_protein_length"),
    )


def uniprot_hla(gene: str) -> dict:
    """Reviewed human UniProt entry for one HLA gene, cached as JSON."""
    UNIPROT_HLA_DIR.mkdir(parents=True, exist_ok=True)
    path = UNIPROT_HLA_DIR / f"{gene}.json"
    # files cached before 27 Oct 2026 lack the sequence; fetch those again
    if not path.exists() or "value" not in json.loads(path.read_text())["results"][0]["sequence"]:
        query = f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true"
        url = ("https://rest.uniprot.org/uniprotkb/search?format=json&fields=accession,gene_primary,"
               "length,sequence,ft_transmem,ft_domain,ft_signal,ft_region,ft_topo_dom&query=" + urllib.parse.quote(query))
        path.write_bytes(urllib.request.urlopen(url, timeout=120).read())
    results = json.loads(path.read_text())["results"]
    assert len(results) == 1, (gene, len(results))
    return results[0]


def hla_part_by_residue(entry: dict, gene: str) -> list[str]:
    """Name of the HLA protein part at each residue (index 0 = residue 1)."""
    length = entry["sequence"]["length"]
    part = ["other"] * length
    groove = {"Alpha-1", "Beta-1"} | ({"Alpha-2"} if gene in CLASS_I else set())
    # later layers win: groove, then Ig-like, then cytoplasmic, signal and transmembrane
    layers = [
        ("peptide-binding groove", lambda f: f["type"] == "Region" and f.get("description") in groove),
        ("Ig-like domain", lambda f: f["type"] == "Domain" and f.get("description", "").startswith("Ig-like")),
        ("cytoplasmic tail", lambda f: f["type"] == "Topological domain" and f.get("description") == "Cytoplasmic"),
        ("signal peptide", lambda f: f["type"] == "Signal"),
        ("transmembrane helix", lambda f: f["type"] == "Transmembrane"),
    ]
    for name, keep in layers:
        for f in entry["features"]:
            if keep(f):
                start, end = f["location"]["start"]["value"], f["location"]["end"]["value"]
                for i in range(start - 1, end):
                    part[i] = name
    return part


class HlaParts:
    """Places kmerseek matches on the parts of the HLA proteins they cover."""

    def __init__(self, gencode_length_by_gene: dict[str, int]):
        self.part_by_residue: dict[str, list[str]] = {}
        self.uniprot_sequence: dict[str, str] = {}
        rows = []
        for gene, gencode_length in sorted(gencode_length_by_gene.items()):
            entry = uniprot_hla(gene)
            self.part_by_residue[gene] = hla_part_by_residue(entry, gene)
            self.uniprot_sequence[gene] = entry["sequence"]["value"]
            rows.append({"human_gene": gene, "uniprot": entry["primaryAccession"],
                         "uniprot_length": entry["sequence"]["length"], "gencode_length": gencode_length})
        self.length_check = pl.DataFrame(rows).with_columns(
            (pl.col("uniprot_length") == pl.col("gencode_length")).alias("same_length"))
        self.same_length = dict(zip(self.length_check["human_gene"], self.length_check["same_length"]))

    def main_part(self, gene: str, target_start: int, target_end: int, target_subseq: str) -> str:
        """HLA part covering most residues of a match (target_start 0-based, target_end exclusive).

        When the GENCODE protein differs in length from UniProt's (another isoform), the
        match is placed where its human residues occur in the UniProt sequence; if they do
        not occur there exactly, the match is not placed.
        """
        if not self.same_length[gene]:
            target_start = self.uniprot_sequence[gene].find(target_subseq)
            if target_start < 0:
                return NOT_MAPPED
            target_end = target_start + len(target_subseq)
        parts = self.part_by_residue[gene][target_start:target_end]
        # a tie goes to the part that comes first along the protein, so reruns agree
        return max(dict.fromkeys(parts), key=parts.count)

    def add_columns(self, matches: pl.DataFrame) -> pl.DataFrame:
        """Add hla_part_matched and hla_class to HLA matches."""
        return matches.with_columns(
            pl.struct("human_gene", "target_start", "target_end", "target_subseq")
            .map_elements(lambda r: self.main_part(r["human_gene"], r["target_start"], r["target_end"],
                                                   r["target_subseq"]),
                          return_dtype=pl.Utf8)
            .alias("hla_part_matched"),
            pl.when(pl.col("human_gene").is_in(list(CLASS_I))).then(pl.lit("class I")).otherwise(pl.lit("class II"))
            .alias("hla_class"),
        )


# UniProtKB/Swiss-Prot release 2026_03, "6.1 Composition in percent for the complete database"
# (https://web.expasy.org/docs/relnotes/relstat.html), the numbers notebook 274 of
# seanome/2024-kmerseek-analysis and 235_plan_arms.py use.
SWISSPROT_PERCENT = {
    "A": 8.25, "Q": 3.93, "L": 9.64, "S": 6.66,
    "R": 5.52, "E": 6.71, "K": 5.79, "T": 5.36,
    "N": 4.06, "G": 7.07, "M": 2.41, "W": 1.10,
    "D": 5.46, "H": 2.27, "F": 3.86, "Y": 2.92,
    "C": 1.38, "I": 5.90, "P": 4.75, "V": 6.85,
}  # fmt: skip
SWISSPROT_SHARE = {aa: pct / sum(SWISSPROT_PERCENT.values()) for aa, pct in SWISSPROT_PERCENT.items()}
# chance that two residues drawn from Swiss-Prot are the same amino acid
SWISSPROT_SAME_RESIDUE = sum(q * q for q in SWISSPROT_SHARE.values())


def alphabet_classes(alphabet: str) -> dict[str, str]:
    """Letter class of each amino acid in a kmerseek alphabet (partitions in 233_plan_arms.py)."""
    spec = importlib.util.spec_from_file_location("plan_233", Path(__file__).with_name("233_plan_arms.py"))
    plan = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(plan)
    return {aa: cls for cls in plan.CLUSTERS[alphabet] for aa in cls}


def match_line(query: str, target: str) -> str:
    """'|' under each position where the two ungapped sequences carry the same amino acid."""
    return "".join("|" if a == b else " " for a, b in zip(query, target))


def _at_least(k: int, probabilities: list[float]) -> float:
    """Chance of k or more successes in independent trials with these probabilities."""
    dist = [1.0]
    for p in probabilities:
        nxt = [0.0] * (len(dist) + 1)
        for i, v in enumerate(dist):
            nxt[i] += v * (1 - p)
            nxt[i + 1] += v * p
        dist = nxt
    return sum(dist[k:])


def same_residue_vs_chance(query: str, target: str) -> dict:
    """Positions with the same amino acid, against unrelated sequences with Swiss-Prot frequencies.

    Returns n_positions, n_same, n_expected (n_positions x the chance two Swiss-Prot residues
    are the same amino acid) and p_at_least (chance of n_same or more).
    """
    n = min(len(query), len(target))
    n_same = sum(a == b for a, b in zip(query, target))
    probs = [SWISSPROT_SAME_RESIDUE] * n
    return {"n_positions": n, "n_same": n_same, "n_expected": sum(probs), "p_at_least": _at_least(n_same, probs)}


def same_residue_vs_chance_in_class(query: str, target: str, alphabet: str) -> dict:
    """As same_residue_vs_chance, but given that each pair is in the same letter class.

    A kmerseek k-mer match puts both residues of a position in the same class of the alphabet,
    so unrelated sequences share more amino acids there than the plain chance level says (in
    gbmr7, G, P, S and T are classes of their own and must be identical). At a position in
    class c, the chance the target carries the query's amino acid a is share(a) / share(c),
    with Swiss-Prot shares. Only positions where both residues are in one class count.
    """
    classes = alphabet_classes(alphabet)
    pairs = [(a, b) for a, b in zip(query, target) if a in classes and b in classes and classes[a] == classes[b]]
    probs = [SWISSPROT_SHARE[a] / sum(SWISSPROT_SHARE[x] for x in classes[a]) for a, _ in pairs]
    n_same = sum(a == b for a, b in pairs)
    return {"n_positions": len(pairs), "n_same": n_same, "n_expected": sum(probs), "p_at_least": _at_least(n_same, probs)}
