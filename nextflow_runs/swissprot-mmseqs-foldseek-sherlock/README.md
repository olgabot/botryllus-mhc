# Botryllus proteome vs Swiss-Prot: MMseqs2 and Foldseek, looking for MHC

Question: does any *Botryllus schlosseri* protein match a vertebrate MHC class I or
class II protein (HLA-A/B/C/E/F/G, HLA-D*, CD1, MR1, FcRn, ZAG, MICA/B, HFE, and
the same proteins in other vertebrates)?

All 72,617 Botryllus proteins (`~/data/botryllus/Bs_proteins.fa`, `*` stop
characters removed) were searched against UniProtKB/Swiss-Prot release 2026_03
(2 September 2026) two ways:

1. **MMseqs2** (sequence against sequence): `easy-search -s 7.5 -e 1e-3 --max-seqs 2000`
   against the 575,748 Swiss-Prot sequences.
2. **Foldseek** (predicted structure against structure). Botryllus has no structures,
   so ProstT5, a protein language model, predicts each protein's 3Di letters from its
   sequence. 3Di is Foldseek's 20-letter alphabet that describes how each residue sits
   against its nearest neighbour in space. These were searched with
   `search -e 1e-3 --max-seqs 2000` against the AlphaFold models of Swiss-Prot
   (`foldseek databases Alphafold/Swiss-Prot`, 590,183 models). The queries have no
   coordinates, so Foldseek scores 3Di and amino-acid letters together and reports no
   TM-score or lDDT.

`--max-seqs 2000` (MMseqs2 default 300, Foldseek 1000) is there because Botryllus
Ig-domain proteins match hundreds of Swiss-Prot antibodies and receptors, and a weaker
MHC match could otherwise be dropped before it is aligned.

## Which Swiss-Prot entries count as MHC

`reference/swissprot_mhc_c1set_pfam.tsv` is a UniProt query (release 2026_03) for
reviewed entries with any of these Pfam families, plus PROCR:

| Pfam | Name | Counted as |
|---|---|---|
| PF00129 | MHC_I | MHC class I or class I-like |
| PF06623 | MHC_I_C | MHC class I or class I-like |
| PF16497 | MHC_I_3 (CD1, EPCR) | MHC class I or class I-like |
| PF00993 | MHC_II_alpha | MHC class II |
| PF00969 | MHC_II_beta | MHC class II |
| PF07654 | C1-set (Ig constant domain) | B2M if the gene is B2M, else "other C1-set" |

That gives 187 MHC entries (MHC class I, CD1, MR1, FcRn, ZAG, MICA/B, HFE, ULBP,
EPCR, class II alpha and beta, from many vertebrates) and 197 B2M and other C1-set
entries. The C1-set entries are context: an MHC chain has one Ig-like C1 domain
(class I alpha3, class II alpha2 or beta2) next to its peptide-binding part, and that
Ig-like domain is shared with antibodies, T-cell receptors and many other proteins.

For every MHC hit the tables give `aligned_aa_in_peptide_binding_region`: how many
aligned residues of the MHC protein fall between the end of its signal peptide and the
start of its first Ig-like domain (UniProt feature coordinates). That stretch is the
peptide-binding part (class I alpha1 and alpha2, class II alpha1 or beta1). A hit
that covers only the Ig-like domain has 0 or a few residues there.

## Results

No Botryllus protein has an MHC class I or class II protein as its best Swiss-Prot
match, by sequence or by predicted structure. Every MHC hit is weak (e-value between
9.6e-6 and 1e-3) and ranks below other Ig-domain proteins in the same protein's hit list.

![Rank of each Botryllus protein's best MHC, B2M and C1-set hit among all its
Swiss-Prot hits](results/rank_of_best_hit_in_class.png)

**Figure.** One dot per Botryllus protein and target class. x: where the protein's
best hit to that class ranks among all its Swiss-Prot hits, by bit score (1 = its
best hit of all; the black line). A Botryllus MHC homolog would put a vermillion or
blue dot on the line. (a) MMseqs2, sequence search. (b) Foldseek with 3Di letters
predicted by ProstT5. Grey dots on the line are Botryllus proteins whose best hit is
a Swiss-Prot protein with a C1-set Ig domain, shown for comparison.

| | MMseqs2 | Foldseek (ProstT5 3Di) |
|---|---|---|
| Botryllus proteins with any hit to MHC class II | 1 | 69 |
| ... to MHC class I or class I-like (incl. CD1, MR1, FcRn, ZAG, MICA/B, HFE) | 0 | 55 |
| ... to B2M | 0 | 15 |
| ... to other C1-set Ig-domain proteins | 8 | 676 |
| Lowest MHC e-value | 7.1e-4 | 9.6e-6 |
| Proteins whose best hit is an MHC protein | 0 | 0 |
| Highest rank of an MHC hit | 28 | 17 |
| Median rank of each protein's best MHC hit | 28 | 1,346 |

All MHC hits are to vertebrate entries. The only non-vertebrate Swiss-Prot entries
in the MHC families are two cytomegalovirus UL18 proteins, a viral copy of MHC
class I, and no Botryllus protein hits them.

**MMseqs2.** One protein, g62973.t1, matches rhesus HLA-DR alpha (DRA_MACMU,
e = 7.1e-4, 46 bits). Its best hit is fibronectin (66 bits), and the DRA hit ranks
28th of its 30 hits. The alignment covers DRA residues 107-208, which is DRA's
Ig-like C1 domain (112-204), and only 5 aligned residues fall in DRA's
peptide-binding alpha1 region (26-111). 33 of 104 aligned positions are identical:

```text
MMseqs2: g62973.t1 115-211 vs DRA_MACMU 107-208, 46 bits, 33 of 104 identical
  g62973.t1   115 IETIPVNIVTLPMAPFRIGHRVLLNCLYQTFPDGTLSVVWKKNGTPISQTRPSERVTIRN
                  |   |     |   |   |    | |    |      | | ||| |   |  || |
  Q30631      107 ITNVPPEVTVLTNSPVELGEPNVLICFIDKFSPPVVKVTWLKNGKPVT-TGVSETVFLPR

  g62973.t1   175 ERQLVIK----EFQDSDVGTYTCSVEY---DLPTALKQGEASIP
                  |  |  |     |  |    | | ||    | |  ||  |   |
  Q30631      166 EDHLFRKFHYLPFLPSTEDIYDCKVEHWCLDAP-LLKHWEFDAP
```

**Foldseek.** 82 Botryllus proteins match at least one MHC protein. The highest-ranked
MHC hit is g22762.t1 vs chimpanzee HLA-DOB (e = 1.5e-5, rank 17 of 365); its best
hit is a titin homolog and the alignment lies in the DOB Ig-like domain only. Only
three proteins have an MHC hit in their top 100.

Some Foldseek alignments do reach into the peptide-binding region: 30 proteins align
40 or more residue pairs there (no gap on either side), up to 139 pairs
(g19880.t1 vs MICB). These are still poor matches. The best-ranked of them,
g64756.t1 vs dog MHC class I (HA19_CANLF), is 189th of the protein's 295 hits, with
17 of its 87 aligned peptide-binding-region pairs identical. mFuHC, the Botryllus
histocompatibility gene, matches mouse H-2 class I (HA23_MOUSE, 62 pairs in the
peptide-binding region) at rank 922 of 1,089, and rat class II beta (HB2D_RAT) at
11% identity, gaps throughout (`results/example_alignments.txt`).

### Completeness check

The full search keeps at most 2,000 candidate targets per protein from its k-mer
step before aligning. 200 Botryllus proteins in Foldseek and 55 in MMseqs2 ended with
exactly 2,000 hits, and more may have been cut before the e-value filter. To make sure no MHC match was cut, every Botryllus protein was also
aligned to every one of the 384 MHC, B2M and C1-set entries with no candidate step
(`make run-mhc-exhaustive`).

| | MMseqs2 | Foldseek |
|---|---|---|
| Pairs (Botryllus protein, MHC/C1-set entry) at e < 1e-3, full search | 18 | 7,135 |
| Found again by the exhaustive search | 18 | 7,135 |
| Found only by the exhaustive search, e < 1e-3 on the full-database scale | 7 | 14,453 |
| Botryllus proteins with such a pair | 4 | 647 |
| E-value scale factor, full / exhaustive (min-max over shared pairs) | 2,024-2,027 | 10.81-10.83 |

The hits the full search missed come mostly from proteins well under the 2,000 limit
(mFuHC has 1,089 hits), so its k-mer candidate step skipped them, not the limit.
They are in the tables with `found_by = exhaustive check only`, their e-values
rescaled by the factor above, and a rank computed from bit scores (`rank_by_bits`).
None of them is a protein's best hit either. The figure and the counts above include
them.

## Files

| File | What it holds |
|---|---|
| `results/mhc_hits.tsv` | every hit to an MHC class I or II protein, both tools |
| `results/mhc_best_hit_per_protein.tsv` | one row per tool and Botryllus protein: its best MHC hit |
| `results/best_hit_per_protein_and_class.tsv` | each protein's best hit to MHC class I, class II, B2M and other C1-set proteins |
| `results/b2m_hits.tsv` | every hit to B2M |
| `results/b2m_and_c1set_summary.tsv` | number of Botryllus proteins and hits per tool and target class |
| `results/exhaustive_check_summary.tsv` | the completeness check table above |
| `results/example_alignments.txt` | aligned residues for g62973.t1 vs DRA_MACMU and mFuHC vs two MHC proteins |
| `results/rank_of_best_hit_in_class.{pdf,png,svg}` | the figure |

Columns: `evalue`, `bits`, `qcov`/`tcov` (fraction of the Botryllus / Swiss-Prot
protein covered by the alignment), `pident` (MMseqs2) or `fident` (Foldseek) identity,
`hit_rank` (position in the full search's list, 1 = best), `rank_by_bits` (1 + the
number of the protein's full-search hits with a higher bit score), `n_hits` (all
full-search hits of that protein at e < 1e-3; 2,000 is the limit), the protein's best
hit of any kind (`best_hit_*`), `target_id` (the Swiss-Prot accession or AlphaFold
model; `target_is_canonical` is false for AlphaFold models of isoforms, whose
residue numbers do not match UniProt's), the peptide-binding-region counts described
above, and the aligned residues (`qaln`, `taln`).

## Running it

On the laptop: `make push` copies the cleaned proteome, scripts and Makefile to
`$SCRATCH/botryllus-mhc-swissprot` on Sherlock. There: `make databases` (login node,
downloads only), `make test` (first 1,000 proteins), `make run`, `make prefilter`.
Then the completeness check: `make mhc-databases run-mhc-exhaustive`,
`make prefilter RUN=mhc_exhaustive`, `make rank-bits`. Back on the laptop:
`make pull`, `make tables`.

Software: MMseqs2 module `mmseqs2/18-8cc5c`; Foldseek 10 static GPU build (commit
463739e) in `$GROUP_HOME/opt/foldseek-10`.

Resources, measured:

| Step | Request | 1,000-protein test | Full proteome |
|---|---|---|---|
| MMseqs2 | 32 CPU, 32 GB, 3 h | 35 s, 3.2 GB (16 CPU) | 4 min 42 s, 7.7 GB |
| ProstT5 3Di prediction | 1 GPU (Ampere or newer), 8 CPU, 16 GB, 4 h | 49 s, 1.0 GB (L40S) | 36 min, 1.0 GB (L40S) |
| Foldseek search | 32 CPU, 32 GB, 8 h | 3 min 31 s, 3.0 GB (16 CPU) | 27 min 46 s |
| Exhaustive check, MMseqs2 / Foldseek | 32 CPU, 32 GB, 2 h | | 34 s / 43 s |

Problems met and fixed:

- `mmseqs databases UniProtKB/Swiss-Prot` downloads from ftp.expasy.org, which
  Sherlock's system curl cannot reach (TLS protocol error). The Makefile takes
  `uniprot_sprot.fasta.gz` from ftp.uniprot.org instead and runs `mmseqs createdb`.
- ProstT5 in this Foldseek build stopped with a CUDA error on a TITAN Xp (compute
  capability 6.1). The job now asks for compute capability 8.0 or higher.
- `foldseek convertalis` cannot write the `prob` column for ProstT5 queries (it needs
  query coordinates), so the output has no `prob`.
- The `dev` partition does not take batch jobs, so the test runs on `hns,normal`.
