## Summary and conclusions

Tool for every number below: kmerseek 0.4.0, the 19 alphabets of notebook 233, each at the k
that carries at least 24 bits per k-mer (protein20 at k = 6), scaled = 1, low-complexity
k-mers removed. Query: 45,339 FUN Botryllus proteins. Target: 19,429 human GENCODE v49
canonical proteins.

1. **No FUN protein is an MHC-like candidate under the rule written at the top.** Across the
   19 alphabets there are 4,981 HLA matches with E ≤ 10, from 1,584 FUN proteins. The 17
   matches with E ≤ 0.001 come from 9 proteins, all from funcgroups8 (k = 9), and all cover
   the HLA transmembrane helix. The strongest is FUN018309 (`VIVLAVIVAVAIISIILLLIIVWK`)
   against the HLA-G helix (`GIVAGLVVLAAVVTGAAVAAVLWR`), E = 0.000075.

   At E ≤ 1 there are 1,247 HLA matches from 450 FUN proteins, and 1,224 of them come from
   funcgroups8. Only 10 cover the groove or the Ig-like domain. They come from 7 FUN records
   (6 genes, since FUN028357 has two copies, -T1 and -T2), all from HP alphabets, with E
   between 0.21 and 0.79. They have the same residue at 0 to 12 positions out of 55 to 94,
   about what two unrelated sequences share. Most of the Ig-like matches run on into the HLA
   membrane helix (human residues 245 to 349). The FUN names of these proteins are not
   MHC-like: flrt2 (a leucine-rich repeat membrane protein), PIK3AP1, PRR11, FADF, and three
   with no gene name.

2. **The new gene models change little.** At E ≤ 0.001, summed over the 19 alphabets (a
   protein found by several alphabets counts once for each), 16,212 FUN proteins have an
   immune match against 16,040 old proteins (notebook 233), although the FUN proteome has
   38% fewer proteins. HLA matches at E ≤ 0.001 are
   funcgroups8 helix matches in both: 9 proteins each, 17 matches against 14. At E ≤ 1,
   polarity4 and gbmr4 lose their few HLA matches. funcgroups8 rises from 285 to 438 FUN
   proteins, all of them helix matches.

3. **Well-known immune genes are found strongly through domains every cell has, as in
   notebook 233.** 59 genes from the hand-picked immune families have E ≤ 0.001 in at least
   one alphabet:
   * PSMB8: best E = 5e-46, 12 alphabets.
   * TRAF3: E = 1e-42, 13 alphabets.
   * PSMB10, LCK, BTK, NFKB1, SYK and ZAP70: E between 2e-13 and 3e-32.

   These are proteasome, kinase, TRAF and Rel domains. The Rel family (NFKB1, RELA, REL,
   RELB) is found by 8 to 10 alphabets each.

4. **Unextended matches again look like shared composition.** gbmr7, polarity4, gbmr4,
   hp_lehninger_hpc3 and hp_kyte_doolittle2 have 90,055, 65,480, 36,262, 5,958 and 71
   unextended HLA matches (E = inf). 1,755 of them pass the k-mer p-value corrected for the
   0.88 billion protein pairs each alphabet compares, and 215 have p = 0 exactly. The
   strongest is the same charged stretch notebook 233 found, now named FUN002657_TAF1
   (`REKRRIQEQLRRIKRNQEKEK`), matched to the HLA-G α2 helix in gbmr4.

   gbmr7 leaves all but 38,483 of its 523.9 million matches unextended, as in notebook 233,
   so it has no HLA match at E ≤ 10 and only 744 immune matches. Its 13
   unextended HLA matches with p = 0 have the same residue at 4 to 11 of 21 to 54
   positions. The one with 11 is FUN024501 against the HLA-DOA Ig-like domain:

   ```
   FUN024501 549  LRLKCTTWGDDVIDVTWLKDGRTLTAGQAR
   HLA-DOA   130  LICIVDNIFPPVINITWLRNGQTVTEGVAQ
   ```

   The shared `TWLxxG` stretch is the tryptophan found in nearly every Ig domain, so it
   marks FUN024501 as an Ig-superfamily protein, not as MHC-like.

5. **kmerseek's best immune-gene match has the gene name funannotate gave for 28% to 71% of
   named FUN proteins**, depending on the alphabet. protein20 is lowest (535 of 1,891).
   gbmr7 (1 of 1) and dayhoff6 (none named) have too few matches to count.
   The differences are mostly paralogues: five tubulins match TUBB4B, DNM1 matches DNM2,
   and PPP2CB matches PPP2CA. Only matches to human immune genes were kept, so a protein
   whose own human gene is not on the immune list cannot agree. This is a check of the
   run, not a test of funannotate.

**Verdict against the rule written at the top:** none of the 19 alphabets gives a FUN
protein matching an HLA groove or Ig-like domain at E ≤ 0.001. The verdict of notebooks 232
and 233, no MHC-like candidate from kmerseek, holds for the new gene models.

**What this does not settle.** An MHC fold that has lost sequence similarity needs structure
comparison (Foldseek against HLA structures), not k-mers. The FUN proteome has whole-protein
models, so that comparison could now start from it.

Output tables are in `notebooks/234_outputs/`. They hold the same tables as notebook 233,
plus:
* FUN names against the best immune-gene match;
* old and new proteome counts per alphabet;
* the header changes made before the search.

Figures are in `figures/234/`, `234a` to `234f` as in notebook 233, plus:
* `234g`: FUN name agreement;
* `234h`: old against new proteome.

Full immune-match files: `~/data/botryllus/kmerseek-0.4-alphabets-fun/` and
`$SCRATCH/234-kmerseek-alphabets-fun/` on Sherlock.
