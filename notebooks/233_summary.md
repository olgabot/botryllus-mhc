## Summary and conclusions

Tool for every number below: kmerseek 0.4.0, 19 alphabets each at the k that carries at least
24 bits per k-mer, scaled = 1, low-complexity k-mers removed. Query: 72,617 Botryllus proteins.
Target: 19,429 human GENCODE v49 canonical proteins.

1. **No alphabet finds an MHC-like candidate under the rule written at the top.** Across the 19
   alphabets there are 3,950 HLA matches with E ≤ 10, from 1,363 Botryllus proteins. The
   strongest match to the peptide-binding groove or the Ig-like domain, in any alphabet, is the
   one notebook 232 found: g564.t1 to the HLA-DRA Ig-like domain, E = 0.13 (hp_lehninger2,
   k = 24). No other groove or Ig-like match has E below 0.17.

   At E ≤ 1 there are 842 HLA matches from 302 Botryllus proteins, and 832 of them are to the
   transmembrane helix (805 from funcgroups8). Only 10 cover the groove or the Ig-like domain:
   6 Botryllus proteins, all with E between 0.13 and 0.92, and all from HP alphabets except one
   gbmr4 match. Six of the ten share no more amino acids than two unrelated sequences would, with Swiss-Prot
   residue frequencies (2 to 5 of 33 to 78 positions). The other four are one 94-residue
   stretch, carried by g65873.t1 and g56468.t1 and found by two alphabets each: 12 of 94
   positions carry the same amino acid where 5.6 are expected (chance 0.010, about 0.1 once
   the ten matches are counted). That stretch runs from the HLA-B Ig-like domain into its
   membrane helix (human residues 245-338). The matches come as near-copies: g42928.t1,
   g59107.t1 and g67578.t1 carry almost the same 62 residues and all match the same HLA-DOA
   stretch. With 72,617 Botryllus proteins searched, E ≤ 1 lets
   each one have about one match by chance, so 10 weak matches are what chance alone gives.

2. **Only funcgroups8 (k = 9) gives HLA matches below E = 0.001, and all 14 are to the HLA
   transmembrane helix.** They come from 9 Botryllus proteins, for example g7000.t1
   (`TSDTGLIIGVVAALLGIAALASLAV`) against the HLA-DPA1 helix (`TTETVLCALGLVLGLVGIIVGTVLI`),
   E = 0.000087. Of funcgroups8's 2,710 HLA matches, 2,645 cover the helix. The same alphabet
   gives HLA, Fc receptor and LILR genes, all single-pass membrane proteins, near-identical
   best E-values around 0.0001 (figure 233e). funcgroups8 is also the one alphabet with no
   measured extension penalty (it runs on the fallback 2 / 8), so its E-values are the least
   trustworthy of the 19.

3. **Five alphabets leave almost every match unextended, so they say little either way.** For
   dayhoff6 (k = 11), gbmr7 (k = 13), gbmr4 (k = 16), hp_kyte_doolittle2 (k = 26) and
   hp_thomas_dill_no_c2 (k = 26), over 99.7% of kept matches have E = inf; dayhoff6 has 41
   matches with a finite E-value out of 5.1 million. "No HLA match at E ≤ 10" in those
   alphabets means no E-value was computed, not that nothing matched.

4. **The unextended HLA matches look like shared composition, not shared ancestry.** gbmr7,
   gbmr4, polarity4 and hp_lehninger_hpc3 have 76,412, 31,044, 33,977 and 3,751 unextended
   HLA matches. 1,535 of them pass a k-mer p-value corrected for the 1.4 billion protein pairs
   each alphabet compares, and 204 have p = 0 exactly. The strongest pair a charged,
   repetitive Botryllus stretch (`REKRRIQEQLRRIKRNQEKEK`) with the HLA α2 helix: in gbmr4 every
   polar residue falls into one class, so any polar run matches any other. A p-value of exactly
   0 for residues this unlike each other points to a limit of kmerseek's unextended-region
   statistic, not to homology.

5. **Well-known immune genes are found strongly by many alphabets, but through domains every
   cell has.** PSMB8 (best E = 8e-46, 12 alphabets at E ≤ 0.001), LCK (13 alphabets), ZAP70
   (12), SYK (11), PSMB10 (11) and TAP2 (10) are kinase, proteasome and ABC-transporter domains.
   The NF-κB family (NFKB1, RELA, REL, RELB) is found by 7 or 8 alphabets each, which fits a
   Rel-family protein in Botryllus. The alphabets with 12 or more classes and protein20 find
   these best; most HP alphabets find them less strongly.

**Verdict against the rule written at the top:** none of the 19 alphabets gives a Botryllus
protein matching an HLA groove or Ig-like domain at E ≤ 0.001, so no protein goes forward as an
MHC-like candidate from kmerseek. The "no HLA match" of notebook 232 does not depend on the
alphabet. The one alphabet with strong HLA matches (funcgroups8) finds membrane helices.

**What this does not settle.** Five alphabets have almost no E-values, so their silence is not
evidence; a kmerseek follow-up should find why so many of their matches are not extended (the
share ranges from 0% for funcgroups8 to 99.999% for dayhoff6 at similar bits per k-mer), and why
unrelated residues can give a region p-value of 0. A search for an MHC fold that has lost
sequence similarity needs structure comparison (Foldseek against HLA structures), not k-mers.

Output tables, in `notebooks/233_outputs/`: runs per alphabet, immune and HLA counts per E-value
cut-off, every HLA match with E ≤ 10 (with residues and the HLA part), HLA matches shared by
alphabets, best E-value per well-known immune gene per alphabet, the unextended totals, and
the 1,535 unextended HLA matches that pass the corrected p-value (all 145,241 are in
`~/data/botryllus/kmerseek-0.4-alphabets/` and on Sherlock). Figures, in `figures/233/`: `233a` runs, `233b` counts per alphabet, `233c` every
HLA match, `233d` matches shared by alphabets, `233e` well-known immune genes, `233f` unextended
matches.
