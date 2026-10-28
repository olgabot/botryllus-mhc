## Summary and conclusions

Tool for every number below: kmerseek 0.4.0, the 19 alphabets of notebook 233, each at the k
that carries at least the 28.45 bits of protein20 at k = 7 (Swiss-Prot residue frequencies),
scaled = 1, low-complexity k-mers removed. Query: 45,339 FUN Botryllus proteins. Target:
19,429 human GENCODE v49 canonical proteins.

1. **No FUN protein is an MHC-like candidate, under either rule written before the run.**
   * **16 alphabets with E-values:** 0 HLA matches at E ≤ 0.001 (notebook 234 had 17, all
     funcgroups8 matches to the HLA membrane helix). At E ≤ 1 there are 31 HLA matches from
     19 FUN proteins; 17 cover the groove or the Ig-like domain, all with E between 0.09 and
     0.99. The one with the most shared amino acids beyond chance is FUN000975_XIAP against
     the HLA-A groove (hp_lehninger2, E = 0.94): 9 of 50 positions the same, chance 0.0026
     over all positions and 0.044 over the 45 positions where both letters are in one
     class. With 17 matches tested, neither passes. FUN016327_flrt2, the 12-of-94 match of
     notebook 234, is back at E = 0.12; counted within letter classes its chance is 0.2.
   * **3 alphabets without E-values** (hp_kyte_doolittle2 k = 31, gbmr7 k = 23, hsdm17
     k = 9): their fit failed even with 19,000 calibration proteins. 263 unextended HLA
     matches; 92 have p below 0.05 / 880,891,431, and 74 of these cover the groove or
     Ig-like domain (70 from hp_kyte_doolittle2, 4 from gbmr7). None shares more amino
     acids than chance within letter classes: the smallest chance is 0.099, where step 3
     needs below 0.00071 (hp_kyte_doolittle2) or 0.0125 (gbmr7). hsdm17's 15 matches below
     the threshold are all membrane helix or signal peptide.

2. **The longer k removes most chance matches, and the funcgroups8 helix matches with them.**
   Summed over the 16 alphabets with E-values in both runs (a protein found by several
   alphabets counts once for each), immune matches at E ≤ 0.001 fall from 14,637 FUN
   proteins in notebook 234 to 6,722, and HLA matches at E ≤ 1 from 452 to 25. funcgroups8
   at k = 13 keeps one HLA match at E ≤ 1, where at k = 9 it had 438 proteins. hp_lehninger2
   and hp_lehninger_c_nonpolar2 rise at E ≤ 1 (13 and 9 proteins, from 7 and 6).

3. **Well-known immune genes are still found through domains every cell has.** 23 genes
   from the hand-picked immune families have E ≤ 0.001 in at least one alphabet. The
   strongest: PSMB8 (E = 3.5e-24, 12 alphabets), TRAF3 (1.7e-22, 12),
   PSMB10 (5.6e-17, 6) and LCK (1.1e-16, 12): proteasome, TRAF and kinase domains.

4. **Unextended matches in the alphabets with E-values add no candidate.** 4,196
   unextended HLA matches come from polarity4, gbmr4 and hp_lehninger_hpc3; 131 have p
   below the corrected threshold, 40 of them on the groove or Ig-like domain. None shares
   more amino acids than chance within letter classes (smallest chance 0.022).

5. **dayhoff6 at k = 14 gives no usable E-values.** Its fit slope is 0.014 where the closed
   form gives 0.464, and only 16 of its kept matches have E ≤ 10. It has no HLA match,
   extended or not, so it cannot hide a candidate, but its row in sections 2 to 4 is
   empty for this reason, not because Botryllus lacks immune matches.

6. **The run is sound where it can be checked.** For FUN proteins named by funannotate, the
   best immune match has the same human gene name in 39% to 54% of named proteins, over
   the 13 alphabets with at least 100 of them (protein20: 429 of 1,100).

**Verdict.** The "no MHC-like candidate from kmerseek" verdict of notebooks 232 to 234 holds
at the k of protein20 k = 7, in all 19 alphabets. Longer k removed chance matches without
bringing up a groove or Ig-like match, so further kmerseek settings are unlikely to change
this; a structure search (FoldSeek on predicted structures) is the remaining test.
