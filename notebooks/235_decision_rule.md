# Notebook 235: decision rules, written 28 Oct 2026 before any notebook 235 result was read

Run: kmerseek 0.4.0, FUN Botryllus proteome (45,339 proteins) against human GENCODE v49
canonical (19,429 proteins), each of 19 alphabets at the k carrying at least the bits of
protein20 k = 7 (28.45 bits, Swiss-Prot 2026_03 residue frequencies; `235_arms.tsv`).

## Alphabets with E-values (the run extends matches and fits E-values)

The rule of notebooks 232 to 234, unchanged. A FUN protein is an MHC-like candidate if, in
any alphabet, it matches an HLA gene at E ≤ 0.001 and the matched part of the HLA protein
is the peptide-binding groove (class I α1/α2, class II α1/β1) or the Ig-like domain. A
match only to the transmembrane helix or the signal peptide is not MHC evidence. Agreement
between alphabets is not evidence on its own.

## Alphabets without E-values (fallback)

If an alphabet's E-value fit fails even with 19,000 calibration sequences (job 49956388),
it is searched without extension (`--extend-mismatch-penalty 0`), so its matches have no
E-value, only kmerseek's region p-value (the chance of that many shared k-mers). These
alphabets are reported in their own section, never pooled with the E-value alphabets.

The region p-value alone is not enough: in notebooks 233 and 234, unrelated charged
stretches got p = 0. So a FUN protein is an MHC-like candidate from such an alphabet only
if one match passes all three:

1. region p-value below 0.05 / 880,891,431 (the protein pairs one alphabet compares);
2. the matched part of the HLA protein is the groove or the Ig-like domain (as above);
3. the two stretches carry the same amino acid at more positions than unrelated
   sequences would, with the chance level computed within each letter class
   (`hla_utils.same_residue_vs_chance_in_class`), at a chance below 0.05 divided by the
   number of matches passing 1 and 2 in that alphabet.

If no match passes, the verdict of notebooks 232 to 234 (no MHC-like candidate from
kmerseek) holds for that alphabet too.
