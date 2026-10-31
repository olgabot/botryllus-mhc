# Floor-rounded pairs requested on 31 October 2026

The added comparison uses hp_lehninger2 k = 24 (23.89096060 bits) and k = 27
(26.87733067 bits) as its reference. For each alphabet, divide the reference bits
by that alphabet's bits per letter and round down. Protein20 therefore uses k = 5
and k = 6. This is a new rule requested after the original run; the original
protein20 k6/k7 comparison used ceiling rounding and is preserved.

All 38 pairs are within notebook 274's k ranges. Eight pairs already completed
with identical query, target, kmerseek version, alphabet, k and search settings.
Their 184 best-match tables are reused. The additional work is 30 indexes and
690 searches (23 query chunks per pair).

| Alphabet | hp_lehninger2 k24 threshold | hp_lehninger2 k27 threshold |
|---|---:|---:|
| protein20 | 5 | 6 |
| dayhoff6 | 11 | 12 |
| hp_lehninger2 | 24 | 27 |
| hp_thomas_dill2 | 24 | 27 |
| hp_kyte_doolittle2 | 25 | 29 |
| hp_thomas_dill_no_c2 | 25 | 28 |
| hp_lehninger_c_nonpolar2 | 24 | 27 |
| hp_lehninger_hpc3 | 23 | 26 |
| hp_pbotc_1st_ed2 | 24 | 27 |
| gbmr4 | 18 | 20 |
| polarity4 | 15 | 17 |
| wwmj5 | 11 | 12 |
| gbmr7 | 19 | 21 |
| funcgroups8 | 10 | 11 |
| sdm12 | 8 | 9 |
| mmseqs12 | 7 | 8 |
| wass14 | 6 | 7 |
| hsdm17 | 7 | 8 |
| uniprot18 | 6 | 7 |

## Resource prediction

New searches request 72–168 GB and two hours on their first attempt. The request
is 1.5 times the larger observed peak from the completed Botryllus run and the
earlier dark-set run. A measurement at a larger k is doubled before that margin;
where available, the historical measurement uses the same or a smaller k.
These are estimates for the new pairs, not measurements of this run.

Indexes retain 128 GB and six hours. The original indexes peaked at 53.9 GB and
finished within six hours. Twice that largest peak, with rounding to the existing
128 GB request, covers the chosen allowance for smaller k; this is an estimate.
The head requests 8 GB, two CPUs and two days. Retries are capped at 500 GB.

The original run finished in 2 h 10 min. A few hours after scheduling is a rough
estimate for this addition; the larger memory requests may wait longer. Eight
pairs already have usable results. No new result is claimed before completion.

Each launch regenerates resource requests from the current and historical traces,
then incorporates failures from any earlier floor-run attempts. Requests are
written to the run directory, leaving tracked assets unchanged.

Outputs: `$SCRATCH/botryllus-kmerseek-swissprot-floor-lehninger/results`.
Indexes: the original run's shared index directory. The original results are retained.
Launch target on Sherlock: `make sherlock-floor-run`. Pull target on the laptop:
`make pull-floor`. No kmerseek or Nextflow runs on the laptop.

The pair table and new summaries preserve reference_alphabet, reference_k,
reference_bits and rounding. The legacy protein20_k field records the floor-rounded
protein20 size (5 or 6); it is not the new comparison's reference label.

Historical trace SHA256: `aba305259fe1c39e1ce868947dd525fad1e479421ac5acd866890447963e25be`.
