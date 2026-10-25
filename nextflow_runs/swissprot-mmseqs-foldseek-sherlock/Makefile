# Search the Botryllus schlosseri proteome against UniProtKB/Swiss-Prot with
# MMseqs2 (sequence) and Foldseek (predicted 3Di structure letters from
# ProstT5, against AlphaFold/Swiss-Prot), then pull out every hit to an MHC
# class I or class II protein.
#
# Two places run this Makefile:
#
#   On the laptop, from this directory:
#     make push        copy the cleaned proteome, scripts and this Makefile to Sherlock
#     make pull        copy the pre-filtered hit tables back
#     make tables      build the MHC hit tables in results/ (polars)
#
#   On a Sherlock login node, from $(SHERLOCK_DIR):
#     make databases   download Swiss-Prot (MMseqs2), AlphaFold/Swiss-Prot and
#                      ProstT5 weights (Foldseek). Login node only: compute
#                      nodes have no internet.
#     make test        1,000-protein test of both searches, to size the full run
#     make run         full proteome, 72,617 proteins
#     make prefilter   cut the full tables down to MHC/C1-set targets and the
#                      best hit per protein, so `make pull` copies little
#     make mhc-databases run-mhc-exhaustive, then make prefilter RUN=mhc_exhaustive
#                      and make rank-bits: every protein against every MHC/C1-set
#                      entry, no candidate step (see below)
#
# Never run the searches on the laptop.

SHERLOCK_HOST ?= sherlock
SHERLOCK_DIR  ?= /scratch/users/olgabot/botryllus-mhc-swissprot

FASTA_FILE  ?= $(HOME)/data/botryllus/Bs_proteins.fa
CLEAN_FASTA ?= input/botryllus_clean.fa
TEST_FASTA  ?= input/botryllus_test1000.fa
TEST_N      ?= 1000

PYTHON ?= $(HOME)/anaconda3/envs/2025-kmerseek-analysis/bin/python
FOLDSEEK ?= $(GROUP_HOME)/opt/foldseek-10/foldseek/bin/foldseek
SP_MMSEQS_DB   = db/swissprot/swissprot
AFDB_SP_DB     = db/afdb_swissprot/afdb_swissprot
PROSTT5        = db/prostt5_weights

# Slurm requests. Sized from the 1,000-protein test (see README.md); the
# full-run values scale that runtime by 72.6 and add a margin.
PARTITION      ?= hns,normal
GPU_PARTITION  ?= gpu,hns
# ProstT5 in this Foldseek build crashed with a CUDA error on a TITAN Xp
# (compute capability 6.1) and ran on an L40S (8.9), so ask for Ampere or newer.
GPU_CONSTRAINT ?= GPU_CC:8.0|GPU_CC:8.6|GPU_CC:8.9|GPU_CC:9.0
TEST_PARTITION ?= hns,normal

MMSEQS_CPUS ?= 32
MMSEQS_MEM  ?= 32G
MMSEQS_TIME ?= 03:00:00

PROSTT5_CPUS ?= 8
PROSTT5_MEM  ?= 16G
PROSTT5_TIME ?= 04:00:00

FOLDSEEK_CPUS ?= 32
FOLDSEEK_MEM  ?= 32G
FOLDSEEK_TIME ?= 08:00:00

.PHONY: help push pull tables databases test run prefilter

help:
	@sed -n '1,25p' Makefile

# ---------------------------------------------------------------- laptop ---

# Strip `*` stop characters the same way ../interproscan6-local/Makefile does.
$(CLEAN_FASTA): $(wildcard $(FASTA_FILE))
	mkdir -p input
	awk '/^>/{print; next} {gsub(/\*/, ""); if (length > 0) print}' $< > $@
	@echo "$$(grep -c '^>' $@) proteins in $@"

push: $(CLEAN_FASTA)
	ssh $(SHERLOCK_HOST) 'mkdir -p $(SHERLOCK_DIR)/input $(SHERLOCK_DIR)/scripts $(SHERLOCK_DIR)/reference'
	rsync -a Makefile $(SHERLOCK_HOST):$(SHERLOCK_DIR)/
	rsync -a scripts/ $(SHERLOCK_HOST):$(SHERLOCK_DIR)/scripts/
	rsync -a reference/ $(SHERLOCK_HOST):$(SHERLOCK_DIR)/reference/
	rsync -a $(CLEAN_FASTA) $(SHERLOCK_HOST):$(SHERLOCK_DIR)/input/

pull:
	mkdir -p prefiltered
	rsync -a $(SHERLOCK_HOST):$(SHERLOCK_DIR)/prefiltered/ prefiltered/
	rsync -a $(SHERLOCK_HOST):$(SHERLOCK_DIR)/db/swissprot/swissprot_headers.tsv.gz prefiltered/

tables:
	$(PYTHON) scripts/mhc_hit_tables.py --prefiltered prefiltered --reference reference/swissprot_mhc_c1set_pfam.tsv --outdir results
	$(PYTHON) scripts/plot_mhc_hits.py results/best_hit_per_protein_and_class.tsv results/rank_of_best_hit_in_class

# -------------------------------------------------------------- Sherlock ---

databases: $(SP_MMSEQS_DB) $(AFDB_SP_DB) $(PROSTT5)/prostt5-f16.gguf db/swissprot/swissprot_headers.tsv.gz

# `mmseqs databases UniProtKB/Swiss-Prot` fetches from ftp.expasy.org, which
# Sherlock's system curl cannot reach (TLS protocol error), so take the FASTA
# from ftp.uniprot.org and build the database here.
UNIPROT_FTP = https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete
$(SP_MMSEQS_DB):
	mkdir -p db/swissprot
	curl -sSf -o db/swissprot/reldate.txt $(UNIPROT_FTP)/reldate.txt
	curl -sSf -o db/swissprot/uniprot_sprot.fasta.gz $(UNIPROT_FTP)/uniprot_sprot.fasta.gz
	bash -lc 'ml biology mmseqs2/18-8cc5c && mmseqs createdb db/swissprot/uniprot_sprot.fasta.gz $@'
	cat db/swissprot/reldate.txt

# accession -> full Swiss-Prot header (protein name, OS=, OX=, GN=), for the
# Foldseek hits, whose AlphaFold targets carry only a model ID.
db/swissprot/swissprot_headers.tsv.gz: $(SP_MMSEQS_DB)
	tr -d '\0' < $(SP_MMSEQS_DB)_h | awk '{acc=$$1; sub(/^sp\|/,"",acc); sub(/\|.*/,"",acc); print acc "\t" $$0}' | gzip > $@

$(AFDB_SP_DB):
	mkdir -p db/afdb_swissprot
	$(FOLDSEEK) databases Alphafold/Swiss-Prot $@ db/tmp_afsp --threads 2

$(PROSTT5)/prostt5-f16.gguf:
	$(FOLDSEEK) databases ProstT5 $(PROSTT5) db/tmp_prostt5

$(TEST_FASTA): $(CLEAN_FASTA)
	awk -v n=$(TEST_N) '/^>/{c++} c>n{exit} {print}' $< > $@

logs:
	mkdir -p logs

# Each search writes out/<run>/<tool>.tsv. Foldseek's ProstT5 step runs first
# on a GPU and the search waits for it with --dependency=afterok.
# Arguments: run name, query FASTA, CPU partition, mmseqs time | prostt5 time, foldseek time.
define submit_mmseqs
	mkdir -p out/$(1)
	sbatch --parsable -p $(3) --time=$(4) -c $(MMSEQS_CPUS) --mem=$(MMSEQS_MEM) \
	    scripts/mmseqs_search.sbatch $(2) $(SP_MMSEQS_DB) out/$(1)/mmseqs_swissprot.tsv
endef
define submit_foldseek
	mkdir -p out/$(1)
	jid=$$(sbatch --parsable -p $(GPU_PARTITION) -G 1 -C '$(GPU_CONSTRAINT)' --time=$(4) \
	    -c $(PROSTT5_CPUS) --mem=$(PROSTT5_MEM) --export=ALL,FOLDSEEK=$(FOLDSEEK) \
	    scripts/foldseek_prostt5_createdb.sbatch $(2) $(PROSTT5) out/$(1)/qdb/botryllus) && \
	echo "prostt5 job $$jid" && \
	sbatch --parsable -p $(3) --time=$(5) -c $(FOLDSEEK_CPUS) --mem=$(FOLDSEEK_MEM) \
	    --dependency=afterok:$$jid --export=ALL,FOLDSEEK=$(FOLDSEEK) \
	    scripts/foldseek_search.sbatch out/$(1)/qdb/botryllus $(AFDB_SP_DB) out/$(1)/foldseek_afdb_swissprot.tsv
endef

.PHONY: run-mmseqs run-foldseek
test: $(TEST_FASTA) | logs
	$(call submit_mmseqs,test1000,$(TEST_FASTA),$(TEST_PARTITION),01:00:00)
	$(call submit_foldseek,test1000,$(TEST_FASTA),$(TEST_PARTITION),01:00:00,01:00:00)

run: run-mmseqs run-foldseek
run-mmseqs: $(CLEAN_FASTA) | logs
	$(call submit_mmseqs,full,$(CLEAN_FASTA),$(PARTITION),$(MMSEQS_TIME))
run-foldseek: $(CLEAN_FASTA) | logs
	$(call submit_foldseek,full,$(CLEAN_FASTA),$(PARTITION),$(PROSTT5_TIME),$(FOLDSEEK_TIME))

# Completeness check for the 2,000-hit limit. A Botryllus Ig-domain protein can
# have more than 2,000 Swiss-Prot candidates, and the search keeps only the top
# 2,000 before aligning, so a weak MHC match could be cut. Here every Botryllus
# protein is aligned to every one of the 384 entries in
# reference/swissprot_mhc_c1set_pfam.tsv with no candidate step at all
# (MMseqs2 --prefilter-mode 2, Foldseek --exhaustive-search 1). Bit scores do
# not depend on database size; e-values are rescaled to the full database in
# scripts/mhc_hit_tables.py. Results go to out/mhc_exhaustive/, then
# `make prefilter RUN=mhc_exhaustive`.
SP_MHC_DB   = db/mhc_c1set/swissprot_mhc_c1set
AFDB_MHC_DB = db/mhc_c1set/afdb_mhc_c1set
.PHONY: mhc-databases run-mhc-exhaustive
mhc-databases: $(SP_MHC_DB) $(AFDB_MHC_DB)
$(SP_MHC_DB): reference/swissprot_mhc_c1set_pfam.tsv
	mkdir -p db/mhc_c1set
	awk -F'\t' 'NR == FNR {if (FNR > 1) k[$$1] = 1; next} ($$2 in k) {print $$1}' $< $(SP_MMSEQS_DB).lookup > db/mhc_c1set/swissprot_keys.txt
	bash -lc 'ml biology mmseqs2/18-8cc5c && mmseqs createsubdb db/mhc_c1set/swissprot_keys.txt $(SP_MMSEQS_DB) $@'
	@echo "$$(wc -l < $@.index) Swiss-Prot entries in $@"
$(AFDB_MHC_DB): reference/swissprot_mhc_c1set_pfam.tsv
	mkdir -p db/mhc_c1set
	awk -F'\t' 'NR == FNR {if (FNR > 1) k[$$1] = 1; next} {split($$2, a, "-")} (a[2] in k) {print $$1}' $< $(AFDB_SP_DB).lookup > db/mhc_c1set/afdb_keys.txt
	$(FOLDSEEK) createsubdb db/mhc_c1set/afdb_keys.txt $(AFDB_SP_DB) $@
	@echo "$$(wc -l < $@.index) AlphaFold models in $@"

run-mhc-exhaustive: $(SP_MHC_DB) $(AFDB_MHC_DB) | logs
	mkdir -p out/mhc_exhaustive
	sbatch --parsable -p $(PARTITION) --time=02:00:00 -c $(MMSEQS_CPUS) --mem=$(MMSEQS_MEM) \
	    --export=ALL,EXTRA_ARGS="--prefilter-mode 2",EXTRA_COLUMNS="qaln taln" \
	    scripts/mmseqs_search.sbatch $(CLEAN_FASTA) $(SP_MHC_DB) out/mhc_exhaustive/mmseqs_swissprot.tsv
	sbatch --parsable -p $(PARTITION) --time=02:00:00 -c $(FOLDSEEK_CPUS) --mem=$(FOLDSEEK_MEM) \
	    --export=ALL,FOLDSEEK=$(FOLDSEEK),EXTRA_ARGS="--exhaustive-search 1 -a",EXTRA_COLUMNS="qaln taln" \
	    scripts/foldseek_search.sbatch out/full/qdb/botryllus $(AFDB_MHC_DB) out/mhc_exhaustive/foldseek_afdb_swissprot.tsv

# Bit scores of every full-search hit of the proteins that have a hit in the
# exhaustive check, so that a hit only the exhaustive check found can still be
# ranked among that protein's hits. Run after both prefilter runs.
.PHONY: rank-bits
rank-bits:
	for tool in mmseqs_swissprot foldseek_afdb_swissprot; do \
	    awk -F'\t' 'NR == FNR {q[$$1] = 1; next} FNR > 1 && ($$1 in q) {print $$1 "\t" $$13}' \
	        prefiltered/mhc_exhaustive.$$tool.n_hits.tsv out/full/$$tool.tsv \
	        > prefiltered/full.$$tool.bits_for_ranking.tsv; \
	done
	wc -l prefiltered/*bits_for_ranking.tsv

# Keep (a) every hit whose target is in reference/swissprot_mhc_c1set_pfam.tsv,
# with its rank among that protein's hits, and (b) each protein's best hit.
RUN ?= full
prefilter:
	mkdir -p prefiltered
	for tool in mmseqs_swissprot foldseek_afdb_swissprot; do \
	    [ -f out/$(RUN)/$$tool.tsv ] || { echo "no out/$(RUN)/$$tool.tsv yet"; continue; }; \
	    awk -F'\t' -f scripts/prefilter.awk reference/swissprot_mhc_c1set_pfam.tsv out/$(RUN)/$$tool.tsv \
	        > prefiltered/$(RUN).$$tool.mhc_c1set_hits.tsv; \
	    awk -F'\t' 'NR==1 || $$1 != prev {print} {prev=$$1}' out/$(RUN)/$$tool.tsv \
	        > prefiltered/$(RUN).$$tool.best_hit.tsv; \
	    awk -F'\t' 'NR>1 {n[$$1]++} END {for (q in n) print q "\t" n[q]}' out/$(RUN)/$$tool.tsv \
	        > prefiltered/$(RUN).$$tool.n_hits.tsv; \
	done
	wc -l prefiltered/*
