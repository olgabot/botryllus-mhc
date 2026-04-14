# Makefile for nf-core/proteinannotator (PR #79: interproscan batches)
# Branch: olgabot/interproscan-batches
# Run from this directory: nextflow_runs/nf-core-proteinannotator/

FASTA_FILE ?= $(HOME)/data/botryllus/Bs_proteins.fa
CLEAN_FASTA ?= botryllus_clean.fa
SAMPLE_ID ?= botryllus
OUTPUT_DIR ?= results
TEST_FASTA ?= test_botryllus.fa
TEST_SAMPLE_ID ?= botryllus_test
TEST_OUTPUT_DIR ?= test_results

PROFILE ?= docker
PIPELINE_BRANCH ?= olgabot/interproscan-batches
WORK_DIR ?= work

# InterProScan database (downloaded once and reused)
INTERPROSCAN_DB ?= $(HOME)/data/interproscan

# InterProScan batch size (PR #79 feature: groups proteins into batches)
INTERPROSCAN_BATCH_SIZE ?= 1000

# InterProScan applications
# Full set (slow - PANTHER alone can take >12h at --cpu 1): Hamap,PANTHER,PIRSF,TIGRFAM,sfld
# Fast subset for testing (Hamap only - profile alignment, no HMMER search; TIGRFAM/PIRSF take >1h even for 10 proteins):
TEST_INTERPROSCAN_APPS ?= Hamap
# Pfam uses HMMER internally in InterProScan
INTERPROSCAN_APPS ?= Hamap,PANTHER,PIRSF,TIGRFAM,sfld,Pfam

# No skip flags: run s4pred (secondary structure), pfam (standalone HMMER/Pfam), and funfam (FunFAM)
SKIP_FLAGS =

.PHONY: help test run resume clean check-fasta

help:
	@echo "nf-core/proteinannotator - Botryllus MHC annotation"
	@echo "Branch: $(PIPELINE_BRANCH)"
	@echo ""
	@echo "Targets:"
	@echo "  make test     Run on first 10 proteins (InterProScan only)"
	@echo "  make run      Run on full Botryllus proteins"
	@echo "  make resume   Resume a previous run"
	@echo "  make clean    Remove work dir and logs (keeps results)"
	@echo ""
	@echo "Variables:"
	@echo "  FASTA_FILE=$(FASTA_FILE)"
	@echo "  SAMPLE_ID=$(SAMPLE_ID)"
	@echo "  PROFILE=$(PROFILE)"
	@echo "  INTERPROSCAN_BATCH_SIZE=$(INTERPROSCAN_BATCH_SIZE)"

# Extract first 10 proteins and run pipeline
test:
	@echo "Extracting first 10 proteins from $(FASTA_FILE) (stripping stop codons)..."
	@awk '/^>/{n++} n>10{exit} /^>/{print; next} {gsub(/\*/, ""); if (length > 0) print}' $(FASTA_FILE) > $(TEST_FASTA)
	@echo "Created $(TEST_FASTA) with $$(grep -c '^>' $(TEST_FASTA)) sequences"
	@echo "id,fasta" > samplesheet_$(TEST_SAMPLE_ID).csv
	@echo "$(TEST_SAMPLE_ID),$$(realpath $(TEST_FASTA))" >> samplesheet_$(TEST_SAMPLE_ID).csv
	nextflow run nf-core/proteinannotator \
	    -r $(PIPELINE_BRANCH) \
	    -profile $(PROFILE) \
	    -c nextflow.config \
	    --input samplesheet_$(TEST_SAMPLE_ID).csv \
	    --outdir $(TEST_OUTPUT_DIR) \
	    -work-dir $(WORK_DIR) \
	    $(SKIP_FLAGS) \
	    --interproscan_db $(INTERPROSCAN_DB) \
	    --interproscan_applications $(TEST_INTERPROSCAN_APPS) \
	    --interproscan_batch_size $(INTERPROSCAN_BATCH_SIZE)

# Strip stop codons (*) from the full FASTA — InterProScan rejects them
$(CLEAN_FASTA): check-fasta
	@echo "Stripping stop codons from $(FASTA_FILE) -> $(CLEAN_FASTA)..."
	@awk '/^>/{print; next} {gsub(/\*/, ""); if (length > 0) print}' $(FASTA_FILE) > $(CLEAN_FASTA)
	@echo "Done: $$(grep -c '^>' $(CLEAN_FASTA)) sequences in $(CLEAN_FASTA)"

# Run on full Botryllus proteome
run: $(CLEAN_FASTA)
	@echo "id,fasta" > samplesheet_$(SAMPLE_ID).csv
	@echo "$(SAMPLE_ID),$$(realpath $(CLEAN_FASTA))" >> samplesheet_$(SAMPLE_ID).csv
	nextflow run nf-core/proteinannotator \
	    -r $(PIPELINE_BRANCH) \
	    -profile $(PROFILE) \
	    -c nextflow.config \
	    --input samplesheet_$(SAMPLE_ID).csv \
	    --outdir $(OUTPUT_DIR) \
	    -work-dir $(WORK_DIR) \
	    $(SKIP_FLAGS) \
	    --interproscan_db $(INTERPROSCAN_DB) \
	    --interproscan_applications $(INTERPROSCAN_APPS) \
	    --interproscan_batch_size $(INTERPROSCAN_BATCH_SIZE)

resume: $(CLEAN_FASTA)
	@echo "id,fasta" > samplesheet_$(SAMPLE_ID).csv
	@echo "$(SAMPLE_ID),$$(realpath $(CLEAN_FASTA))" >> samplesheet_$(SAMPLE_ID).csv
	nextflow run nf-core/proteinannotator \
	    -r $(PIPELINE_BRANCH) \
	    -profile $(PROFILE) \
	    -c nextflow.config \
	    --input samplesheet_$(SAMPLE_ID).csv \
	    --outdir $(OUTPUT_DIR) \
	    -work-dir $(WORK_DIR) \
	    $(SKIP_FLAGS) \
	    --interproscan_db $(INTERPROSCAN_DB) \
	    --interproscan_applications $(INTERPROSCAN_APPS) \
	    --interproscan_batch_size $(INTERPROSCAN_BATCH_SIZE) \
	    -resume

check-fasta:
	@if [ ! -f "$(FASTA_FILE)" ]; then \
	    echo "Error: FASTA file '$(FASTA_FILE)' not found"; \
	    exit 1; \
	fi
	@echo "Found FASTA: $(FASTA_FILE)"

clean:
	rm -rf $(WORK_DIR)/
	rm -f .nextflow.log*
	rm -f samplesheet_*.csv
	rm -f $(TEST_FASTA) $(CLEAN_FASTA)
	@echo "Cleaned (results preserved)"

clean-all: clean
	rm -rf $(OUTPUT_DIR)/ $(TEST_OUTPUT_DIR)/
