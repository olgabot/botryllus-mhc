# Makefile to run ebi-pf-team/interproscan6 locally with Docker.
# https://github.com/ebi-pf-team/interproscan6
#
# Use this instead of ../interproscan6-aws/ -- the AWS/Seqera path is
# currently blocked by an upstream bug (interproscan6#340): its sequence-
# loading processes run natively in the Nextflow head job and can't handle
# an S3 work directory. See memory: project-interproscan6-aws-s3-workdir-incompatible.
#
# Requires Nextflow >=25.10.4 and Docker. The nf-core-v2 conda env already
# has a compatible Nextflow (25.10.4); NEXTFLOW defaults to it below.
#
#   make test    # smoke test on the pipeline's bundled data
#   make mini    # first 10 real Botryllus sequences
#   make run     # full Botryllus proteome run
#
# Run from this directory: nextflow_runs/interproscan6-local/

NEXTFLOW ?= /Users/olga/anaconda3/envs/nf-core-v2/bin/nextflow

PIPELINE ?= https://github.com/ebi-pf-team/interproscan6
REVISION ?= 6.0.1

WORK_DIR ?= work

# Persistent InterPro/member-DB cache (multi-GB), reused across runs and
# across test/mini/run so the data isn't re-downloaded every time.
DATA_DIR ?= $(HOME)/data/interproscan6

# InterPro data release, pinned per the README's reproducibility guidance.
INTERPRO_VERSION ?= 109.0

# Comma-separated member DBs to run, e.g. Hamap,PANTHER,Pfam. Leave empty
# to run the pipeline's full default application set.
APPLICATIONS ?=
APPLICATIONS_FLAG = $(if $(APPLICATIONS),--applications $(APPLICATIONS),)

# --- Input ---
FASTA_FILE ?= $(HOME)/data/botryllus/Bs_proteins.fa
CLEAN_FASTA ?= botryllus_clean.fa
MINI_FASTA ?= botryllus_mini10.fa

OUTDIR ?= results
MINI_OUTDIR ?= mini_results
TEST_OUTDIR ?= test_results

.PHONY: help check-nextflow install-nextflow test mini run resume clean clean-all

help:
	@echo "ebi-pf-team/interproscan6, local Docker execution"
	@echo ""
	@echo "Targets:"
	@echo "  make install-nextflow  Install Nextflow to ./nextflow (if not using the conda env)"
	@echo "  make test              Smoke test on the pipeline's bundled data"
	@echo "  make mini              First 10 real Botryllus sequences"
	@echo "  make run               Full run on the Botryllus proteome"
	@echo "  make resume            Resume the full run with -resume"
	@echo "  make clean             Remove work dir and logs (keeps results)"
	@echo ""
	@echo "Variables:"
	@echo "  NEXTFLOW=$(NEXTFLOW)"
	@echo "  FASTA_FILE=$(FASTA_FILE)"
	@echo "  DATA_DIR=$(DATA_DIR)"
	@echo "  INTERPRO_VERSION=$(INTERPRO_VERSION)"
	@echo "  APPLICATIONS=$(APPLICATIONS)"

check-nextflow:
	@command -v $(NEXTFLOW) >/dev/null 2>&1 || { echo "nextflow not found at '$(NEXTFLOW)' -- run 'make install-nextflow' or set NEXTFLOW=/path/to/nextflow"; exit 1; }

install-nextflow:
	curl -s https://get.nextflow.io | bash
	chmod +x nextflow
	@echo "Installed ./nextflow. Use NEXTFLOW=$$(pwd)/nextflow, or move it onto your PATH."

# Strip stop codons (*) -- InterProScan rejects them
$(CLEAN_FASTA):
	@if [ ! -f "$(FASTA_FILE)" ]; then echo "Error: FASTA file '$(FASTA_FILE)' not found"; exit 1; fi
	@echo "Stripping stop codons from $(FASTA_FILE) -> $(CLEAN_FASTA)..."
	@awk '/^>/{print; next} {gsub(/\*/, ""); if (length > 0) print}' $(FASTA_FILE) > $(CLEAN_FASTA)
	@echo "Done: $$(grep -c '^>' $(CLEAN_FASTA)) sequences in $(CLEAN_FASTA)"

# First 10 sequences, stop codons stripped
$(MINI_FASTA):
	@if [ ! -f "$(FASTA_FILE)" ]; then echo "Error: FASTA file '$(FASTA_FILE)' not found"; exit 1; fi
	@echo "Extracting first 10 proteins from $(FASTA_FILE) -> $(MINI_FASTA)..."
	@awk '/^>/{n++} n>10{exit} /^>/{print; next} {gsub(/\*/, ""); if (length > 0) print}' $(FASTA_FILE) > $(MINI_FASTA)
	@echo "Done: $$(grep -c '^>' $(MINI_FASTA)) sequences in $(MINI_FASTA)"

# Smoke test on the pipeline's bundled test FASTA (assets/test.faa, 103 seqs)
test: check-nextflow
	$(NEXTFLOW) run $(PIPELINE) \
	    -r $(REVISION) \
	    -c nextflow.config \
	    -profile docker,test \
	    --datadir $(DATA_DIR) \
	    --interpro $(INTERPRO_VERSION) \
	    --outdir $(TEST_OUTDIR) \
	    -work-dir $(WORK_DIR)

mini: check-nextflow $(MINI_FASTA)
	$(NEXTFLOW) run $(PIPELINE) \
	    -r $(REVISION) \
	    -c nextflow.config \
	    -profile docker \
	    --input $(MINI_FASTA) \
	    --datadir $(DATA_DIR) \
	    --interpro $(INTERPRO_VERSION) \
	    --outdir $(MINI_OUTDIR) \
	    -work-dir $(WORK_DIR) \
	    $(APPLICATIONS_FLAG)

run: check-nextflow $(CLEAN_FASTA)
	$(NEXTFLOW) run $(PIPELINE) \
	    -r $(REVISION) \
	    -c nextflow.config \
	    -profile docker \
	    --input $(CLEAN_FASTA) \
	    --datadir $(DATA_DIR) \
	    --interpro $(INTERPRO_VERSION) \
	    --outdir $(OUTDIR) \
	    -work-dir $(WORK_DIR) \
	    $(APPLICATIONS_FLAG)

resume: check-nextflow $(CLEAN_FASTA)
	$(NEXTFLOW) run $(PIPELINE) \
	    -r $(REVISION) \
	    -c nextflow.config \
	    -profile docker \
	    --input $(CLEAN_FASTA) \
	    --datadir $(DATA_DIR) \
	    --interpro $(INTERPRO_VERSION) \
	    --outdir $(OUTDIR) \
	    -work-dir $(WORK_DIR) \
	    $(APPLICATIONS_FLAG) \
	    -resume

clean:
	rm -rf $(WORK_DIR)/
	rm -f .nextflow.log*
	rm -f $(CLEAN_FASTA) $(MINI_FASTA)
	@echo "Cleaned (results preserved)"

clean-all: clean
	rm -rf $(OUTDIR)/ $(MINI_OUTDIR)/ $(TEST_OUTDIR)/
