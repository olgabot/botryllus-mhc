# Makefile to launch ebi-pf-team/interproscan6 on AWS Batch via Seqera Platform.
# https://github.com/ebi-pf-team/interproscan6
#
# Requires the Seqera CLI (`tw`) and your workspace/compute-env/S3 config.
# None of that config lives in this repo -- source your private env file
# first (see ~/.config/seqera/seanome-kmerseek.env):
#
#   source ~/.config/seqera/seanome-kmerseek.env
#   make install-tw   # first time only
#   make test         # smoke test on the pipeline's bundled data
#   make run          # full Botryllus proteome run
#
# Run from this directory: nextflow_runs/interproscan6-aws/

TW ?= tw

PIPELINE ?= https://github.com/ebi-pf-team/interproscan6
REVISION ?= 6.0.1

# --- Seqera Platform / AWS Batch (from the sourced private env file) ---
WORKSPACE ?= $(TOWER_WORKSPACE_ID)
COMPUTE_ENV ?= $(TOWER_COMPUTE_ENV)
S3_BUCKET ?= $(TOWER_S3_BUCKET)

S3_PREFIX ?= interproscan6
WORK_DIR ?= s3://$(S3_BUCKET)/$(S3_PREFIX)/work
# Persistent InterPro/member-DB cache, reused across runs so the multi-GB
# reference data isn't re-downloaded on every launch.
DATA_DIR ?= s3://$(S3_BUCKET)/$(S3_PREFIX)/data

# Nextflow config profiles to apply. Seqera sets the awsbatch executor
# automatically from the selected compute environment, so `-profile
# awsbatch` is not needed here. Leave empty unless you want e.g. `wave`.
PROFILE ?=
PROFILE_FLAG = $(if $(PROFILE),--profile $(PROFILE),)

# --- Input ---
FASTA_FILE ?= $(HOME)/data/botryllus/Bs_proteins.fa
CLEAN_FASTA ?= botryllus_clean.fa
SAMPLE_ID ?= botryllus
S3_INPUT ?= s3://$(S3_BUCKET)/$(S3_PREFIX)/input/$(CLEAN_FASTA)
OUTDIR ?= s3://$(S3_BUCKET)/$(S3_PREFIX)/results/$(SAMPLE_ID)

# --- Mini run: first 10 real Botryllus sequences, on AWS (unlike `test`,
# which uses the pipeline's own bundled data) ---
MINI_FASTA ?= botryllus_mini10.fa
MINI_SAMPLE_ID ?= botryllus_mini10
S3_MINI_INPUT ?= s3://$(S3_BUCKET)/$(S3_PREFIX)/input/$(MINI_FASTA)
MINI_OUTDIR ?= s3://$(S3_BUCKET)/$(S3_PREFIX)/results/$(MINI_SAMPLE_ID)

# --- Smoke test: pipeline's own bundled data (assets/test.faa) ---
TEST_OUTDIR ?= s3://$(S3_BUCKET)/$(S3_PREFIX)/results/test

# InterPro data release. Pinned per the README's reproducibility guidance
# (checked https://ftp.ebi.ac.uk/pub/databases/interpro/releases/ on 2026-07-28).
INTERPRO_VERSION ?= 109.0

# Comma-separated member DBs to run, e.g. Hamap,PANTHER,Pfam. Leave empty
# to run the pipeline's full default application set.
APPLICATIONS ?=

RUN_NAME ?= botryllus-interproscan6
LAST_RUN_ID_FILE = .last_run_id
PARAMS_FILE = params.$(SAMPLE_ID).yml
MINI_PARAMS_FILE = params.$(MINI_SAMPLE_ID).yml
TEST_PARAMS_FILE = params.test.yml

.PHONY: help check-env install-tw upload-input upload-mini-input test mini run resume status view clean clean-all

help:
	@echo "ebi-pf-team/interproscan6 on AWS Batch via Seqera Platform"
	@echo ""
	@echo "Targets:"
	@echo "  make install-tw  Install the tw (Seqera) CLI to ~/.local/bin"
	@echo "  make test        Smoke test on the pipeline's bundled data"
	@echo "  make mini        First 10 real Botryllus sequences, on AWS"
	@echo "  make run         Full run on the Botryllus proteome"
	@echo "  make resume      Relaunch the last run with -resume"
	@echo "  make status      Show status of the last launched run"
	@echo "  make view        Show full details of the last launched run"
	@echo "  make clean       Remove generated params file and run-id cache"
	@echo ""
	@echo "Variables:"
	@echo "  WORKSPACE=$(WORKSPACE)"
	@echo "  COMPUTE_ENV=$(COMPUTE_ENV)"
	@echo "  S3_BUCKET=$(S3_BUCKET)"
	@echo "  WORK_DIR=$(WORK_DIR)"
	@echo "  DATA_DIR=$(DATA_DIR)"
	@echo "  INTERPRO_VERSION=$(INTERPRO_VERSION)"
	@echo "  APPLICATIONS=$(APPLICATIONS)"

check-env:
	@command -v $(TW) >/dev/null 2>&1 || { echo "tw CLI not found -- run 'make install-tw'"; exit 1; }
	@if [ -z "$(TOWER_ACCESS_TOKEN)" ]; then echo "TOWER_ACCESS_TOKEN is not set -- source your private Seqera env file first"; exit 1; fi
	@if [ -z "$(WORKSPACE)" ]; then echo "WORKSPACE is not set -- source your private Seqera env file first"; exit 1; fi
	@if [ -z "$(COMPUTE_ENV)" ]; then echo "COMPUTE_ENV is not set -- source your private Seqera env file first"; exit 1; fi
	@if [ -z "$(S3_BUCKET)" ]; then echo "S3_BUCKET is not set -- source your private Seqera env file first"; exit 1; fi

install-tw:
	@echo "Installing tw CLI to $(HOME)/.local/bin/tw ..."
	@mkdir -p $(HOME)/.local/bin
	@os=$$(uname -s); arch=$$(uname -m); \
	if [ "$$os" = "Darwin" ] && [ "$$arch" = "arm64" ]; then asset=tw-osx-arm64; \
	elif [ "$$os" = "Darwin" ]; then asset=tw-osx-x86_64; \
	else asset=tw-linux-x86_64; fi; \
	curl -fsSL -o $(HOME)/.local/bin/tw "https://github.com/seqeralabs/tower-cli/releases/latest/download/$$asset"
	@chmod +x $(HOME)/.local/bin/tw
	@echo "Installed. Make sure $(HOME)/.local/bin is on your PATH."
	@$(HOME)/.local/bin/tw --version

# Strip stop codons (*) -- InterProScan rejects them
$(CLEAN_FASTA):
	@if [ ! -f "$(FASTA_FILE)" ]; then echo "Error: FASTA file '$(FASTA_FILE)' not found"; exit 1; fi
	@echo "Stripping stop codons from $(FASTA_FILE) -> $(CLEAN_FASTA)..."
	@awk '/^>/{print; next} {gsub(/\*/, ""); if (length > 0) print}' $(FASTA_FILE) > $(CLEAN_FASTA)
	@echo "Done: $$(grep -c '^>' $(CLEAN_FASTA)) sequences in $(CLEAN_FASTA)"

upload-input: check-env $(CLEAN_FASTA)
	aws s3 cp $(CLEAN_FASTA) $(S3_INPUT)

$(PARAMS_FILE): upload-input
	@{ \
	  echo 'input: "$(S3_INPUT)"'; \
	  echo 'datadir: "$(DATA_DIR)"'; \
	  echo 'outdir: "$(OUTDIR)"'; \
	  echo 'interpro: "$(INTERPRO_VERSION)"'; \
	  if [ -n "$(APPLICATIONS)" ]; then echo 'applications: "$(APPLICATIONS)"'; fi; \
	} > $(PARAMS_FILE)
	@echo "Wrote $(PARAMS_FILE):"
	@cat $(PARAMS_FILE)

# First 10 sequences, stop codons stripped
$(MINI_FASTA):
	@if [ ! -f "$(FASTA_FILE)" ]; then echo "Error: FASTA file '$(FASTA_FILE)' not found"; exit 1; fi
	@echo "Extracting first 10 proteins from $(FASTA_FILE) -> $(MINI_FASTA)..."
	@awk '/^>/{n++} n>10{exit} /^>/{print; next} {gsub(/\*/, ""); if (length > 0) print}' $(FASTA_FILE) > $(MINI_FASTA)
	@echo "Done: $$(grep -c '^>' $(MINI_FASTA)) sequences in $(MINI_FASTA)"

upload-mini-input: check-env $(MINI_FASTA)
	aws s3 cp $(MINI_FASTA) $(S3_MINI_INPUT)

$(MINI_PARAMS_FILE): upload-mini-input
	@{ \
	  echo 'input: "$(S3_MINI_INPUT)"'; \
	  echo 'datadir: "$(DATA_DIR)"'; \
	  echo 'outdir: "$(MINI_OUTDIR)"'; \
	  echo 'interpro: "$(INTERPRO_VERSION)"'; \
	  if [ -n "$(APPLICATIONS)" ]; then echo 'applications: "$(APPLICATIONS)"'; fi; \
	} > $(MINI_PARAMS_FILE)
	@echo "Wrote $(MINI_PARAMS_FILE):"
	@cat $(MINI_PARAMS_FILE)

# `-profile test` sets --input to the pipeline's bundled test.faa, but
# --datadir is still required, so we still need a params file.
$(TEST_PARAMS_FILE): check-env
	@{ \
	  echo 'datadir: "$(DATA_DIR)"'; \
	  echo 'outdir: "$(TEST_OUTDIR)"'; \
	  echo 'interpro: "$(INTERPRO_VERSION)"'; \
	} > $(TEST_PARAMS_FILE)
	@echo "Wrote $(TEST_PARAMS_FILE):"
	@cat $(TEST_PARAMS_FILE)

# Smoke test on the pipeline's bundled test FASTA -- no S3 upload needed
test: check-env $(TEST_PARAMS_FILE)
	$(TW) launch $(PIPELINE) \
	    -r $(REVISION) \
	    --workspace $(WORKSPACE) \
	    --compute-env "$(COMPUTE_ENV)" \
	    --work-dir $(WORK_DIR) \
	    --params-file $(TEST_PARAMS_FILE) \
	    --profile test \
	    --name $(RUN_NAME)-test \
	    --wait SUBMITTED \
	    2>&1 | tee /tmp/tw_launch_test.log
	@grep -oE 'Workflow [A-Za-z0-9]+' /tmp/tw_launch_test.log | awk '{print $$2}' > $(LAST_RUN_ID_FILE)
	@echo "Run ID: $$(cat $(LAST_RUN_ID_FILE))"

# First 10 real Botryllus sequences, on AWS -- quick sanity check before
# committing to a full run (uploads to S3, runs the real pipeline params)
mini: check-env $(MINI_PARAMS_FILE)
	$(TW) launch $(PIPELINE) \
	    -r $(REVISION) \
	    --workspace $(WORKSPACE) \
	    --compute-env "$(COMPUTE_ENV)" \
	    --work-dir $(WORK_DIR) \
	    --params-file $(MINI_PARAMS_FILE) \
	    $(PROFILE_FLAG) \
	    --name $(RUN_NAME)-mini \
	    --wait SUBMITTED \
	    2>&1 | tee /tmp/tw_launch_mini.log
	@grep -oE 'Workflow [A-Za-z0-9]+' /tmp/tw_launch_mini.log | awk '{print $$2}' > $(LAST_RUN_ID_FILE)
	@echo "Run ID: $$(cat $(LAST_RUN_ID_FILE))"

run: check-env $(PARAMS_FILE)
	$(TW) launch $(PIPELINE) \
	    -r $(REVISION) \
	    --workspace $(WORKSPACE) \
	    --compute-env "$(COMPUTE_ENV)" \
	    --work-dir $(WORK_DIR) \
	    --params-file $(PARAMS_FILE) \
	    $(PROFILE_FLAG) \
	    --name $(RUN_NAME) \
	    --wait SUBMITTED \
	    2>&1 | tee /tmp/tw_launch_run.log
	@grep -oE 'Workflow [A-Za-z0-9]+' /tmp/tw_launch_run.log | awk '{print $$2}' > $(LAST_RUN_ID_FILE)
	@echo "Run ID: $$(cat $(LAST_RUN_ID_FILE))"

resume: check-env
	@if [ ! -f "$(LAST_RUN_ID_FILE)" ]; then echo "No previous run recorded -- run 'make run' first"; exit 1; fi
	$(TW) runs relaunch -i $$(cat $(LAST_RUN_ID_FILE)) --workspace $(WORKSPACE) --wait SUBMITTED \
	    2>&1 | tee /tmp/tw_launch_resume.log
	@grep -oE 'Workflow [A-Za-z0-9]+' /tmp/tw_launch_resume.log | awk '{print $$2}' > $(LAST_RUN_ID_FILE)
	@echo "Run ID: $$(cat $(LAST_RUN_ID_FILE))"

status: check-env
	@if [ ! -f "$(LAST_RUN_ID_FILE)" ]; then echo "No previous run recorded"; exit 1; fi
	$(TW) runs view -i $$(cat $(LAST_RUN_ID_FILE)) --workspace $(WORKSPACE) --status

view: check-env
	@if [ ! -f "$(LAST_RUN_ID_FILE)" ]; then echo "No previous run recorded"; exit 1; fi
	$(TW) runs view -i $$(cat $(LAST_RUN_ID_FILE)) --workspace $(WORKSPACE)

clean:
	rm -f $(CLEAN_FASTA) $(MINI_FASTA) $(PARAMS_FILE) $(MINI_PARAMS_FILE) $(TEST_PARAMS_FILE) $(LAST_RUN_ID_FILE)
	@echo "Cleaned (S3 results and run history preserved)"

clean-all: clean
	@echo "Note: results live in $(OUTDIR) on S3 -- delete manually with 'aws s3 rm --recursive' if needed"
