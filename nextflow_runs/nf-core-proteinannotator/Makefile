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

# InterProScan database (downloaded once and reused).
# If the path doesn't exist, omit the flag so the pipeline downloads it automatically.
INTERPROSCAN_DB ?= $(HOME)/data/interproscan
INTERPROSCAN_DB_FLAG = $(if $(wildcard $(INTERPROSCAN_DB)),--interproscan_db $(INTERPROSCAN_DB),)

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
	    $(INTERPROSCAN_DB_FLAG) \
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
	    $(INTERPROSCAN_DB_FLAG) \
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
	    $(INTERPROSCAN_DB_FLAG) \
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

# =======================================================================================
# Sherlock (Stanford HPC)
#
# Compute nodes have no internet, so everything the pipeline would download is staged on
# a login node first, and the head job runs with NXF_OFFLINE=true. Login-node targets only
# download files and submit SLURM jobs; unpacking and all compute run as jobs.
#
#   make sherlock-sync            # LAPTOP: push this branch, clone/pull it under $SCRATCH
#   make sherlock-upload-fasta    # LAPTOP: copy Bs_proteins.fa to Sherlock
#   ssh sherlock
#   cd $SCRATCH/botryllus-mhc-<branch>/nextflow_runs/nf-core-proteinannotator
#   make sherlock-stage           # LOGIN NODE: Nextflow, plugin, pipeline, HMM DBs, images
#   make sherlock-untar-ips       # SLURM job: unpack + index the InterProScan data
#   make sherlock-mini            # SLURM head job, first 10 proteins
#   make sherlock-run             # SLURM head job, full proteome
# =======================================================================================

SHERLOCK_HOST   ?= sherlock
REPO_URL        ?= https://github.com/olgabot/botryllus-mhc
BRANCH          ?= $(shell git rev-parse --abbrev-ref HEAD)
SHERLOCK_CLONE  ?= $$SCRATCH/botryllus-mhc-$(notdir $(BRANCH))

# Sherlock's nextflow module (25.04.7) is a self-contained build that ignores NXF_VER, and
# the pipeline needs >=25.10.0, so the 25.10.4 release file goes in $HOME/software.
NXF_VERSION     ?= 25.10.4
SH_NEXTFLOW     ?= $(HOME)/software/nextflow/$(NXF_VERSION)/nextflow
# Nextflow 25.10 needs Java 17 or newer. The path of Sherlock's java/18.0.2 module, set
# directly because `sbatch --wrap` runs /bin/sh, where Lmod's `ml` function may not exist.
SH_JAVA_HOME    ?= /share/software/user/open/java/18.0.2
SH_LOAD         := export JAVA_HOME=$(SH_JAVA_HOME) PATH=$(SH_JAVA_HOME)/bin:$$PATH;

# Pinned commit of the PR 79 branch (olgabot/interproscan-batches), so a rerun is the same code
PIPELINE_COMMIT ?= c70cd5849c80b5b1b7477a451a03ad108eb4191a
NF_SCHEMA       ?= nf-schema@2.6.1

SH_BASE         ?= $(SCRATCH)/botryllus-mhc-proteinannotator
SH_DB_DIR       := $(SH_BASE)/dbs
SH_INPUT_DIR    := $(SH_BASE)/input
SH_IMG_DIR      := $(SCRATCH)/singularity-cache/proteinannotator

PFAM_VERSION    ?= 38.2
PFAM_URL        := https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam$(PFAM_VERSION)/Pfam-A.hmm.gz
FUNFAM_URL      := https://download.cathdb.info/cath/releases/all-releases/v4_3_0/sequence-data/funfam-hmm3-v4_3_0.lib.gz
NMPFAMS_URL     := https://pavlopoulos-lab.org/envofams/databases/hmmer/nmpfamsdb.hmm.gz
METAGROOT_URL   := https://pavlopoulos-lab.org/envofams/databases/hmmer/metagroot.hmm.gz
# Must match the module's container (biocontainers interproscan:5.59_91.0). Its
# interproscan.properties names each member database by version (pfam/35.0, panther/17.0,
# tigrfam/15.0, ...), so the pipeline's default 5.72-103.0 data fails: on the Sherlock mini
# run every application but PIRSF and SFLD pointed at a folder that did not exist.
IPS_VERSION     ?= 5.59-91.0
IPS_URL         := https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/$(IPS_VERSION)/interproscan-$(IPS_VERSION)-64-bit.tar.gz
SH_IPS_DIR      := $(SH_DB_DIR)/interproscan-$(IPS_VERSION)

# Images the singularity profile asks for (grep of modules/**/main.nf at PIPELINE_COMMIT).
# Cached under the file name Nextflow looks for: the URL without https://, '/' and ':' -> '-'.
IMAGES := \
  https://depot.galaxyproject.org/singularity/hmmer:3.4--hdbdd923_1 \
  https://depot.galaxyproject.org/singularity/interproscan:5.59_91.0--hec16e2b_1 \
  https://depot.galaxyproject.org/singularity/s4pred:1.2.1--pyhdfd78af_1 \
  https://depot.galaxyproject.org/singularity/seqfu:1.22.3--hfd12232_2 \
  https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0 \
  https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1b/1bef8af6be88c5733461959c46ac8ef73d18f65277f62a1695d0e1633054f9c2/data \
  https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data \
  https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/95/95c0d3d867f5bc805b926b08ee761a993b24062739743eb82cc56363e0f7817d/data

# Speed floor + resume: a bare curl on a multi-GB file can stall for hours without failing
CURL := curl -L --fail --continue-at - --speed-limit 10240 --speed-time 60 --retry 10
# $(1)=file in SH_DB_DIR  $(2)=url. Downloads to .part and renames when complete, so a
# file that exists is whole and a rerun of sherlock-stage skips it.
GET = cd $(SH_DB_DIR) && { [ -s $(1) ] || { $(CURL) -o $(1).part $(2) && mv -f $(1).part $(1); }; }

# Head job: holds the Nextflow JVM for the whole run and only submits and polls jobs
HEAD_TIME ?= 2-00:00:00
HEAD_MEM  ?= 8G
SH_PROFILE_ARGS = -profile singularity -c $(CURDIR)/sherlock.config

DB_ARGS = --pfam_db $(SH_DB_DIR)/Pfam-A.$(PFAM_VERSION).hmm.gz \
          --funfam_db $(SH_DB_DIR)/funfam-hmm3-v4_3_0.lib.gz \
          --nmpfams_db $(SH_DB_DIR)/nmpfamsdb.hmm.gz \
          --metagroot_db $(SH_DB_DIR)/metagroot.hmm.gz \
          --interproscan_db $(SH_IPS_DIR)/data

.PHONY: sherlock-sync sherlock-upload-fasta sherlock-check-login sherlock-stage sherlock-untar-ips sherlock-mini sherlock-run

sherlock-sync:
	git push origin $(BRANCH)
	ssh $(SHERLOCK_HOST) 'set -e; \
	  if [ ! -d $(SHERLOCK_CLONE)/.git ]; then git clone --branch $(BRANCH) $(REPO_URL) $(SHERLOCK_CLONE); \
	  else cd $(SHERLOCK_CLONE) && git fetch origin $(BRANCH) && git checkout $(BRANCH) && git pull --ff-only; fi'

sherlock-upload-fasta:
	ssh $(SHERLOCK_HOST) 'mkdir -p $$SCRATCH/botryllus-mhc-proteinannotator/input'
	rsync -avP $(FASTA_FILE) $(SHERLOCK_HOST):'$$SCRATCH/botryllus-mhc-proteinannotator/input/'

sherlock-check-login:
	@case "$$(hostname)" in *-ln*) ;; *) echo "run this on a Sherlock login node"; exit 1;; esac

# LOGIN NODE. Downloads only (about 11 GB); nothing here computes.
sherlock-stage: sherlock-check-login
	mkdir -p $(dir $(SH_NEXTFLOW)) $(SH_DB_DIR) $(SH_IMG_DIR)
	[ -x $(SH_NEXTFLOW) ] || { $(CURL) -o $(SH_NEXTFLOW) \
	  https://github.com/nextflow-io/nextflow/releases/download/v$(NXF_VERSION)/nextflow-$(NXF_VERSION)-dist \
	  && chmod +x $(SH_NEXTFLOW); }
	$(SH_LOAD) $(SH_NEXTFLOW) -version | grep -q $(NXF_VERSION)
	$(SH_LOAD) $(SH_NEXTFLOW) plugin install $(NF_SCHEMA)
	$(SH_LOAD) $(SH_NEXTFLOW) pull nf-core/proteinannotator -r $(PIPELINE_COMMIT)
	$(call GET,Pfam-A.$(PFAM_VERSION).hmm.gz,$(PFAM_URL))
	$(call GET,funfam-hmm3-v4_3_0.lib.gz,$(FUNFAM_URL))
	$(call GET,nmpfamsdb.hmm.gz,$(NMPFAMS_URL))
	$(call GET,metagroot.hmm.gz,$(METAGROOT_URL))
	$(call GET,interproscan-$(IPS_VERSION)-64-bit.tar.gz,$(IPS_URL))
	$(call GET,interproscan-$(IPS_VERSION)-64-bit.tar.gz.md5,$(IPS_URL).md5)
	cd $(SH_DB_DIR) && md5sum -c interproscan-$(IPS_VERSION)-64-bit.tar.gz.md5
	for url in $(IMAGES); do \
	  name=$$(echo "$$url" | sed -E 's#^https?://##; s#[:/]#-#g').img; \
	  if [ ! -s $(SH_IMG_DIR)/$$name ]; then \
	    $(CURL) -o $(SH_IMG_DIR)/$$name.part "$$url" && mv -f $(SH_IMG_DIR)/$$name.part $(SH_IMG_DIR)/$$name || exit 1; \
	  fi; \
	done
	ls -la $(SH_IMG_DIR)
	@echo "staged. next: make sherlock-untar-ips"

# SLURM job: unpacking 6 GB into ~30 GB is heavy I/O, so it does not run on the login node.
# setup.py runs hmmpress on the member databases; without it PANTHER/Pfam/PIRSF/TIGRFAM fail.
sherlock-untar-ips:
	mkdir -p $(SH_BASE)/logs
	sbatch --wait --job-name=ips-untar --partition=hns,normal --account=ayeletv \
	  --cpus-per-task=2 --mem=8G --time=03:00:00 --output=$(SH_BASE)/logs/ips-untar-%j.out \
	  --wrap='set -eo pipefail; cd $(SH_DB_DIR); \
	    if [ ! -e $(SH_IPS_DIR)/.unpacked ]; then \
	      rm -rf $(SH_IPS_DIR) $(SH_IPS_DIR).tmp; mkdir $(SH_IPS_DIR).tmp; \
	      tar -pxzf interproscan-$(IPS_VERSION)-64-bit.tar.gz -C $(SH_IPS_DIR).tmp; \
	      mv $(SH_IPS_DIR).tmp/interproscan-$(IPS_VERSION) $(SH_IPS_DIR); rmdir $(SH_IPS_DIR).tmp; \
	      touch $(SH_IPS_DIR)/.unpacked; \
	    fi; \
	    cd $(SH_IPS_DIR); \
	    if [ -f setup.py ]; then python3 setup.py -f interproscan.properties; \
	    else python3 initial_setup.py; fi'
	ls $(SH_IPS_DIR)/data
	ls $(SH_IPS_DIR)/data/pfam/*/

# The module stages --interproscan_db as `data` and reads data/<db>/<version>/..., so it
# takes the tarball's data/ subfolder, not the tarball's top folder.

# First 10 proteins, stop codons stripped; the same command as the full run
$(SH_INPUT_DIR)/botryllus_clean.fa: $(SH_INPUT_DIR)/Bs_proteins.fa
	awk '/^>/{print; next} {gsub(/\*/, ""); if (length > 0) print}' $< > $@
$(SH_INPUT_DIR)/botryllus_mini10.fa: $(SH_INPUT_DIR)/Bs_proteins.fa
	awk '/^>/{n++} n>10{exit} /^>/{print; next} {gsub(/\*/, ""); if (length > 0) print}' $< > $@

# $(1)=id  $(2)=fasta  $(3)=outdir  $(4)=extra args
define SH_LAUNCH
	mkdir -p $(SH_BASE)/$(1) $(SH_BASE)/logs
	printf 'id,fasta\n$(1),$(2)\n' > $(SH_BASE)/$(1)/samplesheet.csv
	sbatch --job-name=pa-$(1) --partition=hns,normal --account=ayeletv \
	  --cpus-per-task=2 --mem=$(HEAD_MEM) --time=$(HEAD_TIME) \
	  --output=$(SH_BASE)/logs/head-$(1)-%j.out \
	  --wrap='$(SH_LOAD) set -euo pipefail; cd $(SH_BASE)/$(1); \
	    export NXF_OFFLINE=true NXF_OPTS="-Xms1g -Xmx6g" NXF_SINGULARITY_CACHEDIR=$(SH_IMG_DIR); \
	    $(SH_NEXTFLOW) run nf-core/proteinannotator -r $(PIPELINE_COMMIT) \
	      $(SH_PROFILE_ARGS) \
	      --input $(SH_BASE)/$(1)/samplesheet.csv --outdir $(3) -work-dir $(SH_BASE)/work \
	      $(DB_ARGS) --interproscan_applications $(INTERPROSCAN_APPS) \
	      --interproscan_batch_size $(INTERPROSCAN_BATCH_SIZE) $(4) -resume'
	@echo "head job submitted. log: $(SH_BASE)/logs/head-$(1)-<jobid>.out"
endef

sherlock-mini: sherlock-check-login $(SH_INPUT_DIR)/botryllus_mini10.fa
	$(call SH_LAUNCH,botryllus_mini10,$(SH_INPUT_DIR)/botryllus_mini10.fa,$(SH_BASE)/botryllus_mini10/results,)

sherlock-run: sherlock-check-login $(SH_INPUT_DIR)/botryllus_clean.fa
	$(call SH_LAUNCH,botryllus,$(SH_INPUT_DIR)/botryllus_clean.fa,$(SH_BASE)/botryllus/results,)
