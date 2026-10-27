# Timing pilot: how long would nf-core/proteinannotator, nf-core/proteinfold (ESMFold),
# InterProScan with every member database, and UniFIRE take on the 45,339-protein
# Botryllus proteome? Each runs on 500 proteins drawn to match the proteome's length
# spread, on Sherlock, and bin/scale_timings.py scales the measured times up.
#
#   make sample                   # LAPTOP: draw the 500 proteins (committed, so done once)
#   make sherlock-sync            # LAPTOP: push this branch, clone or pull it under $SCRATCH
#   ssh sherlock; cd $SCRATCH/botryllus-mhc-annotation-timing-pilot/nextflow_runs/timing-pilot
#   make stage                    # LOGIN NODE: ESMFold weights and image, UniFIRE lite image
#   make run                      # LOGIN NODE: submit every job, chained by dependencies
#   make running                  # what is queued or running
#   make collect                  # LAPTOP: copy the traces back and write the timing table
#
# Uses what ../nf-core-proteinannotator staged on Sherlock (`make sherlock-stage` there):
# Nextflow 25.10.4, the pipeline at its pinned commit, the HMM databases, the images, and
# the InterProScan 5.59-91.0 tarball.

SHELL := /bin/bash

SHERLOCK_HOST  ?= sherlock
REPO_URL       ?= https://github.com/olgabot/botryllus-mhc
BRANCH         ?= $(shell git rev-parse --abbrev-ref HEAD)
SHERLOCK_CLONE ?= $$SCRATCH/botryllus-mhc-$(notdir $(BRANCH))

PROTEOME   ?= $(HOME)/data/botryllus/botryllus_new_genome_protein.fasta
N_PROTEINS ?= 500
SAMPLE     := botryllus_pilot500.fasta
SAMPLE_OX  := botryllus_pilot500_unifire.fasta

# Sherlock side. PILOT holds every output; PA_* is what the proteinannotator folder staged.
PILOT        ?= $(SCRATCH)/botryllus-annotation-timing-pilot
PA_BASE      ?= $(SCRATCH)/botryllus-mhc-proteinannotator
PA_DB_DIR    := $(PA_BASE)/dbs
IPS_VERSION  ?= 5.59-91.0
IPS_DIR      := $(PA_DB_DIR)/interproscan-$(IPS_VERSION)
PA_IMG_DIR   ?= $(SCRATCH)/singularity-cache/proteinannotator
IPS_IMAGE    := $(PA_IMG_DIR)/depot.galaxyproject.org-singularity-interproscan-5.59_91.0--hec16e2b_1.img
NXF_VERSION  ?= 25.10.4
SH_NEXTFLOW  ?= $(HOME)/software/nextflow/$(NXF_VERSION)/nextflow
SH_JAVA_HOME ?= /share/software/user/open/java/18.0.2
SH_LOAD      := export JAVA_HOME=$(SH_JAVA_HOME) PATH=$(SH_JAVA_HOME)/bin:$$PATH;
PIPELINE_COMMIT ?= c70cd5849c80b5b1b7477a451a03ad108eb4191a
PA_CONFIG    := $(CURDIR)/../nf-core-proteinannotator/sherlock.config
# The six databases ../nf-core-proteinannotator runs; the full set is interproscan_full.sbatch
INTERPROSCAN_APPS ?= Hamap,PANTHER,PIRSF,TIGRFAM,sfld,Pfam
# 100 per batch gives 5 batches, so per-batch times can be compared
INTERPROSCAN_BATCH_SIZE ?= 100

ESMFOLD_IMAGE := docker://quay.io/nf-core/proteinfold_esmfold:2.0.0
ESM_URL       := https://dl.fbaipublicfiles.com/fair-esm
UNIFIRE_LITE  := docker://dockerhub.ebi.ac.uk/uniprot-public/unifire:2026.3-lite

CURL := curl -L --fail --continue-at - --speed-limit 10240 --speed-time 60 --retry 10
# $(1)=destination  $(2)=url. Downloads to .part and renames when whole.
GET = [ -s $(1) ] || { $(CURL) -o $(1).part $(2) && mv -f $(1).part $(1); }

.PHONY: help sample sherlock-sync check-login stage run running collect

help:
	@sed -n '1,20p' Makefile

sample:
	python3 -I bin/sample_by_length.py $(PROTEOME) $(N_PROTEINS) $(SAMPLE) botryllus_pilot500_lengths.tsv
	awk '/^>/{id=substr($$1,2); print ">" id "|" id " OS=Botryllus schlosseri OX=30301"; next} {print}' \
	  $(SAMPLE) > $(SAMPLE_OX)

sherlock-sync:
	git push origin $(BRANCH)
	ssh $(SHERLOCK_HOST) 'set -e; \
	  if [ ! -d $(SHERLOCK_CLONE)/.git ]; then git clone --branch $(BRANCH) $(REPO_URL) $(SHERLOCK_CLONE); \
	  else cd $(SHERLOCK_CLONE) && git fetch origin $(BRANCH) && git checkout $(BRANCH) && git pull --ff-only; fi'

check-login:
	@case "$$(hostname)" in *-ln*) ;; *) echo "run this on a Sherlock login node"; exit 1;; esac

# LOGIN NODE: compute nodes have no internet. About 8 GB of weights and 14 GB of images.
stage: check-login
	mkdir -p $(PILOT)/esmfold_params/checkpoints $(PILOT)/containers $(PILOT)/logs
	$(call GET,$(PILOT)/esmfold_params/checkpoints/esmfold_3B_v1.pt,$(ESM_URL)/models/esmfold_3B_v1.pt)
	$(call GET,$(PILOT)/esmfold_params/checkpoints/esm2_t36_3B_UR50D.pt,$(ESM_URL)/models/esm2_t36_3B_UR50D.pt)
	$(call GET,$(PILOT)/esmfold_params/checkpoints/esm2_t36_3B_UR50D-contact-regression.pt,$(ESM_URL)/regression/esm2_t36_3B_UR50D-contact-regression.pt)
	[ -s $(PILOT)/containers/proteinfold_esmfold_2.0.0.sif ] || \
	  apptainer pull $(PILOT)/containers/proteinfold_esmfold_2.0.0.sif $(ESMFOLD_IMAGE)
	[ -s $(PILOT)/containers/unifire-2026.3-lite.sif ] || \
	  apptainer pull $(PILOT)/containers/unifire-2026.3-lite.sif $(UNIFIRE_LITE)
	[ -s $(IPS_IMAGE) ] || { echo "missing $(IPS_IMAGE): run make sherlock-stage in ../nf-core-proteinannotator"; exit 1; }
	[ -s $(PA_DB_DIR)/interproscan-$(IPS_VERSION)-64-bit.tar.gz ] || { echo "missing the InterProScan tarball"; exit 1; }
	ls -la $(PILOT)/esmfold_params/checkpoints $(PILOT)/containers

# LOGIN NODE: submits only. Job ids go to $(PILOT)/jobs.tsv.
#   untar InterProScan 5.59 (skipped if done) -> proteinannotator head job
#                                             -> full InterProScan -> UniFIRE lite
#   ESMFold, independent
run: check-login
	mkdir -p $(PILOT)/logs $(PILOT)/proteinannotator
	printf 'id,fasta\npilot500,$(CURDIR)/$(SAMPLE)\n' > $(PILOT)/proteinannotator/samplesheet.csv
	set -euo pipefail; \
	dep=""; \
	if [ ! -e $(IPS_DIR)/.setup_done ]; then \
	  untar=$$(sbatch --parsable --job-name=ips-untar-5.59 --partition=hns,normal \
	    --cpus-per-task=2 --mem=8G --time=03:00:00 --output=$(PILOT)/logs/ips-untar-%j.out \
	    --wrap='set -eo pipefail; cd $(PA_DB_DIR); md5sum -c interproscan-$(IPS_VERSION)-64-bit.tar.gz.md5; \
	      [ -d $(IPS_DIR)/data ] || tar -pxzf interproscan-$(IPS_VERSION)-64-bit.tar.gz; \
	      cd $(IPS_DIR) && { python3 setup.py -f interproscan.properties || python3 setup.py interproscan.properties; } && touch .setup_done'); \
	  dep="--dependency=afterok:$$untar"; echo -e "ips-untar\t$$untar" >> $(PILOT)/jobs.tsv; \
	fi; \
	head=$$(sbatch --parsable $$dep --job-name=pilot-proteinannotator --partition=hns,normal \
	  --cpus-per-task=2 --mem=8G --time=2-00:00:00 --output=$(PILOT)/logs/head-%j.out \
	  --wrap='$(SH_LOAD) set -euo pipefail; cd $(PILOT)/proteinannotator; \
	    export NXF_OFFLINE=true NXF_OPTS="-Xms1g -Xmx6g" NXF_SINGULARITY_CACHEDIR=$(PA_IMG_DIR); \
	    $(SH_NEXTFLOW) run nf-core/proteinannotator -r $(PIPELINE_COMMIT) \
	      -profile singularity -c $(PA_CONFIG) -c $(CURDIR)/pilot.config \
	      --input samplesheet.csv --outdir $(PILOT)/proteinannotator/results \
	      -work-dir $(PILOT)/proteinannotator/work \
	      --pfam_db $(PA_DB_DIR)/Pfam-A.38.2.hmm.gz --funfam_db $(PA_DB_DIR)/funfam-hmm3-v4_3_0.lib.gz \
	      --nmpfams_db $(PA_DB_DIR)/nmpfamsdb.hmm.gz --metagroot_db $(PA_DB_DIR)/metagroot.hmm.gz \
	      --interproscan_db $(IPS_DIR)/data --interproscan_applications $(INTERPROSCAN_APPS) \
	      --interproscan_batch_size $(INTERPROSCAN_BATCH_SIZE)'); \
	echo -e "proteinannotator\t$$head" >> $(PILOT)/jobs.tsv; \
	ips=$$(sbatch --parsable $$dep --output=$(PILOT)/logs/ips-full-%j.out \
	  interproscan_full.sbatch $(PILOT) $(IPS_DIR)/data $(IPS_IMAGE) $(CURDIR)/$(SAMPLE_OX)); \
	echo -e "interproscan_full\t$$ips" >> $(PILOT)/jobs.tsv; \
	uf=$$(sbatch --parsable --dependency=afterok:$$ips --output=$(PILOT)/logs/unifire-%j.out \
	  unifire_lite.sbatch $(PILOT)); \
	echo -e "unifire\t$$uf" >> $(PILOT)/jobs.tsv; \
	esm=$$(sbatch --parsable --output=$(PILOT)/logs/esmfold-%j.out \
	  esmfold.sbatch $(PILOT) $(CURDIR)/$(SAMPLE)); \
	echo -e "esmfold\t$$esm" >> $(PILOT)/jobs.tsv; \
	cat $(PILOT)/jobs.tsv

running:
	@squeue --me -o "%.10i %.28j %.8T %.10M %.12l %R" | awk 'NR==1 || /pilot|ips-untar/'

# LAPTOP: the small files the scaling needs, then the table and figure
collect:
	mkdir -p results
	rsync -av --prune-empty-dirs --exclude='work/' --exclude='pdb/' --exclude='containers/' \
	  --exclude='esmfold_params/' --include='*/' \
	  --include='execution_trace_*.txt' --include='esm-fold.log' --include='gpu.txt' \
	  --include='*_epoch.txt' --include='jobs.tsv' --include='pilot500_full.tsv' \
	  --include='predictions*' --exclude='*' \
	  $(SHERLOCK_HOST):'$$SCRATCH/botryllus-annotation-timing-pilot/' results/raw/
	python3 -I bin/scale_timings.py results/raw botryllus_pilot500_lengths.tsv $(PROTEOME) results
