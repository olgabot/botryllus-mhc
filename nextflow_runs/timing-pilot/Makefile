# Timing pilot: how long would nf-core/proteinannotator (with InterProScan's every member
# database), UniFIRE, nf-core/proteinfold (ESMFold) and Foldseek take on the 45,339-protein
# Botryllus proteome? Each runs on 500 proteins drawn to match the proteome's length
# spread, on Sherlock, and bin/scale_timings.py scales the measured times up.
#
#   make sample                   # LAPTOP: draw the 500 proteins (committed, so done once)
#   make sherlock-sync            # LAPTOP: push this branch, clone or pull it under $SCRATCH
#   ssh sherlock; cd $SCRATCH/botryllus-mhc-annotation-timing-pilot/nextflow_runs/timing-pilot
#   make stage                    # LOGIN NODE: ESMFold weights and image, UniFIRE lite image,
#                                 #   Foldseek image and AlphaFold DB Swiss-Prot
#   make run                      # LOGIN NODE: one Nextflow run (main.nf); head job in tmux
#   make watch                    # the live Nextflow task table (Ctrl-b d to leave it)
#   make pa-progress              # proteinannotator's own progress (a run inside the run)
#   make running                  # what is queued or running
#   make collect                  # LAPTOP: copy the traces back and write the timing table
#   make full-input               # LAPTOP: the renamed proteome to Sherlock
#   make run-full                 # LOGIN NODE: proteinannotator + UniFIRE on all 45,339 proteins
#   make watch FULL=1             # (likewise stop, pa-progress) for that run
#
# Uses what ../nf-core-proteinannotator staged on Sherlock (`make sherlock-stage` there):
# Nextflow 25.10.4, the pipeline at its pinned commit, the HMM databases, the images, and
# the InterProScan 5.59-91.0 tarball.

SHELL := /bin/bash
PYTHON ?= /Users/olga/anaconda3/envs/2025-kmerseek-analysis/bin/python3

SHERLOCK_HOST  ?= sherlock
REPO_URL       ?= https://github.com/olgabot/botryllus-mhc
BRANCH         ?= $(shell git rev-parse --abbrev-ref HEAD)
SHERLOCK_CLONE ?= $$SCRATCH/botryllus-mhc-$(notdir $(BRANCH))

PROTEOME   ?= $(HOME)/data/botryllus/botryllus_new_genome_protein.fasta
PROTEOME_UNIQUE ?= $(HOME)/data/botryllus/botryllus_new_genome_protein_unique_names.fasta
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
# The pipeline commit, databases, batch size and image paths main.nf uses are in nextflow.config

ESMFOLD_IMAGE := docker://quay.io/nf-core/proteinfold_esmfold:2.0.0
ESM_URL       := https://dl.fbaipublicfiles.com/fair-esm
UNIFIRE_LITE  := docker://dockerhub.ebi.ac.uk/uniprot-public/unifire:2026.3-lite
# The Foldseek nf-core/proteinfold 2.0.0 pins; a ready .sif, so no build step
FOLDSEEK_SIF  := https://depot.galaxyproject.org/singularity/foldseek:9.427df8a--pl5321hb365157_0

CURL := curl -L --fail --continue-at - --speed-limit 10240 --speed-time 60 --retry 10
# $(1)=destination  $(2)=url. Downloads to .part and renames when whole.
GET = [ -s $(1) ] || { $(CURL) -o $(1).part $(2) && mv -f $(1).part $(1); }

.PHONY: help sample sherlock-sync check-login stage stub sif run watch pa-progress stop running collect full-input run-full _run-full

help:
	@sed -n '1,20p' Makefile

# The proteome reuses names; proteinannotator would keep one protein per name. So every
# step reads this copy, where each protein has its own name. The full run uses it too.
$(PROTEOME_UNIQUE): $(PROTEOME) bin/unique_names.py
	python3 -I bin/unique_names.py $(PROTEOME) $@ botryllus_renamed_proteins.tsv

sample: $(PROTEOME_UNIQUE)
	python3 -I bin/sample_by_length.py $(PROTEOME_UNIQUE) $(N_PROTEINS) $(SAMPLE) botryllus_pilot500_lengths.tsv
	awk '/^>/{id=substr($$1,2); print ">" id "|" id " OS=Botryllus schlosseri OX=30301"; next} {print}' \
	  $(SAMPLE) > $(SAMPLE_OX)

sherlock-sync:
	git push origin $(BRANCH)
	ssh $(SHERLOCK_HOST) 'set -e; \
	  if [ ! -d $(SHERLOCK_CLONE)/.git ]; then git clone --branch $(BRANCH) $(REPO_URL) $(SHERLOCK_CLONE); \
	  else cd $(SHERLOCK_CLONE) && git fetch origin $(BRANCH) && git checkout $(BRANCH) && git pull --ff-only; fi'
	# The repo ignores FASTA files (.gitignore), so the two sample FASTAs go over by rsync
	rsync -av $(SAMPLE) $(SAMPLE_OX) $(SHERLOCK_HOST):'$(SHERLOCK_CLONE)/nextflow_runs/timing-pilot/'

check-login:
	@case "$$(hostname)" in *-ln*) ;; *) echo "run this on a Sherlock login node"; exit 1;; esac

# LOGIN NODE: about 8 GB of weights, 14 GB of images and a 3.8 GB Foldseek database.
stage: check-login
	mkdir -p $(PILOT)/esmfold_params/checkpoints $(PILOT)/containers $(PILOT)/logs
	$(call GET,$(PILOT)/esmfold_params/checkpoints/esmfold_3B_v1.pt,$(ESM_URL)/models/esmfold_3B_v1.pt)
	$(call GET,$(PILOT)/esmfold_params/checkpoints/esm2_t36_3B_UR50D.pt,$(ESM_URL)/models/esm2_t36_3B_UR50D.pt)
	$(call GET,$(PILOT)/esmfold_params/checkpoints/esm2_t36_3B_UR50D-contact-regression.pt,$(ESM_URL)/regression/esm2_t36_3B_UR50D-contact-regression.pt)
	# `apptainer pull` was killed on the login node (exit 137 in mksquashfs, 2026-10-27), so
	# the login node only unpacks each image into a folder; `make run` squashes it in a job.
	[ -s $(PILOT)/containers/proteinfold_esmfold_2.0.0.sif ] || [ -d $(PILOT)/containers/proteinfold_esmfold_2.0.0 ] || \
	  apptainer build --sandbox $(PILOT)/containers/proteinfold_esmfold_2.0.0 $(ESMFOLD_IMAGE)
	[ -s $(PILOT)/containers/unifire-2026.3-lite.sif ] || [ -d $(PILOT)/containers/unifire-2026.3-lite ] || \
	  apptainer build --sandbox $(PILOT)/containers/unifire-2026.3-lite $(UNIFIRE_LITE)
	$(call GET,$(PILOT)/containers/foldseek-9.427df8a.sif,$(FOLDSEEK_SIF))
	# AlphaFold DB Swiss-Prot as a Foldseek database: a download and an untar, one thread
	[ -s $(PILOT)/foldseek_db/afdb_swissprot.dbtype ] || { mkdir -p $(PILOT)/foldseek_db && \
	  apptainer exec --bind $(SCRATCH) $(PILOT)/containers/foldseek-9.427df8a.sif foldseek databases \
	    Alphafold/Swiss-Prot $(PILOT)/foldseek_db/afdb_swissprot $(PILOT)/foldseek_db/tmp --threads 1; }
	[ -s $(IPS_IMAGE) ] || { echo "missing $(IPS_IMAGE): run make sherlock-stage in ../nf-core-proteinannotator"; exit 1; }
	[ -s $(PA_DB_DIR)/interproscan-$(IPS_VERSION)-64-bit.tar.gz ] || { echo "missing the InterProScan tarball"; exit 1; }
	ls -la $(PILOT)/esmfold_params/checkpoints $(PILOT)/containers $(PILOT)/foldseek_db

# Head job: runs `nextflow run main.nf` in tmux on a compute node (bin/nf-head, copied from
# 2024-kmerseek-analysis/scripts/nf-head), so the run survives logging out and `make watch`
# attaches to the live task table. 2 days is the normal QOS cap (a 3-day ask was refused
# 2026-10-27, QOSMaxWallDurationPerJobLimit); PROTEINANNOTATOR gets 46 h inside it.
SLURM_ACCOUNT ?= ayeletv
HEAD_PARTITION ?= hns,normal
HEAD_QOS     ?= normal
HEAD_TIME    ?= 2-00:00:00
HEAD_CPUS    ?= 2
HEAD_MEM     ?= 8G
HEAD_SESSION ?= timing-pilot
HEAD_JOBNAME ?= nfdrv-timing-pilot
# `make run-full` and FULL=1 for watch, stop and pa-progress: the whole proteome. Without
# s4pred it fits the same 2-day head job as the pilot (see full.config).
FULL_DIR := $(SCRATCH)/botryllus-annotation-full
ifeq ($(FULL),1)
PILOT        := $(FULL_DIR)
HEAD_SESSION := annotation-full
HEAD_JOBNAME := nfdrv-annotation-full
endif
HEAD_SUBMIT = HEAD_ACCOUNT=$(SLURM_ACCOUNT) HEAD_PARTITION=$(HEAD_PARTITION) \
              HEAD_QOS=$(HEAD_QOS) HEAD_TIME=$(HEAD_TIME) \
              HEAD_CPUS=$(HEAD_CPUS) HEAD_MEM=$(HEAD_MEM) \
              HEAD_SESSION=$(HEAD_SESSION) HEAD_JOBNAME=$(HEAD_JOBNAME) \
              $(CURDIR)/bin/nf-head

# LOGIN NODE: one Nextflow run. The tasks and their order are in main.nf; sizes are there too.
run: check-login
	@squeue --me -h -n $(HEAD_JOBNAME) | grep -q . && { echo "$(HEAD_JOBNAME) is already queued"; exit 1; } || true
	[ -s $(PILOT)/containers/proteinfold_esmfold_2.0.0.sif ] && [ -s $(PILOT)/containers/unifire-2026.3-lite.sif ] \
	  && [ -s $(PILOT)/containers/foldseek-9.427df8a.sif ] && [ -s $(PILOT)/foldseek_db/afdb_swissprot.dbtype ] \
	  || { echo "images missing: run make stage, then make sif"; exit 1; }
	mkdir -p $(PILOT)/launch
	cd $(PILOT)/launch && export JAVA_HOME=$(SH_JAVA_HOME) PATH=$(SH_JAVA_HOME)/bin:$$PATH && \
	  $(HEAD_SUBMIT) $(SH_NEXTFLOW) run $(CURDIR)/main.nf -resume

# LAPTOP: the renamed proteome, plain and with OX=30301 headers for UniFIRE, to Sherlock
full-input: $(PROTEOME_UNIQUE)
	awk '/^>/{id=substr($$1,2); print ">" id "|" id " OS=Botryllus schlosseri OX=30301"; next} {print}' \
	  $(PROTEOME_UNIQUE) > $(basename $(PROTEOME_UNIQUE))_unifire.fasta
	ssh $(SHERLOCK_HOST) 'mkdir -p $$SCRATCH/botryllus-annotation-full/input'
	rsync -av $(PROTEOME_UNIQUE) $(SHERLOCK_HOST):'$$SCRATCH/botryllus-annotation-full/input/botryllus_unique_names.fasta'
	rsync -av $(basename $(PROTEOME_UNIQUE))_unifire.fasta \
	  $(SHERLOCK_HOST):'$$SCRATCH/botryllus-annotation-full/input/botryllus_unique_names_unifire.fasta'

# LOGIN NODE: proteinannotator and UniFIRE on all 45,339 proteins, no ESMFold (full.config)
run-full: check-login
	$(MAKE) --no-print-directory _run-full FULL=1
_run-full:
	@squeue --me -h -n $(HEAD_JOBNAME) | grep -q . && { echo "$(HEAD_JOBNAME) is already queued"; exit 1; } || true
	[ -s $(FULL_DIR)/input/botryllus_unique_names_unifire.fasta ] || { echo "run make full-input on the laptop"; exit 1; }
	mkdir -p $(FULL_DIR)/launch
	cd $(FULL_DIR)/launch && export JAVA_HOME=$(SH_JAVA_HOME) PATH=$(SH_JAVA_HOME)/bin:$$PATH && \
	  $(HEAD_SUBMIT) $(SH_NEXTFLOW) run $(CURDIR)/main.nf -c $(CURDIR)/full.config -resume

# LOGIN NODE: main.nf's stub blocks in a 30-minute job on one node (stub_test.config).
# Checks the wiring, the input paths and that each image starts, without real work.
stub: check-login
	mkdir -p $(PILOT)/stub-test/launch $(PILOT)/logs
	cd $(PILOT)/stub-test/launch && srun --job-name=pilot-stub --partition=hns,normal --cpus-per-task=4 \
	  --mem=16G --time=00:30:00 bash -c 'export JAVA_HOME=$(SH_JAVA_HOME) PATH=$(SH_JAVA_HOME)/bin:$$PATH; \
	  $(SH_NEXTFLOW) run $(CURDIR)/main.nf -stub -ansi-log false -c $(CURDIR)/stub_test.config'

# LOGIN NODE: squash the unpacked image folders from `make stage` into .sif files, in a job
# (mksquashfs is killed on login nodes). Already done on 2026-10-27; kept for a fresh setup.
sif: check-login
	sbatch --wait --job-name=pilot-sif-build --partition=hns,normal --cpus-per-task=4 --mem=16G \
	  --time=02:00:00 --output=$(PILOT)/logs/sif-build-%j.out \
	  --wrap='set -euo pipefail; cd $(PILOT)/containers; \
	    for img in proteinfold_esmfold_2.0.0 unifire-2026.3-lite; do \
	      [ -s $$img.sif ] || apptainer build $$img.sif $$img/; done; ls -la'

# The live Nextflow table, read-only. Leave with Ctrl-b d; the run keeps going.
watch:
	@jid=$$(squeue -u $$USER -h -n $(HEAD_JOBNAME) -o "%i %T" | awk '$$2=="RUNNING"{print $$1; exit}'); \
	 [ -n "$$jid" ] || { echo "no $(HEAD_JOBNAME) job is running"; squeue -u $$USER -h -n $(HEAD_JOBNAME); exit 1; }; \
	 exec srun --jobid=$$jid --overlap --pty tmux -S /tmp/nf-head-$$jid attach -r -t $(HEAD_SESSION)

# proteinannotator is a Nextflow run of its own inside the PROTEINANNOTATOR task, so its
# tasks are not in the table above. This prints its latest progress lines once.
pa-progress:
	@f=$(PILOT)/proteinannotator/nextflow.log; [ -s $$f ] || { echo "proteinannotator has not started"; exit 0; }; \
	 grep -E "Submitted process|Task completed|ERROR|WARN" $$f | sed 's/^.*\] //' | tail -n 25

# Ctrl-C the run, so Nextflow cancels its own jobs
stop:
	@jid=$$(squeue -u $$USER -h -n $(HEAD_JOBNAME) -o "%i %T" | awk '$$2=="RUNNING"{print $$1; exit}'); \
	 [ -n "$$jid" ] || { echo "no $(HEAD_JOBNAME) job is running"; exit 1; }; \
	 scancel --signal=TERM --batch $$jid && echo "asked job $$jid to shut down"

running:
	@squeue --me -o "%.10i %.28j %.8T %.10M %.12l %R" | awk 'NR==1 || /nf-|nfdrv|pilot/'

# LAPTOP: the small files the scaling needs, then the table and figure
collect:
	mkdir -p results
	rsync -av --prune-empty-dirs --include='*/' \
	  --include='pilot_trace.txt' --include='execution_trace_*.txt' --include='esm-fold.log' \
	  --include='gpu.txt' --include='predictions*' --include='unifire.log' \
	  --include='*.m8' --include='foldseek.log' \
	  --exclude='*' \
	  $(SHERLOCK_HOST):'$$SCRATCH/botryllus-annotation-timing-pilot/results/' results/raw/
	rsync -av $(SHERLOCK_HOST):'$$SCRATCH/botryllus-annotation-timing-pilot/pipeline_info/pilot_trace.txt' results/raw/
	$(PYTHON) -I bin/scale_timings.py results/raw botryllus_pilot500_lengths.tsv $(PROTEOME_UNIQUE) results
