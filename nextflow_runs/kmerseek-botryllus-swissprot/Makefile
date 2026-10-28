# kmerseek search of the whole Botryllus proteome against UniProtKB/Swiss-Prot 2026_03,
# 19 alphabets x 2 k-mer sizes (the k carrying the bits of protein20 k = 6 and k = 7).
# Never run on the laptop: kmerseek and Nextflow run on Sherlock only.
#
#   make sherlock-sync     LAPTOP: push this branch, clone or pull it under $SCRATCH
#   ssh sherlock
#   cd $SCRATCH/botryllus-mhc-kmerseek-swissprot/nextflow_runs/kmerseek-botryllus-swissprot
#   make sherlock-inputs   LOGIN NODE: copy the query and Swiss-Prot FASTA into the run folder
#   make sherlock-stub     SLURM job: the workflow wiring, no kmerseek (minutes)
#   make sherlock-run      SLURM head job: the full run
#   make sherlock-status   LOGIN NODE: one look at the head job and its tasks (wait >= 60 s between looks)
#   make pull              LAPTOP: copy the summary tables and trace into results/

SHERLOCK_HOST   ?= sherlock
BRANCH          ?= $(shell git rev-parse --abbrev-ref HEAD)
REPO_URL        ?= https://github.com/olgabot/botryllus-mhc
SHERLOCK_CLONE  ?= $$SCRATCH/botryllus-mhc-kmerseek-swissprot
PIPELINE_SUBDIR := nextflow_runs/kmerseek-botryllus-swissprot

# Run folder: inputs, indexes, work and results. Outside the clone, so a fresh clone keeps them.
SH_BASE         ?= $(or $(SCRATCH),/scratch/users/olgabot)/botryllus-kmerseek-swissprot
QUERY_SRC       ?= $(SCRATCH)/botryllus-annotation-full/input/botryllus_unique_names.fasta
QUERY_MD5       := 414f4a938497854ef52a69d35e924a73
SWISSPROT_SRC   ?= $(SCRATCH)/botryllus-mhc-swissprot/db/swissprot/uniprot_sprot.fasta.gz
SWISSPROT_REL   ?= $(SCRATCH)/botryllus-mhc-swissprot/db/swissprot/reldate.txt

NXF_VERSION     ?= 25.10.4
SH_NEXTFLOW     ?= $(HOME)/software/nextflow/$(NXF_VERSION)/nextflow
SH_JAVA_HOME    ?= /share/software/user/open/java/18.0.2
SH_LOAD         := export JAVA_HOME=$(SH_JAVA_HOME) PATH=$(SH_JAVA_HOME)/bin:$$PATH NXF_OPTS="-Xms1g -Xmx6g";

# The head job only submits and watches tasks. 2 days is the normal QOS limit; the
# predicted run is under a day (README).
HEAD_TIME       ?= 2-00:00:00
HEAD_MEM        ?= 8G
# Extra sbatch flags, e.g. AFTER='--dependency=afterok:<stub job id>' to start only if the stub passed
AFTER           ?=

NF_ARGS = -profile sherlock \
	--query $(SH_BASE)/input/botryllus_unique_names.fasta \
	--swissprot_fasta $(SH_BASE)/input/uniprot_sprot_2026_03.fasta.gz \
	--outdir $(SH_BASE)/results --index_dir $(SH_BASE)/indexes \
	-work-dir $(SH_BASE)/work

.PHONY: sherlock-sync sherlock-inputs sherlock-stub sherlock-run sherlock-status pull plan

plan:
	python3 bin/plan_pairs.py

sherlock-sync:
	git push -u origin $(BRANCH)
	ssh $(SHERLOCK_HOST) 'set -e; if [ -d $(SHERLOCK_CLONE)/.git ]; then \
	    cd $(SHERLOCK_CLONE) && git fetch origin && git checkout $(BRANCH) && git reset --hard origin/$(BRANCH); \
	  else git clone -b $(BRANCH) $(REPO_URL) $(SHERLOCK_CLONE); fi; \
	  cd $(SHERLOCK_CLONE) && git log --oneline -1'

sherlock-inputs:
	@[ -z "$$SLURM_JOB_ID" ] || { echo "run on a login node"; exit 1; }
	mkdir -p $(SH_BASE)/input $(SH_BASE)/logs
	echo "$(QUERY_MD5)  $(QUERY_SRC)" | md5sum -c -
	install -m 444 $(QUERY_SRC) $(SH_BASE)/input/botryllus_unique_names.fasta
	grep -q 'Swiss-Prot Release 2026_03' $(SWISSPROT_REL)
	install -m 444 $(SWISSPROT_SRC) $(SH_BASE)/input/uniprot_sprot_2026_03.fasta.gz
	install -m 444 $(SWISSPROT_REL) $(SH_BASE)/input/uniprot_reldate.txt
	ls -la $(SH_BASE)/input

sherlock-stub:
	sbatch --job-name=kss-stub --partition=hns,normal --account=ayeletv \
	  --cpus-per-task=2 --mem=8G --time=02:00:00 --output=$(SH_BASE)/logs/stub-%j.out \
	  --wrap='$(SH_LOAD) set -euo pipefail; mkdir -p $(SH_BASE)/stub; cd $(SH_BASE)/stub; \
	    $(SH_NEXTFLOW) run $(CURDIR)/main.nf -stub -profile stub \
	      --query $(SH_BASE)/input/botryllus_unique_names.fasta \
	      --swissprot_fasta $(SH_BASE)/input/uniprot_sprot_2026_03.fasta.gz \
	      --outdir $(SH_BASE)/stub/results --index_dir $(SH_BASE)/stub/indexes -work-dir $(SH_BASE)/stub/work'

sherlock-run:
	sbatch --job-name=kss-head --partition=hns,normal --account=ayeletv $(AFTER) \
	  --cpus-per-task=2 --mem=$(HEAD_MEM) --time=$(HEAD_TIME) --output=$(SH_BASE)/logs/head-%j.out \
	  --wrap='$(SH_LOAD) set -euo pipefail; cd $(SH_BASE); \
	    $(SH_NEXTFLOW) run $(CURDIR)/main.nf $(NF_ARGS) -resume'

sherlock-status:
	squeue --me --name=kss-head -o '%i %j %T %M %l'
	@t=$(SH_BASE)/results/pipeline_info/trace.txt; [ -f $$t ] && \
	  awk -F'\t' 'NR>1 { split($$3, a, " "); n[a[1] " " $$5]++ } END { for (k in n) print n[k], k }' $$t | sort -k2 || true

pull:
	mkdir -p results
	rsync -av $(SHERLOCK_HOST):$(SH_BASE)/results/summary/ results/summary/
	rsync -av $(SHERLOCK_HOST):$(SH_BASE)/results/query/ results/query/
	rsync -av $(SHERLOCK_HOST):$(SH_BASE)/results/pipeline_info/trace.txt results/trace.txt
