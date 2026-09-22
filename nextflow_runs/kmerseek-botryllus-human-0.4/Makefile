# Botryllus schlosseri proteins (query) vs human GENCODE v49 canonical
# translations (target), searched with kmerseek 0.4.0.
#
#   make run       # hp_lehninger2, k=24, ungapped extension, E-values
#   make resume    # same, with -resume
#
# Run from this directory.

NEXTFLOW ?= /Users/olga/anaconda3/envs/nf-core-v2/bin/nextflow
SCRIPT   ?= kmerseek_botryllus_human_0.4.nf

# nf-core-v2 keeps its JVM under lib/jvm, which the plain shell does not find.
export NXF_JAVA_HOME ?= /Users/olga/anaconda3/envs/nf-core-v2/lib/jvm

# kmerseek 0.4.0: the release-0.4.0 branch. 0.4 refuses to open a pre-0.4 index,
# so this run builds its own and shares nothing with the 0.3 pipeline.
KMERSEEK_SRC      ?= $(HOME)/code/kmerseek/.claude/worktrees/release-0.4.0
KMERSEEK_MAIN_SRC ?= $(HOME)/code/kmerseek/.claude/worktrees/release-0.4.0

ALPHABET ?= hp_lehninger2
KSIZE    ?= 24

# Mismatch penalty C and X-drop for the ungapped extension, and how many database
# sequences the Karlin-Altschul fit behind region_evalue is read off.
EXTEND_MISMATCH_PENALTY ?= 1.63
EXTEND_XDROP            ?= 6.52
KA_QUERIES              ?= 200

# A hit is kept if either the whole-query p-value or the best region's score passes.
MAX_QUERY_PVALUE ?= 0.05
MIN_REGION_SCORE ?= 1.3

INPUT  ?= $(HOME)/data/gencode/human/v49/gencode.v49.pc_translations.canonical.fa
QUERY  ?= $(HOME)/data/botryllus/Bs_proteins.fa
OUTDIR ?= $(HOME)/data/botryllus/kmerseek-0.4

VISUALIZE_QUERY_NAME ?= BHF
VISUALIZE_MAX_HITS   ?= 20

WORK_DIR ?= work

NF_ARGS = -c nextflow.config \
	    --input $(INPUT) \
	    --query $(QUERY) \
	    --outdir $(OUTDIR) \
	    --alphabet $(ALPHABET) \
	    --ksize $(KSIZE) \
	    --kmerseek_src $(KMERSEEK_SRC) \
	    --kmerseek_main_src $(KMERSEEK_MAIN_SRC) \
	    --extend_mismatch_penalty $(EXTEND_MISMATCH_PENALTY) \
	    --extend_xdrop $(EXTEND_XDROP) \
	    --ka_queries $(KA_QUERIES) \
	    --max_query_pvalue $(MAX_QUERY_PVALUE) \
	    --min_region_score $(MIN_REGION_SCORE) \
	    --visualize_query_name $(VISUALIZE_QUERY_NAME) \
	    --visualize_max_hits $(VISUALIZE_MAX_HITS) \
	    -work-dir $(WORK_DIR)

.PHONY: help check-nextflow run resume clean

help:
	@echo "kmerseek 0.4.0: Botryllus vs human GENCODE"
	@echo ""
	@echo "  make run      Build, index (with the Karlin-Altschul fit), search, filter, plot"
	@echo "  make resume   Same, with -resume"
	@echo "  make clean    Remove the work dir and logs (published outputs are kept)"
	@echo ""
	@echo "  ALPHABET=$(ALPHABET)  KSIZE=$(KSIZE)"
	@echo "  EXTEND_MISMATCH_PENALTY=$(EXTEND_MISMATCH_PENALTY)  EXTEND_XDROP=$(EXTEND_XDROP)"
	@echo "  QUERY=$(QUERY)"
	@echo "  INPUT=$(INPUT)"
	@echo "  OUTDIR=$(OUTDIR)"

check-nextflow:
	@command -v $(NEXTFLOW) >/dev/null 2>&1 || { echo "nextflow not found at '$(NEXTFLOW)'"; exit 1; }

run: check-nextflow
	$(NEXTFLOW) run $(SCRIPT) $(NF_ARGS)

resume: check-nextflow
	$(NEXTFLOW) run $(SCRIPT) $(NF_ARGS) -resume

clean:
	/bin/rm -rf $(WORK_DIR)/
	/bin/rm -f .nextflow.log*
	@echo "Cleaned (published outputs in $(OUTDIR) are kept)"
