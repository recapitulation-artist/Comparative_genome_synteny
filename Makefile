RSCRIPT ?= $(if $(wildcard .tools/r/bin/Rscript),.tools/r/bin/Rscript,Rscript)
PYTHON ?= python3
BEDTOOLS ?=
BLASTN ?=
MAKEBLASTDB ?=
MAFFT ?=
IQTREE ?=
THREADS ?= 56
MODE ?= default
TARGET ?= verify
TARGET_GENE ?= Nramp
NEIGHBORS ?= 10
GENOME_DIR ?= data/genomes
GENOME_METADATA ?= data/genomes.tsv
GENE_ANNOTATIONS ?= data/genes.tsv
GENE_PALETTE ?= main-muted
LINK_PALETTE ?= blue
LINK_COLORS ?=
SLURM_PARTITION ?=
SLURM_NODE ?=
SLURM_TIME ?=
SLURM_MEMORY ?=
SUBMIT_DRY_RUN ?= 0

# Quote values as single shell arguments, including paths with spaces or quotes.
shquote = '$(subst ','"'"',$(1))'
TREE_DIRECTORY = results/$(1)/$(if $(filter replay,$(MODE)),phylogeny_replay,phylogeny)
PREPARE_OPTIONS = --genome-dir=$(call shquote,$(GENOME_DIR)) --genome-metadata=$(call shquote,$(GENOME_METADATA)) --annotations=$(call shquote,$(GENE_ANNOTATIONS)) --target-gene=$(call shquote,$(TARGET_GENE)) --neighbors=$(call shquote,$(NEIGHBORS)) --bedtools=$(call shquote,$(BEDTOOLS)) --threads=$(call shquote,$(THREADS))
PLOT_OPTIONS = --gene-palette=$(call shquote,$(GENE_PALETTE)) --link-palette=$(call shquote,$(LINK_PALETTE)) --link-colors=$(call shquote,$(LINK_COLORS))

# Every stage owns the complete CPU budget, even when the caller uses make -j.
.NOTPARALLEL:
.PHONY: all setup setup-bedtools setup-blast setup-phylogeny setup-r check-data check-tools check-blast check-phylogeny check-r check run verify submit plot-main verify-main
all: run

setup: setup-bedtools setup-blast setup-phylogeny setup-r
setup-bedtools:
	$(call shquote,$(PYTHON)) scripts/setup_bedtools.py --bedtools=$(call shquote,$(BEDTOOLS))
setup-blast:
	$(call shquote,$(PYTHON)) scripts/setup_blast.py --blastn=$(call shquote,$(BLASTN)) --makeblastdb=$(call shquote,$(MAKEBLASTDB))
setup-phylogeny:
	$(call shquote,$(PYTHON)) scripts/setup_phylogeny.py --mafft=$(call shquote,$(MAFFT)) --iqtree=$(call shquote,$(IQTREE))
setup-r:
	$(call shquote,$(RSCRIPT)) --vanilla scripts/environment.R --install
check-data:
	sha256sum --check --quiet data/checksums.sha256
check-tools:
	$(call shquote,$(PYTHON)) scripts/setup_bedtools.py --check --bedtools=$(call shquote,$(BEDTOOLS))
check-blast:
	$(call shquote,$(PYTHON)) scripts/run_blast.py --check-tools --blastn=$(call shquote,$(BLASTN)) --makeblastdb=$(call shquote,$(MAKEBLASTDB))
check-phylogeny:
	$(call shquote,$(PYTHON)) scripts/run_phylogeny.py --check-tools --mafft=$(call shquote,$(MAFFT)) --iqtree=$(call shquote,$(IQTREE))
check-r:
	$(call shquote,$(RSCRIPT)) --vanilla scripts/environment.R
check: check-data check-tools check-blast check-phylogeny check-r

# The main figure needs only R and the bundled plotting inputs.
plot-main: check-r
	$(call shquote,$(RSCRIPT)) --vanilla scripts/plot_synteny.R --entry main $(PLOT_OPTIONS)
verify-main: plot-main
	$(call shquote,$(RSCRIPT)) --vanilla scripts/verify_results.R --entry main $(PLOT_OPTIONS)

# Reconstructed entries keep their inputs, intermediates and trees separate.
define reconstructed_entry
.PHONY: prepare-$(1) blast-$(1) phylogeny-$(1) plot-$(1) verify-$(1)
prepare-$(1): check-tools
	$$(call shquote,$$(PYTHON)) scripts/prepare_entry.py --entry $(1) $$(PREPARE_OPTIONS) --output results/$(1)
blast-$(1): prepare-$(1) check-blast
	$$(call shquote,$$(PYTHON)) scripts/run_blast.py --input-dir results/$(1) --work-dir results/$(1)/blast --output results/$(1)/blast_hits.csv --threads=$$(call shquote,$$(THREADS)) --blastn=$$(call shquote,$$(BLASTN)) --makeblastdb=$$(call shquote,$$(MAKEBLASTDB))
phylogeny-$(1): prepare-$(1) check-phylogeny
	$$(call shquote,$$(PYTHON)) scripts/infer_entry.py --entry $(1) --input-dir results/$(1) --mode=$$(call shquote,$$(MODE)) --output-root=$$(call shquote,$$(call TREE_DIRECTORY,$(1))) --mafft=$$(call shquote,$$(MAFFT)) --iqtree=$$(call shquote,$$(IQTREE)) --threads=$$(call shquote,$$(THREADS))
plot-$(1): blast-$(1) phylogeny-$(1) check-r
	$$(call shquote,$$(RSCRIPT)) --vanilla scripts/plot_synteny.R --entry $(1) --target-gene=$$(call shquote,$$(TARGET_GENE)) --phylogeny-dir=$$(call shquote,$$(call TREE_DIRECTORY,$(1))) $$(PLOT_OPTIONS)
verify-$(1): plot-$(1)
	$$(call shquote,$$(PYTHON)) scripts/prepare_entry.py --entry $(1) $$(PREPARE_OPTIONS) --output results/$(1) --verify-only
	$$(call shquote,$$(RSCRIPT)) --vanilla scripts/verify_results.R --entry $(1) --target-gene=$$(call shquote,$$(TARGET_GENE)) --phylogeny-dir=$$(call shquote,$$(call TREE_DIRECTORY,$(1))) $$(PLOT_OPTIONS)
endef
$(eval $(call reconstructed_entry,curated))
$(eval $(call reconstructed_entry,auto))

run: plot-main plot-curated plot-auto
verify: verify-main verify-curated verify-auto

submit:
	$(call shquote,$(PYTHON)) scripts/submit_workflow.py --mode=$(call shquote,$(MODE)) --target=$(call shquote,$(TARGET)) --threads=$(call shquote,$(THREADS)) --python=$(call shquote,$(PYTHON)) --target-gene=$(call shquote,$(TARGET_GENE)) --neighbors=$(call shquote,$(NEIGHBORS)) --genome-dir=$(call shquote,$(GENOME_DIR)) --genome-metadata=$(call shquote,$(GENOME_METADATA)) --annotations=$(call shquote,$(GENE_ANNOTATIONS)) --rscript=$(call shquote,$(RSCRIPT)) --bedtools=$(call shquote,$(BEDTOOLS)) --blastn=$(call shquote,$(BLASTN)) --makeblastdb=$(call shquote,$(MAKEBLASTDB)) --mafft=$(call shquote,$(MAFFT)) --iqtree=$(call shquote,$(IQTREE)) --partition=$(call shquote,$(SLURM_PARTITION)) --node=$(call shquote,$(SLURM_NODE)) --time=$(call shquote,$(SLURM_TIME)) --memory=$(call shquote,$(SLURM_MEMORY)) $(PLOT_OPTIONS) $(if $(filter 1,$(SUBMIT_DRY_RUN)),--dry-run,)
