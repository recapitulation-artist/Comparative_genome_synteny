RSCRIPT ?= Rscript
PYTHON ?= python3
BEDTOOLS ?= .tools/bedtools
BLASTN ?=
MAKEBLASTDB ?=
MAFFT ?=
IQTREE ?=
THREADS ?= 96
MODE ?= default
TARGET ?= verify
TARGET_GENE ?= Nramp
NEIGHBORS ?= 10
GENOME_DIR ?= data/genomes
GENOME_METADATA ?= data/genomes.tsv
GENE_ANNOTATIONS ?= data/genes.tsv
TREE_DIRECTORY = results/$(1)/$(if $(filter replay,$(MODE)),phylogeny_replay,phylogeny)
PREPARE_OPTIONS = --genome-dir "$(GENOME_DIR)" --genome-metadata "$(GENOME_METADATA)" --annotations "$(GENE_ANNOTATIONS)" --target-gene "$(TARGET_GENE)" --neighbors "$(NEIGHBORS)" --bedtools "$(BEDTOOLS)"

.PHONY: all setup setup-bedtools setup-blast setup-phylogeny setup-r check-data check-tools check-blast check-phylogeny check-r check-node check run verify submit plot plot-main verify-main extract extract-raw blast blast-raw phylogeny phylogeny-raw plot-reference plot-inferred plot-raw-inferred verify-raw
all: run

setup: setup-bedtools setup-blast setup-phylogeny setup-r
setup-bedtools:
	$(PYTHON) scripts/setup_bedtools.py --bedtools "$(BEDTOOLS)"
setup-blast:
	$(PYTHON) scripts/setup_blast.py
setup-phylogeny:
	$(PYTHON) scripts/setup_phylogeny.py
setup-r:
	$(RSCRIPT) --vanilla scripts/environment.R --install
check-data:
	sha256sum --check --quiet data/checksums.sha256
check-tools:
	$(PYTHON) scripts/setup_bedtools.py --check --bedtools "$(BEDTOOLS)"
check-blast:
	$(PYTHON) scripts/run_blast.py --check-tools --blastn "$(BLASTN)" --makeblastdb "$(MAKEBLASTDB)"
check-phylogeny:
	$(PYTHON) scripts/run_phylogeny.py --check-tools --mafft "$(MAFFT)" --iqtree "$(IQTREE)"
check-r:
	$(RSCRIPT) --vanilla scripts/environment.R
check-node:
	$(PYTHON) scripts/submit_workflow.py --check-node
check: check-data check-tools check-blast check-phylogeny check-r

# Main figure reads only the supplied plotting tables, tree and metadata.
plot-main: check-r
	$(RSCRIPT) --vanilla scripts/plot_synteny.R --entry main
verify-main: plot-main
	$(RSCRIPT) --vanilla scripts/verify_results.R --entry main

# Each reconstructed entry uses its own normalized inputs and recorded tree.
define reconstructed_entry
.PHONY: prepare-$(1) blast-$(1) phylogeny-$(1) plot-$(1) verify-$(1)
prepare-$(1): check-tools
	$$(PYTHON) scripts/prepare_entry.py --entry $(1) $$(PREPARE_OPTIONS) --output results/$(1)
blast-$(1): check-node prepare-$(1) check-blast
	$$(PYTHON) scripts/run_blast.py --input-dir results/$(1) --work-dir results/$(1)/blast --output results/$(1)/blast_hits.csv --threads "$$(THREADS)" --blastn "$$(BLASTN)" --makeblastdb "$$(MAKEBLASTDB)"
phylogeny-$(1): check-node prepare-$(1) check-phylogeny
	$$(PYTHON) scripts/infer_entry.py --entry $(1) --input-dir results/$(1) --mode "$$(MODE)" --output-root "$$(call TREE_DIRECTORY,$(1))" --mafft "$$(MAFFT)" --iqtree "$$(IQTREE)"
plot-$(1): blast-$(1) phylogeny-$(1) check-r
	$$(RSCRIPT) --vanilla scripts/plot_synteny.R --entry $(1) --target-gene "$$(TARGET_GENE)" --phylogeny-dir "$$(call TREE_DIRECTORY,$(1))"
verify-$(1): plot-$(1)
	$$(PYTHON) scripts/prepare_entry.py --entry $(1) $$(PREPARE_OPTIONS) --output results/$(1) --verify-only
	$$(RSCRIPT) --vanilla scripts/verify_results.R --entry $(1) --target-gene "$$(TARGET_GENE)" --phylogeny-dir "$$(call TREE_DIRECTORY,$(1))"
endef
$(eval $(call reconstructed_entry,curated))
$(eval $(call reconstructed_entry,auto))

run: plot-main plot-curated plot-auto
verify: verify-main verify-curated verify-auto

# Compatibility aliases resolve to the three supported entries.
plot plot-reference: plot-main
plot-inferred: plot-curated
plot-raw-inferred: plot-auto
verify-raw: verify-auto
extract: prepare-curated
extract-raw: prepare-auto
blast: blast-curated
blast-raw: blast-auto
phylogeny: phylogeny-curated
phylogeny-raw: phylogeny-auto

submit:
	$(PYTHON) scripts/submit_workflow.py --mode "$(MODE)" --target "$(TARGET)" --target-gene "$(TARGET_GENE)" --neighbors "$(NEIGHBORS)" --genome-dir "$(GENOME_DIR)" --genome-metadata "$(GENOME_METADATA)" --annotations "$(GENE_ANNOTATIONS)" --rscript "$(RSCRIPT)" --bedtools "$(BEDTOOLS)" --blastn "$(BLASTN)" --makeblastdb "$(MAKEBLASTDB)" --mafft "$(MAFFT)" --iqtree "$(IQTREE)"
