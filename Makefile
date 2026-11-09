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
PHYLOGENY_DIR ?= $(if $(filter replay,$(MODE)),results/phylogeny_replay,results/phylogeny)
RAW_PHYLOGENY_DIR ?= $(if $(filter replay,$(MODE)),results/raw/phylogeny_replay,results/raw/phylogeny)

.PHONY: all setup setup-bedtools setup-blast setup-phylogeny setup-r check-data check-tools check-blast check-phylogeny check-r check-node check extract extract-raw blast blast-raw phylogeny phylogeny-raw record-phylogeny plot plot-reference plot-inferred plot-raw-inferred run verify verify-raw submit

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

extract: check-data check-tools
	$(PYTHON) scripts/extract_regions.py --bedtools "$(BEDTOOLS)"

extract-raw: check-tools
	$(PYTHON) scripts/extract_raw_regions.py --check-inputs
	$(PYTHON) scripts/extract_raw_regions.py --bedtools "$(BEDTOOLS)"

blast: check-node extract check-blast
	$(PYTHON) scripts/run_blast.py --threads "$(THREADS)" --blastn "$(BLASTN)" --makeblastdb "$(MAKEBLASTDB)"

blast-raw: check-node extract-raw check-blast
	$(PYTHON) scripts/run_blast.py --input-dir results/raw --work-dir results/raw/blast --output results/raw/blast_hits.csv --threads "$(THREADS)" --blastn "$(BLASTN)" --makeblastdb "$(MAKEBLASTDB)"

phylogeny: check-node extract check-phylogeny check-r
	$(PYTHON) scripts/run_phylogeny.py --mode "$(MODE)" --output-root "$(PHYLOGENY_DIR)" --mafft "$(MAFFT)" --iqtree "$(IQTREE)"
	$(if $(filter replay,$(MODE)),$(PYTHON) scripts/record_phylogeny.py --prepare-replay --phylogeny-dir "$(PHYLOGENY_DIR)")
	$(RSCRIPT) --vanilla scripts/report_phylogeny.R --phylogeny-dir "$(PHYLOGENY_DIR)"

phylogeny-raw: check-node extract-raw check-phylogeny
	$(PYTHON) scripts/run_raw_phylogeny.py --mode "$(MODE)" --output-root "$(RAW_PHYLOGENY_DIR)" --mafft "$(MAFFT)" --iqtree "$(IQTREE)"

record-phylogeny:
	$(PYTHON) scripts/record_phylogeny.py --phylogeny-dir "$(PHYLOGENY_DIR)"

plot-reference: blast check-r
	$(RSCRIPT) --vanilla scripts/plot_synteny.R --tree-source reference

plot-inferred: blast phylogeny check-r
	$(RSCRIPT) --vanilla scripts/plot_synteny.R --tree-source inferred --phylogeny-dir "$(PHYLOGENY_DIR)"

plot-raw-inferred: blast-raw phylogeny-raw check-r
	$(RSCRIPT) --vanilla scripts/plot_synteny.R --tree-source inferred --annotation-source raw --phylogeny-dir "$(RAW_PHYLOGENY_DIR)"

plot: plot-reference

run: plot-reference plot-inferred plot-raw-inferred

verify: run
	$(PYTHON) scripts/extract_raw_regions.py --bedtools "$(BEDTOOLS)" --verify-only
	$(RSCRIPT) --vanilla scripts/verify_results.R --tree-source reference
	$(RSCRIPT) --vanilla scripts/verify_results.R --tree-source inferred --phylogeny-dir "$(PHYLOGENY_DIR)"
	$(RSCRIPT) --vanilla scripts/verify_results.R --tree-source inferred --annotation-source raw --phylogeny-dir "$(RAW_PHYLOGENY_DIR)"

verify-raw: plot-raw-inferred
	$(PYTHON) scripts/extract_raw_regions.py --bedtools "$(BEDTOOLS)" --verify-only
	$(RSCRIPT) --vanilla scripts/verify_results.R --tree-source inferred --annotation-source raw --phylogeny-dir "$(RAW_PHYLOGENY_DIR)"

submit:
	$(PYTHON) scripts/submit_workflow.py --mode "$(MODE)" --target "$(TARGET)" --rscript "$(RSCRIPT)" --bedtools "$(BEDTOOLS)" --blastn "$(BLASTN)" --makeblastdb "$(MAKEBLASTDB)" --mafft "$(MAFFT)" --iqtree "$(IQTREE)" --raw-phylogeny-dir "$(RAW_PHYLOGENY_DIR)"
