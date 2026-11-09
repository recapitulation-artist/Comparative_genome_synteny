RSCRIPT ?= Rscript
PYTHON ?= python3
BEDTOOLS ?= .tools/bedtools
BLASTN ?=
MAKEBLASTDB ?=

.PHONY: all setup setup-bedtools setup-blast setup-r check-data check-tools check-blast check-r check extract blast plot run verify

all: run

setup: setup-bedtools setup-blast setup-r

setup-bedtools:
	$(PYTHON) scripts/setup_bedtools.py --bedtools "$(BEDTOOLS)"

setup-blast:
	$(PYTHON) scripts/setup_blast.py

setup-r:
	$(RSCRIPT) --vanilla scripts/environment.R --install

check-data:
	sha256sum --check data/checksums.sha256

check-tools:
	$(PYTHON) scripts/setup_bedtools.py --check --bedtools "$(BEDTOOLS)"

check-blast:
	$(PYTHON) scripts/run_blast.py --check-tools --blastn "$(BLASTN)" --makeblastdb "$(MAKEBLASTDB)"

check-r:
	$(RSCRIPT) --vanilla scripts/environment.R

check: check-data check-tools check-blast check-r

extract: check-data check-tools
	$(PYTHON) scripts/extract_regions.py --bedtools "$(BEDTOOLS)"

blast: extract check-blast
	$(PYTHON) scripts/run_blast.py --blastn "$(BLASTN)" --makeblastdb "$(MAKEBLASTDB)"

plot: blast check-r
	$(RSCRIPT) --vanilla scripts/plot_synteny.R

run: plot

verify: run
	$(RSCRIPT) --vanilla scripts/verify_results.R
