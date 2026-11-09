RSCRIPT ?= Rscript
PYTHON ?= python3
BEDTOOLS ?= .tools/bedtools

.PHONY: all setup setup-bedtools setup-r check-data check-tools check-r check extract plot run verify

all: run

setup: setup-bedtools setup-r

setup-bedtools:
	$(PYTHON) scripts/setup_bedtools.py --bedtools "$(BEDTOOLS)"

setup-r:
	$(RSCRIPT) --vanilla scripts/environment.R --install

check-data:
	sha256sum --check data/checksums.sha256

check-tools:
	$(PYTHON) scripts/setup_bedtools.py --check --bedtools "$(BEDTOOLS)"

check-r:
	$(RSCRIPT) --vanilla scripts/environment.R

check: check-data check-tools check-r

extract: check-data check-tools
	$(PYTHON) scripts/extract_regions.py --bedtools "$(BEDTOOLS)"

plot: check-data check-r
	$(RSCRIPT) --vanilla scripts/plot_synteny.R

run: extract plot

verify: run
	$(RSCRIPT) --vanilla scripts/verify_results.R
