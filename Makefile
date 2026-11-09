RSCRIPT ?= Rscript

.PHONY: all setup check run verify

all: run

setup:
	$(RSCRIPT) --vanilla scripts/environment.R --install

check:
	sha256sum --check data/checksums.sha256
	$(RSCRIPT) --vanilla scripts/environment.R

run: check
	$(RSCRIPT) --vanilla scripts/plot_synteny.R

verify: run
	$(RSCRIPT) --vanilla scripts/verify_results.R
