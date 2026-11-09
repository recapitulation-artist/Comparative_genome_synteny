# Cd Uptake Synteny

Reproduce an aligned phylogenetic topology and gene-synteny figure for 29
sequences. The workflow validates inputs and the R environment, renders the
figure, and checks the results against the bundled dataset.

## Run

Linux x86_64 with glibc >=2.28 and <3.0 is required for the Conda lock.
Conda and `sha256sum` must be available. Initial setup requires internet access.

```sh
conda create --name cd_uptake_synteny --file environment/conda_linux_64.lock
conda activate cd_uptake_synteny
make setup
make verify
```

`make setup` installs exact CRAN package versions into `.r_library/`, using
current or archived source releases in dependency order. It does not modify
the user R library. The Conda lock includes R 4.6.0, compilers, Cairo and fonts;
`environment/r_packages.lock.tsv` pins all 23 non-base R dependencies.

With an existing environment matching the locks, `make verify` can run
directly. Override the executable with `make verify RSCRIPT=/path/to/Rscript`.
Use `make check` to validate inputs and dependencies, or `make run` to render.
All Make targets invoke R with `--vanilla` to skip user startup files.

## Files

```text
config/       Fixed colors and expected dataset metrics
data/         Standardized input copies, checksums and schemas
environment/  Conda runtime lock and R package lock
scripts/      Environment setup, plotting and result verification
results/      Generated figures, tables and session information (Git-ignored)
```

All entry points resolve paths relative to the repository. No original file or
external project directory is needed. Input provenance is documented in
[`data/README.md`](data/README.md). The plotting script is a cleaned adaptation
of the 29-sequence section of `gggenome_xishou.R`.

## Outputs and methods

- `results/synteny.pdf`: one-page vector figure, 28 by 30 inches.
- `results/synteny.png`: 4,200 by 4,500 pixels at 150 dpi.
- `results/sequence_order.csv`: top-to-bottom sequence order and annotation counts.
- `results/plotted_blast_hits.csv`: displayed hits with original directed coordinates.
- `results/metrics.tsv`: dataset and filtering counts.
- `results/session_info.txt`: actual R and package environment.

The supplied topology is ladderized without branch lengths or an inferred
biological root. Tree tips and gene tracks share the same vertical coordinates.
Only hits between adjacent tracks with identity >=70% are drawn, preserving
reverse endpoints to show inversions. Identity bins are >=90%, 80 to <90%, and
70 to <80%. Gene arrows follow annotation strands; red stars mark `Uptake-Gene`.
Colors and legend order are fixed in `config/colors.csv`.

## Verification

`make verify` checks input SHA256 hashes, every locked package version, identifier
and coordinate consistency, expected metrics, selected BLAST coordinates, all
28 adjacent pairs, PDF signature, and PNG signature and dimensions. Expected
results include 579 genes, 29 uptake markers and 135 displayed BLAST hits.

The workflow has been tested with the matching local R environment and in a
separate copy of the tracked repository. Fresh network restoration has not
been tested. Fonts and PDF timestamps can vary across machines; file hashes of
rendered figures are not used as a reproducibility criterion.
