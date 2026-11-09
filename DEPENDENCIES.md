# Dependencies

Runtime and installation requirements. See [README.md](README.md) for setup and
usage; `environment/` contains exact builds, URLs and hashes.

## Platform and host

| Dependency | Requirement / role |
| --- | --- |
| Platform | Linux x86_64 / Windows WSL2 x86_64, glibc >=2.28 |
| Python | >=3.9; standard library only |
| GNU Make | Tested with 4.4.1; workflow entry points |
| Bash | Slurm runner and installation shell |
| GNU coreutils | `sha256sum`, `cmp` and filesystem commands |
| curl | BEDTools download; included in R Conda lock |
| Conda / Miniforge | Restore explicit environments |
| Git | Obtain and version the repository |
| Slurm | Optional `sbatch`/`srun` for submission |

Unix file locking is used. Memory depends on genome count and alignment size.

## Analysis software

| Software | Required version | Role / lock |
| --- | --- | --- |
| BEDTools | 2.31.0 | Extract and intersect; `environment/bedtools.lock.json` |
| NCBI BLAST+ | 2.17.0+ | `blastn` / `makeblastdb`; `environment/blast_linux_64.lock` |
| MAFFT | 7.526 | Alignment; `environment/phylogeny_linux_64.lock` |
| IQ-TREE | 3.1.4 stable | Inference; `environment/phylogeny_linux_64.lock` |
| R | 4.6.0 with Cairo | Figures/verification; `environment/conda_linux_64.lock` |

Main plotting needs R and basic commands; reconstructed entries need all tools.
Official BEDTools downloads require the locked SHA256; existing installations
must match the required version.

## R packages

`scripts/environment.R` restores `environment/r_packages.lock.tsv` into
`.r_library/`.

| Package | Version | Role |
| --- | --- | --- |
| ape | 5.8-1 | Phylogenetic tree input/layout |
| digest | 0.6.39 | ape dependency |
| lattice | 0.23-1 | ape/nlme dependency |
| nlme | 3.1-171 | ape dependency |
| Rcpp | 1.1.2 | Compiled ape dependency |

`grid`, `grDevices`, `graphics`, `methods`, `stats`, `utils` and `parallel` ship
with R. The R Conda lock includes compilers, Fortran, Make, Cairo and fontconfig.

## Fonts

Liberation Sans regular, italic and bold are bundled in `environment/fonts/`
with their license and SHA256 manifest. Preserve all locks and fonts when
distributing the repository; restore environments at the destination.
