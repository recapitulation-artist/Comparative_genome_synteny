# Dependencies

Central inventory of runtime, installation and package dependencies. See
[README.md](README.md) for installation and entry commands. Machine-readable
locks in `environment/` retain exact Conda builds, package URLs and hashes.

## Platform and host

| Dependency | Requirement / role |
| --- | --- |
| Platform | Linux x86_64 / Windows WSL2 x86_64, glibc >=2.28 |
| Python | >=3.9; standard library only |
| GNU Make | Tested with 4.4.1; workflow entry points |
| Bash | Slurm runner and installation shell |
| GNU coreutils | `sha256sum`, `cmp` and filesystem commands |
| curl | Verified BEDTools download; included in R Conda lock |
| Conda / Miniforge | On PATH to restore explicit environments |
| Git | Obtain and version the repository |
| Slurm | Optional `sbatch`/`srun`, only for submission |

Unix file locking is used. Native macOS, Windows and ARM are outside this
version's supported scope. Default CPU budget: 56, overridable with `THREADS`.
Memory depends on genomes, region count and alignment size; no specific node or
fixed memory request is required locally.

## Analysis software

| Software | Required version | Role / lock |
| --- | --- | --- |
| BEDTools | 2.31.0 | Extract and intersect; `environment/bedtools.lock.json` |
| NCBI BLAST+ | 2.17.0+ | `blastn` / `makeblastdb`; `environment/blast_linux_64.lock` |
| MAFFT | 7.526 | Alignment; `environment/phylogeny_linux_64.lock` |
| IQ-TREE | 3.1.4 stable | Inference; `environment/phylogeny_linux_64.lock` |
| R | 4.6.0 with Cairo | Figures/verification; `environment/conda_linux_64.lock` |

Main plotting needs R and basic commands only. Curated/automatic additionally
need Python and all sequence-analysis tools above. Executable versions and hashes
are recorded. Official BEDTools downloads require the locked SHA256 before
execution; other installations of the required version are accepted.

## Complete non-base R package set

`environment/r_packages.lock.tsv` is restored by `scripts/environment.R` into
`.r_library/`.

| Package | Version | Role |
| --- | --- | --- |
| ape | 5.8-1 | Phylogenetic tree input/layout |
| digest | 0.6.39 | ape dependency |
| lattice | 0.23-1 | ape/nlme dependency |
| nlme | 3.1-171 | ape dependency |
| Rcpp | 1.1.2 | Compiled ape dependency |

`grid`, `grDevices`, `graphics`, `methods`, `stats`, `utils` and `parallel` ship
with R. ggplot2/gridExtra are not needed. The R Conda lock also fixes compilers,
Fortran, Make, Cairo, fontconfig and transitive libraries for source installation
and consistent graphics.

## Fonts and environment records

Liberation Sans regular, italic and bold are bundled in `environment/fonts/`
with license and SHA256 manifest. Relative fontconfig settings avoid dependence
on personal fonts.

Preserve every lock when distributing tracked files. Conda prefixes and compiled
R libraries are installation artifacts and should be restored at the destination.
Initial installation downloads software from locked URLs and CRAN/CRAN Archive;
analysis itself needs no genome downloads.
