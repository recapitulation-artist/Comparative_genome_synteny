#!/usr/bin/env Rscript
# Verify manual annotations, tree-derived row order, links and output formats.
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
source(file.path(dirname(script), "environment.R"))
check_environment()
args <- commandArgs(trailingOnly = TRUE)
tree_source <- "reference"
phylogeny_dir <- "results/phylogeny"
usage <- "Usage: Rscript scripts/verify_results.R [--tree-source reference|inferred] [--phylogeny-dir PATH]"
while (length(args)) {
  name <- sub("=.*$", "", args[1])
  if (!name %in% c("--tree-source", "--phylogeny-dir")) stop(usage)
  if (grepl("=", args[1], fixed = TRUE)) {
    value <- sub("^[^=]+=", "", args[1])
    args <- args[-1]
  } else {
    if (length(args) < 2L) stop(usage)
    value <- args[2]
    args <- args[-c(1, 2)]
  }
  if (!nzchar(value)) stop(usage)
  if (name == "--tree-source") tree_source <- value else phylogeny_dir <- value
}
if (!tree_source %in% c("reference", "inferred")) stop(usage)
resolve_path <- function(path) if (startsWith(path, "/")) path else file.path(project_root, path)
input <- function(name) file.path(project_root, "data", name)
output <- function(name) file.path(project_root, "results", tree_source, name)
tree_path <- if (tree_source == "reference") "data/phylogeny.nwk" else file.path(phylogeny_dir, "default", "tree.treefile")
tree_file <- normalizePath(resolve_path(tree_path), mustWork = TRUE)
root_prefix <- paste0(normalizePath(project_root), "/")
tree_path <- if (startsWith(tree_file, root_prefix)) substring(tree_file, nchar(root_prefix) + 1L) else tree_file
tree <- ape::read.tree(tree_file)
tree$tip.label <- trimws(tree$tip.label)
tree <- ape::ladderize(tree, right = FALSE)

seqs <- read.csv(input("sequences.csv"), stringsAsFactors = FALSE)
genes <- read.delim(input("genes.tsv"), stringsAsFactors = FALSE, quote = "", comment.char = "")
genomes <- read.delim(input("genomes.tsv"), stringsAsFactors = FALSE, quote = "", comment.char = "")
blast <- read.csv(input("blast_hits.csv"), stringsAsFactors = FALSE)
seqs$seq_id <- trimws(seqs$seq_id)
for (name in c("seq_id", "strand", "note", "organisms")) genes[[name]] <- trimws(genes[[name]])
blast$seq_id <- trimws(blast$seq_id)
blast$seq_id2 <- trimws(blast$seq_id2)
stars <- genes[genes$note == "Uptake-Gene", ]
genus_map <- unique(genes[, c("seq_id", "organisms")])
stopifnot(nrow(seqs) == 29L, nrow(genes) == 579L, nrow(stars) == 29L,
          !anyDuplicated(seqs$seq_id), !anyDuplicated(tree$tip.label),
          !anyDuplicated(genomes$seq_id), !anyDuplicated(genus_map$seq_id),
          setequal(tree$tip.label, seqs$seq_id), setequal(genes$seq_id, seqs$seq_id),
          setequal(genomes$seq_id, seqs$seq_id), !anyNA(genomes$organism_name),
          all(nzchar(genomes$organism_name)), all(table(stars$seq_id) == 1L),
          all(is.finite(seqs$length)), all(seqs$length > 0),
          !anyNA(genes[, c("start", "end", "strand", "note", "organisms")]),
          all(genes$strand %in% c("+", "-")), all(genes$start >= 1),
          all(genes$end >= genes$start),
          all(genes$end <= seqs$length[match(genes$seq_id, seqs$seq_id)]),
          all(blast$seq_id %in% seqs$seq_id), all(blast$seq_id2 %in% seqs$seq_id),
          all(is.finite(as.matrix(blast[, c("pident", "start", "end", "start2", "end2")]))),
          all(blast$pident >= 0 & blast$pident <= 100),
          all(pmin(blast$start, blast$end) >= 1), all(pmin(blast$start2, blast$end2) >= 1),
          all(pmax(blast$start, blast$end) <= seqs$length[match(blast$seq_id, seqs$seq_id)]),
          all(pmax(blast$start2, blast$end2) <= seqs$length[match(blast$seq_id2, seqs$seq_id)]))

rows <- read.csv(output("sequence_order.csv"), stringsAsFactors = FALSE)
tree_rows <- data.frame(seq_id = tree$tip.label, y_pos = ape::node.height(tree)[seq_len(ape::Ntip(tree))])
tree_rows <- tree_rows[order(tree_rows$y_pos, decreasing = TRUE), ]
stopifnot(identical(rows$seq_id, tree_rows$seq_id), !anyDuplicated(rows$seq_id),
          identical(rows$y_pos, rev(seq_len(nrow(rows)))),
          all(rows$y_pos == tree_rows$y_pos),
          all(rows$length == seqs$length[match(rows$seq_id, seqs$seq_id)]),
          all(rows$genes == as.integer(table(genes$seq_id)[rows$seq_id])),
          all(rows$uptake_genes == as.integer(table(stars$seq_id)[rows$seq_id])),
          all(rows$uptake_genes == 1L),
          identical(rows$organisms, genus_map$organisms[match(rows$seq_id, genus_map$seq_id)]),
          identical(rows$organism_name, genomes$organism_name[match(rows$seq_id, genomes$seq_id)]))
y1 <- rows$y_pos[match(blast$seq_id, rows$seq_id)]
y2 <- rows$y_pos[match(blast$seq_id2, rows$seq_id)]
chosen <- abs(y1 - y2) == 1 & blast$pident >= 70
selected <- blast[chosen, ]
hits <- read.csv(output("plotted_blast_hits.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(hits) == nrow(selected), all(names(blast) %in% names(hits)),
          all(hits$y1 == y1[chosen]), all(hits$y2 == y2[chosen]))
if (nrow(selected)) {
  stopifnot(isTRUE(all.equal(selected, hits[, names(blast)], check.attributes = FALSE)))
  colors <- read.csv(file.path(project_root, "config", "colors.csv"), stringsAsFactors = FALSE)
  levels <- colors$label[colors$group == "identity"]
  stopifnot(identical(hits$grad, ifelse(hits$pident >= 90, levels[1],
                                      ifelse(hits$pident >= 80, levels[2], levels[3]))))
}
calculated <- c(
  sequences = nrow(seqs), tree_tips = ape::Ntip(tree), genes = nrow(genes),
  gene_categories = length(unique(genes$note)), genera = length(unique(rows$organisms)),
  uptake_genes = nrow(stars), input_blast_hits = nrow(blast),
  plotted_blast_hits = nrow(hits), adjacent_pairs = length(unique(pmin(hits$y1, hits$y2))),
  reverse_hits = sum((hits$end - hits$start) * (hits$end2 - hits$start2) < 0),
  identity_ge90 = sum(hits$pident >= 90),
  identity_80to90 = sum(hits$pident >= 80 & hits$pident < 90),
  identity_70to80 = sum(hits$pident >= 70 & hits$pident < 80))
expected <- read.delim(file.path(project_root, "config", "expected_metrics.tsv"))
actual <- read.delim(output("metrics.tsv"))
stopifnot(identical(actual$metric, expected$metric), identical(actual$metric, names(calculated)),
          isTRUE(all.equal(actual$value, unname(calculated), check.attributes = FALSE)),
          identical(actual$value[1:7], expected$value[1:7]))
if (tree_source == "reference") {
  stopifnot(identical(actual$value, expected$value), calculated[["adjacent_pairs"]] == nrow(rows) - 1L)
}

metadata_paths <- c(tree = tree_path, annotations = "data/genes.tsv", sequences = "data/sequences.csv",
                    blast_hits = "data/blast_hits.csv", genomes = "data/genomes.tsv", colors = "config/colors.csv")
sha256 <- function(path) {
  value <- system2("sha256sum", shQuote(resolve_path(path)), stdout = TRUE)
  hash <- substr(value[1], 1, 64)
  if (length(hash) != 1L || !grepl("^[0-9a-f]{64}$", hash)) stop("Cannot hash input: ", path)
  hash
}
metadata <- read.delim(output("plot_metadata.tsv"), stringsAsFactors = FALSE)
stopifnot(identical(metadata$input, names(metadata_paths)), identical(metadata$path, unname(metadata_paths)),
          identical(metadata$sha256, unname(vapply(metadata_paths, sha256, character(1)))),
          all(metadata$tree_source == tree_source))

pdf <- readBin(output("synteny.pdf"), "raw", n = 5)
png <- readBin(output("synteny.png"), "raw", n = 24)
stopifnot(identical(rawToChar(pdf), "%PDF-"),
          identical(as.integer(png[1:8]), c(137L, 80L, 78L, 71L, 13L, 10L, 26L, 10L)))
dimensions <- vapply(list(png[17:20], png[21:24]), function(bytes) sum(as.integer(bytes) * 256^(3:0)), numeric(1))
stopifnot(identical(dimensions, c(4200, 4500)),
          all(file.info(c(output("synteny.pdf"), output("synteny.png")))$size > 10000))
cat(sprintf("Verification passed (%s): 29 sequences, 579 manual annotations, %d tree-adjacent links, input hashes, PDF and PNG.\n",
             tree_source, nrow(hits)))
