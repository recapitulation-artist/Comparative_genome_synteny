#!/usr/bin/env Rscript
# Annotate inferred region trees and compare unrooted splits with the reference.
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
source(file.path(dirname(script), "environment.R"))
check_environment()
suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
phylogeny_dir <- file.path(project_root, "results", "phylogeny")
if (length(args)) {
  if (length(args) != 2L || args[1] != "--phylogeny-dir") {
    stop("Usage: Rscript scripts/report_phylogeny.R [--phylogeny-dir PATH]")
  }
  phylogeny_dir <- normalizePath(args[2], mustWork = TRUE)
}
write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}
require_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}
species <- read.delim(file.path(phylogeny_dir, "species.tsv"), stringsAsFactors = FALSE,
                      check.names = FALSE, quote = "", comment.char = "", na.strings = character())
columns <- c("seq_id", "organism_name", "strain", "assembly_accession", "contig", "target_strand")
require_true(all(columns %in% names(species)), "species.tsv lacks required metadata columns")
expected_ids <- read.csv(file.path(project_root, "data", "sequences.csv"),
                         stringsAsFactors = FALSE)$seq_id
require_true(nrow(species) == 29L && !anyDuplicated(species$seq_id) &&
               length(expected_ids) == 29L && !anyDuplicated(expected_ids) &&
               setequal(species$seq_id, expected_ids), "Expected matching metadata for 29 unique sequence IDs")
require_true(!anyNA(species[, columns]) && all(nzchar(species$organism_name)) &&
               all(species$target_strand %in% c("+", "-")), "Invalid species metadata")
species <- species[match(expected_ids, species$seq_id), , drop = FALSE]
ids <- sort(expected_ids)
safe_label <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)
species$annotated_tip <- paste(safe_label(species$organism_name), safe_label(species$strain),
                              species$seq_id, sep = "|")
species$display_label <- paste0(species$organism_name, " / ", species$strain, " [", species$seq_id, "]")
require_true(!anyDuplicated(species$annotated_tip), "Annotated tip labels are not unique")

tree_paths <- c(reference = file.path(project_root, "data", "phylogeny.nwk"),
                default = file.path(phylogeny_dir, "default", "tree.treefile"),
                oriented = file.path(phylogeny_dir, "oriented", "tree.treefile"))
recorded_paths <- c(recorded_default = file.path(phylogeny_dir, "recorded_default.nwk"),
                    recorded_oriented = file.path(phylogeny_dir, "recorded_oriented.nwk"))
has_recorded <- file.exists(recorded_paths)
require_true(all(has_recorded) || !any(has_recorded), "Both recorded baseline trees must be supplied together")
if (all(has_recorded)) tree_paths <- c(tree_paths, recorded_paths)
trees <- lapply(names(tree_paths), function(name) {
  path <- tree_paths[[name]]
  tree <- ape::read.tree(path)
  require_true(inherits(tree, "phylo") && length(tree$tip.label) == 29L &&
                 !anyDuplicated(tree$tip.label) && setequal(tree$tip.label, ids),
               paste("Tree tip IDs differ from metadata:", name))
  require_true(length(tree$edge.length) == nrow(tree$edge) &&
                 all(is.finite(tree$edge.length)) && all(tree$edge.length >= 0),
               paste("Missing or invalid branch lengths:", name))
  tree
})
names(trees) <- names(tree_paths)

# Canonicalize each unrooted bipartition by its smaller side, then lexical order.
# Complementary root edges collapse to one split when an input tree is rooted.
canonical_side <- function(side) {
  side <- sort(side)
  other <- setdiff(ids, side)
  if (length(side) > length(other) ||
      (length(side) == length(other) && paste(side, collapse = "|") > paste(other, collapse = "|"))) other else side
}
split_key <- function(side) paste(canonical_side(side), collapse = "|")
tree_splits <- function(tree) {
  ntip <- length(tree$tip.label)
  children <- split(tree$edge[, 2], tree$edge[, 1])
  descendants <- function(node) {
    if (node <= ntip) return(tree$tip.label[node])
    unlist(lapply(children[[as.character(node)]], descendants), use.names = FALSE)
  }
  rows <- list()
  for (i in seq_len(nrow(tree$edge))) {
    node <- tree$edge[i, 2]
    side <- canonical_side(descendants(node))
    if (length(side) < 2L || length(side) > ntip - 2L) next
    label <- if (node > ntip && length(tree$node.label)) tree$node.label[node - ntip] else ""
    if (is.na(label)) label <- ""
    rows[[length(rows) + 1L]] <- data.frame(key = paste(side, collapse = "|"),
      side_ids = paste(side, collapse = ";"), complement_ids = paste(setdiff(ids, side), collapse = ";"),
      node_label = label, stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(data.frame(key = character(), side_ids = character(),
    complement_ids = character(), node_label = character(), stringsAsFactors = FALSE))
  result <- do.call(rbind, rows)
  result <- result[order(result$key, !nzchar(result$node_label)), , drop = FALSE]
  result[!duplicated(result$key), , drop = FALSE]
}
splits <- lapply(trees, tree_splits)
reference_labels <- suppressWarnings(as.numeric(splits$reference$node_label))
high_reference <- splits$reference$key[!is.na(reference_labels) & reference_labels >= 95]
positive_ids <- species$seq_id[species$target_strand == "+"]
negative_ids <- species$seq_id[species$target_strand == "-"]
require_true(length(positive_ids) == 17L && length(negative_ids) == 12L,
             "Expected the historical 17-positive / 12-negative strand groups")
strand_key <- split_key(positive_ids)
distances <- lapply(trees, function(tree) {
  matrix <- ape::cophenetic.phylo(tree)[ids, ids]
  matrix[upper.tri(matrix)]
})
require_true(all(vapply(distances, length, integer(1)) == choose(29L, 2L)),
             "Expected 406 pairwise tip distances")

tree_metadata <- do.call(rbind, lapply(names(trees), function(name) {
  tree <- trees[[name]]
  data.frame(tree = name, tips = length(tree$tip.label), unrooted_splits = nrow(splits[[name]]),
    total_branch_length = sum(tree$edge.length), maximum_branch_length = max(tree$edge.length),
    reference_raw_labels_ge95_total = length(high_reference),
    reference_raw_labels_ge95_retained = sum(high_reference %in% splits[[name]]$key),
    strand_17_12_split_present = strand_key %in% splits[[name]]$key, stringsAsFactors = FALSE)
}))
write_tsv(tree_metadata, file.path(phylogeny_dir, "tree_metadata.tsv"))

read_bytes <- function(path) readBin(path, "raw", n = file.info(path)$size)
comparisons <- list(c("reference", "default"), c("reference", "oriented"), c("default", "oriented"))
if (all(has_recorded)) comparisons <- c(comparisons,
  list(c("recorded_default", "default"), c("recorded_oriented", "oriented")))
difference_rows <- list()
comparison <- do.call(rbind, lapply(comparisons, function(pair) {
  a <- pair[1]; b <- pair[2]
  sa <- splits[[a]]; sb <- splits[[b]]
  only_a <- setdiff(sa$key, sb$key); only_b <- setdiff(sb$key, sa$key)
  for (key in c(only_a, only_b)) {
    in_a <- key %in% sa$key; in_b <- key %in% sb$key
    row <- if (in_a) sa[match(key, sa$key), ] else sb[match(key, sb$key), ]
    difference_rows[[length(difference_rows) + 1L]] <<- data.frame(tree_a = a, tree_b = b,
      present_in_a = in_a, present_in_b = in_b, split_side_ids = row$side_ids,
      split_complement_ids = row$complement_ids,
      node_label_a = if (in_a) sa$node_label[match(key, sa$key)] else "",
      node_label_b = if (in_b) sb$node_label[match(key, sb$key)] else "", stringsAsFactors = FALSE)
  }
  rf <- length(only_a) + length(only_b)
  denominator <- nrow(sa) + nrow(sb)
  pearson <- if (sd(distances[[a]]) > 0 && sd(distances[[b]]) > 0)
    cor(distances[[a]], distances[[b]], method = "pearson") else NA_real_
  data.frame(tree_a = a, tree_b = b, matching_unique_tips = 29L,
    splits_a = nrow(sa), splits_b = nrow(sb), shared_splits = length(intersect(sa$key, sb$key)),
    unique_splits_a = length(only_a), unique_splits_b = length(only_b), rf_distance = rf,
    normalized_rf = if (denominator) rf / denominator else 0,
    patristic_pairs = length(distances[[a]]),
    patristic_rmse = sqrt(mean((distances[[a]] - distances[[b]])^2)), patristic_pearson = pearson,
    reference_raw_labels_ge95_total = length(high_reference),
    reference_raw_labels_ge95_retained_a = sum(high_reference %in% sa$key),
    reference_raw_labels_ge95_retained_b = sum(high_reference %in% sb$key),
    strand_17_12_split_a = strand_key %in% sa$key, strand_17_12_split_b = strand_key %in% sb$key,
    exact_file_bytes_equal = identical(read_bytes(tree_paths[[a]]), read_bytes(tree_paths[[b]])),
    stringsAsFactors = FALSE)
}))
write_tsv(comparison, file.path(phylogeny_dir, "comparison.tsv"))
differences <- if (length(difference_rows)) do.call(rbind, difference_rows) else data.frame(
  tree_a = character(), tree_b = character(), present_in_a = logical(), present_in_b = logical(),
  split_side_ids = character(), split_complement_ids = character(), node_label_a = character(),
  node_label_b = character(), stringsAsFactors = FALSE)
write_tsv(differences, file.path(phylogeny_dir, "split_differences.tsv"))

draw_tree <- function(tree, title) {
  tree$tip.label <- species$display_label[match(tree$tip.label, species$seq_id)]
  tree <- ape::ladderize(tree)
  depth <- max(ape::node.depth.edgelength(tree))
  if (depth == 0) depth <- 1
  par(mar = c(3, 1, 4, 1), family = "sans")
  plot(tree, type = "phylogram", use.edge.length = TRUE, cex = 0.68, font = 1, underscore = TRUE,
       label.offset = depth * 0.012, x.lim = c(0, depth * 1.9),
       no.margin = FALSE, edge.width = 1, main = title)
  ticks <- pretty(c(0, depth)); ticks <- ticks[ticks >= 0 & ticks <= depth]
  axis(1, at = ticks, labels = format(ticks, digits = 3), cex.axis = 0.8)
  mtext("Branch length (estimated substitutions/site)", side = 1, line = 1.8, cex = 0.8)
  mtext("Region phylogeny; labels retain unique sequence IDs", side = 3, line = 0.4, cex = 0.8)
}
for (name in c("default", "oriented")) {
  directory <- file.path(phylogeny_dir, name)
  annotated <- trees[[name]]
  annotated$tip.label <- species$annotated_tip[match(annotated$tip.label, species$seq_id)]
  ape::write.tree(annotated, file = file.path(directory, "annotated_tree.nwk"), digits = 10)
  reread <- ape::read.tree(file.path(directory, "annotated_tree.nwk"))
  require_true(setequal(reread$tip.label, species$annotated_tip), "Annotated Newick tip mapping changed")
  write_tsv(species, file.path(directory, "metadata.tsv"))
  title <- if (name == "default") "Default region alignment" else "Target-oriented region alignment"
  grDevices::cairo_pdf(file.path(directory, "tree.pdf"), width = 14, height = 11)
  tryCatch(draw_tree(trees[[name]], title), finally = grDevices::dev.off())
  grDevices::png(file.path(directory, "tree.png"), width = 2240, height = 1760, res = 160, type = "cairo")
  tryCatch(draw_tree(trees[[name]], title), finally = grDevices::dev.off())
}

yes_no <- function(x) if (x) "present" else "absent"
summary <- c("# Region phylogeny comparison", "",
  sprintf("All compared trees contain the same 29 unique sequence IDs, representing %d recorded species.",
          length(unique(species$organism_name))),
  "These are local region phylogenies; they are not inferred species trees.", "")
for (i in seq_len(nrow(comparison))) {
  row <- comparison[i, ]
  summary <- c(summary, sprintf(
    "- %s vs %s: unrooted RF %d (normalized %.4f), %d shared splits; 406-distance RMSE %.6g and Pearson %.6g; exact file bytes equal: %s.",
    row$tree_a, row$tree_b, row$rf_distance, row$normalized_rf, row$shared_splits,
    row$patristic_rmse, row$patristic_pearson, tolower(row$exact_file_bytes_equal)))
}
summary <- c(summary, "", sprintf(
  "The reference has %d nontrivial unrooted splits and %d raw numeric node labels >=95. The historical support method is unknown; these labels are not assumed to be a particular bootstrap statistic.",
  nrow(splits$reference), length(high_reference)))
for (i in seq_len(nrow(tree_metadata))) {
  row <- tree_metadata[i, ]
  summary <- c(summary, sprintf(
    "- %s: %d/%d reference high-label splits retained; total branch length %.6g, maximum %.6g; the 17-positive/12-negative target-strand split is %s.",
    row$tree, row$reference_raw_labels_ge95_retained, row$reference_raw_labels_ge95_total,
    row$total_branch_length, row$maximum_branch_length, yes_no(row$strand_17_12_split_present)))
}
summary <- c(summary, "",
  "The strand groups refer to the original GFF target directions. A split matching those groups can reflect sequence orientation; this check alone does not establish its cause.",
  "RF comparisons ignore rooting, node rotations and Newick serialization. Normalized RF divides by the total number of nontrivial splits in the two compared trees.",
  "Patristic distances use the same 406 unordered tip pairs. Branch-length comparisons are descriptive when alignments or fitted models differ.",
  "Species and strain labels come from species.tsv; stable sequence IDs preserve distinct tips when species or strain labels repeat.")
writeLines(summary, file.path(phylogeny_dir, "summary.md"))
cat(sprintf("Annotated two region trees and compared %d pairs for 29 tips: %s\n", length(comparisons), phylogeny_dir))
