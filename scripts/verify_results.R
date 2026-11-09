#!/usr/bin/env Rscript
# Verify the inputs, palette assignments, topology and exports of one plotting entry.
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
source(file.path(dirname(script), "environment.R"))
check_environment()
source(file.path(dirname(script), "palettes.R"))
invisible(Sys.setlocale("LC_COLLATE", "C"))
args <- commandArgs(trailingOnly = TRUE)
entry <- "main"
entry_explicit <- FALSE
target_gene <- "Nramp"
input_dir <- NULL
phylogeny_dir <- NULL
gene_palette <- "main-muted"
link_palette <- "blue"
link_custom <- ""
legacy_tree <- "reference"
legacy_annotations <- "curated"
usage <- paste("Usage: Rscript scripts/verify_results.R [--entry main|curated|auto]",
               "[--target-gene NAME] [--input-dir PATH] [--phylogeny-dir PATH]",
               "[--gene-palette NAME] [--link-palette NAME] [--link-colors HEX,HEX,HEX]")
while (length(args)) {
  name <- sub("=.*$", "", args[1])
  if (!name %in% c("--entry", "--target-gene", "--input-dir", "--phylogeny-dir",
                  "--tree-source", "--annotation-source", "--gene-palette", "--link-palette",
                  "--link-colors")) stop(usage)
  if (grepl("=", args[1], fixed = TRUE)) {
    value <- sub("^[^=]+=", "", args[1]); args <- args[-1]
  } else {
    if (length(args) < 2L) stop(usage)
    value <- args[2]; args <- args[-c(1, 2)]
  }
  if (!nzchar(value) && name != "--link-colors") stop(usage)
  if (name == "--entry") { entry <- value; entry_explicit <- TRUE }
  if (name == "--target-gene") target_gene <- value
  if (name == "--input-dir") input_dir <- value
  if (name == "--phylogeny-dir") phylogeny_dir <- value
  if (name == "--gene-palette") gene_palette <- value
  if (name == "--link-palette") link_palette <- value
  if (name == "--link-colors") link_custom <- value
  if (name == "--tree-source") legacy_tree <- value
  if (name == "--annotation-source") legacy_annotations <- value
}
if (!legacy_tree %in% c("reference", "inferred") ||
    !legacy_annotations %in% c("curated", "raw")) stop(usage)
if (!entry_explicit) entry <- if (legacy_annotations == "raw") "auto" else {
  if (legacy_tree == "inferred") "curated" else "main"
}
if (!entry %in% c("main", "curated", "auto")) stop(usage)
if (!gene_palette %in% gene_palette_names || !link_palette %in% link_palette_names) stop(usage)
main_entry <- entry == "main"
tree_source <- if (main_entry) "reference" else "inferred"
resolve_path <- function(path) if (startsWith(path, "/")) path else file.path(project_root, path)
relative_path <- function(path) {
  path <- normalizePath(resolve_path(path), mustWork = TRUE)
  prefix <- paste0(normalizePath(project_root), "/")
  if (startsWith(path, prefix)) substring(path, nchar(prefix) + 1L) else path
}
read_tsv <- function(path, character = FALSE) read.delim(path, stringsAsFactors = FALSE,
  quote = "", comment.char = "", check.names = FALSE, na.strings = character(),
  colClasses = if (character) "character" else NA)
output <- function(name) file.path(project_root, "results", entry, name)
if (main_entry) {
  sequence_path <- "data/sequences.csv"
  annotation_path <- "data/genes.tsv"
  blast_path <- "data/reference_blast_hits.csv"
  genome_path <- "data/genomes.tsv"
  tree_path <- "data/phylogeny.nwk"
} else {
  if (is.null(input_dir)) input_dir <- file.path("results", entry)
  if (is.null(phylogeny_dir)) phylogeny_dir <- file.path(input_dir, "phylogeny")
  sequence_path <- file.path(input_dir, "sequences.csv")
  annotation_path <- file.path(input_dir, "genes.tsv")
  blast_path <- file.path(input_dir, "blast_hits.csv")
  genome_path <- file.path(input_dir, "genomes.tsv")
  tree_path <- file.path(phylogeny_dir, "default", "tree.treefile")
}
metadata_paths <- c(tree = tree_path, annotations = annotation_path, sequences = sequence_path,
                    blast_hits = blast_path, genomes = genome_path, colors = "config/colors.csv",
                    figure_style = "config/figure_style.tsv", gene_palette = "config/gene_palette.csv",
                    palette_code = "scripts/palettes.R")
metadata_paths <- vapply(metadata_paths, relative_path, character(1))
sha256 <- function(path) {
  value <- system2("sha256sum", shQuote(resolve_path(path)), stdout = TRUE)
  hash <- substr(value[1], 1, 64)
  if (length(hash) != 1L || !grepl("^[0-9a-f]{64}$", hash)) stop("Cannot hash input: ", path)
  hash
}
seqs <- read.csv(resolve_path(sequence_path), stringsAsFactors = FALSE)
genes_original <- read_tsv(resolve_path(annotation_path), TRUE)
genes <- read_tsv(resolve_path(annotation_path))
genomes <- read_tsv(resolve_path(genome_path))
blast <- read.csv(resolve_path(blast_path), stringsAsFactors = FALSE)
tree <- ape::read.tree(resolve_path(tree_path))
tree$tip.label <- trimws(tree$tip.label)
tree <- ape::ladderize(tree, right = FALSE)
seqs$seq_id <- trimws(seqs$seq_id)
for (name in c("seq_id", "strand", "note", "organisms")) genes[[name]] <- trimws(genes[[name]])
blast$seq_id <- trimws(blast$seq_id)
blast$seq_id2 <- trimws(blast$seq_id2)
is_true <- function(x) tolower(as.character(x)) %in% c("true", "1", "yes")
is_target <- if (main_entry) genes$note == "Uptake-Gene" else {
  stopifnot("is_target" %in% names(genes),
    all(tolower(as.character(genes$is_target)) %in% c("true", "false", "0", "1", "yes", "no")))
  is_true(genes$is_target)
}
stars <- genes[is_target, , drop = FALSE]
genus_map <- unique(genes[, c("seq_id", "organisms")])
stopifnot(nrow(seqs) >= if (main_entry) 2L else 4L, nrow(genes) > 0L,
          !anyDuplicated(seqs$seq_id), !anyDuplicated(tree$tip.label),
          !anyDuplicated(genomes$seq_id), !anyDuplicated(genus_map$seq_id),
          setequal(tree$tip.label, seqs$seq_id), setequal(genes$seq_id, seqs$seq_id),
          setequal(stars$seq_id, seqs$seq_id), setequal(genomes$seq_id, seqs$seq_id),
          !anyNA(genomes$organism_name), all(nzchar(genomes$organism_name)),
          nrow(stars) == nrow(seqs), all(table(stars$seq_id) == 1L),
          all(is.finite(seqs$length)), all(seqs$length > 0),
          !anyNA(genes[, c("start", "end", "strand", "note", "organisms")]),
          all(genes$strand %in% c("+", "-")), all(genes$start >= 1), all(genes$end >= genes$start),
          all(genes$end <= seqs$length[match(genes$seq_id, seqs$seq_id)]),
          all(blast$seq_id %in% seqs$seq_id), all(blast$seq_id2 %in% seqs$seq_id))
blast_numeric <- c("pident", "length", "mismatch", "gapopen", "start", "end", "start2", "end2", "evalue", "bitscore")
stopifnot(all(blast_numeric %in% names(blast)), all(is.finite(as.matrix(blast[, blast_numeric]))),
          all(blast$pident >= 0 & blast$pident <= 100), all(blast$length > 0),
          all(blast$mismatch >= 0 & blast$mismatch <= blast$length),
          all(blast$gapopen >= 0 & blast$gapopen <= blast$length),
          all(blast$evalue >= 0), all(blast$bitscore >= 0),
          all(pmin(blast$start, blast$end) >= 1), all(pmin(blast$start2, blast$end2) >= 1),
          all(pmax(blast$start, blast$end) <= seqs$length[match(blast$seq_id, seqs$seq_id)]),
          all(pmax(blast$start2, blast$end2) <= seqs$length[match(blast$seq_id2, seqs$seq_id)]),
          all(abs(blast$end - blast$start) + 1 <= blast$length),
          all(abs(blast$end2 - blast$start2) + 1 <= blast$length))

if (!main_entry) {
  regions <- read_tsv(resolve_path(file.path(input_dir, "regions.tsv")))
  neighbors <- read_tsv(resolve_path(file.path(input_dir, "neighbors.tsv")))
  stopifnot(!anyDuplicated(regions$seq_id), setequal(regions$seq_id, seqs$seq_id),
            !anyDuplicated(paste(neighbors$seq_id, neighbors$gene_id, sep = "\t")),
            setequal(neighbors$seq_id, seqs$seq_id),
            all(regions$region_start0 >= 0), all(regions$region_end0 <= regions$contig_length),
            all(regions$length == regions$region_end0 - regions$region_start0),
            all(regions$length == seqs$length[match(regions$seq_id, seqs$seq_id)]))
  region <- regions[match(neighbors$seq_id, regions$seq_id), ]
  stopifnot(all(neighbors$local_start_1based == pmax(neighbors$start_1based - 1L, region$region_start0) - region$region_start0 + 1L),
            all(neighbors$local_end_1based == pmin(neighbors$end_1based, region$region_end0) - region$region_start0),
            all(neighbors$local_start0 == neighbors$local_start_1based - 1L),
            all(neighbors$local_end0 == neighbors$local_end_1based),
            all(is_true(neighbors$clipped_left) == (neighbors$start_1based - 1L < region$region_start0)),
            all(is_true(neighbors$clipped_right) == (neighbors$end_1based > region$region_end0)),
            identical(neighbors$assembly_accession, region$assembly_accession),
            identical(neighbors$contig, region$contig))
  selected_targets <- neighbors[is_true(neighbors$is_target), ]
  target_regions <- regions[match(selected_targets$seq_id, regions$seq_id), ]
  stopifnot(nrow(selected_targets) == nrow(seqs), !anyDuplicated(selected_targets$seq_id),
            identical(selected_targets$locus_tag, target_regions$target_locus_tag),
            identical(selected_targets$strand, target_regions$target_strand),
            all(vapply(seq_len(nrow(selected_targets)), function(i) {
              target_regions$target_protein_id[i] %in% strsplit(selected_targets$protein_ids[i], ";", fixed = TRUE)[[1]]
            }, logical(1))))
  if (entry == "auto") {
    selected <- neighbors[is_true(neighbors$selected_for_plot), , drop = FALSE]
    original_selected <- read_tsv(resolve_path(file.path(input_dir, "neighbors.tsv")), TRUE)
    original_selected <- original_selected[is_true(original_selected$selected_for_plot), , drop = FALSE]
    rownames(original_selected) <- NULL
    original_plotted <- genes_original[, names(original_selected), drop = FALSE]
    rownames(original_plotted) <- NULL
    stopifnot(identical(original_plotted, original_selected), all(is_true(genes$is_protein_coding)),
              all(is_true(genes$selected_for_plot)), all(genes$feature_type == "gene"),
              all(genes$biotype %in% c("", "protein_coding")),
              identical(genes$start, genes$local_start_1based), identical(genes$end, genes$local_end_1based),
              !any(is_true(genes$clipped_left)), !any(is_true(genes$clipped_right)))
    for (i in seq_len(nrow(regions))) {
      g <- genes[genes$seq_id == regions$seq_id[i], , drop = FALSE]
      g <- g[order(g$start_1based, g$end_1based, g$gene_id, method = "radix"), ]
      target <- which(is_true(g$is_target))
      stopifnot(length(target) == 1L, g$gene_id[target] == regions$target_gene_id[i],
        target - 1L == regions$upstream_genes_selected[i],
        nrow(g) - target == regions$downstream_genes_selected[i],
        regions$upstream_genes_selected[i] <= regions$requested_neighbors_per_side[i],
        regions$downstream_genes_selected[i] <= regions$requested_neighbors_per_side[i],
        nrow(g) <= 2L * regions$requested_neighbors_per_side[i] + 1L)
    }
    matches <- read_tsv(resolve_path(file.path(input_dir, "target_matches.tsv")))
    matches <- matches[is_true(matches$selected), , drop = FALSE]
    stopifnot(!anyDuplicated(matches$seq_id), setequal(matches$seq_id, seqs$seq_id),
              identical(matches$gene_id[match(regions$seq_id, matches$seq_id)], regions$target_gene_id))
  }
  source_hashes <- read_tsv(resolve_path(file.path(input_dir, "input_hashes.tsv")), TRUE)
  stopifnot(!anyDuplicated(source_hashes$file),
            identical(source_hashes$sha256, unname(vapply(source_hashes$file, sha256, character(1)))))
  if (entry == "auto") stopifnot(!any(source_hashes$file %in%
    c("data/genes.tsv", "data/region_manifest.tsv", "data/phylogeny.nwk", "data/reference_sequences.fasta")))
}

rows <- read.csv(output("sequence_order.csv"), stringsAsFactors = FALSE)
tree_rows <- data.frame(seq_id = tree$tip.label, y_pos = ape::node.height(tree)[seq_len(ape::Ntip(tree))])
tree_rows <- tree_rows[order(tree_rows$y_pos, decreasing = TRUE), ]
stopifnot(identical(rows$seq_id, tree_rows$seq_id), !anyDuplicated(rows$seq_id),
          identical(rows$y_pos, rev(seq_len(nrow(rows)))), all(rows$y_pos == tree_rows$y_pos),
          all(rows$length == seqs$length[match(rows$seq_id, seqs$seq_id)]),
          all(rows$genes == as.integer(table(genes$seq_id)[rows$seq_id])),
          all(rows$uptake_genes == 1L),
          identical(rows$organisms, genus_map$organisms[match(rows$seq_id, genus_map$seq_id)]),
          identical(rows$organism_name, genomes$organism_name[match(rows$seq_id, genomes$seq_id)]))
plotted_genes <- read_tsv(output("plotted_genes.tsv"), TRUE)
stopifnot(identical(names(plotted_genes), c(names(genes_original), "y_pos", "color_category", "color")),
          identical(plotted_genes[, names(genes_original)], genes_original),
          identical(as.numeric(plotted_genes$y_pos), as.numeric(rows$y_pos[match(genes$seq_id, rows$seq_id)])))

style_config <- read_tsv(resolve_path("config/figure_style.tsv"), TRUE)
style_report <- read_tsv(output("figure_style.tsv"), TRUE)
required_style <- c("width_mm", "height_mm", "png_dpi", "font_family", "text_size_pt",
  "species_size_pt", "sequence_size_pt", "line_width_pt", "gene_line_width_pt",
  "genus_dot_mm", "target_dot_mm", "main_saturation_factor")
derived_style <- c("entry", "tree_representation", "row_pitch_mm", "minimum_font_pt", "maximum_font_pt",
                   "legend_columns", "legend_rows", "gene_palette", "gene_palette_capacity",
                   "link_palette", "link_colors")
stopifnot(identical(names(style_config), c("parameter", "value")),
          identical(names(style_report), c("parameter", "value")),
          !anyDuplicated(style_config$parameter), !anyDuplicated(style_report$parameter),
          setequal(style_config$parameter, required_style),
          setequal(style_report$parameter, c(required_style, derived_style)))
configured <- setNames(style_config$value, style_config$parameter)
reported <- setNames(style_report$value, style_report$parameter)
stopifnot(identical(reported[setdiff(required_style, "height_mm")], configured[setdiff(required_style, "height_mm")]),
          configured[["font_family"]] == "Liberation Sans", as.numeric(configured[["genus_dot_mm"]]) == 2.2,
          as.numeric(configured[["main_saturation_factor"]]) == 0.70,
          reported[["entry"]] == entry, reported[["tree_representation"]] == "cladogram")
colors <- read.csv(resolve_path("config/colors.csv"), stringsAsFactors = FALSE)
color_group <- function(group) {
  x <- colors[colors$group == group, ]; stopifnot(!anyDuplicated(x$label)); setNames(x$color, x$label)
}
hex <- function(x) grDevices::rgb(t(grDevices::col2rgb(x)) / 255)
desaturate <- function(x) {
  hsv <- grDevices::rgb2hsv(grDevices::col2rgb(x))
  hsv["s", ] <- hsv["s", ] * as.numeric(configured[["main_saturation_factor"]])
  grDevices::hsv(hsv["h", ], hsv["s", ], hsv["v", ])
}
palette_config <- read.csv(resolve_path("config/gene_palette.csv"), stringsAsFactors = FALSE)
stopifnot(identical(names(palette_config), c("palette", "rank", "color")),
          nrow(palette_config) == 42L * length(gene_palette_names),
          setequal(palette_config$palette, gene_palette_names),
          all(grepl("^#[0-9A-F]{6}$", palette_config$color)))
for (preset in gene_palette_names) {
  p <- palette_config[palette_config$palette == preset, , drop = FALSE]
  stopifnot(identical(p$rank, seq_len(42L)), !anyDuplicated(p$color))
}
ordinary_main <- setdiff(names(color_group("gene")), c("Uptake-Gene", "Unknown"))
main_anchor <- unname(hex(desaturate(color_group("gene")[ordinary_main])))
stopifnot(length(main_anchor) == 26L)
darkened <- grDevices::rgb(t(round(grDevices::col2rgb(main_anchor[seq_len(16L)]) * 0.85)),
                           maxColorValue = 255)
stopifnot(identical(palette_config$color[palette_config$palette == "main-muted"],
                    c(main_anchor, darkened)))
mix_white <- function(base, fractions) {
  unname(head(unlist(lapply(fractions, function(fraction) {
    grDevices::rgb(t(round(grDevices::col2rgb(base) * (1 - fraction) + 255 * fraction)),
                  maxColorValue = 255)
  })), 42L))
}
expected_palettes <- list(
  automatic = c("#B84C4C", "#BE762F", "#B39B24", "#43845A", "#328D91", "#4B70A9", "#8056A1",
    "#C46363", "#CA8949", "#C1AB40", "#5B966E", "#4DA0A3", "#6484B7", "#946DAC",
    "#CD7B7B", "#D49C64", "#CDBA5F", "#73A783", "#68B2B5", "#7D99C5", "#A785BD",
    "#D79393", "#DEAE7F", "#D9C97F", "#8CB998", "#84C4C6", "#96AED3", "#BA9DCD",
    "#E1ABAB", "#E8C19A", "#E5D89F", "#A5CAAD", "#A0D6D7", "#AFC3E1", "#CDB6DD",
    "#EBC3C3", "#F1D4B5", "#F0E8BF", "#BEDCC3", "#BDE8E8", "#C8D8EE", "#E0CEEB"),
  nature = mix_white(c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4",
                      "#91D1C2", "#DC0000", "#7E6148", "#B09C85"), c(0, 0.15, 0.30, 0.45, 0.60)),
  cell = mix_white(c("#3B6FB6", "#D8736F", "#4A9576", "#BD9A46", "#61A6B1", "#8C79B4", "#AD6C92"),
                   c(0, 0.12, 0.24, 0.36, 0.48, 0.60)),
  science = mix_white(c("#3B4992", "#EE0000", "#008B45", "#631879", "#008280", "#BB0021",
                       "#5F559B", "#A20056", "#808180", "#1B1919"), c(0, 0.15, 0.30, 0.45, 0.60)))
for (preset in names(expected_palettes)) stopifnot(identical(
  palette_config$color[palette_config$palette == preset], expected_palettes[[preset]]))
palette <- read_gene_palette(resolve_path("config/gene_palette.csv"), gene_palette)
stopifnot(identical(palette, palette_config$color[palette_config$palette == gene_palette]))
expected_links <- list(blue = c("#AECDE3", "#D0E1EE", "#EBF3F8"),
  green = c("#B4D5BD", "#D5E8DA", "#EEF6F0"), purple = c("#C7BDDE", "#E0DAED", "#F3F0F8"),
  orange = c("#E8C4AE", "#F2DDCF", "#FAF0E9"), original = c("#A6A6A6", "#CECECE", "#EBEBEB"))
for (preset in names(expected_links)) {
  group <- color_group(if (preset == "original") "identity" else paste0("identity_", preset))
  stopifnot(identical(names(group), names(color_group("identity"))),
            identical(unname(group), expected_links[[preset]]))
}
link_selection <- resolve_link_colors(colors, link_palette, link_custom)
link_colors <- link_selection$colors
expected_link <- if (nzchar(link_custom)) {
  values <- toupper(trimws(strsplit(link_custom, ",", fixed = TRUE)[[1]]))
  stopifnot(length(values) == 3L, all(grepl("^#[0-9A-F]{6}$", values)), !anyDuplicated(values))
  setNames(values, names(color_group("identity")))
} else color_group(if (link_palette == "original") "identity" else paste0("identity_", link_palette))
stopifnot(identical(link_selection$palette, if (nzchar(link_custom)) "custom" else link_palette),
          identical(link_colors, expected_link), length(link_colors) == 3L,
          identical(names(link_colors), names(color_group("identity"))),
          reported[["gene_palette"]] == gene_palette,
          as.numeric(reported[["gene_palette_capacity"]]) == length(palette),
          reported[["link_palette"]] == link_selection$palette,
          reported[["link_colors"]] == paste(link_colors, collapse = ","))
link_report <- read_tsv(output("link_colors.tsv"))
stopifnot(identical(names(link_report), c("category", "color", "lower_bound", "upper_bound", "upper_inclusive")),
          identical(link_report$category, names(link_colors)),
          identical(link_report$color, unname(link_colors)),
          identical(link_report$lower_bound, c(90L, 80L, 70L)),
          identical(link_report$upper_bound, c(100L, 90L, 80L)),
          identical(link_report$upper_inclusive, c(TRUE, FALSE, FALSE)))
target_category <- if (main_entry) "Uptake-Gene" else paste("Target", target_gene)
if (main_entry) {
  gene_colors <- hex(desaturate(color_group("gene")))
  names(gene_colors) <- names(color_group("gene"))
  if (gene_palette != "main-muted") {
    counts <- table(genes$note[genes$note %in% ordinary_main])
    ranked <- names(counts)[order(as.integer(counts), names(counts), decreasing = TRUE, method = "radix")]
    gene_colors[ranked] <- palette[seq_along(ranked)]
  }
  categories <- names(gene_colors)
  assigned <- genes$note
  category_ranks <- seq_along(categories)
  labels <- categories
  labels[labels == "Uptake-Gene"] <- paste("Target", target_gene)
  labels[labels == "Unknown"] <- "Unassigned"
} else {
  candidate <- genes$note
  if (entry == "auto") {
    candidate <- trimws(genes$gene)
    product <- if ("products" %in% names(genes)) trimws(genes$products) else rep("", nrow(genes))
    if (!any(nzchar(product)) && "product" %in% names(genes)) product <- trimws(genes$product)
    missing <- !nzchar(candidate) | candidate %in% c("-", ".", "NA")
    candidate[missing] <- product[missing]
  }
  unknown <- !nzchar(candidate) | tolower(candidate) %in% c("unknown", "other", "na", "-", ".") |
    grepl("hypothetical|uncharacterized|uncharacterised", candidate, ignore.case = TRUE)
  candidate[unknown] <- "Other"
  counts <- table(candidate[!is_target & candidate != "Other"])
  retained <- head(names(counts)[order(as.integer(counts), names(counts), decreasing = TRUE, method = "radix")], length(palette))
  categories <- c(target_category, retained, "Other")
  gene_colors <- setNames(c("#FFFFFF", palette[seq_along(retained)], "#BDBDBD"), categories)
  assigned <- ifelse(is_target, target_category, ifelse(candidate %in% retained, candidate, "Other"))
  category_ranks <- c(0L, seq_along(retained), length(palette) + 1L)
  labels <- c(target_category, sprintf("%02d %s", seq_along(retained), retained), "Other")
}
color_table <- read_tsv(output("gene_colors.tsv"))
stopifnot(identical(names(color_table), c("category", "label", "count", "rank", "color", "display_label")),
          identical(color_table$category, categories), identical(color_table$label, labels),
          identical(color_table$count, as.integer(table(factor(assigned, levels = categories)))),
          identical(color_table$rank, category_ranks), identical(color_table$color, unname(gene_colors)),
          all(nzchar(color_table$display_label)),
          identical(plotted_genes$color_category, assigned),
          identical(plotted_genes$color, unname(gene_colors[assigned])))
if (!main_entry) stopifnot(file.exists(output("legend.pdf")), file.info(output("legend.pdf"))$size > 1000)
genus_colors <- color_group("genus")
if (main_entry) genus_colors[] <- desaturate(genus_colors)
extra_genera <- sort(setdiff(unique(genus_map$organisms), names(genus_colors)))
if (length(extra_genera)) genus_colors <- c(genus_colors, setNames(grDevices::hcl.colors(length(extra_genera), "Dark 3"), extra_genera))
genus_colors <- genus_colors[names(genus_colors) %in% genus_map$organisms]
genus_colors[] <- hex(genus_colors)
genus_report <- read_tsv(output("genus_colors.tsv"))
stopifnot(identical(names(genus_report), c("category", "color", "count", "diameter_mm")),
          identical(genus_report$category, names(genus_colors)), identical(genus_report$color, unname(genus_colors)),
          identical(genus_report$count, as.integer(table(factor(rows$organisms, levels = names(genus_colors))))),
          all(genus_report$diameter_mm == 2.2))

y1 <- rows$y_pos[match(blast$seq_id, rows$seq_id)]
y2 <- rows$y_pos[match(blast$seq_id2, rows$seq_id)]
chosen <- abs(y1 - y2) == 1 & blast$pident >= 70
selected <- blast[chosen, , drop = FALSE]
hits <- read.csv(output("plotted_blast_hits.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(hits) == nrow(selected), all(names(blast) %in% names(hits)),
          all(hits$y1 == y1[chosen]), all(hits$y2 == y2[chosen]),
          isTRUE(all.equal(selected, hits[, names(blast), drop = FALSE], check.attributes = FALSE)))
levels <- colors$label[colors$group == "identity"]
stopifnot(identical(hits$grad, ifelse(hits$pident >= 90, levels[1], ifelse(hits$pident >= 80, levels[2], levels[3]))))
calculated <- c(sequences = nrow(seqs), tree_tips = ape::Ntip(tree), genes = nrow(genes),
  gene_categories = length(unique(assigned)), genera = length(unique(rows$organisms)), uptake_genes = nrow(stars),
  input_blast_hits = nrow(blast), plotted_blast_hits = nrow(hits),
  adjacent_pairs = length(unique(pmin(hits$y1, hits$y2))),
  reverse_hits = sum((hits$end - hits$start) * (hits$end2 - hits$start2) < 0),
  identity_ge90 = sum(hits$pident >= 90), identity_80to90 = sum(hits$pident >= 80 & hits$pident < 90),
  identity_70to80 = sum(hits$pident >= 70 & hits$pident < 80),
  assigned_color_categories = sum(!color_table$category %in% c(target_category, "Other") & color_table$count > 0),
  other_genes = sum(assigned == "Other"))
actual <- read_tsv(output("metrics.tsv"))
stopifnot(identical(actual$metric, names(calculated)),
          isTRUE(all.equal(actual$value, unname(calculated), check.attributes = FALSE)))
if (main_entry) {
  expected <- read_tsv(resolve_path("config/expected_metrics.tsv"))
  stopifnot(identical(actual$metric[seq_len(nrow(expected))], expected$metric),
            identical(actual$value[seq_len(nrow(expected))], expected$value),
            calculated[["adjacent_pairs"]] == nrow(rows) - 1L)
}
metadata <- read_tsv(output("plot_metadata.tsv"))
stopifnot(identical(names(metadata), c("input", "path", "sha256", "tree_source", "entry")),
          identical(metadata$input, names(metadata_paths)), identical(metadata$path, unname(metadata_paths)),
          identical(metadata$sha256, unname(vapply(metadata_paths, sha256, character(1)))),
          all(metadata$tree_source == tree_source), all(metadata$entry == entry))
font_hashes <- read.table(resolve_path("environment/fonts/checksums.sha256"), col.names = c("sha256", "path"), colClasses = "character")
stopifnot(identical(unname(vapply(font_hashes$path, sha256, character(1))), font_hashes$sha256))
legend_rows <- ceiling(nrow(color_table) / 7)
bottom_margin <- 30 + 4 * legend_rows
height_mm <- if (main_entry) as.numeric(configured[["height_mm"]]) else max(as.numeric(configured[["height_mm"]]), bottom_margin + 10 + 3 * nrow(seqs))
font_sizes <- as.numeric(configured[c("text_size_pt", "species_size_pt", "sequence_size_pt")])
stopifnot(all(is.finite(as.numeric(configured[setdiff(required_style, "font_family")]))),
          all(font_sizes >= 5), all(font_sizes <= 7),
          as.numeric(reported[["height_mm"]]) == height_mm,
          as.numeric(reported[["legend_columns"]]) == 7,
          as.numeric(reported[["legend_rows"]]) == legend_rows,
          as.numeric(reported[["minimum_font_pt"]]) == min(font_sizes),
          as.numeric(reported[["maximum_font_pt"]]) == 7,
          isTRUE(all.equal(as.numeric(reported[["row_pitch_mm"]]), (height_mm - bottom_margin - 10) / nrow(rows))))
caption <- paste(readLines(output("caption.txt"), warn = FALSE), collapse = " ")
stopifnot(grepl(paste0("of ", nrow(seqs), " target-associated regions"), caption, fixed = TRUE),
          grepl(if (main_entry) "supplied reference" else "newly inferred", caption, fixed = TRUE),
          grepl(gene_palette, caption, fixed = TRUE),
          grepl(link_selection$palette, caption, fixed = TRUE))
if (entry == "auto") stopifnot(grepl("original protein-coding", caption, fixed = TRUE), !grepl("manually curated", caption, fixed = TRUE))
if (!main_entry) stopifnot(grepl("top 42", caption, fixed = TRUE), grepl("Other", caption, fixed = TRUE))
pdf <- readBin(output("synteny.pdf"), "raw", n = 5)
png <- readBin(output("synteny.png"), "raw", n = 24)
stopifnot(identical(rawToChar(pdf), "%PDF-"), identical(as.integer(png[1:8]), c(137L, 80L, 78L, 71L, 13L, 10L, 26L, 10L)))
dimensions <- vapply(list(png[17:20], png[21:24]), function(bytes) sum(as.integer(bytes) * 256^(3:0)), numeric(1))
expected_dimensions <- trunc(as.numeric(reported[c("width_mm", "height_mm")]) / 25.4 * as.numeric(reported[["png_dpi"]]))
stopifnot(identical(dimensions, expected_dimensions), all(file.info(c(output("synteny.pdf"), output("synteny.png")))$size > 10000))

# Inspect the uncompressed font dictionaries emitted by this Cairo exporter.
# This deliberately supports these generated PDFs, rather than arbitrary PDF files.
pdf_raw <- readBin(output("synteny.pdf"), "raw", n = file.info(output("synteny.pdf"))$size)
pdf_raw[as.integer(pdf_raw) < 32L | as.integer(pdf_raw) > 126L] <- as.raw(32)
pdf_text <- rawToChar(pdf_raw)
objects <- regmatches(pdf_text, gregexpr("[[:digit:]]+[[:space:]]+[[:digit:]]+[[:space:]]+obj\\b.*?\\bendobj\\b",
                                       pdf_text, perl = TRUE))[[1]]
object_ids <- gsub("[[:space:]]+", " ", sub("[[:space:]]+obj.*$", "", objects))
if (!length(objects) || anyDuplicated(object_ids)) stop("Cannot inspect Cairo PDF objects.")
names(objects) <- object_ids
pdf_field <- function(object, pattern, field) {
  matches <- regmatches(object, regexec(pattern, object, perl = TRUE))[[1]]
  if (length(matches) != 2L) stop("Missing or unsupported Cairo PDF ", field, ".")
  matches[2]
}
pdf_object <- function(reference) {
  id <- gsub("[[:space:]]+", " ", reference)
  if (!id %in% names(objects)) stop("Missing Cairo PDF object: ", id)
  objects[[id]]
}
font_objects <- objects[grepl("/Type[[:space:]]+/Font[[:space:]/]", objects)]
expected_fonts <- c("LiberationSans", "LiberationSans-Italic", "LiberationSans-Bold")
if (length(font_objects) != 3L || !identical(configured[["font_family"]], "Liberation Sans")) {
  stop("Cairo PDF must contain exactly three Liberation Sans font dictionaries.")
}
font_names <- vapply(font_objects, function(object) {
  if (!grepl("/Subtype[[:space:]]+/TrueType[[:space:]/]", object)) {
    stop("Cairo PDF font is not TrueType.")
  }
  name <- pdf_field(object, "/BaseFont[[:space:]]+/([^[:space:]/<>]+)", "BaseFont")
  descriptor_ref <- pdf_field(object, "/FontDescriptor[[:space:]]+([[:digit:]]+[[:space:]]+[[:digit:]]+)[[:space:]]+R\\b", "FontDescriptor")
  unicode_ref <- pdf_field(object, "/ToUnicode[[:space:]]+([[:digit:]]+[[:space:]]+[[:digit:]]+)[[:space:]]+R\\b", "ToUnicode")
  descriptor <- pdf_object(descriptor_ref)
  descriptor_name <- pdf_field(descriptor, "/FontName[[:space:]]+/([^[:space:]/<>]+)", "FontName")
  if (!grepl("/Type[[:space:]]+/FontDescriptor[[:space:]/]", descriptor) || descriptor_name != name) {
    stop("Cairo PDF font descriptor does not match its font.")
  }
  file_ref <- pdf_field(descriptor, "/FontFile2[[:space:]]+([[:digit:]]+[[:space:]]+[[:digit:]]+)[[:space:]]+R\\b", "FontFile2")
  font_file <- pdf_object(file_ref)
  unicode_map <- pdf_object(unicode_ref)
  if (!grepl("\\bstream\\b.*\\bendstream\\b", font_file, perl = TRUE) ||
      !grepl("/Length1[[:space:]]+[1-9][[:digit:]]*\\b", font_file, perl = TRUE) ||
      !grepl("\\bstream\\b.*\\bendstream\\b", unicode_map, perl = TRUE)) {
    stop("Cairo PDF lacks an embedded TrueType font stream or Unicode mapping stream.")
  }
  sub("^[A-Z]{6}\\+", "", name)
}, character(1), USE.NAMES = FALSE)
if (!setequal(font_names, expected_fonts)) stop("Unexpected Cairo PDF font substitution: ",
                                              paste(font_names, collapse = ", "))

# Poppler provides an independent audit where available; compute nodes need only R.
poppler <- Sys.which(c("pdfinfo", "pdffonts"))
if (all(nzchar(poppler))) {
  run_poppler <- function(tool, args) {
    result <- suppressWarnings(system2(poppler[[tool]], args, stdout = TRUE, stderr = TRUE,
                                      env = "LC_ALL=C"))
    status <- attr(result, "status")
    if (!is.null(status) && status != 0L) stop(tool, " failed: ", paste(result, collapse = "\n"))
    result
  }
  pdf_info <- run_poppler("pdfinfo", c("-box", shQuote(output("synteny.pdf"))))
  page_count <- grep("^Pages:", pdf_info, value = TRUE)
  media_box <- grep("^MediaBox:", pdf_info, value = TRUE)
  if (length(page_count) != 1L || length(media_box) != 1L) stop("Cannot inspect PDF page geometry.")
  stopifnot(as.integer(trimws(sub("^Pages:", "", page_count))) == 1L)
  box <- as.numeric(strsplit(trimws(sub("^MediaBox:", "", media_box)), "[[:space:]]+")[[1]])
  stopifnot(length(box) == 4L, all(is.finite(box)))
  page_mm <- c(box[3] - box[1], box[4] - box[2]) * 25.4 / 72
  expected_mm <- as.numeric(reported[c("width_mm", "height_mm")])
  if (any(abs(page_mm - expected_mm) >= 0.4)) stop("PDF MediaBox differs from the configured size by >=0.4 mm.")

  font_output <- run_poppler("pdffonts", shQuote(output("synteny.pdf")))
  font_rows <- font_output[seq_along(font_output) > 2L & nzchar(trimws(font_output))]
  font_pattern <- paste0("^([^[:space:]]+)[[:space:]]+(TrueType)[[:space:]]+",
                         "([^[:space:]]+)[[:space:]]+(yes|no)[[:space:]]+(yes|no)",
                         "[[:space:]]+(yes|no)[[:space:]]+[[:digit:]]+[[:space:]]+[[:digit:]]+[[:space:]]*$")
  font_fields <- regmatches(font_rows, regexec(font_pattern, font_rows))
  if (!length(font_fields) || any(lengths(font_fields) != 7L)) {
    stop("PDF must contain only embedded Unicode TrueType fonts; unexpected pdffonts output.")
  }
  font_fields <- do.call(rbind, font_fields)
  font_names <- sub("^[A-Z]{6}\\+", "", font_fields[, 2])
  if (length(font_names) != length(expected_fonts) || !setequal(font_names, expected_fonts) ||
      any(font_fields[, 5] != "yes") || any(font_fields[, 7] != "yes")) {
    stop("PDF fonts must be Liberation Sans regular, italic and bold, embedded with Unicode mappings; found: ",
         paste(font_names, collapse = ", "))
  }
}
cat(sprintf("Verification passed (%s): %d sequences, %d annotations, %d tree-adjacent links; palette ranks, source fields, hashes, 2.2 mm genus markers and embedded-font exports.\n",
             entry, nrow(seqs), nrow(genes), nrow(hits)))
