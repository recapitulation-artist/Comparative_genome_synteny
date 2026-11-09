#!/usr/bin/env Rscript
# Render a topology-aligned gene neighbourhood figure for one reproducible entry.
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
source(file.path(dirname(script), "environment.R"))
check_environment()
source(file.path(dirname(script), "palettes.R"))
base_dir <- project_root
invisible(Sys.setlocale("LC_COLLATE", "C"))
args <- commandArgs(trailingOnly = TRUE)
entry <- "main"
entry_explicit <- FALSE
target_gene <- "Nramp"
gene_palette <- "main-muted"
link_palette <- "blue"
custom_link_colors <- ""
input_dir <- NULL
phylogeny_dir <- NULL
legacy_tree <- "reference"
legacy_annotations <- "curated"
usage <- paste("Usage: Rscript scripts/plot_synteny.R [--entry main|curated|auto]",
               "[--target-gene NAME] [--input-dir PATH] [--phylogeny-dir PATH]",
               "[--gene-palette main-muted|automatic|nature|cell|science]",
               "[--link-palette blue|green|purple|orange|original] [--link-colors HEX,HEX,HEX]")
while (length(args)) {
  name <- sub("=.*$", "", args[1])
  if (!name %in% c("--entry", "--target-gene", "--input-dir", "--phylogeny-dir",
                  "--gene-palette", "--link-palette", "--link-colors",
                  "--tree-source", "--annotation-source")) stop(usage)
  if (grepl("=", args[1], fixed = TRUE)) {
    value <- sub("^[^=]+=", "", args[1]); args <- args[-1]
  } else {
    if (length(args) < 2L) stop(usage)
    value <- args[2]; args <- args[-c(1, 2)]
  }
  if (!nzchar(value) && name != "--link-colors") stop(usage)
  if (name == "--entry") { entry <- value; entry_explicit <- TRUE }
  if (name == "--target-gene") target_gene <- value
  if (name == "--gene-palette") gene_palette <- value
  if (name == "--link-palette") link_palette <- value
  if (name == "--link-colors") custom_link_colors <- value
  if (name == "--input-dir") input_dir <- value
  if (name == "--phylogeny-dir") phylogeny_dir <- value
  if (name == "--tree-source") legacy_tree <- value
  if (name == "--annotation-source") legacy_annotations <- value
}
if (!legacy_tree %in% c("reference", "inferred") ||
    !legacy_annotations %in% c("curated", "raw")) stop(usage)
if (!entry_explicit) entry <- if (legacy_annotations == "raw") "auto" else {
  if (legacy_tree == "inferred") "curated" else "main"
}
if (!entry %in% c("main", "curated", "auto")) stop(usage)
main_entry <- entry == "main"
tree_source <- if (main_entry) "reference" else "inferred"
resolve_path <- function(path) if (startsWith(path, "/")) path else file.path(base_dir, path)
relative_path <- function(path) {
  path <- normalizePath(resolve_path(path), mustWork = TRUE)
  prefix <- paste0(normalizePath(base_dir), "/")
  if (startsWith(path, prefix)) substring(path, nchar(prefix) + 1L) else path
}
sha256 <- function(path) {
  value <- system2("sha256sum", shQuote(resolve_path(path)), stdout = TRUE)
  hash <- substr(value[1], 1, 64)
  if (length(hash) != 1L || !grepl("^[0-9a-f]{64}$", hash)) stop("Cannot hash input: ", path)
  hash
}
font_hashes <- read.table(resolve_path("environment/fonts/checksums.sha256"),
                          col.names = c("sha256", "path"), colClasses = "character")
stopifnot(identical(unname(vapply(font_hashes$path, sha256, character(1))), font_hashes$sha256))
Sys.setenv(FONTCONFIG_FILE = normalizePath(resolve_path("config/fonts.conf"), mustWork = TRUE))
out_dir <- file.path(base_dir, "results", entry)
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
if (any(!file.exists(vapply(metadata_paths, resolve_path, character(1))))) {
  stop("Missing plotting input: ", paste(metadata_paths[!file.exists(
    vapply(metadata_paths, resolve_path, character(1)))], collapse = ", "))
}
metadata_paths <- vapply(metadata_paths, relative_path, character(1))
seq_data <- read.csv(resolve_path(sequence_path), stringsAsFactors = FALSE, check.names = FALSE)
gene_input <- read.delim(resolve_path(annotation_path), stringsAsFactors = FALSE,
                         quote = "", comment.char = "", check.names = FALSE, na.strings = character())
gene_data <- gene_input
blast_data <- read.csv(resolve_path(blast_path), stringsAsFactors = FALSE)
genomes <- read.delim(resolve_path(genome_path), stringsAsFactors = FALSE,
                      quote = "", comment.char = "", na.strings = character())
tree <- ape::read.tree(resolve_path(tree_path))
required <- c("seq_id", "start", "end", "strand", "note", "organisms")
if (!all(required %in% names(gene_data))) stop("Gene table lacks required columns: ",
  paste(setdiff(required, names(gene_data)), collapse = ", "))
for (nm in c("seq_id", "strand", "note", "organisms")) gene_data[[nm]] <- trimws(gene_data[[nm]])
seq_data$seq_id <- trimws(seq_data$seq_id)
blast_data$seq_id <- trimws(blast_data$seq_id)
blast_data$seq_id2 <- trimws(blast_data$seq_id2)
tree$tip.label <- trimws(tree$tip.label)
stopifnot(nrow(seq_data) > 1, nrow(gene_data) > 0,
          all(is.finite(seq_data$length)), all(seq_data$length > 0),
          !anyNA(seq_data$seq_id), !anyNA(tree$tip.label),
          !anyDuplicated(seq_data$seq_id), !anyDuplicated(tree$tip.label),
          setequal(tree$tip.label, seq_data$seq_id),
          setequal(gene_data$seq_id, seq_data$seq_id),
          all(gene_data$strand %in% c("+", "-")),
          !anyNA(gene_data[, c("start", "end", "note", "organisms")]),
          all(gene_data$start >= 1), all(gene_data$end >= gene_data$start),
          all(gene_data$end <= seq_data$length[match(gene_data$seq_id, seq_data$seq_id)]),
          all(blast_data$seq_id %in% seq_data$seq_id),
          all(blast_data$seq_id2 %in% seq_data$seq_id),
          all(is.finite(as.matrix(blast_data[, c("pident", "start", "end", "start2", "end2")]))),
          all(blast_data$pident >= 0 & blast_data$pident <= 100),
          all(pmin(blast_data$start, blast_data$end) >= 1),
          all(pmin(blast_data$start2, blast_data$end2) >= 1),
          all(pmax(blast_data$start, blast_data$end) <= seq_data$length[match(blast_data$seq_id, seq_data$seq_id)]),
          all(pmax(blast_data$start2, blast_data$end2) <= seq_data$length[match(blast_data$seq_id2, seq_data$seq_id)]))
genus_map <- unique(gene_data[, c("seq_id", "organisms")])
stopifnot(!anyDuplicated(genus_map$seq_id), !anyDuplicated(genomes$seq_id),
          setequal(genomes$seq_id, seq_data$seq_id), !anyNA(genomes$organism_name),
          all(nzchar(genomes$organism_name)))

style_table <- read.delim(resolve_path("config/figure_style.tsv"), stringsAsFactors = FALSE)
style <- setNames(style_table$value, style_table$parameter)
number <- function(name) as.numeric(style[[name]])
width_mm <- number("width_mm")
font_family <- style[["font_family"]]
text_pt <- number("text_size_pt")
species_pt <- number("species_size_pt")
sequence_pt <- number("sequence_size_pt")
stopifnot(width_mm == 183, number("height_mm") == 170, font_family == "Liberation Sans",
          all(c(text_pt, species_pt, sequence_pt) >= 5),
          all(c(text_pt, species_pt, sequence_pt) <= 7), number("genus_dot_mm") == 2.2)
colors <- read.csv(resolve_path("config/colors.csv"), stringsAsFactors = FALSE)
color_group <- function(group) {
  x <- colors[colors$group == group, ]
  stopifnot(!anyDuplicated(x$label)); setNames(x$color, x$label)
}
hex_colors <- function(x) grDevices::rgb(t(grDevices::col2rgb(x)) / 255)
desaturate <- function(x) {
  hsv <- grDevices::rgb2hsv(grDevices::col2rgb(x))
  hsv["s", ] <- hsv["s", ] * number("main_saturation_factor")
  grDevices::hsv(hsv["h", ], hsv["s", ], hsv["v", ])
}
palette <- read_gene_palette(resolve_path("config/gene_palette.csv"), gene_palette)
palette_capacity <- length(palette)
resolved_links <- resolve_link_colors(colors, link_palette, custom_link_colors)
link_colors <- resolved_links$colors
genus_colors <- color_group("genus")
if (main_entry) genus_colors[] <- desaturate(genus_colors)
# Extra genera receive deterministic colours without changing the historical mapping.
extra_genera <- sort(setdiff(unique(genus_map$organisms), names(genus_colors)))
if (length(extra_genera)) genus_colors <- c(genus_colors, setNames(
  grDevices::hcl.colors(length(extra_genera), "Dark 3"), extra_genera))
genus_colors <- genus_colors[names(genus_colors) %in% genus_map$organisms]
genus_colors[] <- hex_colors(genus_colors)
link_levels <- names(link_colors)
is_true <- function(x) tolower(as.character(x)) %in% c("true", "1", "yes")
if (main_entry) {
  is_target <- gene_data$note == "Uptake-Gene"
} else {
  if (!"is_target" %in% names(gene_data)) stop("Prepared annotations require is_target.")
  if (any(!tolower(as.character(gene_data$is_target)) %in% c("true", "false", "0", "1", "yes", "no"))) {
    stop("Invalid is_target value.")
  }
  is_target <- is_true(gene_data$is_target)
}
if (!any(is_target)) stop("No target genes are present.")
target_category <- if (main_entry) "Uptake-Gene" else paste("Target", target_gene)
if (main_entry) {
  gene_colors <- color_group("gene")
  stopifnot(setequal(gene_data$note, names(gene_colors)))
  gene_colors[] <- hex_colors(desaturate(gene_colors))
  ordinary <- setdiff(names(gene_colors), c("Uptake-Gene", "Unknown"))
  stopifnot(identical(unname(gene_colors[ordinary]),
    head(read_gene_palette(resolve_path("config/gene_palette.csv"), "main-muted"), length(ordinary))))
  if (gene_palette != "main-muted") {
    counts <- table(gene_data$note[gene_data$note %in% ordinary])
    ranked <- names(counts)[order(as.integer(counts), names(counts), decreasing = TRUE, method = "radix")]
    stopifnot(length(ranked) <= palette_capacity)
    gene_colors[ranked] <- palette[seq_along(ranked)]
  }
  gene_data$color_category <- gene_data$note
  labels <- names(gene_colors)
  labels[labels == "Uptake-Gene"] <- paste("Target", target_gene)
  labels[labels == "Unknown"] <- "Unassigned"
  color_table <- data.frame(category = names(gene_colors), label = labels,
    count = as.integer(table(factor(gene_data$color_category, levels = names(gene_colors)))),
    rank = seq_along(gene_colors), color = unname(gene_colors), display_label = labels,
    stringsAsFactors = FALSE)
} else {
  candidate <- gene_data$note
  if (entry == "auto") {
    column <- function(name) if (name %in% names(gene_data)) trimws(gene_data[[name]]) else rep("", nrow(gene_data))
    candidate <- column("gene")
    product <- column("products")
    if (!any(nzchar(product)) && "product" %in% names(gene_data)) product <- column("product")
    missing <- !nzchar(candidate) | candidate %in% c("-", ".", "NA")
    candidate[missing] <- product[missing]
  }
  unknown <- !nzchar(candidate) | tolower(candidate) %in% c("unknown", "other", "na", "-", ".") |
    grepl("hypothetical|uncharacterized|uncharacterised", candidate, ignore.case = TRUE)
  candidate[unknown] <- "Other"
  known <- candidate[!is_target & candidate != "Other"]
  counts <- table(known)
  # Reverse bytewise name order breaks equal-count ties reproducibly.
  sorted_names <- sort(names(counts), decreasing = TRUE)
  retained <- head(sorted_names[order(-as.integer(counts[sorted_names]), method = "radix")], palette_capacity)
  gene_data$color_category <- ifelse(is_target, target_category,
                                    ifelse(candidate %in% retained, candidate, "Other"))
  categories <- c(target_category, retained, "Other")
  gene_colors <- setNames(c("#FFFFFF", palette[seq_along(retained)], "#BDBDBD"), categories)
  labels <- c(target_category, sprintf("%02d %s", seq_along(retained), retained), "Other")
  color_table <- data.frame(category = categories, label = labels,
    count = as.integer(table(factor(gene_data$color_category, levels = categories))),
    rank = c(0L, seq_along(retained), palette_capacity + 1L), color = unname(gene_colors), display_label = labels,
    stringsAsFactors = FALSE)
}
gene_data$color <- unname(gene_colors[gene_data$color_category])
stopifnot(!anyNA(gene_data$color))

# Display the input topology without branch lengths or an inferred biological root.
tree <- ape::ladderize(tree, right = FALSE)
n_seqs <- ape::Ntip(tree)
node_y <- ape::node.height(tree)
depth <- ape::node.depth(tree, method = 2)
node_x <- max(depth) - depth
max_x_tree <- max(node_x)
leaf_pos <- data.frame(seq_id = tree$tip.label, x = node_x[seq_len(n_seqs)],
                       y_pos = node_y[seq_len(n_seqs)])
leaf_pos$organisms <- genus_map$organisms[match(leaf_pos$seq_id, genus_map$seq_id)]
leaf_pos$organism_name <- genomes$organism_name[match(leaf_pos$seq_id, genomes$seq_id)]
stopifnot(setequal(leaf_pos$y_pos, seq_len(n_seqs)), all(leaf_pos$x == max_x_tree),
          !anyNA(leaf_pos$organisms))
edge <- tree$edge
tree_horizontal <- data.frame(x = node_x[edge[, 1]], xend = node_x[edge[, 2]],
                              y = node_y[edge[, 2]], yend = node_y[edge[, 2]])
tree_vertical <- do.call(rbind, lapply(unique(edge[, 1]), function(parent) {
  kids <- edge[edge[, 1] == parent, 2]
  data.frame(x = node_x[parent], xend = node_x[parent],
             y = min(node_y[kids]), yend = max(node_y[kids]))
}))
seq_data$y_pos <- leaf_pos$y_pos[match(seq_data$seq_id, leaf_pos$seq_id)]
gene_data$y_pos <- leaf_pos$y_pos[match(gene_data$seq_id, leaf_pos$seq_id)]
gene_polys <- do.call(rbind, lapply(seq_len(nrow(gene_data)), function(i) {
  g <- gene_data[i, ]
  aw <- min((g$end - g$start) * 0.35, 600)
  xs <- if (g$strand == "+") c(g$start, g$end - aw, g$end, g$end - aw, g$start) else {
    c(g$end, g$start + aw, g$start, g$start + aw, g$end)
  }
  data.frame(x = xs, y = g$y_pos + c(-0.25, -0.25, 0, 0.25, 0.25), group = i)
}))
blast_data$y1 <- leaf_pos$y_pos[match(blast_data$seq_id, leaf_pos$seq_id)]
blast_data$y2 <- leaf_pos$y_pos[match(blast_data$seq_id2, leaf_pos$seq_id)]
hits <- blast_data[abs(blast_data$y1 - blast_data$y2) == 1 & blast_data$pident >= 70, ]
hits$grad <- ifelse(hits$pident >= 90, link_levels[1],
                     ifelse(hits$pident >= 80, link_levels[2], link_levels[3]))
link_df <- if (!nrow(hits)) data.frame(x = numeric(), y = numeric(), group = integer(), grad = character()) else {
  do.call(rbind, lapply(seq_len(nrow(hits)), function(j) {
    h <- hits[j, ]
    xs <- if (h$y1 > h$y2) c(h$start, h$end, h$end2, h$start2) else c(h$start2, h$end2, h$end, h$start)
    data.frame(x = xs, y = c(rep(max(h$y1, h$y2) - 0.25, 2),
      rep(min(h$y1, h$y2) + 0.25, 2)), group = j, grad = h$grad)
  }))
}
stars <- gene_data[is_target, ]
stars$center_x <- (stars$start + stars$end) / 2
legend_columns <- 7L
legend_rows <- ceiling(nrow(color_table) / legend_columns)
bottom_margin <- 30 + 4 * legend_rows
height_mm <- if (main_entry) number("height_mm") else max(number("height_mm"), bottom_margin + 10 + 3 * n_seqs)
row_pitch <- (height_mm - bottom_margin - 10) / n_seqs
row_y <- function(y) bottom_margin + (y - 0.5) * row_pitch
max_region <- max(seq_data$length)
stroke <- function(points) points * 96 / 72
text_gp <- function(size = text_pt, face = "plain") {
  grid::gpar(fontfamily = font_family, fontsize = size, fontface = face, col = "black")
}
text_width <- function(label, size = text_pt, face = "plain") {
  vapply(label, function(x) grid::convertWidth(grid::grobWidth(
    grid::textGrob(x, gp = text_gp(size, face))), "mm", valueOnly = TRUE), numeric(1))
}
fitting_prefix <- function(label, limit, size = species_pt, face = "plain", suffix = "") {
  lower <- 0L; upper <- nchar(label)
  while (lower < upper) {
    middle <- ceiling((lower + upper) / 2)
    if (text_width(paste0(substr(label, 1, middle), suffix), size, face) <= limit) lower <- middle else {
      upper <- middle - 1L
    }
  }
  lower
}
shorten <- function(label, limit, size = species_pt, face = "plain") {
  if (text_width(label, size, face) <= limit) return(label)
  paste0(substr(label, 1, fitting_prefix(label, limit, size, face, "…")), "…")
}
legend_faces <- if (main_entry) ifelse(color_table$category %in% c("Uptake-Gene", "Unknown", "nc/tRNA"),
                                     "plain", "italic") else rep("plain", nrow(color_table))
cell_width <- (width_mm - 4) / legend_columns
make_layout <- function() {
  tree_left <- 2; tree_right <- 24; species_x <- tree_right + 3.5
  species_width <- max(text_width(leaf_pos$organism_name, species_pt, "italic"))
  id_x <- species_x + species_width + 2.5
  id_width <- max(text_width(leaf_pos$seq_id, sequence_pt))
  genome_left <- id_x + id_width + 2.5
  genome_right <- width_mm - 2
  if (genome_right - genome_left < 75) stop("Species names and sequence IDs leave insufficient plot width.")
  list(tree_left = tree_left, tree_right = tree_right, species_x = species_x,
       species_width = species_width, id_x = id_x, id_width = id_width,
       genome_left = genome_left, genome_right = genome_right)
}
label <- function(text, x, y, size = text_pt, face = "plain", just = "left") {
  grid::grid.text(text, x, y, default.units = "native", just = just, gp = text_gp(size, face))
}
line <- function(x, xend, y, yend, color = "#262626", points = number("line_width_pt")) {
  grid::grid.segments(x, y, xend, yend, default.units = "native",
    gp = grid::gpar(col = color, lwd = stroke(points), lineend = "butt"))
}
symbol <- function(x, y, fill) {
  grid::grid.points(rep_len(x, max(length(x), length(y))), rep_len(y, max(length(x), length(y))),
    default.units = "native", pch = 21, size = grid::unit(number("genus_dot_mm"), "mm"),
    gp = grid::gpar(col = "#333333", fill = fill, lwd = stroke(0.25)))
}
arrow_key <- function(x, y, color, target = FALSE) {
  grid::grid.polygon(x + c(0, 2.4, 3.3, 2.4, 0), y + c(-0.8, -0.8, 0, 0.8, 0.8),
    default.units = "native", gp = grid::gpar(fill = color, col = "#333333",
      lwd = stroke(number("gene_line_width_pt"))))
  if (target) grid::grid.points(x + 1.5, y, default.units = "native", pch = 16,
    size = grid::unit(number("target_dot_mm"), "mm"), gp = grid::gpar(col = "black"))
}
draw_figure <- function(layout) {
  grid::pushViewport(grid::viewport(xscale = c(0, width_mm), yscale = c(0, height_mm)))
  tree_x <- function(x) layout$tree_left + x / max_x_tree * (layout$tree_right - layout$tree_left)
  genome_x <- function(x) layout$genome_left + x / max_region * (layout$genome_right - layout$genome_left)
  heading_y <- height_mm - 4
  label(if (main_entry) "Reference tree" else "Inferred tree", layout$tree_left, heading_y, 6.5, "bold")
  label("Species", layout$species_x, heading_y, 6.5, "bold")
  label("Sequence ID", layout$id_x, heading_y, 6.5, "bold")
  label("Gene neighbourhood", layout$genome_left, heading_y, 6.5, "bold")
  line(2, width_mm - 2, height_mm - 7, height_mm - 7, "#8C8C8C", 0.3)
  line(tree_x(tree_horizontal$x), tree_x(tree_horizontal$xend), row_y(tree_horizontal$y), row_y(tree_horizontal$yend))
  line(tree_x(tree_vertical$x), tree_x(tree_vertical$xend), row_y(tree_vertical$y), row_y(tree_vertical$yend))
  symbol(layout$tree_right + 1.2, row_y(leaf_pos$y_pos), genus_colors[leaf_pos$organisms])
  label(leaf_pos$organism_name, layout$species_x, row_y(leaf_pos$y_pos), species_pt, "italic")
  label(leaf_pos$seq_id, layout$id_x, row_y(leaf_pos$y_pos), sequence_pt)
  if (nrow(link_df)) grid::grid.path(genome_x(link_df$x), row_y(link_df$y),
    id = link_df$group, pathId = link_df$group, default.units = "native", rule = "evenodd",
    gp = grid::gpar(fill = link_colors[link_df$grad[!duplicated(link_df$group)]], col = "#B0B0B0", lwd = stroke(0.25)))
  line(genome_x(0), genome_x(seq_data$length), row_y(seq_data$y_pos), row_y(seq_data$y_pos), "#5A5A5A", 0.3)
  grid::grid.polygon(genome_x(gene_polys$x), row_y(gene_polys$y), id = gene_polys$group,
    default.units = "native", gp = grid::gpar(fill = gene_data$color, col = "#333333",
      lwd = stroke(number("gene_line_width_pt")), linejoin = "mitre"))
  grid::grid.points(genome_x(stars$center_x), row_y(stars$y_pos), default.units = "native",
    pch = 16, size = grid::unit(number("target_dot_mm"), "mm"), gp = grid::gpar(col = "black"))
  tick_step <- max(1000, ceiling(max_region / 8 / 1000) * 1000)
  ticks <- seq(0, max_region, by = tick_step)
  axis_y <- bottom_margin - 1.8
  line(genome_x(0), genome_x(max_region), axis_y, axis_y)
  line(genome_x(ticks), genome_x(ticks), axis_y, axis_y - 1)
  label(format(ticks / 1000, trim = TRUE), genome_x(ticks), axis_y - 2.3, just = "centre")
  label("Position in extracted region (kb)", (layout$genome_left + layout$genome_right) / 2,
        axis_y - 6.6, 6.5, just = "centre")
  label("Topology only; branch lengths not shown", 2, axis_y - 5.2, sequence_pt)
  line(2, width_mm - 2, axis_y - 10.4, axis_y - 10.4, "#8C8C8C", 0.3)
  genus_y <- axis_y - 13.7
  label("Genus", 2, genus_y, text_pt, "bold")
  cursor <- 15
  for (genus in names(genus_colors)) {
    symbol(cursor, genus_y, genus_colors[[genus]])
    label(genus, cursor + 1.8, genus_y, species_pt, "italic")
    cursor <- cursor + 1.8 + text_width(genus, species_pt, "italic") + 3.2
  }
  if (cursor > width_mm - 1) stop("Genus legend exceeds figure width.")
  identity_y <- genus_y - 5.4
  label("BLAST identity", 2, identity_y, text_pt, "bold")
  cursor <- 23
  for (i in seq_along(link_levels)) {
    grid::grid.rect(cursor + 1.5, identity_y, width = 3, height = 1.8, default.units = "native",
      gp = grid::gpar(fill = link_colors[[i]], col = "#666666", lwd = stroke(0.25)))
    identity_text <- c("90–100%", "80–<90%", "70–<80%")[i]
    label(identity_text, cursor + 4.2, identity_y, species_pt)
    cursor <- cursor + 4.2 + text_width(identity_text, species_pt) + 4
  }
  title_y <- identity_y - 5
  label(if (main_entry) "Curated annotation" else "Gene categories (ranked by count)",
        2, title_y, text_pt, "bold")
  if (!main_entry) label("Full names: legend.pdf", width_mm - 2, title_y, sequence_pt, just = "right")
  for (i in seq_len(nrow(color_table))) {
    x <- 2 + ((i - 1L) %% legend_columns) * cell_width
    y <- title_y - 4.1 - ((i - 1L) %/% legend_columns) * 4
    arrow_key(x, y, color_table$color[i], color_table$category[i] == target_category)
    label(color_table$display_label[i], x + 4.5, y, species_pt, legend_faces[i])
  }
  grid::popViewport()
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
pdf_file <- file.path(out_dir, "synteny.pdf")
grDevices::cairo_pdf(pdf_file, width = width_mm / 25.4, height = height_mm / 25.4,
                     family = font_family, pointsize = text_pt, bg = "white")
grid::grid.newpage()
layout <- make_layout()
for (i in seq_len(nrow(color_table))) color_table$display_label[i] <- shorten(color_table$label[i],
  cell_width - 5.0, species_pt, legend_faces[i])
tryCatch(draw_figure(layout), finally = grDevices::dev.off())
grDevices::png(file.path(out_dir, "synteny.png"), width = width_mm, height = height_mm, units = "mm",
               res = number("png_dpi"), type = "cairo", family = font_family, pointsize = text_pt, bg = "white")
grid::grid.newpage()
tryCatch(draw_figure(layout), finally = grDevices::dev.off())

if (!main_entry) {
  # Full category names remain searchable vector text, wrapping across pages as needed.
  grDevices::cairo_pdf(file.path(out_dir, "legend.pdf"), width = width_mm / 25.4, height = 250 / 25.4,
    family = font_family, pointsize = text_pt, bg = "white", onefile = TRUE)
  wrap_label <- function(text, limit) {
    result <- character(); rest <- text
    while (text_width(rest, 6) > limit) {
      n <- max(1L, fitting_prefix(rest, limit, 6))
      head <- substr(rest, 1, n)
      spaces <- gregexpr(" ", head, fixed = TRUE)[[1]]
      if (any(spaces > 0) && max(spaces) > n * 0.5) n <- max(spaces)
      result <- c(result, trimws(substr(rest, 1, n)))
      rest <- trimws(substr(rest, n + 1, nchar(rest)))
    }
    c(result, rest)
  }
  new_legend_page <- function() {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(xscale = c(0, width_mm), yscale = c(0, 250)))
    label(paste("Gene colours:", entry, "entry"), 5, 244, 7, "bold")
    label("Categories ranked by count; ties use descending name order.", 5, 239, 6)
  }
  new_legend_page(); y <- 233
  for (i in seq_len(nrow(color_table))) {
    wrapped <- wrap_label(sprintf("%s (n = %d)", color_table$label[i], color_table$count[i]), width_mm - 18)
    required_height <- length(wrapped) * 3 + 2
    if (y - required_height < 5) { grid::popViewport(); new_legend_page(); y <- 233 }
    arrow_key(5, y, color_table$color[i], color_table$category[i] == target_category)
    for (j in seq_along(wrapped)) label(wrapped[j], 10, y - (j - 1L) * 3, 6)
    y <- y - required_height
  }
  grid::popViewport(); grDevices::dev.off()
}
style_table$value[style_table$parameter == "height_mm"] <- as.character(height_mm)
extra_style <- data.frame(parameter = c("entry", "tree_representation", "row_pitch_mm", "minimum_font_pt",
  "maximum_font_pt", "legend_columns", "legend_rows", "gene_palette", "gene_palette_capacity",
  "link_palette", "link_colors"),
  value = c(entry, "cladogram", row_pitch, min(text_pt, species_pt, sequence_pt), 7, legend_columns, legend_rows,
    gene_palette, palette_capacity, resolved_links$palette, paste(unname(link_colors), collapse = ",")))
write.table(rbind(style_table, extra_style), file.path(out_dir, "figure_style.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
caption <- paste0("Gene neighbourhoods of ", n_seqs, " target-associated regions. Rows follow the ",
  if (main_entry) "supplied reference" else "newly inferred default-parameter", " region-tree topology. ",
  "Branch lengths and node supports are not shown; the display root does not specify an outgroup. ",
  "Coloured tip dots indicate genus. Species names are italicized and unique sequence IDs are retained. ",
  if (entry == "auto") "Arrows show selected original protein-coding gene annotations and strand direction. " else
    "Arrows show manually curated annotations and strand direction. ",
  "A black dot marks each target ", target_gene, " gene. ",
  if (main_entry) {
    if (gene_palette == "main-muted") "Gene and genus colours retain their historical mapping with reduced saturation. " else
      "Non-target gene colours are assigned by count, with descending names breaking ties; all curated categories are retained. "
  } else
    paste0("Non-target categories are ranked by gene count, with descending names breaking ties. ",
      "The top ", palette_capacity, " categories receive distinct palette colours; remaining and unassigned genes are grey Other. ",
      "Numbered keys map to complete category names in legend.pdf and gene_colors.tsv. "),
  "Gene palette: ", gene_palette, ". ",
  "Ribbons connect adjacent regions with BLAST nucleotide identity of at least 70%; ribbon colours encode ",
  "the indicated identity intervals (link palette: ", resolved_links$palette, "). ",
  "Coordinates are relative to each extracted region; original genomic ",
  "orientations are retained. Species annotations come from assembly metadata.")
writeLines(caption, file.path(out_dir, "caption.txt"))
row_order <- leaf_pos[order(leaf_pos$y_pos, decreasing = TRUE), c("seq_id", "organisms", "organism_name", "y_pos")]
row_order$length <- seq_data$length[match(row_order$seq_id, seq_data$seq_id)]
row_order$genes <- as.integer(table(gene_data$seq_id)[row_order$seq_id])
row_order$uptake_genes <- as.integer(table(stars$seq_id)[row_order$seq_id])
stopifnot(all(row_order$uptake_genes == 1), nrow(gene_polys) == 5 * nrow(gene_data),
          nrow(link_df) == 4 * nrow(hits))
write.csv(row_order, file.path(out_dir, "sequence_order.csv"), row.names = FALSE)
plotted_genes <- gene_input
plotted_genes$y_pos <- gene_data$y_pos
plotted_genes$color_category <- gene_data$color_category
plotted_genes$color <- gene_data$color
write.table(plotted_genes, file.path(out_dir, "plotted_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(color_table, file.path(out_dir, "gene_colors.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(category = names(link_colors), color = unname(link_colors),
  lower_bound = c(90L, 80L, 70L), upper_bound = c(100L, 90L, 80L),
  upper_inclusive = c(TRUE, FALSE, FALSE)), file.path(out_dir, "link_colors.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(category = names(genus_colors), color = unname(genus_colors),
  count = as.integer(table(factor(leaf_pos$organisms, levels = names(genus_colors)))),
  diameter_mm = number("genus_dot_mm")), file.path(out_dir, "genus_colors.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(hits, file.path(out_dir, "plotted_blast_hits.csv"), row.names = FALSE)
metrics <- c(sequences = nrow(seq_data), tree_tips = n_seqs, genes = nrow(gene_data),
  gene_categories = length(unique(gene_data$color_category)), genera = length(unique(leaf_pos$organisms)),
  uptake_genes = nrow(stars), input_blast_hits = nrow(blast_data), plotted_blast_hits = nrow(hits),
  adjacent_pairs = length(unique(pmin(hits$y1, hits$y2))),
  reverse_hits = sum((hits$end - hits$start) * (hits$end2 - hits$start2) < 0),
  identity_ge90 = sum(hits$pident >= 90), identity_80to90 = sum(hits$pident >= 80 & hits$pident < 90),
  identity_70to80 = sum(hits$pident >= 70 & hits$pident < 80),
  assigned_color_categories = sum(!color_table$category %in% c(target_category, "Other") & color_table$count > 0),
  other_genes = sum(gene_data$color_category == "Other"))
write.table(data.frame(metric = names(metrics), value = unname(metrics)), file.path(out_dir, "metrics.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(input = names(metadata_paths), path = unname(metadata_paths),
  sha256 = unname(vapply(metadata_paths, sha256, character(1))), tree_source = tree_source, entry = entry),
  file.path(out_dir, "plot_metadata.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
cat(sprintf("Rendered %s: %d sequences, %d genes and %d adjacent BLAST hits.\n", entry, n_seqs, nrow(gene_data), nrow(hits)))
cat(sprintf("Results: results/%s/synteny.pdf and results/%s/synteny.png\n", entry, entry))
