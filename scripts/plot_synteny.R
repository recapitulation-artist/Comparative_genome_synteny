#!/usr/bin/env Rscript
# Render an aligned topology and gene-synteny figure from validated inputs.
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
source(file.path(dirname(script), "environment.R"))
check_environment()
base_dir <- project_root
args <- commandArgs(trailingOnly = TRUE)
tree_source <- "reference"
phylogeny_dir <- "results/phylogeny"
usage <- "Usage: Rscript scripts/plot_synteny.R [--tree-source reference|inferred] [--phylogeny-dir PATH]"
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
resolve_path <- function(path) if (startsWith(path, "/")) path else file.path(base_dir, path)
sha256 <- function(path) {
  value <- system2("sha256sum", shQuote(resolve_path(path)), stdout = TRUE)
  hash <- substr(value[1], 1, 64)
  if (length(hash) != 1L || !grepl("^[0-9a-f]{64}$", hash)) stop("Cannot hash input: ", path)
  hash
}
# Use the same bundled TrueType faces on every node, without system substitution.
font_hashes <- read.table(resolve_path("environment/fonts/checksums.sha256"),
                           col.names = c("sha256", "path"), colClasses = "character")
stopifnot(identical(unname(vapply(font_hashes$path, sha256, character(1))), font_hashes$sha256))
Sys.setenv(FONTCONFIG_FILE = normalizePath(resolve_path("config/fonts.conf"), mustWork = TRUE))
out_dir <- file.path(base_dir, "results", tree_source)
inputs <- file.path(base_dir, "data", c("sequences.csv", "genes.tsv", "blast_hits.csv"))
tree_path <- if (tree_source == "reference") "data/phylogeny.nwk" else file.path(phylogeny_dir, "default", "tree.treefile")
tree_file <- resolve_path(tree_path)
if (!file.exists(tree_file)) stop("Missing ", tree_source, " tree: ", tree_path)
tree_file <- normalizePath(tree_file)
root_prefix <- paste0(normalizePath(base_dir), "/")
tree_path <- if (startsWith(tree_file, root_prefix)) substring(tree_file, nchar(root_prefix) + 1L) else tree_file
genomes <- read.delim(file.path(base_dir, "data", "genomes.tsv"), stringsAsFactors = FALSE,
                      quote = "", comment.char = "")
colors <- read.csv(file.path(base_dir, "config", "colors.csv"), stringsAsFactors = FALSE)
color_group <- function(group) {
  x <- colors[colors$group == group, ]
  stopifnot(!anyDuplicated(x$label))
  setNames(x$color, x$label)
}

seq_data <- read.csv(inputs[1], stringsAsFactors = FALSE)
gene_data <- read.delim(inputs[2], stringsAsFactors = FALSE, quote = "", comment.char = "")
blast_data <- read.csv(inputs[3], stringsAsFactors = FALSE)
tree <- ape::read.tree(tree_file)
for (nm in c("seq_id", "strand", "note", "organisms")) {
  gene_data[[nm]] <- trimws(gene_data[[nm]])
}
seq_data$seq_id <- trimws(seq_data$seq_id)
blast_data$seq_id <- trimws(blast_data$seq_id)
blast_data$seq_id2 <- trimws(blast_data$seq_id2)
tree$tip.label <- trimws(tree$tip.label)
stopifnot(nrow(seq_data) > 1, nrow(gene_data) > 0, nrow(blast_data) > 0,
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

# Display topology without branch lengths or an inferred biological root.
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
stopifnot(setequal(leaf_pos$y_pos, seq_len(n_seqs)),
          all(leaf_pos$x == max_x_tree), !anyNA(leaf_pos$organisms))
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

# Preserve strand direction in sharp gene arrows.
gene_h <- 0.25
link_offset <- 0.25
gene_polys <- do.call(rbind, lapply(seq_len(nrow(gene_data)), function(i) {
  g <- gene_data[i, ]
  aw <- min((g$end - g$start) * 0.35, 600)
  xs <- if (g$strand == "+") c(g$start, g$end - aw, g$end, g$end - aw, g$start) else {
    c(g$end, g$start + aw, g$start, g$start + aw, g$end)
  }
  data.frame(x = xs, y = g$y_pos + c(-gene_h, -gene_h, 0, gene_h, gene_h),
             group = i, note = g$note)
}))
gene_colors <- color_group("gene")
link_colors <- color_group("identity")
genus_colors <- color_group("genus")
legend_breaks <- names(gene_colors)
link_levels <- names(link_colors)
stopifnot(setequal(gene_data$note, names(gene_colors)),
          setequal(leaf_pos$organisms, names(genus_colors)))

# Keep only adjacent rows, preserving directed endpoints of reverse BLAST hits.
blast_data$y1 <- leaf_pos$y_pos[match(blast_data$seq_id, leaf_pos$seq_id)]
blast_data$y2 <- leaf_pos$y_pos[match(blast_data$seq_id2, leaf_pos$seq_id)]
hits <- blast_data[abs(blast_data$y1 - blast_data$y2) == 1 & blast_data$pident >= 70, ]
hits$grad <- ifelse(hits$pident >= 90, link_levels[1],
                     ifelse(hits$pident >= 80, link_levels[2], link_levels[3]))
link_df <- if (!nrow(hits)) data.frame(x = numeric(), y = numeric(), group = integer(), grad = character()) else {
  do.call(rbind, lapply(seq_len(nrow(hits)), function(j) {
  h <- hits[j, ]
  if (h$y1 > h$y2) {
    xs <- c(h$start, h$end, h$end2, h$start2)
  } else xs <- c(h$start2, h$end2, h$end, h$start)
  data.frame(x = xs, y = c(rep(max(h$y1, h$y2) - link_offset, 2),
                           rep(min(h$y1, h$y2) + link_offset, 2)), group = j, grad = h$grad)
  }))
}
stars <- gene_data[gene_data$note == "Uptake-Gene", ]
stars$center_x <- (stars$start + stars$end) / 2

# Use physical dimensions and point-sized type at the final publication size.
style_path <- file.path(base_dir, "config", "figure_style.tsv")
style_table <- read.delim(style_path, stringsAsFactors = FALSE)
style <- setNames(style_table$value, style_table$parameter)
number <- function(name) as.numeric(style[[name]])
width_mm <- number("width_mm")
height_mm <- number("height_mm")
font_family <- style[["font_family"]]
text_pt <- number("text_size_pt")
species_pt <- number("species_size_pt")
sequence_pt <- number("sequence_size_pt")
stopifnot(width_mm == 183, height_mm == 170, font_family == "Liberation Sans",
          all(c(text_pt, species_pt, sequence_pt) >= 5),
          all(c(text_pt, species_pt, sequence_pt) <= 7))
# grid line widths use 1/96 inch, whereas the style is specified in points.
stroke <- function(points) points * 96 / 72
text_gp <- function(size = text_pt, face = "plain") {
  grid::gpar(fontfamily = font_family, fontsize = size, fontface = face, col = "black")
}
text_width <- function(label, size = text_pt, face = "plain") {
  vapply(label, function(x) grid::convertWidth(grid::grobWidth(
    grid::textGrob(x, gp = text_gp(size, face))), "mm", valueOnly = TRUE), numeric(1))
}
# Reserve 114 mm for aligned rows, 10 mm for headings and 46 mm for axes/keys.
row_pitch <- (height_mm - 56) / n_seqs
row_y <- function(y) 46 + (y - 0.5) * row_pitch
max_region <- max(seq_data$length)

make_layout <- function() {
  tree_left <- 2
  tree_right <- 24
  species_x <- tree_right + 3
  species_width <- max(text_width(leaf_pos$organism_name, species_pt, "italic"))
  id_x <- species_x + species_width + 2.5
  id_width <- max(text_width(leaf_pos$seq_id, sequence_pt))
  genome_left <- id_x + id_width + 2.5
  genome_right <- width_mm - 2
  stopifnot(genome_right - genome_left >= 100)
  list(tree_left = tree_left, tree_right = tree_right, species_x = species_x,
       species_width = species_width, id_x = id_x, id_width = id_width,
       genome_left = genome_left, genome_right = genome_right)
}
legend_labels <- legend_breaks
legend_labels[legend_labels == "Uptake-Gene"] <- "Target Nramp"
legend_labels[legend_labels == "Unknown"] <- "Unassigned"
legend_faces <- ifelse(legend_breaks %in% c("Uptake-Gene", "Unknown", "nc/tRNA"),
                       "plain", "italic")

draw_figure <- function(layout) {
  grid::pushViewport(grid::viewport(xscale = c(0, width_mm), yscale = c(0, height_mm)))
  label <- function(text, x, y, size = text_pt, face = "plain", just = "left") {
    grid::grid.text(text, x, y, default.units = "native", just = just,
                    gp = text_gp(size, face))
  }
  line <- function(x, xend, y, yend, color = "#262626", points = number("line_width_pt")) {
    grid::grid.segments(x, y, xend, yend, default.units = "native",
                        gp = grid::gpar(col = color, lwd = stroke(points), lineend = "butt"))
  }
  symbol <- function(x, y, fill, size = 1.1) {
    grid::grid.points(rep_len(x, max(length(x), length(y))),
                      rep_len(y, max(length(x), length(y))),
                      default.units = "native", pch = 21, size = grid::unit(size, "mm"),
                      gp = grid::gpar(col = "#333333", fill = fill, lwd = stroke(0.25)))
  }
  tree_x <- function(x) layout$tree_left + x / max_x_tree * (layout$tree_right - layout$tree_left)
  genome_x <- function(x) layout$genome_left + x / max_region *
    (layout$genome_right - layout$genome_left)

  heading_y <- height_mm - 4
  label(if (tree_source == "reference") "Reference tree" else "Inferred tree",
        layout$tree_left, heading_y, 6.5, "bold")
  label("Species", layout$species_x, heading_y, 6.5, "bold")
  label("Sequence ID", layout$id_x, heading_y, 6.5, "bold")
  label("Gene neighbourhood", layout$genome_left, heading_y, 6.5, "bold")
  line(2, width_mm - 2, height_mm - 7, height_mm - 7, "#8C8C8C", 0.3)

  # A cladogram preserves the input display root, without implying an outgroup.
  line(tree_x(tree_horizontal$x), tree_x(tree_horizontal$xend),
       row_y(tree_horizontal$y), row_y(tree_horizontal$yend))
  line(tree_x(tree_vertical$x), tree_x(tree_vertical$xend),
       row_y(tree_vertical$y), row_y(tree_vertical$yend))
  symbol(layout$tree_right + 1.1, row_y(leaf_pos$y_pos), genus_colors[leaf_pos$organisms])
  label(leaf_pos$organism_name, layout$species_x, row_y(leaf_pos$y_pos), species_pt, "italic")
  label(leaf_pos$seq_id, layout$id_x, row_y(leaf_pos$y_pos), sequence_pt)

  if (nrow(link_df)) {
    fills <- link_colors[link_df$grad[!duplicated(link_df$group)]]
    grid::grid.path(genome_x(link_df$x), row_y(link_df$y),
                    id = link_df$group, pathId = link_df$group,
                       default.units = "native", rule = "evenodd",
                       gp = grid::gpar(fill = fills, col = "#B0B0B0", lwd = stroke(0.25)))
  }
  line(genome_x(0), genome_x(seq_data$length), row_y(seq_data$y_pos), row_y(seq_data$y_pos),
       "#5A5A5A", 0.3)
  fills <- gene_colors[gene_polys$note[!duplicated(gene_polys$group)]]
  grid::grid.polygon(genome_x(gene_polys$x), row_y(gene_polys$y), id = gene_polys$group,
                     default.units = "native", gp = grid::gpar(fill = fills, col = "#333333",
                     lwd = stroke(number("gene_line_width_pt")), linejoin = "mitre"))
  # A small black dot identifies the white target arrow without covering its outline.
  grid::grid.points(genome_x(stars$center_x), row_y(stars$y_pos), default.units = "native",
                    pch = 16, size = grid::unit(0.8, "mm"), gp = grid::gpar(col = "black"))

  ticks <- seq(0, max_region, by = 5000)
  axis_y <- 44.2
  line(genome_x(0), genome_x(max_region), axis_y, axis_y)
  line(genome_x(ticks), genome_x(ticks), axis_y, axis_y - 1.0)
  label(format(ticks / 1000, trim = TRUE), genome_x(ticks), axis_y - 2.3, just = "centre")
  label("Position in extracted region (kb)", (layout$genome_left + layout$genome_right) / 2,
        37.6, 6.5, just = "centre")
  label("Topology only; branch lengths not shown", 2, 39.0, sequence_pt)
  line(2, width_mm - 2, 33.8, 33.8, "#8C8C8C", 0.3)

  # Place keys explicitly to prevent clipping or crowding after journal-size export.
  label("Genus", 2, 30.5, text_pt, "bold")
  cursor <- 14
  for (genus in names(genus_colors)) {
    symbol(cursor, 30.5, genus_colors[[genus]])
    label(genus, cursor + 1.5, 30.5, species_pt, "italic")
    cursor <- cursor + 1.5 + text_width(genus, species_pt, "italic") + 3
  }
  stopifnot(cursor <= width_mm - 1)
  label("BLAST identity", 2, 25.1, text_pt, "bold")
  cursor <- 23
  for (i in seq_along(link_levels)) {
    grid::grid.rect(cursor + 1.5, 25.1, width = 3, height = 1.8, default.units = "native",
                    gp = grid::gpar(fill = link_colors[[i]], col = "#666666", lwd = stroke(0.25)))
    identity_text <- c("90–100%", "80–<90%", "70–<80%")[i]
    label(identity_text, cursor + 4.2, 25.1, species_pt)
    cursor <- cursor + 4.2 + text_width(identity_text, species_pt) + 4
  }
  label("Curated annotation", 2, 20.1, text_pt, "bold")
  cell_width <- (width_mm - 4) / 7
  for (i in seq_along(legend_breaks)) {
    col <- (i - 1L) %% 7L
    row <- (i - 1L) %/% 7L
    x <- 2 + col * cell_width
    y <- 16.0 - row * 4
    # Arrow-shaped keys match the geometry in the gene tracks.
    grid::grid.polygon(x + c(0, 2.4, 3.3, 2.4, 0), y + c(-0.8, -0.8, 0, 0.8, 0.8),
                       default.units = "native", gp = grid::gpar(fill = gene_colors[[i]],
                       col = "#333333", lwd = stroke(number("gene_line_width_pt"))))
    if (legend_breaks[i] == "Uptake-Gene") {
      grid::grid.points(x + 1.5, y, default.units = "native", pch = 16,
                        size = grid::unit(0.8, "mm"), gp = grid::gpar(col = "black"))
    }
    label(legend_labels[i], x + 4.5, y, species_pt, legend_faces[i])
    stopifnot(text_width(legend_labels[i], species_pt, legend_faces[i]) + 4.5 < cell_width)
  }
  grid::popViewport()
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
pdf_file <- file.path(out_dir, "synteny.pdf")
png_file <- file.path(out_dir, "synteny.png")
grDevices::cairo_pdf(pdf_file, width = width_mm / 25.4, height = height_mm / 25.4,
                     family = font_family, pointsize = text_pt, bg = "white")
grid::grid.newpage()
layout <- make_layout()
tryCatch(draw_figure(layout), finally = grDevices::dev.off())
grDevices::png(png_file, width = width_mm, height = height_mm, units = "mm",
               res = number("png_dpi"), type = "cairo", family = font_family,
               pointsize = text_pt, bg = "white")
grid::grid.newpage()
tryCatch(draw_figure(layout), finally = grDevices::dev.off())
write.table(rbind(style_table, data.frame(parameter = c("tree_representation", "row_pitch_mm",
                    "minimum_font_pt", "maximum_font_pt"),
                value = c("cladogram", row_pitch, min(text_pt, species_pt, sequence_pt), 6.5))),
            file.path(out_dir, "figure_style.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
caption <- paste0("Gene neighbourhoods of 29 uptake-associated loci. Rows follow the ",
  if (tree_source == "reference") "supplied reference" else "newly inferred default-parameter",
  " region-tree topology. Branch lengths and node supports are not shown; the display root does not ",
  "specify an outgroup. Coloured tip dots indicate genus. Species names are italicized and unique ",
  "sequence IDs are retained. ",
  "Arrows show manually curated annotations and strand direction; a black dot marks each target ",
  "Nramp divalent-metal transporter. Ribbons connect adjacent regions with BLAST nucleotide identity ",
  "of at least 70%; grey levels encode the indicated identity intervals. Coordinates are relative ",
  "to each extracted region; original genomic orientations ",
  "are retained. Species annotations come from assembly metadata.")
writeLines(caption, file.path(out_dir, "caption.txt"))

row_order <- leaf_pos[order(leaf_pos$y_pos, decreasing = TRUE),
                       c("seq_id", "organisms", "organism_name", "y_pos")]
row_order$length <- seq_data$length[match(row_order$seq_id, seq_data$seq_id)]
row_order$genes <- as.integer(table(gene_data$seq_id)[row_order$seq_id])
row_order$uptake_genes <- as.integer(table(stars$seq_id)[row_order$seq_id])
stopifnot(all(row_order$uptake_genes == 1), nrow(gene_polys) == 5 * nrow(gene_data),
          nrow(link_df) == 4 * nrow(hits))
write.csv(row_order, file.path(out_dir, "sequence_order.csv"), row.names = FALSE)
write.csv(hits, file.path(out_dir, "plotted_blast_hits.csv"), row.names = FALSE)
metrics <- c(
  sequences = nrow(seq_data), tree_tips = n_seqs, genes = nrow(gene_data),
  gene_categories = length(unique(gene_data$note)), genera = length(unique(leaf_pos$organisms)),
  uptake_genes = nrow(stars), input_blast_hits = nrow(blast_data),
  plotted_blast_hits = nrow(hits), adjacent_pairs = length(unique(pmin(hits$y1, hits$y2))),
  reverse_hits = sum((hits$end - hits$start) * (hits$end2 - hits$start2) < 0),
  identity_ge90 = sum(hits$pident >= 90),
  identity_80to90 = sum(hits$pident >= 80 & hits$pident < 90),
  identity_70to80 = sum(hits$pident >= 70 & hits$pident < 80))
write.table(data.frame(metric = names(metrics), value = unname(metrics)),
            file.path(out_dir, "metrics.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
metadata_paths <- c(tree = tree_path, annotations = "data/genes.tsv", sequences = "data/sequences.csv",
                    blast_hits = "data/blast_hits.csv", genomes = "data/genomes.tsv", colors = "config/colors.csv",
                    figure_style = "config/figure_style.tsv")
write.table(data.frame(input = names(metadata_paths), path = unname(metadata_paths),
                        sha256 = unname(vapply(metadata_paths, sha256, character(1))), tree_source = tree_source),
             file.path(out_dir, "plot_metadata.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
cat(sprintf("Rendered %d sequences, %d genes and %d adjacent BLAST hits.\n", n_seqs, nrow(gene_data), nrow(hits)))
cat(sprintf("Results: results/%s/synteny.pdf and results/%s/synteny.png\n", tree_source, tree_source))
