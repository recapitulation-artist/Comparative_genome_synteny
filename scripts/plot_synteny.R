#!/usr/bin/env Rscript
# Render an aligned topology and gene-synteny figure from validated inputs.
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
source(file.path(dirname(script), "environment.R"))
check_environment()
suppressPackageStartupMessages(library(ggplot2))
options(device = function(...) grDevices::pdf(NULL))
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
leaf_pos$display_label <- paste(leaf_pos$organism_name, leaf_pos$seq_id, sep = " | ")
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

# Use independent legend grobs with a single fill scale.
y_scale <- function() scale_y_continuous(limits = c(0.5, n_seqs + 0.5), expand = c(0, 0))
base_theme <- theme_minimal(base_size = 14) + theme(
  panel.grid = element_blank(), axis.title.y = element_blank(),
  axis.text.y = element_blank(), axis.ticks.y = element_blank(),
  plot.margin = margin(10, 10, 10, 0), legend.position = "none")
p_tree <- ggplot() +
  geom_segment(data = tree_horizontal, aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.35) +
  geom_segment(data = tree_vertical, aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.35) +
  geom_point(data = leaf_pos, aes(x = x, y = y_pos, color = organisms), size = 3) +
  geom_text(data = leaf_pos, aes(x = x + max_x_tree * 0.10, y = y_pos, label = display_label),
            hjust = 0, size = 3.5) +
  scale_color_manual(values = genus_colors) +
  scale_x_continuous(limits = c(0, max_x_tree * 3.5), expand = c(0, 0)) +
  y_scale() + labs(x = " ") + base_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p_genom <- ggplot() +
  geom_polygon(data = link_df, aes(x = x, y = y, group = group, fill = grad),
               alpha = 0.5, color = "black", linewidth = 0.1, rule = "evenodd") +
  geom_segment(data = seq_data, aes(x = 0, xend = length, y = y_pos, yend = y_pos),
               color = "black", linewidth = 0.3) +
  geom_polygon(data = gene_polys, aes(x = x, y = y, group = group, fill = note),
               color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c(link_colors, gene_colors)) +
  geom_text(data = stars, aes(x = center_x, y = y_pos - 0.2, label = "*"),
            color = "red", size = 20, fontface = "bold", vjust = 0.4, hjust = 0.5) +
  scale_x_continuous(expand = expansion(mult = c(0.003, 0.005)), labels = scales::comma) +
  y_scale() + labs(x = "Genomic position (bp)") + base_theme

extract_legend <- function(p) {
  g <- ggplotGrob(p + theme(legend.position = "bottom", legend.text = element_text(size = 12),
                          legend.title = element_text(size = 13, face = "bold"),
                          legend.key.size = grid::unit(0.55, "cm")))
  idx <- which(g$layout$name == "guide-box-bottom")
  if (length(idx) != 1 || inherits(g$grobs[[idx]], "zeroGrob")) stop("Legend extraction failed.")
  g$grobs[[idx]]
}
genus_legend <- extract_legend(ggplot(leaf_pos, aes(x, y_pos, color = organisms)) +
  geom_point(size = 3) + scale_color_manual(values = genus_colors, name = "Genus") +
  guides(color = guide_legend(nrow = 1)))
identity_legend <- extract_legend(ggplot(data.frame(grad = factor(link_levels, levels = link_levels), x = 1:3),
  aes(x, x, fill = grad)) + geom_tile(alpha = 0.5, color = "black", linewidth = 0.1) +
  scale_fill_manual(values = link_colors, name = "BLAST identity", breaks = link_levels,
                    labels = c(">=90%", "80 to <90%", "70 to <80%"), drop = FALSE) +
  guides(fill = guide_legend(nrow = 1)))
gene_legend <- extract_legend(ggplot(data.frame(note = factor(legend_breaks, levels = legend_breaks),
                                              x = seq_along(legend_breaks)), aes(x, x, fill = note)) +
  geom_tile(color = "black", linewidth = 0.2) +
  scale_fill_manual(values = gene_colors, breaks = legend_breaks, name = "Gene note") +
  guides(fill = guide_legend(ncol = 7, byrow = TRUE)))

tree_grob <- ggplotGrob(p_tree)
genome_grob <- ggplotGrob(p_genom)
stopifnot(length(tree_grob$heights) == length(genome_grob$heights))
aligned_heights <- grid::unit.pmax(tree_grob$heights, genome_grob$heights)
tree_grob$heights <- genome_grob$heights <- aligned_heights
main_grob <- gridExtra::arrangeGrob(tree_grob, genome_grob, ncol = 2, widths = c(1.3, 5))
top_legends <- gridExtra::arrangeGrob(genus_legend, identity_legend, ncol = 2, widths = c(2.4, 1))
heading <- grid::textGrob(sprintf("Cd-Uptake Gene Cluster: %d sequences | %s tree", n_seqs, tree_source),
                          gp = grid::gpar(fontsize = 19, fontface = "bold"))
final_plot <- gridExtra::arrangeGrob(heading, main_grob, top_legends, gene_legend, ncol = 1,
  heights = grid::unit.c(grid::unit(0.5, "in"), grid::unit(1, "null"),
                         grid::unit(0.65, "in"), grid::unit(1.25, "in")))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
pdf_file <- file.path(out_dir, "synteny.pdf")
png_file <- file.path(out_dir, "synteny.png")
ggsave(pdf_file, final_plot, device = grDevices::cairo_pdf, width = 28, height = 30, limitsize = FALSE)
ggsave(png_file, final_plot, device = "png", type = "cairo", bg = "white",
       width = 28, height = 30, dpi = 150, limitsize = FALSE)

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
                    blast_hits = "data/blast_hits.csv", genomes = "data/genomes.tsv", colors = "config/colors.csv")
sha256 <- function(path) {
  value <- system2("sha256sum", shQuote(resolve_path(path)), stdout = TRUE)
  hash <- substr(value[1], 1, 64)
  if (length(hash) != 1L || !grepl("^[0-9a-f]{64}$", hash)) stop("Cannot hash input: ", path)
  hash
}
write.table(data.frame(input = names(metadata_paths), path = unname(metadata_paths),
                        sha256 = unname(vapply(metadata_paths, sha256, character(1))), tree_source = tree_source),
             file.path(out_dir, "plot_metadata.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
cat(sprintf("Rendered %d sequences, %d genes and %d adjacent BLAST hits.\n", n_seqs, nrow(gene_data), nrow(hits)))
cat(sprintf("Results: results/%s/synteny.pdf and results/%s/synteny.png\n", tree_source, tree_source))
