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
                    blast_hits = "data/blast_hits.csv", genomes = "data/genomes.tsv", colors = "config/colors.csv",
                    figure_style = "config/figure_style.tsv")
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

style_config <- read.delim(resolve_path("config/figure_style.tsv"), colClasses = "character")
style_report <- read.delim(output("figure_style.tsv"), colClasses = "character")
required_style <- c("width_mm", "height_mm", "png_dpi", "font_family", "text_size_pt",
                    "species_size_pt", "sequence_size_pt", "line_width_pt", "gene_line_width_pt")
derived_style <- c("tree_representation", "row_pitch_mm", "minimum_font_pt", "maximum_font_pt")
stopifnot(identical(names(style_config), c("parameter", "value")),
          identical(names(style_report), c("parameter", "value")),
          !anyDuplicated(style_config$parameter), !anyDuplicated(style_report$parameter),
          setequal(style_config$parameter, required_style),
          setequal(style_report$parameter, c(required_style, derived_style)),
          !anyNA(style_config), !anyNA(style_report))
configured <- setNames(style_config$value, style_config$parameter)
reported <- setNames(style_report$value, style_report$parameter)
numeric_parameters <- setdiff(required_style, "font_family")
configured_numbers <- as.numeric(configured[numeric_parameters])
reported_numbers <- as.numeric(reported[numeric_parameters])
stopifnot(all(is.finite(configured_numbers)), all(configured_numbers > 0),
          isTRUE(all.equal(reported_numbers, configured_numbers)),
          identical(reported[["font_family"]], configured[["font_family"]]),
          identical(reported[["tree_representation"]], "cladogram"))
minimum_font <- as.numeric(reported[["minimum_font_pt"]])
maximum_font <- as.numeric(reported[["maximum_font_pt"]])
font_sizes <- as.numeric(configured[c("text_size_pt", "species_size_pt", "sequence_size_pt")])
stopifnot(is.finite(minimum_font), is.finite(maximum_font),
          minimum_font >= 5, maximum_font <= 7, minimum_font <= maximum_font,
          all(font_sizes >= minimum_font), all(font_sizes <= maximum_font),
          isTRUE(all.equal(as.numeric(reported[["row_pitch_mm"]]),
                           (as.numeric(configured[["height_mm"]]) - 56) / nrow(rows))))

pdf <- readBin(output("synteny.pdf"), "raw", n = 5)
png <- readBin(output("synteny.png"), "raw", n = 24)
stopifnot(identical(rawToChar(pdf), "%PDF-"),
          identical(as.integer(png[1:8]), c(137L, 80L, 78L, 71L, 13L, 10L, 26L, 10L)))
dimensions <- vapply(list(png[17:20], png[21:24]), function(bytes) sum(as.integer(bytes) * 256^(3:0)), numeric(1))
expected_dimensions <- trunc(as.numeric(configured[c("width_mm", "height_mm")]) /
                               25.4 * as.numeric(configured[["png_dpi"]]))
stopifnot(identical(dimensions, expected_dimensions),
          all(file.info(c(output("synteny.pdf"), output("synteny.png")))$size > 10000))

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
  expected_mm <- as.numeric(configured[c("width_mm", "height_mm")])
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
cat(sprintf("Verification passed (%s): 29 sequences, 579 manual annotations, %d tree-adjacent links, input hashes, figure style, PDF with embedded Liberation Sans and PNG.\n",
             tree_source, nrow(hits)))
