# Fixed, local palettes shared by the renderer and independent output verifier.
gene_palette_names <- c("main-muted", "automatic", "nature", "cell", "science")
link_palette_names <- c("blue", "green", "purple", "orange", "original")

read_gene_palette <- function(path, palette) {
  if (length(palette) != 1L || !palette %in% gene_palette_names) {
    stop("Unknown gene palette: ", palette, ". Choose ", paste(gene_palette_names, collapse = ", "), ".")
  }
  table <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!identical(names(table), c("palette", "rank", "color")) ||
      !setequal(unique(table$palette), gene_palette_names) || anyNA(table) ||
      any(!grepl("^#[0-9A-Fa-f]{6}$", table$color))) stop("Invalid gene palette table: ", path)
  for (name in gene_palette_names) {
    part <- table[table$palette == name, , drop = FALSE]
    if (nrow(part) != 42L || !identical(part$rank, seq_len(42L)) || anyDuplicated(toupper(part$color))) {
      stop("Gene palette must contain 42 ordered, distinct colours: ", name)
    }
  }
  unname(toupper(table$color[table$palette == palette]))
}

resolve_link_colors <- function(colors, palette, custom = "") {
  if (length(palette) != 1L || !palette %in% link_palette_names) {
    stop("Unknown link palette: ", palette, ". Choose ", paste(link_palette_names, collapse = ", "), ".")
  }
  levels <- colors$label[colors$group == "identity"]
  if (!identical(levels, c(">=90%", "80-90%", "70-80%"))) stop("Invalid identity intervals.")
  group <- if (palette == "original") "identity" else paste0("identity_", palette)
  part <- colors[colors$group == group, , drop = FALSE]
  if (!identical(part$label, levels) || any(!grepl("^#[0-9A-Fa-f]{6}$", part$color)) ||
      anyDuplicated(toupper(part$color))) stop("Invalid link palette table: ", palette)
  values <- toupper(part$color)
  if (nzchar(custom)) {
    # Appending a sentinel preserves an empty trailing field for validation.
    values <- trimws(head(strsplit(paste0(custom, ",END"), ",", fixed = TRUE)[[1]], -1L))
    if (length(values) != 3L || any(!grepl("^#[0-9A-Fa-f]{6}$", values)) ||
        anyDuplicated(toupper(values))) {
      stop("--link-colors requires three distinct #RRGGBB colours, ordered high, medium, low identity.")
    }
    values <- toupper(values)
    palette <- "custom"
  }
  list(palette = palette, colors = setNames(values, levels))
}
