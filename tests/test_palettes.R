#!/usr/bin/env Rscript
# Freeze historical colours and validate user-facing palette selection.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
root <- dirname(dirname(normalizePath(script)))
source(file.path(root, "scripts", "palettes.R"))
path <- file.path(root, "config", "gene_palette.csv")
colors <- read.csv(file.path(root, "config", "colors.csv"), stringsAsFactors = FALSE)
checks <- 0L
check <- function(value) {
  stopifnot(isTRUE(value))
  checks <<- checks + 1L
}
fails <- function(expression) {
  check(inherits(tryCatch(force(expression), error = identity), "error"))
}

for (name in c("main-muted", "automatic", "nature", "cell", "science")) {
  palette <- read_gene_palette(path, name)
  check(length(palette) == 42L && !anyDuplicated(palette) &&
          all(grepl("^#[0-9A-F]{6}$", palette)))
}

# These 42 colours are the original automatic-entry palette, in rank order.
automatic <- c(
  "#B84C4C", "#BE762F", "#B39B24", "#43845A", "#328D91", "#4B70A9", "#8056A1",
  "#C46363", "#CA8949", "#C1AB40", "#5B966E", "#4DA0A3", "#6484B7", "#946DAC",
  "#CD7B7B", "#D49C64", "#CDBA5F", "#73A783", "#68B2B5", "#7D99C5", "#A785BD",
  "#D79393", "#DEAE7F", "#D9C97F", "#8CB998", "#84C4C6", "#96AED3", "#BA9DCD",
  "#E1ABAB", "#E8C19A", "#E5D89F", "#A5CAAD", "#A0D6D7", "#AFC3E1", "#CDB6DD",
  "#EBC3C3", "#F1D4B5", "#F0E8BF", "#BEDCC3", "#BDE8E8", "#C8D8EE", "#E0CEEB")
check(identical(read_gene_palette(path, "automatic"), automatic))

# Preserve the old main-entry gene mapping after its single HSV desaturation.
historical <- c(
  "#FF9999", "#CC0000", "#FFBB66", "#CC7A00", "#EF6C00", "#FFFF99", "#FFFF33",
  "#999900", "#99FF99", "#7CFC00", "#00CD00", "#006400", "#66B2FF", "#0033CC",
  "#0000FF", "#C999FF", "#9933FF", "#5C0099", "#00CED1", "#33FFFF", "#008080",
  "#EE6DB3", "#FF1493", "#804000", "#E0E0E0", "#616161")
hsv <- grDevices::rgb2hsv(grDevices::col2rgb(historical))
hsv["s", ] <- hsv["s", ] * 0.70
muted <- grDevices::hsv(hsv["h", ], hsv["s", ], hsv["v", ])
check(identical(read_gene_palette(path, "main-muted")[seq_len(26L)], muted))

presets <- list(
  blue = c("#AECDE3", "#D0E1EE", "#EBF3F8"),
  green = c("#B4D5BD", "#D5E8DA", "#EEF6F0"),
  purple = c("#C7BDDE", "#E0DAED", "#F3F0F8"),
  orange = c("#E8C4AE", "#F2DDCF", "#FAF0E9"),
  original = c("#A6A6A6", "#CECECE", "#EBEBEB"))
intervals <- c(">=90%", "80-90%", "70-80%")
for (name in names(presets)) {
  actual <- resolve_link_colors(colors, name)
  check(identical(actual$palette, name) &&
          identical(actual$colors, setNames(presets[[name]], intervals)))
  check(identical(resolve_link_colors(colors, name, ""), actual))
  custom <- resolve_link_colors(colors, name, " #aBcDeF, #123456 ,#fFeEdD ")
  check(identical(custom$palette, "custom") && identical(custom$colors,
    setNames(c("#ABCDEF", "#123456", "#FFEEDD"), intervals)))
}

fails(read_gene_palette(path, "missing"))
fails(resolve_link_colors(colors, "missing"))
fails(resolve_link_colors(colors, "missing", "#ABCDEF,#123456,#FFEEDD"))
for (custom in c("#123456,#ABCDEF", "#123456,#ABCDEF,#FFEEDD,#010203",
                 "#123456,#ABCDEF,#FFEEDD,", "#123456,,#FFEEDD", ",#ABCDEF,#FFEEDD",
                 "#12345G,#ABCDEF,#FFEEDD", "123456,#ABCDEF,#FFEEDD",
                 "#123,#ABCDEF,#FFEEDD", "#abcdef,#ABCDEF,#FFEEDD", " ")) {
  fails(resolve_link_colors(colors, "blue", custom))
}
cat(sprintf("Palette checks passed: %d.\n", checks))
