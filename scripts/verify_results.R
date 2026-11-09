#!/usr/bin/env Rscript
# Verify scientific counts, selected links and output formats.
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
source(file.path(dirname(script), "environment.R"))
check_environment()
input <- function(name) file.path(project_root, "data", name)
output <- function(name) file.path(project_root, "results", name)
expected <- read.delim(file.path(project_root, "config", "expected_metrics.tsv"))
actual <- read.delim(output("metrics.tsv"))
stopifnot(identical(actual$metric, expected$metric), identical(actual$value, expected$value))

seqs <- read.csv(input("sequences.csv"))
genes <- read.delim(input("genes.tsv"))
rows <- read.csv(output("sequence_order.csv"))
stopifnot(setequal(rows$seq_id, seqs$seq_id), !anyDuplicated(rows$seq_id),
          identical(rows$y_pos, rev(seq_len(nrow(rows)))),
          all(rows$length == seqs$length[match(rows$seq_id, seqs$seq_id)]),
          sum(rows$genes) == nrow(genes), all(rows$uptake_genes == 1))
blast <- read.csv(input("blast_hits.csv"))
y1 <- rows$y_pos[match(blast$seq_id, rows$seq_id)]
y2 <- rows$y_pos[match(blast$seq_id2, rows$seq_id)]
selected <- blast[abs(y1 - y2) == 1 & blast$pident >= 70, ]
hits <- read.csv(output("plotted_blast_hits.csv"))
stopifnot(isTRUE(all.equal(selected, hits[, names(blast)], check.attributes = FALSE)),
          length(unique(pmin(hits$y1, hits$y2))) == nrow(rows) - 1)

pdf <- readBin(output("synteny.pdf"), "raw", n = 5)
png <- readBin(output("synteny.png"), "raw", n = 24)
stopifnot(identical(rawToChar(pdf), "%PDF-"),
          identical(as.integer(png[1:8]), c(137L, 80L, 78L, 71L, 13L, 10L, 26L, 10L)))
dimensions <- vapply(list(png[17:20], png[21:24]), function(bytes) sum(as.integer(bytes) * 256^(3:0)), numeric(1))
stopifnot(identical(dimensions, c(4200, 4500)),
          all(file.info(c(output("synteny.pdf"), output("synteny.png")))$size > 10000))
cat("Verification passed: expected counts, all adjacent pairs, directed BLAST coordinates, PDF and PNG.\n")
