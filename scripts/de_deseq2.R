#!/usr/bin/env Rscript
# Differential expression of ObP1 against ObN1 with DESeq2, for one library
# preparation, on two quantification routes:
#   - featureCounts read-pair counts (step 06), with DESeqDataSetFromMatrix;
#   - kallisto estimated counts (step 07), imported with tximport and
#     DESeqDataSetFromTximport, which also uses the effective lengths as offsets.
#
# - The obelisk sense and antisense strands are added up into one gene
#   (Obelisk_000003), as in Zheludev et al. (2024).
# - The 12 rRNA genes are left out: their level reflects the rRNA depletion, not
#   the strain, and featureCounts leaves them almost empty (multimapping reads)
#   while kallisto does not.
# - Design ~ strain, contrast ObP1 vs ObN1: a positive log2 fold change means
#   higher in ObP1. A gene is called DE when padj < alpha; alpha is also the
#   target of the DESeq2 independent filtering.
# - The PCA uses host genes only: the obelisk, absent from ObN1, would make up
#   the first component by itself.
# - size_factors.tsv recomputes median-of-ratios size factors on the counts of
#   both routes, with and without the obelisk, next to plain library-size
#   factors, to see whether the obelisk shifts the normalisation. (For the
#   kallisto route DESeq2 itself uses normalisation factors that also include
#   the transcript lengths.)
#
# Outputs (in --outdir):
#   deseq2_featurecounts.tsv, deseq2_kallisto.tsv   DESeq2 results per gene
#   size_factors.tsv                                 size factors per route and run
#   summary.tsv                                      one row per route
#   pca.png, ma.png, volcano.png, lfc_concordance.png
#
# Usage (inside the in-silico-task environment):
#   Rscript scripts/de_deseq2.R --samples config/samples_A.tsv --prep rRNAdep \
#       --featurecounts results/counts/A/featurecounts.tsv \
#       --kallisto-dir results/quant/A --genes data/ref/genes.tsv \
#       --outdir results/de/A/rRNAdep [--alpha 0.05]

suppressPackageStartupMessages({
    library(DESeq2)
    library(tximport)
    library(ggplot2)
})

OBELISK <- "Obelisk_000003"
OBELISK_ANTISENSE <- "Obelisk_000003_antisense"
STRAINS <- c("ObN1", "ObP1")  # reference level first
MIN_BASE_MEAN <- 10           # genes used to set axis limits and for the concordance

# Figure style, shared with scripts/plot_obelisk_coverage.py.
INK_PRIMARY <- "#0b0b0b"
INK_SECONDARY <- "#52514e"
INK_MUTED <- "#898781"
GRID <- "#e1e0d9"
BASELINE <- "#c3c2b7"
SERIES <- c("#2a78d6", "#eb6834")  # categorical slots 1 and 2
STRAIN_COLORS <- setNames(SERIES, STRAINS)
STATUS_COLORS <- c("Not DE" = BASELINE, "DE host gene" = SERIES[1], "Obelisk-S.s" = SERIES[2])
CONCORDANCE_COLORS <- c("Not DE" = BASELINE, "DE in one route" = SERIES[2], "DE in both routes" = SERIES[1])
SHAPES <- c(inside = 16, above = 2, below = 6)

# ---------------------------------------------------------------------------
# Input
# ---------------------------------------------------------------------------

parse_args <- function(argv) {
    required <- c("samples", "prep", "featurecounts", "kallisto-dir", "genes", "outdir")
    if (length(argv) %% 2 != 0) stop("Arguments must come as --name value pairs")
    keys <- sub("^--", "", argv[c(TRUE, FALSE)])
    unknown <- setdiff(keys, c(required, "alpha"))
    if (length(unknown) > 0) stop("Unknown arguments: ", paste0("--", unknown, collapse = ", "))
    args <- modifyList(list(alpha = "0.05"), as.list(setNames(argv[c(FALSE, TRUE)], keys)))
    missing <- setdiff(required, names(args))
    if (length(missing) > 0) stop("Missing arguments: ", paste0("--", missing, collapse = ", "))
    args$alpha <- as.numeric(args$alpha)
    args
}

read_samples <- function(path, prep) {
    samples <- read.delim(path, colClasses = "character")
    samples <- samples[samples$prep == prep & samples$strain %in% STRAINS, ]
    if (nrow(samples) == 0) stop("No ObP1/ObN1 runs for preparation ", prep)
    samples$strain <- factor(samples$strain, levels = STRAINS)
    rownames(samples) <- samples$run
    samples
}

read_genes <- function(path) {
    read.delim(path, colClasses = "character", quote = "")
}

# Host genes without rRNA, plus the obelisk as a single gene.
analysed_genes <- function(genes) {
    genes$gene_id[genes$biotype != "rRNA" & genes$gene_id != OBELISK_ANTISENSE]
}

merge_obelisk_strands <- function(counts) {
    counts[OBELISK, ] <- counts[OBELISK, ] + counts[OBELISK_ANTISENSE, ]
    counts[rownames(counts) != OBELISK_ANTISENSE, , drop = FALSE]
}

read_featurecounts <- function(path, runs, keep) {
    table <- read.delim(path, comment.char = "#", check.names = FALSE)
    counts <- as.matrix(table[, -(1:6)])
    rownames(counts) <- table$Geneid
    colnames(counts) <- sub("\\.bam$", "", basename(colnames(counts)))
    counts <- merge_obelisk_strands(counts[, runs, drop = FALSE])
    counts[keep, , drop = FALSE]
}

read_kallisto <- function(dir, runs, genes, keep) {
    files <- setNames(file.path(dir, runs, "abundance.tsv"), runs)
    # Both obelisk strands point to one gene, so tximport adds up their counts.
    tx2gene <- data.frame(
        transcript = genes$gene_id,
        gene = ifelse(genes$gene_id == OBELISK_ANTISENSE, OBELISK, genes$gene_id)
    )
    txi <- suppressMessages(tximport(files, type = "kallisto", tx2gene = tx2gene))
    for (assay in c("abundance", "counts", "length")) {
        txi[[assay]] <- txi[[assay]][keep, , drop = FALSE]
    }
    txi
}

# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------

write_tsv <- function(table, path) {
    write.table(table, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

results_table <- function(res, genes) {
    res$gene_id <- rownames(res)
    table <- merge(genes[, c("gene_id", "biotype", "name", "product")], res, by = "gene_id")
    table[order(table$padj, table$pvalue, na.last = TRUE), ]
}

library_factors <- function(counts) {
    totals <- colSums(counts)
    totals / exp(mean(log(totals)))
}

size_factor_table <- function(counts_by_route, samples) {
    rows <- lapply(names(counts_by_route), function(route) {
        counts <- counts_by_route[[route]]
        host <- counts[rownames(counts) != OBELISK, , drop = FALSE]
        data.frame(
            route = route,
            run = colnames(counts),
            sample = samples[colnames(counts), "sample"],
            strain = samples[colnames(counts), "strain"],
            size_factor = estimateSizeFactorsForMatrix(counts),
            size_factor_no_obelisk = estimateSizeFactorsForMatrix(host),
            library_factor = library_factors(counts),
            library_factor_no_obelisk = library_factors(host),
            row.names = NULL
        )
    })
    do.call(rbind, rows)
}

is_de <- function(padj, alpha) !is.na(padj) & padj < alpha

# Host genes with enough counts in both routes, for the log2FC comparison.
concordance_data <- function(results, alpha) {
    fc <- results$featureCounts
    ka <- results$kallisto
    shared <- setdiff(intersect(rownames(fc), rownames(ka)), OBELISK)
    fc <- fc[shared, ]
    ka <- ka[shared, ]
    keep <- fc$baseMean >= MIN_BASE_MEAN & ka$baseMean >= MIN_BASE_MEAN &
        is.finite(fc$log2FoldChange) & is.finite(ka$log2FoldChange)
    n_de <- is_de(fc$padj, alpha) + is_de(ka$padj, alpha)
    status <- c("Not DE", "DE in one route", "DE in both routes")[n_de + 1]
    data.frame(gene_id = shared, featurecounts = fc$log2FoldChange,
               kallisto = ka$log2FoldChange, status = status)[keep, ]
}

summary_table <- function(results, prep, alpha, concordance) {
    rows <- lapply(names(results), function(route) {
        res <- results[[route]]
        host <- res[rownames(res) != OBELISK, ]
        de <- is_de(host$padj, alpha)
        data.frame(
            prep = prep,
            route = route,
            genes = nrow(res),
            host_tested = sum(!is.na(host$padj)),
            host_de = sum(de),
            host_de_up = sum(de & host$log2FoldChange > 0),
            host_de_down = sum(de & host$log2FoldChange < 0),
            obelisk_base_mean = res[OBELISK, "baseMean"],
            obelisk_log2fc = res[OBELISK, "log2FoldChange"],
            obelisk_lfc_se = res[OBELISK, "lfcSE"],
            obelisk_padj = res[OBELISK, "padj"],
            host_de_both_routes = sum(concordance$status == "DE in both routes"),
            lfc_genes = nrow(concordance),
            lfc_spearman = cor(concordance$featurecounts, concordance$kallisto, method = "spearman")
        )
    })
    do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

theme_report <- function() {
    theme_minimal(base_size = 9) +
        theme(
            text = element_text(colour = INK_SECONDARY),
            plot.title = element_text(colour = INK_PRIMARY, size = 10),
            plot.subtitle = element_text(colour = INK_MUTED, size = 8),
            axis.text = element_text(colour = INK_MUTED),
            panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
            panel.grid.minor = element_blank(),
            strip.text = element_text(colour = INK_PRIMARY, hjust = 0),
            legend.position = "top",
            legend.justification = "left",
            plot.background = element_rect(fill = "white", colour = NA)
        )
}

save_plot <- function(plot, path, width = 7.5, height = 3.8) {
    ggsave(path, plot, width = width, height = height, dpi = 200, bg = "white")
}

clamp <- function(x, limit) pmax(pmin(x, limit), -limit)

position_flag <- function(x, limit) {
    ifelse(x > limit, "above", ifelse(x < -limit, "below", "inside"))
}

plot_pca <- function(datasets, samples, prep, path) {
    parts <- lapply(names(datasets), function(route) {
        vsd <- varianceStabilizingTransformation(datasets[[route]], blind = TRUE)
        vsd <- vsd[rownames(vsd) != OBELISK, ]
        pca <- plotPCA(vsd, intgroup = "strain", ntop = 500, returnData = TRUE)
        variance <- round(100 * attr(pca, "percentVar"))
        data.frame(
            panel = sprintf("%s (PC1 %d %%, PC2 %d %%)", route, variance[1], variance[2]),
            PC1 = pca$PC1, PC2 = pca$PC2, strain = pca$strain,
            replicate = samples[as.character(pca$name), "replicate"]
        )
    })
    table <- do.call(rbind, parts)
    plot <- ggplot(table, aes(PC1, PC2, colour = strain)) +
        geom_point(size = 2.5) +
        geom_text(aes(label = replicate), vjust = -1.1, size = 2.8,
                  colour = INK_SECONDARY, show.legend = FALSE) +
        scale_colour_manual(values = STRAIN_COLORS, name = "Strain") +
        scale_y_continuous(expand = expansion(mult = 0.2)) +
        facet_wrap(~ panel, scales = "free") +
        labs(title = sprintf("%s: PCA of host genes", prep),
             subtitle = "Variance-stabilised counts, 500 most variable host genes; labels are replicates") +
        theme_report()
    save_plot(plot, path)
}

# One row per gene and route, with the status used to colour the points.
gene_points <- function(results, alpha) {
    parts <- lapply(names(results), function(route) {
        res <- results[[route]]
        res <- res[res$baseMean > 0, ]
        data.frame(route = route, gene_id = rownames(res), base_mean = res$baseMean,
                   lfc = res$log2FoldChange, padj = res$padj)
    })
    table <- do.call(rbind, parts)
    status <- ifelse(table$gene_id == OBELISK, "Obelisk-S.s",
                     ifelse(is_de(table$padj, alpha), "DE host gene", "Not DE"))
    table$status <- factor(status, levels = names(STATUS_COLORS))
    table[order(table$status), ]  # highlighted points are drawn last
}

# Symmetric log2FC limit set by the host genes, so the obelisk does not squash the axis.
lfc_limit <- function(table) {
    host <- table[table$status != "Obelisk-S.s" & table$base_mean >= MIN_BASE_MEAN, ]
    max(2, ceiling(max(abs(host$lfc[is.finite(host$lfc)]))))
}

obelisk_label <- function(table) {
    obelisk <- table[table$status == "Obelisk-S.s", ]
    obelisk$label <- sprintf("Obelisk-S.s\nlog2FC = %.1f", obelisk$lfc)
    obelisk
}

plot_ma <- function(results, alpha, prep, path) {
    table <- gene_points(results, alpha)
    limit <- lfc_limit(table)
    table$y <- clamp(table$lfc, limit)
    table$where <- position_flag(table$lfc, limit)
    plot <- ggplot(table, aes(base_mean, y)) +
        geom_hline(yintercept = 0, colour = INK_MUTED, linewidth = 0.3) +
        geom_point(aes(colour = status, shape = where), size = 1.2) +
        geom_text(data = obelisk_label(table), aes(label = label), hjust = 1.1, size = 2.6,
                  lineheight = 0.9, colour = INK_PRIMARY) +
        scale_x_log10() +
        scale_colour_manual(values = STATUS_COLORS, name = NULL, drop = FALSE) +
        scale_shape_manual(values = SHAPES, guide = "none") +
        facet_wrap(~ route) +
        labs(title = sprintf("%s: MA plot, ObP1 vs ObN1", prep),
             subtitle = sprintf("DE: padj < %s. Triangles: |log2FC| > %d, drawn at the edge", alpha, limit),
             x = "Mean of normalised counts (log scale)",
             y = "log2 fold change (ObP1 / ObN1)") +
        theme_report()
    save_plot(plot, path)
}

plot_volcano <- function(results, alpha, prep, path) {
    table <- gene_points(results, alpha)
    table <- table[!is.na(table$padj), ]
    limit <- lfc_limit(table)
    table$neg_log_padj <- -log10(table$padj)
    host_y <- table$neg_log_padj[table$status != "Obelisk-S.s" & is.finite(table$neg_log_padj)]
    y_limit <- ceiling(max(-log10(alpha), host_y)) + 1
    clamped <- abs(table$lfc) > limit | table$neg_log_padj > y_limit
    table$x <- clamp(table$lfc, limit)
    table$y <- pmin(table$neg_log_padj, y_limit)
    table$where <- ifelse(clamped, "above", "inside")
    plot <- ggplot(table, aes(x, y)) +
        geom_hline(yintercept = -log10(alpha), colour = INK_MUTED, linewidth = 0.3, linetype = "dashed") +
        geom_point(aes(colour = status, shape = where), size = 1.2) +
        geom_text(data = obelisk_label(table), aes(label = label), hjust = 1.1, size = 2.6,
                  lineheight = 0.9, colour = INK_PRIMARY) +
        scale_colour_manual(values = STATUS_COLORS, name = NULL, drop = FALSE) +
        scale_shape_manual(values = SHAPES, guide = "none") +
        facet_wrap(~ route) +
        labs(title = sprintf("%s: volcano plot, ObP1 vs ObN1", prep),
             subtitle = sprintf(paste("Dashed line: padj = %s. Triangles: beyond |log2FC| = %d or",
                                      "-log10 padj = %d, drawn at the edge"), alpha, limit, y_limit),
             x = "log2 fold change (ObP1 / ObN1)",
             y = "-log10 adjusted p-value") +
        theme_report()
    save_plot(plot, path)
}

plot_concordance <- function(concordance, prep, path) {
    rho <- cor(concordance$featurecounts, concordance$kallisto, method = "spearman")
    concordance$status <- factor(concordance$status, levels = names(CONCORDANCE_COLORS))
    concordance <- concordance[order(concordance$status), ]
    limit <- ceiling(max(abs(c(concordance$featurecounts, concordance$kallisto))))
    plot <- ggplot(concordance, aes(featurecounts, kallisto)) +
        geom_abline(slope = 1, intercept = 0, colour = INK_MUTED, linewidth = 0.3, linetype = "dashed") +
        geom_point(aes(colour = status), size = 1.2) +
        scale_colour_manual(values = CONCORDANCE_COLORS, name = NULL, drop = FALSE) +
        coord_equal(xlim = c(-limit, limit), ylim = c(-limit, limit)) +
        labs(title = sprintf("%s: log2 fold change, featureCounts vs kallisto", prep),
             subtitle = sprintf("%d host genes with mean count >= %d in both routes; Spearman rho = %.3f",
                                nrow(concordance), MIN_BASE_MEAN, rho),
             x = "log2 fold change, featureCounts",
             y = "log2 fold change, kallisto") +
        theme_report()
    save_plot(plot, path, width = 5, height = 5)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

main <- function() {
    args <- parse_args(commandArgs(trailingOnly = TRUE))
    dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
    output <- function(name) file.path(args$outdir, name)

    samples <- read_samples(args$samples, args$prep)
    genes <- read_genes(args$genes)
    keep <- analysed_genes(genes)

    fc_counts <- read_featurecounts(args$featurecounts, samples$run, keep)
    txi <- read_kallisto(args[["kallisto-dir"]], samples$run, genes, keep)
    datasets <- list(
        featureCounts = DESeqDataSetFromMatrix(fc_counts, samples, design = ~ strain),
        kallisto = suppressMessages(DESeqDataSetFromTximport(txi, samples, design = ~ strain))
    )
    datasets <- lapply(datasets, DESeq, quiet = TRUE)
    results <- lapply(datasets, function(dds) {
        as.data.frame(results(dds, contrast = c("strain", "ObP1", "ObN1"), alpha = args$alpha))
    })

    write_tsv(results_table(results$featureCounts, genes), output("deseq2_featurecounts.tsv"))
    write_tsv(results_table(results$kallisto, genes), output("deseq2_kallisto.tsv"))
    write_tsv(size_factor_table(lapply(datasets, counts), samples), output("size_factors.tsv"))
    concordance <- concordance_data(results, args$alpha)
    write_tsv(summary_table(results, args$prep, args$alpha, concordance), output("summary.tsv"))

    plot_pca(datasets, samples, args$prep, output("pca.png"))
    plot_ma(results, args$alpha, args$prep, output("ma.png"))
    plot_volcano(results, args$alpha, args$prep, output("volcano.png"))
    plot_concordance(concordance, args$prep, output("lfc_concordance.png"))
}

main()
