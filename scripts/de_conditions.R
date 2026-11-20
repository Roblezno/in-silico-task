#!/usr/bin/env Rscript
# Differential expression between culture conditions (oxygen and sugar) of
# dataset B (Maddamsetti & You, PRJNA937727) with DESeq2, on the featureCounts
# read-pair counts (step 06). One strain, so the interest is how the conditions
# change the host genes and the obelisk level.
#
# - The obelisk sense and antisense strands are added up into one gene
#   (Obelisk_000003), as in Zheludev et al. (2024).
# - The 12 rRNA genes are left out: their level reflects the rRNA depletion, not
#   the condition, and featureCounts leaves them almost empty (multimapping
#   reads).
# - Design ~ oxygen + sugar, with no interaction (there are too few replicates
#   per condition to estimate it). Two contrasts are extracted from the same
#   fit:
#     oxygen: oxic vs anoxic    (positive log2 fold change = higher in oxic)
#     sugar:  sucrose vs glucose (positive log2 fold change = higher in sucrose)
#   A gene is called DE when padj < alpha; alpha is also the target of the
#   DESeq2 independent filtering.
# - The PCA uses host genes only (the obelisk is removed before choosing the
#   500 most variable genes), so the obelisk does not drive the components.
# - --exclude drops samples (values of the `sample` column) before the analysis,
#   for sensitivity runs; the number of samples and the excluded ones are shown
#   in the plot subtitles.
# - obelisk_counts.tsv lists the raw obelisk count of each run (both strands
#   added up) and the count divided by the DESeq2 size factor.
#
# Outputs (in --outdir):
#   deseq2_oxygen.tsv, deseq2_sugar.tsv   DESeq2 results per gene, with annotation
#   summary.tsv                            one row per contrast
#   obelisk_counts.tsv                     obelisk raw and normalised counts per run
#   pca.png, ma.png, volcano.png, obelisk_counts.png
#
# Usage (inside the in-silico-task environment):
#   Rscript scripts/de_conditions.R --samples config/samples_B.tsv \
#       --featurecounts results/counts/B/featurecounts.tsv \
#       --genes data/ref/genes.tsv --outdir results/de/B/all \
#       [--alpha 0.05] [--exclude glucose_oxic_3,sucrose_oxic_3]

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
})

OBELISK <- "Obelisk_000003"
OBELISK_ANTISENSE <- "Obelisk_000003_antisense"
OXYGEN <- c("anoxic", "oxic")       # reference level first
SUGAR <- c("glucose", "sucrose")    # reference level first
CONDITIONS <- c("glucose anoxic", "sucrose anoxic", "glucose oxic", "sucrose oxic")
MIN_BASE_MEAN <- 10                 # genes used to set axis limits

# Contrasts as c(factor, numerator level, reference level).
CONTRASTS <- list(
    oxygen = c("oxygen", "oxic", "anoxic"),
    sugar = c("sugar", "sucrose", "glucose")
)
CONTRAST_LABELS <- c(oxygen = "Oxygen: oxic vs anoxic", sugar = "Sugar: sucrose vs glucose")

# Figure style, shared with scripts/plot_obelisk_coverage.py.
INK_PRIMARY <- "#0b0b0b"
INK_SECONDARY <- "#52514e"
INK_MUTED <- "#898781"
GRID <- "#e1e0d9"
BASELINE <- "#c3c2b7"
SERIES <- c("#2a78d6", "#eb6834")  # categorical slots 1 and 2
OXYGEN_COLORS <- setNames(SERIES, OXYGEN)
SUGAR_SHAPES <- c(glucose = 16, sucrose = 17)
STATUS_COLORS <- c("Not DE" = BASELINE, "DE host gene" = SERIES[1], "Obelisk-S.s" = SERIES[2])
SHAPES <- c(inside = 16, above = 2, below = 6)

# ---------------------------------------------------------------------------
# Input
# ---------------------------------------------------------------------------

parse_args <- function(argv) {
    required <- c("samples", "featurecounts", "genes", "outdir")
    if (length(argv) %% 2 != 0) stop("Arguments must come as --name value pairs")
    keys <- sub("^--", "", argv[c(TRUE, FALSE)])
    unknown <- setdiff(keys, c(required, "alpha", "exclude"))
    if (length(unknown) > 0) stop("Unknown arguments: ", paste0("--", unknown, collapse = ", "))
    args <- modifyList(list(alpha = "0.05", exclude = ""), as.list(setNames(argv[c(FALSE, TRUE)], keys)))
    missing <- setdiff(required, names(args))
    if (length(missing) > 0) stop("Missing arguments: ", paste0("--", missing, collapse = ", "))
    args$alpha <- as.numeric(args$alpha)
    args
}

split_names <- function(x) {
    names <- trimws(strsplit(x, ",")[[1]])
    names[names != ""]
}

read_samples <- function(path, exclude) {
    samples <- read.delim(path, colClasses = "character")
    unknown <- setdiff(exclude, samples$sample)
    if (length(unknown) > 0) stop("Unknown samples in --exclude: ", paste(unknown, collapse = ", "))
    samples <- samples[!samples$sample %in% exclude, ]
    samples$oxygen <- factor(samples$oxygen, levels = OXYGEN)
    samples$sugar <- factor(samples$sugar, levels = SUGAR)
    if (anyNA(samples$oxygen) || anyNA(samples$sugar)) stop("Unexpected oxygen or sugar values in ", path)
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

is_de <- function(padj, alpha) !is.na(padj) & padj < alpha

contrast_name <- function(spec) sprintf("%s: %s vs %s", spec[1], spec[2], spec[3])

summary_table <- function(results, n_samples, alpha) {
    rows <- lapply(names(results), function(name) {
        res <- results[[name]]
        host <- res[rownames(res) != OBELISK, ]
        de <- is_de(host$padj, alpha)
        data.frame(
            contrast = contrast_name(CONTRASTS[[name]]),
            samples = n_samples,
            genes = nrow(res),
            host_tested = sum(!is.na(host$padj)),
            host_de = sum(de),
            host_de_up = sum(de & host$log2FoldChange > 0),
            host_de_down = sum(de & host$log2FoldChange < 0),
            obelisk_base_mean = res[OBELISK, "baseMean"],
            obelisk_log2fc = res[OBELISK, "log2FoldChange"],
            obelisk_lfc_se = res[OBELISK, "lfcSE"],
            obelisk_padj = res[OBELISK, "padj"]
        )
    })
    do.call(rbind, rows)
}

# Raw obelisk count (strands added up) and count divided by the size factor.
obelisk_counts_table <- function(dds, samples) {
    runs <- colnames(dds)
    data.frame(
        run = runs,
        sample = samples[runs, "sample"],
        oxygen = samples[runs, "oxygen"],
        sugar = samples[runs, "sugar"],
        replicate = samples[runs, "replicate"],
        size_factor = sizeFactors(dds)[runs],
        obelisk_count = counts(dds)[OBELISK, runs],
        obelisk_normalised = counts(dds, normalized = TRUE)[OBELISK, runs],
        row.names = NULL
    )
}

# Text for the plot subtitles: number of samples and the excluded ones.
sample_note <- function(n_samples, excluded) {
    if (length(excluded) == 0) return(sprintf("%d samples", n_samples))
    sprintf("%d samples (excluded: %s)", n_samples, paste(excluded, collapse = ", "))
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

plot_pca <- function(dds, samples, note, path) {
    vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
    vsd <- vsd[rownames(vsd) != OBELISK, ]
    pca <- plotPCA(vsd, intgroup = c("oxygen", "sugar"), ntop = 500, returnData = TRUE)
    variance <- round(100 * attr(pca, "percentVar"))
    pca$replicate <- samples[as.character(pca$name), "replicate"]
    plot <- ggplot(pca, aes(PC1, PC2, colour = oxygen, shape = sugar)) +
        geom_point(size = 2.5) +
        geom_text(aes(label = replicate), vjust = -1.1, size = 2.8,
                  colour = INK_SECONDARY, show.legend = FALSE) +
        scale_colour_manual(values = OXYGEN_COLORS, name = "Oxygen") +
        scale_shape_manual(values = SUGAR_SHAPES, name = "Sugar") +
        scale_x_continuous(expand = expansion(mult = 0.15)) +
        scale_y_continuous(expand = expansion(mult = 0.2)) +
        labs(title = "Dataset B: PCA of host genes",
             subtitle = sprintf(paste("Variance-stabilised counts, 500 most variable host genes;",
                                      "labels are replicates\n%s"), note),
             x = sprintf("PC1 (%d %% of variance)", variance[1]),
             y = sprintf("PC2 (%d %% of variance)", variance[2])) +
        theme_report()
    save_plot(plot, path, width = 5.5, height = 4.2)
}

# One row per gene and contrast, with the status used to colour the points.
gene_points <- function(results, alpha) {
    parts <- lapply(names(results), function(name) {
        res <- results[[name]]
        res <- res[res$baseMean > 0, ]
        data.frame(contrast = CONTRAST_LABELS[[name]], gene_id = rownames(res),
                   base_mean = res$baseMean, lfc = res$log2FoldChange, padj = res$padj)
    })
    table <- do.call(rbind, parts)
    table$contrast <- factor(table$contrast, levels = CONTRAST_LABELS)
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

plot_ma <- function(results, alpha, note, path) {
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
        facet_wrap(~ contrast) +
        labs(title = "Dataset B: MA plot",
             subtitle = sprintf("DE: padj < %s. Triangles: |log2FC| > %d, drawn at the edge\n%s",
                                alpha, limit, note),
             x = "Mean of normalised counts (log scale)",
             y = "log2 fold change") +
        theme_report()
    save_plot(plot, path)
}

plot_volcano <- function(results, alpha, note, path) {
    table <- gene_points(results, alpha)
    table <- table[!is.na(table$padj), ]
    limit <- lfc_limit(table)
    table$neg_log_padj <- -log10(table$padj)
    # Each panel has its own y limit, set by its host genes, so the very strong
    # sugar genes do not squash the oxygen panel.
    host <- table$status != "Obelisk-S.s" & is.finite(table$neg_log_padj)
    panel_max <- tapply(table$neg_log_padj[host], table$contrast[host], max)
    table$y_limit <- ceiling(pmax(-log10(alpha), panel_max[as.character(table$contrast)])) + 1
    clamped <- abs(table$lfc) > limit | table$neg_log_padj > table$y_limit
    table$x <- clamp(table$lfc, limit)
    table$y <- pmin(table$neg_log_padj, table$y_limit)
    table$where <- ifelse(clamped, "above", "inside")
    # The obelisk usually falls inside the cloud of host genes here, so its label
    # goes a quarter of the panel above (or below, if the point is high) on a
    # white background, joined to the point by a line, and is pushed inwards
    # when the point is near a side of the panel.
    labels <- obelisk_label(table)
    shift <- 0.25 * labels$y_limit
    labels$label_y <- ifelse(labels$y > 0.6 * labels$y_limit, labels$y - shift, labels$y + shift)
    labels$hjust <- ifelse(labels$x > 0.6 * limit, 1, ifelse(labels$x < -0.6 * limit, 0, 0.5))
    plot <- ggplot(table, aes(x, y)) +
        geom_hline(yintercept = -log10(alpha), colour = INK_MUTED, linewidth = 0.3, linetype = "dashed") +
        geom_point(aes(colour = status, shape = where), size = 1.2) +
        geom_segment(data = labels, aes(xend = x, yend = label_y), colour = INK_MUTED, linewidth = 0.3) +
        geom_label(data = labels, aes(y = label_y, label = label, hjust = hjust), size = 2.6,
                   lineheight = 0.9, colour = INK_PRIMARY, fill = alpha("white", 0.85),
                   linewidth = 0) +
        scale_colour_manual(values = STATUS_COLORS, name = NULL, drop = FALSE) +
        scale_shape_manual(values = SHAPES, guide = "none") +
        facet_wrap(~ contrast, scales = "free_y") +
        labs(title = "Dataset B: volcano plot",
             subtitle = sprintf(paste("Dashed line: padj = %s. Triangles: beyond |log2FC| = %d,",
                                      "drawn at the edge. Each panel has its own y axis\n%s"),
                                alpha, limit, note),
             x = "log2 fold change",
             y = "-log10 adjusted p-value") +
        theme_report()
    save_plot(plot, path)
}

plot_obelisk_counts <- function(obelisk_counts, note, path) {
    obelisk_counts$condition <- factor(paste(obelisk_counts$sugar, obelisk_counts$oxygen),
                                       levels = CONDITIONS)
    jitter <- position_jitter(width = 0.15, height = 0, seed = 1)
    plot <- ggplot(obelisk_counts, aes(condition, obelisk_normalised, colour = oxygen)) +
        geom_point(size = 2.5, position = jitter) +
        geom_text(aes(label = replicate), vjust = -1.1, size = 2.8,
                  show.legend = FALSE, position = jitter) +
        scale_colour_manual(values = OXYGEN_COLORS, name = "Oxygen") +
        scale_x_discrete(labels = function(x) sub(" ", "\n", x)) +
        scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.12))) +
        labs(title = "Dataset B: obelisk level per sample",
             subtitle = sprintf("Both strands added up, divided by the DESeq2 size factor; labels are replicates\n%s",
                                note),
             x = NULL,
             y = "Normalised obelisk count") +
        theme_report()
    save_plot(plot, path, width = 5.5, height = 4)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

main <- function() {
    args <- parse_args(commandArgs(trailingOnly = TRUE))
    dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
    output <- function(name) file.path(args$outdir, name)

    excluded <- split_names(args$exclude)
    samples <- read_samples(args$samples, excluded)
    genes <- read_genes(args$genes)
    keep <- analysed_genes(genes)
    note <- sample_note(nrow(samples), excluded)

    counts <- read_featurecounts(args$featurecounts, samples$run, keep)
    dds <- DESeqDataSetFromMatrix(counts, samples, design = ~ oxygen + sugar)
    dds <- DESeq(dds, quiet = TRUE)
    results <- lapply(CONTRASTS, function(contrast) {
        as.data.frame(results(dds, contrast = contrast, alpha = args$alpha))
    })

    for (name in names(results)) {
        write_tsv(results_table(results[[name]], genes), output(sprintf("deseq2_%s.tsv", name)))
    }
    summary <- summary_table(results, ncol(dds), args$alpha)
    write_tsv(summary, output("summary.tsv"))
    obelisk_counts <- obelisk_counts_table(dds, samples)
    write_tsv(obelisk_counts, output("obelisk_counts.tsv"))

    plot_pca(dds, samples, note, output("pca.png"))
    plot_ma(results, args$alpha, note, output("ma.png"))
    plot_volcano(results, args$alpha, note, output("volcano.png"))
    plot_obelisk_counts(obelisk_counts, note, output("obelisk_counts.png"))

    for (i in seq_len(nrow(summary))) {
        row <- summary[i, ]
        cat(sprintf("%s (%s): obelisk log2FC = %.2f, padj = %.3g; host DE genes = %d (%d up, %d down)\n",
                    row$contrast, note, row$obelisk_log2fc, row$obelisk_padj,
                    row$host_de, row$host_de_up, row$host_de_down))
    }
}

main()
