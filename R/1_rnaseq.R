# 1_rnaseq.R

count_rnaseq <- function(){
  dds <- rnaseq.lf.hypoxia.molidustat::lf_hyp_bay_rnaseq

  SummarizedExperiment::colData(dds) <-
    SummarizedExperiment::colData(dds) |>
    tibble::as_tibble() |>
    dplyr::mutate(
      experiment = factor(experiment),
      oxygen = factor(oxygen, levels = c("21%", "0.5%"), labels = c("N", "H")),
      treatment = factor(treatment, labels = c("DMSO", "BAY")),
      group = stringr::str_c(oxygen, treatment, sep = "."),
      group = factor(group, levels = c("N.DMSO", "H.DMSO", "N.BAY", "H.BAY"))
    ) |>
    as("DataFrame")

  rownames(SummarizedExperiment::colData(dds)) <- SummarizedExperiment::colData(dds)$id

  design <- ~ experiment + group

  dds <-
    DESeq2::DESeqDataSet(
      dds,
      design = design
    )

  keep <- rowSums(DESeq2::counts(dds)) > 1
  dds <- dds[keep, ]
  DESeq2::DESeq(dds)
}

vst_rnaseq <- function(dds){
  vsd <- DESeq2::vst(dds, blind = FALSE)
  SummarizedExperiment::assay(vsd) <-
    limma::removeBatchEffect(SummarizedExperiment::assay(vsd), vsd$experiment)

  DESeq2::plotPCA(
    vsd,
    intgroup = "group",
    returnData = TRUE
  ) |>
    dplyr::mutate(
      label = dplyr::case_when(
        group == "N.DMSO" ~ "21%\nDMSO",
        group == "H.DMSO" ~ "0.5%\nDMSO",
        group == "N.BAY" ~ "21%\nBAY",
        group == "H.BAY" ~ "0.5%\nBAY"
      )
    )
}

plot_rnaseq_pca <- function(pca_data){
  percent_variance <- round(100 * attr(pca_data, "percentVar"))

  pair_clrs <- c(
    N.DMSO = "#b2df8a",
    H.DMSO = "#33a02c",
    N.BAY = "#cab2d6",
    H.BAY = "#6a3d9a"
  )

  ggplot2::ggplot(pca_data) +
    ggplot2::aes(
      x = PC1,
      y = PC2,
      color = group
    ) +
    ggforce::geom_mark_ellipse(
      ggplot2::aes(
        color = group,
        label = label
      ),
      expand = ggplot2::unit(2, "mm"),
      label.fontsize = 6,
      label.fontface = "plain",
      label.family = "Calibri",
      label.hjust = 0.5,
      label.buffer = ggplot2::unit(0, "mm"),
      label.margin = ggplot2::margin(-1.5, -1.5, -1.5, -1.5, "mm"),
      con.type = "none",
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = group),
      pch = 21,
      color = "white",
      size = 2,
      stroke = 0.2,
      show.legend = FALSE
    ) +
    ggplot2::labs(
      x = paste0("PC1: ", percent_variance[1], "% variance"),
      y = paste0("PC2: ", percent_variance[2], "% variance")
    ) +
    ggplot2::scale_fill_manual(
      name = NULL,
      values = pair_clrs,
      labels = c("21% | DMSO", "0.5% | DMSO", "21% | BAY", "0.5% | BAY")
    ) +
    ggplot2::scale_color_manual(
      name = NULL,
      values = pair_clrs,
      labels = c("21% | DMSO", "0.5% | DMSO", "21% | BAY", "0.5% | BAY")
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = c(3, 3))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = c(2, 2))) +
    theme_plots() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.size = ggplot2::unit(1, units = "lines"),
      axis.line = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, size = 0.25)
    )
}

identify_deg <- function(dds, comp){
  alpha <- 0.05
  fc <- 1

  mod_mat <- model.matrix(
    DESeq2::design(dds),
    SummarizedExperiment::colData(dds)
  )

  n.dmso <- colMeans(mod_mat[dds$group == "N.DMSO", ])
  n.bay <- colMeans(mod_mat[dds$group == "N.BAY", ])
  h.dmso <- colMeans(mod_mat[dds$group == "H.DMSO", ])
  h.bay <- colMeans(mod_mat[dds$group == "H.BAY", ])

  con <- eval(parse(text = comp))
  # con <- (h.dmso - n.dmso) - (n.bay - n.dmso)

  annots <-
    tibble::as_tibble(SummarizedExperiment::rowData(dds), rownames = "row") |>
    dplyr::select(row, hgnc_symbol, description)

  DESeq2::results(
    dds,
    contrast = con,
    alpha = alpha,
    lfcThreshold = log(fc, base = 2),
    tidy = TRUE,
    parallel = TRUE
  ) |>
    dplyr::left_join(annots, by = "row") |>
    dplyr::rename(symbol = hgnc_symbol) |>
    dplyr::relocate(symbol, description, .after = "row") |>
    dplyr::arrange(padj)
}

plot_rnaseq_volcano <- function(results, gois = NULL, xlab = NULL, nudge = 8){
  df <-
    tibble::as_tibble(results)

  left <-
    df |>
    dplyr::filter(log2FoldChange < 0) |>
    dplyr::slice_min(stat, n = 10)

  right <-
    df |>
    dplyr::filter(log2FoldChange > 0) |>
    dplyr::slice_max(stat, n = 10)

  ggplot2::ggplot(df) +
    ggplot2::aes(
      x = log2FoldChange,
      y = padj
    ) +
    ggrepel::geom_text_repel(
      data = left,
      ggplot2::aes(
        label = symbol,
        color = symbol %in% gois
      ),
      size = 6/ggplot2::.pt,
      max.overlaps = 20,
      segment.size = 0.1,
      nudge_x = -nudge - left$log2FoldChange,
      hjust = 0,
      direction = "y",
      family = "Calibri",
      segment.color = "black",
      show.legend = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = right,
      ggplot2::aes(
        label = symbol,
        color = symbol %in% gois
      ),
      size = 6/ggplot2::.pt,
      max.overlaps = 20,
      segment.size = 0.1,
      nudge_x = nudge - right$log2FoldChange,
      hjust = 1,
      direction = "y",
      segment.color = "black",
      family = "Calibri",
      show.legend = FALSE
    ) +
    ggplot2::geom_hex(
      bins = c(50, 4),
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_viridis_c(trans = "log10") +
    ggplot2::scale_color_manual(values = c("black", "darkred")) +
    ggplot2::scale_y_continuous(
      trans = c("log10", "reverse"),
      labels = scales::label_log()
    ) +
    ggplot2::scale_x_continuous(
      breaks = scales::pretty_breaks(n = 7),
      limits = c(-nudge - 0.25, nudge + 0.25),
      labels = scales::math_format(2^.x)
    ) +
    ggplot2::labs(
      x = xlab,
      y = "Adjusted p-value"
    ) +
    theme_plots() +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, size = 0.25)
    ) +
    NULL
}

get_unique_symbol_ids <- function(dds){
  tibble::as_tibble(SummarizedExperiment::rowData(dds), rownames = "entrezid") |>
    dplyr::filter(hgnc_symbol != "" & !is.na(hgnc_symbol)) |>
    dplyr::group_by(hgnc_symbol) |>
    dplyr::filter(baseMean == max(baseMean)) |>
    dplyr::pull(entrezid)
}

dds_to_symbols <- function(dds, unique_symbol_ids){
  df <- dds[rownames(dds) %in% unique_symbol_ids, ]
  rownames(df) <- SummarizedExperiment::rowData(df)$hgnc_symbol
  df
}

plot_goi <- function(dds, goi){

  goi <- toupper(goi)

  idx <- which(SummarizedExperiment::rowData(dds)$hgnc_symbol %in% goi)
  if (length(goi) != length(idx)){
    rlang::abort(message = "Suspect one-to-many mapping of symbol to Entrez ID")
  }

  df <-
    dds[idx, ] |>
    SummarizedExperiment::assay() |>
    tibble::as_tibble(rownames = "entrez") |>
    tidyr::pivot_longer(
      -entrez,
      names_to = "id",
      values_to = "count"
    ) |>
    dplyr::left_join(
      dplyr::select(tibble::as_tibble(SummarizedExperiment::colData(dds[idx, ])), id:group),
      by = "id",
      copy = TRUE
    ) |>
    dplyr::left_join(
      dplyr::select(
        tibble::as_tibble(SummarizedExperiment::rowData(dds[idx, ]), rownames = "entrez"),
        entrez,
        symbol = hgnc_symbol
      ),
      by = "entrez",
      copy = TRUE
    ) |>
    dplyr::mutate(
      symbol = factor(symbol, levels = goi),
      oxygen = factor(oxygen, levels = c("N", "H"), labels = c("21%", "0.5%"))
    )

  ggplot2::ggplot(df) +
    ggplot2::facet_wrap(~ symbol, scales = "free_y", nrow = 1) +
    ggplot2::aes(
      x = treatment,
      y = count/1000,
      fill = oxygen
    ) +
    ggplot2::stat_summary(
      geom = "col",
      fun = "mean",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.6,
      show.legend = TRUE
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = oxygen),
      geom = "errorbar",
      fun.data = "mean_se",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.2,
      size = 0.25,
      show.legend = FALSE
    ) +
    # ggplot2::scale_fill_manual(values = clrs, limits = force) +
    ggplot2::labs(
      x = "Treatment",
      y = expression(paste("Count (x", 10^3, ")")),
      fill = NULL
    ) +
    ggplot2::ylim(c(0, NA)) +
    theme_plots() +
    ggplot2::theme(
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    ) +
    NULL
}

get_msigdb_pathways <- function(species = 'Homo sapiens', category = NULL, subcategory = NULL){
  x <-
    msigdbr::msigdbr(species = species, category = category, subcategory = subcategory) |>
    dplyr::select(gs_name, symbol = gene_symbol) |>
    dplyr::group_by(gs_name) |>
    dplyr::mutate(
      gs_name = ifelse(
        stringr::str_detect(gs_name, 'TARGET_GENES'),
        stringr::str_c(
          'GTRD_',
          stringr::str_replace(gs_name, '_TARGET_GENES', '')
        ),
        gs_name
      )
    ) |>
    tidyr::nest() |>
    dplyr::summarise(gs_symbol = unlist(data, recursive = FALSE)) |>
    tibble::deframe()
}

run_gsea <- function(results, pathways){
  set.seed(42)

  rnks <-
    results |>
    dplyr::select(symbol, stat) |>
    dplyr::filter(!(is.na(symbol) | symbol == "")) |>
    dplyr::group_by(symbol) |>
    dplyr::arrange(dplyr::desc(abs(stat)), .by_group = TRUE) |>
    dplyr::slice(1) |>
    dplyr::arrange(stat) |>
    tibble::deframe()

  fgsea::fgsea(
    pathways = pathways,
    stats = rnks,
    nPermSimple = 10000,
    eps = 0,
    BPPARAM = BiocParallel::bpparam()
  ) |>
    tibble::as_tibble() |>
    dplyr::arrange(desc(NES)) |>
    tidyr::separate(pathway, c("source", "pathway"), "_", extra = "merge")
}

plot_gsea <- function(rnaseq_gsea, sources, lbls, vals){
  rnaseq_gsea |>
    dplyr::filter(padj < 0.05) |>
    dplyr::filter(source %in% sources) |>
    dplyr::select(source, pathway, NES) |>
    dplyr::mutate(
      pathway = stringr::str_replace_all(pathway, "_", " "),
      pathway = tolower(pathway)
    ) |>
    dplyr::distinct() |>
    ggplot2::ggplot() +
    ggplot2::aes(
      x = NES,
      y = reorder(pathway, NES)
    ) +
    ggplot2::geom_col(
      ggplot2::aes(fill = NES > 0)
    ) +
    ggplot2::labs(
      x = "Normalized Enrichment Score",
      y = NULL
    ) +
    ggplot2::scale_y_discrete(position = "right") +
    ggplot2::scale_fill_manual(
      name = NULL,
      labels = lbls,
      values = vals
    ) +
    theme_plots() +
    ggplot2::theme(
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    ) +
    NULL
}

plot_gsea_table <- function(df, title, clr) {
  x <-
    df |>
    dplyr::filter(source == "HALLMARK" & padj < 0.05) |>
    dplyr::select(pathway, NES) |>
    dplyr::mutate(pathway = stringr::str_replace_all(pathway, "_", " "))

  gt::gt(x) |>
    gt::tab_header(
      title = title
    ) |>
    gt::cols_label(
      pathway = "PATHWAY"
    ) |>
    gt::fmt_number(
      columns = c("NES")
    ) |>
    gt::data_color(
      columns = NES,
      colors = scales::col_numeric(
        palette = colorRamp(c(clr[2], "white", clr[1]), interpolate = "linear"),
        domain = c(-3.3, 3.3)
      )
    ) |>
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = list(
        gt::cells_title(),
        gt::cells_column_labels()
      )
    ) |>
    # gt::tab_style(
    #   style = gt::cell_borders(sides = "bottom"),
    #   locations = list(
    #     gt::cells_body(rows = x$NES == min(x$NES[x$NES > 0]))
    #   )
    # ) |>
    gt::cols_align("center", c(NES)) |>
    gtExtras::gt_theme_538() |>
    gt::opt_table_font(font = "Calibri")
}

run_tfea <- function(dds){
  # vst if using scale as method
  # dds <- DESeq2::vst(dds, blind = FALSE)

  # remove duplicates
  ids <-
    SummarizedExperiment::rowData(dds) |>
    tibble::as_tibble(rownames = "entrez") |>
    dplyr::filter(!(is.na(hgnc_symbol) | hgnc_symbol == "")) |>
    dplyr::group_by(hgnc_symbol) |>
    dplyr::slice_max(baseMean) |>
    dplyr::select(entrez, hgnc_symbol) |>
    tibble::deframe()

  # make matrix
  mat <-
    dds[rownames(dds) %in% names(ids), ] |>
    SummarizedExperiment::assay() %>%
    magrittr::set_rownames(ids[rownames(.)])

  # run viper
  regulons <-
    dorothea::dorothea_hs |>
    dplyr::filter(confidence == "A")

  dorothea::run_viper(
    input = mat,
    regulons = regulons,
    options = list(
      method = "rank",
      minsize = 1,
      nes = TRUE,
      cores = 4,
      verbose = FALSE
    ),
    tidy = FALSE
  )
}

fit_tfea <- function(dds, tfea){
  # phenotype
  pheno <-
    SummarizedExperiment::colData(dds) |>
    tibble::as_tibble() |>
    dplyr::mutate(group = forcats::fct_relabel(group, tolower))

  #limma
  design <-
    model.matrix(~ 0 + group, data = pheno) %>%
    magrittr::set_colnames(stringr::str_extract(colnames(.), "(?<=group).*"))

  corfit <- limma::duplicateCorrelation(tfea, design, block = pheno$experiment)

  cm <-
    limma::makeContrasts(
      hyp = h.dmso - n.dmso,
      bay = n.bay - n.dmso,
      hyp_bay = h.bay - n.bay,
      int = (h.dmso - n.dmso) - (n.bay - n.dmso),
      levels = design
    )

  limma::lmFit(
    tfea,
    design,
    block = pheno$experiment,
    correlation = corfit$consensus
  ) |>
    limma::contrasts.fit(cm) |>
    limma::eBayes()
}

index_tfea <- function(tfea, comp){
  limma::topTable(
    tfea,
    coef = comp,
    number = Inf
  ) |>
    tibble::as_tibble(rownames = "tf")
}

plot_rnaseq_venn <- function(hyp, bay) {
  nm <- hyp$symbol
  hyp_symbol <-
    hyp |>
    dplyr::filter(padj < 0.05) |>
    dplyr::pull(symbol)
  bay_symbol <-
    bay|>
    dplyr::filter(padj < 0.05) |>
    dplyr::pull(symbol)
  bay_deg <- nm %in% bay_symbol
  hyp_deg <- nm %in% hyp_symbol

  tibble::tibble(nm, hyp_deg, bay_deg) |>
    ggplot2::ggplot() +
    ggplot2::aes(A = hyp_deg, B = bay_deg) +
    ggvenn::geom_venn(
      set_names = c("0.5%", "BAY"),
      digits = 0,
      show_percentage = TRUE,
      fill_color = clrs[c(2, 4)],
      fill_alpha = 0.25,
      stroke_size = 0.25,
      set_name_size = 8/ggplot2::.pt,
      text_size = 6/ggplot2::.pt
    ) +
    ggplot2::annotate(
      geom = "text",
      x = 1,
      y = -1.2,
      label = "20609 total",
      size = 6/ggplot2::.pt
    ) +
    theme_plots() +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      title = "Transcripts"
    ) +
    ggplot2::coord_fixed(
      xlim = c(-1.75, 1.75),
      ylim = c(-1.2, 1.2),
      clip = "off"
    ) +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 8),
    )
}

plot_gsea_venn <- function(hyp, bay) {
  nm <- hyp$pathway
  hyp_pathway <-
    hyp |>
    dplyr::filter(padj < 0.05) |>
    dplyr::pull(pathway)
  bay_pathway <-
    bay|>
    dplyr::filter(padj < 0.05) |>
    dplyr::pull(pathway)
  bay_deg <- nm %in% bay_pathway
  hyp_deg <- nm %in% hyp_pathway

  tibble::tibble(nm, hyp_deg, bay_deg) |>
    ggplot2::ggplot() +
    ggplot2::aes(A = hyp_deg, B = bay_deg) +
    ggvenn::geom_venn(
      set_names = c("0.5%", "BAY"),
      digits = 0,
      show_percentage = TRUE,
      fill_color = clrs[c(2, 4)],
      fill_alpha = 0.25,
      stroke_size = 0.25,
      set_name_size = 8/ggplot2::.pt,
      text_size = 6/ggplot2::.pt
    ) +
    ggplot2::annotate(
      geom = "text",
      x = 1,
      y = -1.2,
      label = "50 total",
      size = 6/ggplot2::.pt
    ) +
    theme_plots() +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      title = "Hallmark gene sets"
    ) +
    ggplot2::coord_fixed(
      xlim = c(-1.75, 1.75),
      ylim = c(-1.2, 1.2),
      clip = "off"
    ) +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 8),
    )
}

plot_tfea_venn <- function(hyp, bay) {
  nm <- hyp$tf
  hyp_pathway <-
    hyp |>
    dplyr::filter(adj.P.Val < 0.05) |>
    dplyr::pull(tf)
  bay_pathway <-
    bay|>
    dplyr::filter(adj.P.Val < 0.05) |>
    dplyr::pull(tf)
  bay_deg <- nm %in% bay_pathway
  hyp_deg <- nm %in% hyp_pathway

  tibble::tibble(nm, hyp_deg, bay_deg) |>
    ggplot2::ggplot() +
    ggplot2::aes(A = hyp_deg, B = bay_deg) +
    ggvenn::geom_venn(
      set_names = c("0.5%", "BAY"),
      digits = 0,
      show_percentage = TRUE,
      fill_color = clrs[c(2, 4)],
      fill_alpha = 0.25,
      stroke_size = 0.25,
      set_name_size = 8/ggplot2::.pt,
      text_size = 6/ggplot2::.pt
    ) +
    ggplot2::annotate(
      geom = "text",
      x = 1,
      y = -1.2,
      label = "96 total",
      size = 6/ggplot2::.pt
    ) +
    theme_plots() +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      title = "Transcription factors"
    ) +
    ggplot2::coord_fixed(
      xlim = c(-1.75, 1.75),
      ylim = c(-1.2, 1.2),
      clip = "off"
    ) +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 8),
    )
}

plot_pathway_volcanoes <- function(deg, pathways, sets, title, nudge = 3) {
  targets <-
    pathways[sets] |>
    unlist() |>
    unique()

  results <-
    deg |>
    dplyr::filter(symbol %in% targets)

  left <-
    results |>
    dplyr::filter(log2FoldChange < 0  & padj < 0.05) |>
    dplyr::slice_min(stat, n = 10)

  right <-
    results |>
    dplyr::filter(log2FoldChange > 0 & padj < 0.05) |>
    dplyr::slice_max(stat, n = 10)

  ggplot2::ggplot(results) +
    ggplot2::aes(
      x = log2FoldChange,
      y = padj
    ) +
    ggrepel::geom_text_repel(
      data = left,
      ggplot2::aes(label = symbol),
      size = 6/ggplot2::.pt,
      max.overlaps = 20,
      segment.size = 0.1,
      nudge_x = -nudge - left$log2FoldChange,
      hjust = 0,
      segment.color = "black",
      direction = "y",
      family = "Calibri",
      show.legend = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = right,
      ggplot2::aes(label = symbol),
      size = 6/ggplot2::.pt,
      max.overlaps = 20,
      segment.size = 0.1,
      nudge_x = nudge - right$log2FoldChange,
      hjust = 0,
      segment.color = "black",
      direction = "y",
      family = "Calibri",
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = subset(results, padj > 0.05),
      pch = 21,
      color = "white",
      fill = "grey80",
      stroke = 0.2
    ) +
    ggplot2::geom_point(
      data = subset(results, log2FoldChange > 0 & padj < 0.05),
      pch = 21,
      color = "white",
      fill = clrs[[2]],
      stroke = 0.2
    ) +
    ggplot2::geom_point(
      data = subset(results, log2FoldChange < 0 & padj < 0.05),
      pch = 21,
      color = "white",
      fill = clrs[[1]],
      stroke = 0.2
    ) +
    ggplot2::scale_y_continuous(
      trans = c("log10", "reverse"),
      labels = scales::label_log()
    ) +
    ggplot2::scale_x_continuous(
      breaks = scales::pretty_breaks(n = 7),
      limits = c(-nudge - 0.25, nudge + 0.25),
      labels = scales::math_format(2^.x)
    ) +
    ggplot2::labs(
      x = "Hypoxia/Normoxia in BAY",
      y = "Adjusted p-value",
      title = title
    ) +
    theme_plots() +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, size = 0.25)
    ) +
    NULL
}

plot_tfea_volcanoes <- function(tf, xlab,colors, nudge = 5) {
  left <-
    tf |>
    dplyr::filter(logFC < 0  & adj.P.Val < 0.05) |>
    dplyr::slice_min(t, n = 10)

  right <-
    tf |>
    dplyr::filter(logFC > 0 & adj.P.Val < 0.05) |>
    dplyr::slice_max(t, n = 10)

  ggplot2::ggplot(tf) +
    ggplot2::aes(
      x = logFC,
      y = adj.P.Val
    ) +
    ggrepel::geom_text_repel(
      data = left,
      ggplot2::aes(label = tf),
      size = 6/ggplot2::.pt,
      max.overlaps = 20,
      segment.size = 0.1,
      nudge_x = -nudge - left$logFC,
      hjust = 0,
      segment.color = "black",
      direction = "y",
      family = "Calibri",
      show.legend = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = right,
      ggplot2::aes(label = tf),
      size = 6/ggplot2::.pt,
      max.overlaps = 20,
      segment.size = 0.1,
      nudge_x = nudge - right$logFC,
      hjust = 0,
      segment.color = "black",
      direction = "y",
      family = "Calibri",
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = subset(tf, adj.P.Val > 0.05),
      pch = 21,
      color = "white",
      fill = "grey80",
      stroke = 0.2
    ) +
    ggplot2::geom_point(
      data = subset(tf, logFC > 0 & adj.P.Val < 0.05),
      pch = 21,
      color = "white",
      fill = colors[[1]],
      stroke = 0.2
    ) +
    ggplot2::geom_point(
      data = subset(tf, logFC < 0 & adj.P.Val < 0.05),
      pch = 21,
      color = "white",
      fill = colors[[2]],
      stroke = 0.2
    ) +
    ggplot2::scale_y_continuous(
      trans = c("log10", "reverse"),
      labels = scales::label_log()
    ) +
    ggplot2::scale_x_continuous(
      breaks = scales::pretty_breaks(n = 7),
      limits = c(-nudge - 0.25, nudge + 0.25),
      labels = scales::math_format(2^.x)
    ) +
    ggplot2::labs(
      x = xlab,
      y = "Adjusted p-value"
    ) +
    theme_plots() +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, size = 0.25)
    ) +
    NULL
}
