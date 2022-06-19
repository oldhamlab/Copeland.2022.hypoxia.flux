# 1_metab.R

new_tbl_se <- function(
    tbl,
    a_data,
    f_names,
    f_data = NULL,
    s_names,
    s_data = NULL
) {
  structure(
    tbl,
    class = c("tbl_se", class(tbl)),
    a_data = a_data,
    f_names = f_names,
    s_names = s_names,
    f_data = f_data,
    s_data = s_data
  )
}

tbl_to_se <- function(tbl_se, assay_name){
  assay_data <-
    tbl_se |>
    dplyr::select(
      attr(tbl_se, "f_names"),
      attr(tbl_se, "s_names"),
      attr(tbl_se, "a_data")
    ) |>
    tidyr::pivot_wider(
      names_from = attr(tbl_se, "s_names"),
      values_from = attr(tbl_se, "a_data")
    ) |>
    tibble::column_to_rownames(attr(tbl_se, "f_names"))

  feature_data <-
    tbl_se |>
    dplyr::select(
      attr(tbl_se, "f_names"),
      attr(tbl_se, "f_data")
    ) |>
    dplyr::group_by(!!rlang::sym(attr(tbl_se, "f_names"))) |>
    dplyr::summarise(
      metabolite = unique(metabolite),
      mz = mean(mz, na.rm = TRUE),
      mz_min = min(mz, na.rm = TRUE),
      mz_max = max(mz, na.rm = TRUE),
      rt = mean(rt, na.rm = TRUE),
      rt_min = min(rt, na.rm = TRUE),
      rt_max = max(rt, na.rm = TRUE)
    ) |>
    tibble::column_to_rownames(attr(tbl_se, "f_names")) %>%
    .[match(rownames(assay_data), rownames(.)), ]

  sample_data <-
    tbl_se |>
    dplyr::select(
      attr(tbl_se, "s_names"),
      attr(tbl_se, "s_data")
    ) |>
    dplyr::distinct() |>
    tibble::column_to_rownames(attr(tbl_se, "s_names")) %>%
    .[match(colnames(assay_data), rownames(.)), ]

  SummarizedExperiment::SummarizedExperiment(
    assays = assay_data,
    rowData = feature_data,
    colData = sample_data
  )
}

format_metab_tar <- function(files){
  readxl::read_excel(
    files,
    sheet = 2
  ) |>
    dplyr::select(
      id = `Raw File Name`,
      sample = `Sample ID`,
      metabolite = `Compound Name`,
      mz = `Detected Mass`,
      rt = RT,
      area = `Peak Area`
    ) |>
    dplyr::left_join(wmo::hmdb_mappings, by = "metabolite") |>
    dplyr::mutate(
      metabolite = dplyr::case_when(
        metabolite == "glyceraldehyde 3-phosphate" ~ "GAP",
        TRUE ~ metabolite
      )
    ) |>
    dplyr::filter(!is.na(hmdb)) |>
    dplyr::filter(!stringr::str_detect(hmdb, "ISTD")) |>
    dplyr::mutate(
      id = stringr::str_c("S", sprintf("%02d", id)),
      type = dplyr::case_when(
        sample == "water" ~ "blank",
        stringr::str_detect(sample, "mix") ~ "qc",
        TRUE ~ "sample"
      ),
      oxygen = stringr::str_extract(sample, "hyp|norm"),
      oxygen = factor(oxygen, levels = c("norm", "hyp"), labels = c("N", "H")),
      treatment = stringr::str_extract(sample, "bay|dmso"),
      treatment = factor(treatment, levels = c("dmso", "bay"), labels = c("DMSO", "BAY")),
      group = stringr::str_c(oxygen, treatment, sep = "."),
      group = factor(group, levels = c("N.DMSO", "N.BAY", "H.DMSO", "H.BAY")),
      replicate = stringr::str_extract(sample, "(?<=\\.)\\w{1}$"),
      mz = replace(mz, mz == "N/F", NA_real_),
      mz = as.numeric(mz),
      dplyr::across(c(rt, area), ~replace(., . == 0, NA_real_))
    ) |>
    dplyr::select(-sample) |>
    new_tbl_se(
      a_data = "area",
      f_names = "hmdb",
      f_data = c("metabolite", "mz", "rt"),
      s_names = "id",
      s_data = c("type", "oxygen", "treatment", "replicate", "group")
    ) |>
    tbl_to_se()
}

remove_missing_metab <- function(raw){
  qc <- SummarizedExperiment::assay(raw[, raw$type == "qc"])
  missing <- names(which(apply(qc, 1, function(x) sum(is.na(x))) > 0))
  raw[rownames(raw) %nin% missing, raw$type %nin% c("water", "blank")]
}

prepare_assay_data <- function(se){
  SummarizedExperiment::assay(se) |>
    tibble::rownames_to_column("hmdb") |>
    tidyr::pivot_longer(-hmdb, names_to = "sample", values_to = "value") |>
    dplyr::mutate(
      value = log(value),
      run_order = as.numeric(stringr::str_extract(sample, "\\d{2}"))
    ) |>
    dplyr::group_by(hmdb) |>
    tidyr::nest()
}

correct_drift <- function(missing){
  models <-
    missing[, missing$type == "qc"] |>
    prepare_assay_data() |>
    dplyr::mutate(
      model = map(data, ~smooth.spline(x = .x$run_order, y = .x$value, spar = 0.2)),
      mean = map_dbl(data, ~mean(.x$value))
    ) |>
    dplyr::select(-data)

  corrected <-
    missing |>
    prepare_assay_data() |>
    dplyr::left_join(models, by = "hmdb") |>
    dplyr::mutate(pred = purrr::map2(model, data, ~predict(.x, .y$run_order)$y)) |>
    tidyr::unnest(c(data, pred)) |>
    dplyr::mutate(corr = value + mean - pred) |>
    dplyr::select(hmdb, sample, corr) |>
    tidyr::pivot_wider(names_from = sample, values_from = corr) |>
    tibble::column_to_rownames("hmdb") |>
    exp()

  SummarizedExperiment::assay(missing) <- corrected
  missing
}

quality_control <- function(drift){
  rsd <-
    drift[, drift$type == "qc"] |>
    SummarizedExperiment::assay() |>
    apply(1, function(x) 1.4826 * mad(x, na.rm = TRUE) / median(x, na.rm = TRUE))

  mad_qc <-
    drift[, drift$type == "qc"] |>
    SummarizedExperiment::assay() |>
    apply(1, function(x) mad(x, na.rm = TRUE))

  ref <-
    drift[, drift$type == "qc"] |>
    SummarizedExperiment::assay() |>
    apply(1, function(x) median(x, na.rm = TRUE))

  mad_s <-
    drift[, drift$type == "sample"] |>
    SummarizedExperiment::assay() |>
    apply(1, function(x) mad(x, na.rm = TRUE))

  d_ratio <- mad_qc / mad_s

  SummarizedExperiment::rowData(drift)$rsd <- rsd
  SummarizedExperiment::rowData(drift)$rsd <- d_ratio
  SummarizedExperiment::rowData(drift)$good <- rsd < 0.2 & d_ratio < 0.4
  SummarizedExperiment::rowData(drift)$reference <- ref

  drift[SummarizedExperiment::rowData(drift)$good == TRUE, drift$type != "qc"]
}

impute_missing <- function(qc){
  set.seed(42)
  SummarizedExperiment::assay(qc) <-
    missForest::missForest(
      t(SummarizedExperiment::assay(qc)),
      maxiter = 10
    )$ximp |>
    t()
  qc
}

pqn <- function(imputed){
  mat <- SummarizedExperiment::assay(imputed)
  quotients <- mat / SummarizedExperiment::rowData(imputed)$reference
  quotient_medians <- apply(quotients, 2, median)
  SummarizedExperiment::assay(imputed) <- t(t(mat) / quotient_medians)
  imputed
}

log_transform <- function(pqn){
  SummarizedExperiment::assay(pqn) <- log(SummarizedExperiment::assay(pqn), base = 2)
  pqn
}

annot_metabs <- function(se){
  fdata <-
    SummarizedExperiment::rowData(se) %>%
    data.frame() %>%
    tibble::as_tibble(rownames = "HMDB")

  ah <- AnnotationHub::AnnotationHub()
  df <-
    ah[["AH91792"]] %>%
    dplyr::filter(HMDB %in% fdata$HMDB) %>%
    dplyr::group_by(HMDB) %>%
    dplyr::arrange(HMDB, KEGG, ChEBI, .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::select(-Name)

  annot_fdata <-
    dplyr::left_join(fdata, df, by ="HMDB") %>%
    tibble::column_to_rownames("HMDB")

  SummarizedExperiment::rowData(se) <- annot_fdata
  se
}

calc_metab_pca <- function(clean) {
  input <-
    SummarizedExperiment::assay(clean) |>
    t() |>
    scale() |>
    t()

  pheno <- SummarizedExperiment::colData(clean)

  design <- model.matrix(~ 0 + group, data = pheno)
  colnames(design) <- stringr::str_extract(colnames(design), "(?<=group).*")

  limma::removeBatchEffect(input, batch = pheno$replicate, design = design) |>
    t() |>
    pcaMethods::pca(scale = "none", center = TRUE)
}

plot_metab_pca <- function(clean, df) {
  percent_variance <- round(100 * c(df@R2[[1]], df@R2[[2]]))

  pair_clrs <- c(
    N.DMSO = "#b2df8a",
    H.DMSO = "#33a02c",
    N.BAY = "#cab2d6",
    H.BAY = "#6a3d9a"
  )

  merge(pcaMethods::scores(df), SummarizedExperiment::colData(clean), by = 0) |>
    dplyr::mutate(
      label = dplyr::case_when(
        group == "N.DMSO" ~ "21%\nDMSO",
        group == "H.DMSO" ~ "0.5%\nDMSO",
        group == "N.BAY" ~ "21%\nBAY",
        group == "H.BAY" ~ "0.5%\nBAY"
      )
    ) |>
    ggplot2::ggplot() +
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
      label.fontsize = 5,
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
      legend.key.size = ggplot2::unit(1, units = "lines")
    )
}

fit_metab_limma <- function(clean){
  input <-
    SummarizedExperiment::assay(clean) |>
    log() |>
    t() |>
    scale() |>
    t()

  pheno <- SummarizedExperiment::colData(clean)

  design <- model.matrix(~ 0 + group, data = pheno)
  colnames(design) <- tolower(stringr::str_extract(colnames(design), "(?<=group).*"))

  corfit <- limma::duplicateCorrelation(input, design, block = pheno$replicate)

  cm <-
    limma::makeContrasts(
      hyp = h.dmso - n.dmso,
      bay = n.bay - n.dmso,
      hyp_bay = h.bay - n.bay,
      int = (h.dmso - n.dmso) - (n.bay - n.dmso),
      levels = design
    )

  fit <-
    limma::lmFit(
      input,
      design,
      block = pheno$replicate,
      correlation = corfit$consensus
    ) |>
    limma::contrasts.fit(cm) |>
    limma::eBayes()
}

index_metab_limma <- function(se, res, comp){
  limma::topTable(
    res,
    number = Inf,
    p.value = 1,
    coef = comp
  ) |>
    tibble::as_tibble(rownames = "hmdb") |>
    dplyr::left_join(
      tibble::as_tibble(SummarizedExperiment::rowData(se), rownames = "hmdb"),
      by = "hmdb"
    )
}

plot_metab_volcano <- function(results, mois = NULL, colors = NULL, xlab = NULL){
  left <-
    results |>
    dplyr::filter(logFC < 0) |>
    dplyr::slice_min(t, n = 10)

  right <-
    results |>
    dplyr::filter(logFC > 0) |>
    dplyr::slice_max(t, n = 10)

  ggplot2::ggplot(results) +
    ggplot2::aes(
      x = logFC,
      y = adj.P.Val
    ) +
    ggrepel::geom_text_repel(
      data = left,
      ggplot2::aes(
        label = metabolite,
        color = metabolite %in% mois
      ),
      size = 5/ggplot2::.pt,
      max.overlaps = 20,
      segment.size = 0.1,
      # nudge_x = -4,
      nudge_x = -5.5 - left$logFC,
      hjust = 0,
      segment.color = "black",
      direction = "y",
      family = "Calibri",
      show.legend = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = right,
      ggplot2::aes(
        label = metabolite,
        color = metabolite %in% mois
      ),
      size = 5/ggplot2::.pt,
      max.overlaps = 20,
      segment.size = 0.1,
      segment.color = "black",
      # nudge_x = 4.5,
      nudge_x = 5.5 - right$logFC,
      hjust = 1,
      direction = "y",
      family = "Calibri",
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = subset(results, adj.P.Val > 0.05),
      pch = 21,
      color = "white",
      fill = "grey80"
    ) +
    ggplot2::geom_point(
      data = subset(results, logFC > 0 & adj.P.Val < 0.05),
      pch = 21,
      color = "white",
      fill = colors[[1]]
    ) +
    ggplot2::geom_point(
      data = subset(results, logFC < 0 & adj.P.Val < 0.05),
      pch = 21,
      color = "white",
      fill = colors[[2]]
    ) +
    ggplot2::scale_color_manual(values = c("black", "darkred")) +
    ggplot2::scale_y_continuous(trans = c("log10", "reverse")) +
    ggplot2::scale_x_continuous(
      breaks = seq(-4, 4, 2),
      labels = scales::math_format(2^.x),
      expand = ggplot2::expansion(mult = 1)
    ) +
    ggplot2::labs(
      x = xlab,
      y = "Adjusted p-value"
    ) +
    theme_plots()
}

plot_moi <- function(se, moi){

  idx <- which(SummarizedExperiment::rowData(se)$metabolite %in% moi)

  df <-
    se[idx, ] |>
    SummarizedExperiment::assay() |>
    tibble::as_tibble(rownames = "hmdb") |>
    tidyr::pivot_longer(
      -hmdb,
      names_to = "id",
      values_to = "value"
    ) |>
    dplyr::left_join(
      dplyr::select(tibble::as_tibble(SummarizedExperiment::colData(se[idx, ]), rownames = "id"), id:group),
      by = "id",
      copy = TRUE
    ) |>
    dplyr::left_join(
      dplyr::select(
        tibble::as_tibble(SummarizedExperiment::rowData(se[idx, ]), rownames = "hmdb"),
        hmdb,
        metabolite
      ),
      by = "hmdb",
      copy = TRUE
    ) |>
    dplyr::mutate(
      metabolite = factor(metabolite, levels = moi),
      oxygen = factor(oxygen, levels = c("N", "H"), labels = c("21%", "0.5%")),
      treatment = factor(treatment, levels = c("DMSO", "BAY"))
    )

  ggplot2::ggplot(df) +
    ggplot2::facet_wrap(~ metabolite, scales = "free_y", nrow = 1) +
    ggplot2::aes(
      x = treatment,
      y = 2 ^ value,
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
    ggplot2::scale_fill_manual(values = clrs, limits = force) +
    ggplot2::labs(
      x = "Treatment",
      y = "Peak area\n(normalized)",
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

read_metab_pathways <- function(filename){
  cols <- c("entrez_gene_ids", "metabolites")
  readr::read_tsv(filename, col_types = "cccc") %>%
    tidyr::unite("pathway", source, pathway, sep = " | ") %>%
    dplyr::mutate(
      dplyr::across(
        tidyselect::any_of(cols),
        ~stringr::str_split(.x, pattern = ",")
      )
    ) %>%
    dplyr::select(pathway, tidyselect::any_of(cols)) %>%
    tibble::deframe()
}

run_msea <- function(tt, pathways){
  stats <-
    tt %>%
    dplyr::filter(!is.na(KEGG)) %>%
    dplyr::mutate(KEGG = stringr::str_c("kegg:", KEGG)) %>%
    dplyr::select(KEGG, t) %>%
    tibble::deframe()

  fgsea::fgsea(
    pathways = pathways,
    stats = stats,
    minSize = 3,
    BPPARAM = BiocParallel::bpparam()
  ) %>%
    tibble::as_tibble() %>%
    dplyr::filter(pval < 0.05) %>%
    tidyr::separate(pathway, c("source", "pathway"), sep = " \\| ") %>%
    dplyr::arrange(desc(NES))
}
