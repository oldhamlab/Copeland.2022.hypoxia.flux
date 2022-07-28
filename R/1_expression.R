# 1_expression.R

read_data <- function(data_files){
  data_files[stringr::str_detect(data_files, "\\.csv$")] %>%
    rlang::set_names(stringr::str_extract(., "(lf|pasmc)_(02|05-bay|05|bay|)")) |>
    purrr::map_dfr(read_csv, .id = "experiment") |>
    dplyr::mutate(
      oxygen = factor(oxygen, levels = c("21%", "0.5%", "0.2%"), ordered = TRUE),
      treatment = factor(treatment, levels = c("None", "DMSO", "BAY"), ordered = TRUE)
    )
}

plot_blot <- function(blot_image){
  blot_image <- magick::image_read(blot_image)

  ggplot2::ggplot() +
    ggpubr::background_image(blot_image) +
    ggplot2::coord_fixed() +
    wmo::theme_wmo(
      base_family = "Calibri",
      base_size = 8
    ) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(5, 5, 5, 5),
      panel.border = ggplot2::element_blank(),
      plot.tag = ggplot2::element_text(face = "bold")
    )
}

normalize_densities <- function(blot_raw){
  blot_raw |>
    dplyr::filter(.data$time < 96) |>
    dplyr::filter(!(experiment == "lf_05-bay" & gel %in% c("b", "e", "f"))) |>
    tidyr::pivot_longer(
      blot:tidyselect::last_col(),
      names_to = "protein",
      values_to = "value",
      values_drop_na = TRUE
    ) |>
    dplyr::group_by(.data$experiment, .data$gel, .data$protein) |>
    dplyr::mutate(norm = value / mean(value, na.rm = TRUE)) |>
    dplyr::group_by(dplyr::across(.data$experiment:.data$time)) |>
    dplyr::mutate(density = norm / norm[protein == "blot"]) |>
    dplyr::filter(protein != "blot") |>
    dplyr::group_by(.data$experiment, .data$protein) |>
    dplyr::mutate(
      fold_change = density /
        mean(density[oxygen == min(oxygen) & treatment %in% c("None", "DMSO") & time == min(time)]),
      fold_change = replace(fold_change, experiment == "lf_bay" & protein == "hif1a", sqrt(fold_change)),
      group = dplyr::case_when(
        experiment %in% c("lf_02", "lf_05", "lf_bay", "pasmc_05") & treatment == "None" & oxygen == "21%" ~ "21%",
        experiment %in% c("lf_02", "lf_05", "lf_bay", "pasmc_05") & treatment == "DMSO" ~ "DMSO",
        experiment %in% c("lf_02", "lf_05", "lf_bay", "pasmc_05") & treatment == "BAY" ~ "BAY",
        experiment %in% c("lf_02", "lf_05", "lf_bay", "pasmc_05") & oxygen == "0.5%" ~ "0.5%",
        experiment %in% c("lf_02", "lf_05", "lf_bay", "pasmc_05") & oxygen == "0.2%" ~ "0.2%",
      ),
      group = factor(group, levels = c("21%", "0.5%", "0.2%", "DMSO", "BAY"))
    ) |>
    dplyr::group_by(experiment, oxygen, treatment, group, time, protein) |>
    # wmo::remove_nested_outliers(fold_change, remove = TRUE) |>
    dplyr::relocate(group, .after = treatment)
}

normalize_qpcr <- function(raw_mrna){
  raw_mrna |>
    dplyr::mutate(gene = tolower(gene)) |>
    dplyr::group_by(dplyr::across(c(experiment:gene))) |>
    dplyr::summarize(ct = mean(ct, na.rm = TRUE)) |>
    tidyr::pivot_wider(names_from = gene, values_from = ct) |>
    dplyr::ungroup() %>%
    dplyr::mutate(
      dplyr::across(-c(experiment:time, actin), ~ . - actin)
    ) |>
    dplyr::select(-actin) |>
    tidyr::pivot_longer(
      -c(experiment:time),
      names_to = "protein",
      values_to = "dct"
    ) |>
    dplyr::filter(!is.na(dct)) |>
    dplyr::group_by(.data$experiment, .data$protein) |>
    dplyr::mutate(
      ddct = dct - mean(
        dct[oxygen == min(oxygen) & treatment %in% c("None", "DMSO") & time == 0]
      ),
      fold_change = 2 ^ -ddct,
      group = dplyr::case_when(
        experiment %in% c("lf_02", "lf_05", "lf_bay", "pasmc_05") & treatment == "None" & oxygen == "21%" ~ "21%",
        experiment %in% c("lf_02", "lf_05", "lf_bay", "pasmc_05") & treatment == "DMSO" ~ "DMSO",
        experiment %in% c("lf_02", "lf_05", "lf_bay", "pasmc_05") & treatment == "BAY" ~ "BAY",
        experiment %in% c("lf_02", "lf_05", "lf_bay", "pasmc_05") & oxygen == "0.5%" ~ "0.5%",
        experiment %in% c("lf_02", "lf_05", "lf_bay", "pasmc_05") & oxygen == "0.2%" ~ "0.2%",
      ),
      group = factor(group, levels = c("21%", "0.5%", "0.2%", "DMSO", "BAY"))
    ) |>
    dplyr::group_by(experiment, oxygen, treatment, group, time, protein) |>
    wmo::remove_nested_outliers(fold_change, remove = TRUE) |>
    dplyr::relocate(group, .after = treatment)
}

analyze_hyp_bay_densities <- function(x, prot) {
  df <-
    x |>
    dplyr::filter(experiment == "lf_05-bay") |>
    dplyr::filter(protein == prot) |>
    dplyr::group_by(protein, oxygen, treatment) |>
    # wmo::remove_nested_outliers(fold_change, remove = TRUE) |>
    {\(x) x}()

  annot <-
    df |>
    dplyr::group_by(protein) |>
    tidyr::nest() |>
    dplyr::mutate(
      m = purrr::map(data, ~lmerTest::lmer(fold_change ~ oxygen * treatment + (1 | gel), data = .x)),
      res = purrr::map(m, ~emmeans::emmeans(
        .x,
        "pairwise" ~ oxygen * treatment,
        simple = "each",
        adjust = "mvt",
        combine = TRUE
      )[["contrasts"]]
      ),
      out = purrr::map(res, broom::tidy)
    ) |>
    tidyr::unnest(c(out)) |>
    dplyr::select(protein, oxygen, treatment, adj.p.value) |>
    dplyr::mutate(
      oxygen = replace(oxygen, oxygen == ".", "0.5%"),
      oxygen = factor(oxygen, levels = c("21%", "0.5%")),
      treatment = factor(treatment, levels = c("DMSO", "BAY")),
      y_pos = Inf,
      vjust = 1,
      lab = dplyr::case_when(
        adj.p.value < 0.05 ~ "*",
        TRUE ~ NA_character_
      )
    )

  list(data = df, annot = annot)
}
