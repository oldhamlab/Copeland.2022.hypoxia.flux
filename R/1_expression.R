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
    wmo::remove_nested_outliers(fold_change, remove = TRUE) |>
    dplyr::relocate(group, .after = treatment)
}
