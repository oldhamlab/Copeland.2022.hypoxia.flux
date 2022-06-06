# 1_qbias.R

isotope_library <-
  tibble::tribble(
    ~ metabolite, ~ formula, ~ polarity,
    "2HG", "C5H8O5", "negative",
    "2OG", "C5H6O5", "negative",
    "alanine", "C3H7NO2", "negative",
    "aspartate", "C4H7NO4", "negative",
    "citrate", "C6H8O7", "negative",
    "glutamate", "C5H9NO4", "negative",
    "glutamine", "C5H10N2O3", "negative",
    "lactate", "C3H6O3", "negative",
    "malate", "C4H6O5", "negative",
    "pyruvate", "C3H4O3", "negative",
    "serine", "C3H7NO3", "negative",
    "succinate", "C4H6O4", "negative",
    "3PG", "C3H7O7P", "negative",
    "aconitate", "C6H6O6", "negative",
    "FBP", "C6H14O12P2", "negative",
    "G3P", "C3H9O6P", "negative",
    "palmitate", "C16H32O2", "negative",
    "PEP", "C3H5O6P", "negative",
    "sedoheptulose", "C7H14O7", "negative",
    "DHAP", "C3H7O6P", "negative",
    "GAP", "C3H7O6P", "negative",
    "G1P", "C6H13O9P", "negative",
    "G6P", "C6H13O9P", "negative",
    "R5P", "C5H11O8P", "negative"
  )

calculate_ratios <- function(file_name){
  readr::read_csv(file_name) |>
    dplyr::filter(!is.na(area)) |>
    tidyr::separate(
      .data$filename,
      into = c(NA, "window", "replicate"),
      sep = c(1, 2),
      convert = TRUE
    ) |>
    tidyr::separate(
      .data$ion,
      into = c("metabolite", "isotope"),
      sep = " ",
      convert = TRUE
    ) |>
    dplyr::mutate(carbons = case_when(
      .data$metabolite %in% c("citrate") ~ 6,
      .data$metabolite %in% c("2HG", "2OG", "glutamate", "glutamine") ~ 5,
      .data$metabolite %in% c("aspartate", "malate", "succinate") ~ 4,
      .data$metabolite %in% c("lactate", "pyruvate", "alanine", "serine") ~ 3
    )) |>
    dplyr::filter(.data$window <= .data$carbons) |>
    tidyr::pivot_wider(names_from = .data$isotope, values_from = .data$area) |>
    dplyr::mutate(ratio = .data$M1 / .data$M0) |>
    dplyr::group_by(.data$metabolite)
}

import_qbias <- function(file_list){
  file_list[stringr::str_detect(file_list, "\\.csv$")] %>%
    rlang::set_names(
      stringr::str_extract(basename(.), pattern = "(?<=_)\\w(?=\\.csv)")
    ) |>
    purrr::map_dfr(calculate_ratios, .id = "batch") |>
    dplyr::group_by(batch, metabolite) |>
    dplyr::arrange(metabolite)
}

calculate_predicted_ratios <- function(){
  predicted_ratios <-
    isotope_library |>
    dplyr::mutate(table = purrr::map2(
      formula,
      polarity,
      ~mzrtools::mz_iso_quant(molecule = .x, polarity = .y)[["prob_matrix"]])
    ) |>
    dplyr::mutate(pred_ratio = purrr::map_dbl(table, ~ .x[[2, 1]]/.x[[1, 1]])) |>
    dplyr::select(.data$metabolite, .data$pred_ratio)
}

calculate_correction_factors <- function(qbias_ratios, pred_ratios){
  qbias_ratios |>
    tidyr::nest() |>
    dplyr::mutate(
      model = purrr::map(.data$data, ~MASS::rlm(ratio ~ poly(window, 3), data = .x, maxit = 1000)),
      predict = purrr::map2(.data$model, .data$data, predict)
    ) |>
    tidyr::unnest(c(data, predict)) |>
    dplyr::select(.data$batch, .data$metabolite, .data$window, .data$carbons, .data$predict) |>
    dplyr::distinct() |>
    dplyr::left_join(pred_ratios, by = "metabolite") |>
    dplyr::mutate(cf = .data$predict / .data$pred_ratio) |>
    dplyr::select(.data$metabolite, "M" = .data$window, .data$carbons, .data$cf) |>
    dplyr::filter(.data$M < .data$carbons) |>
    dplyr::ungroup() |>
    dplyr::mutate(M = M + 1) |>
    tidyr::pivot_wider(names_from = .data$M, values_from = .data$cf) |>
    dplyr:: mutate(
      M0 = 1,
      M1 = 1 / .data$`1` * .data$M0,
      M2 = 1 / .data$`2` * .data$M1,
      M3 = 1 / .data$`3` * .data$M2,
      M4 = 1 / .data$`4` * .data$M3,
      M5 = 1 / .data$`5` * .data$M4,
      M6 = 1 / .data$`6` * .data$M5
    ) |>
    dplyr::select(.data$batch, .data$metabolite, tidyselect::matches("M[0-9]+")) |>
    tidyr::pivot_longer(
      cols = matches("M[0-9]+"),
      names_to = "M",
      values_to = "cf",
      values_drop_na = TRUE
    ) |>
    dplyr::arrange(.data$batch, .data$metabolite)
}
