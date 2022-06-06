# 1_mids.R

clean_mids <- function(mid_files){
  read_function <- function(filename){
    readr::read_csv(filename) %>%
      replace(is.na(.), 0)
  }

  mid_files %>%
    rlang::set_names(sub("\\..*$", "", basename(.))) |>
    purrr::map_dfr(read_function, .id = "experiment") |>
    tidyr::pivot_longer(
      cols = .data$`2HG M0`:tidyselect::last_col(),
      names_to = "ion",
      values_to = "peak_area",
      values_drop_na = TRUE
    ) |>
    dplyr::mutate(
      batch = stringr::str_extract(.data$experiment, "a|b|c"),
      method = stringr::str_extract(.data$experiment, "fs|sim"),
      cell_type = stringr::str_extract(.data$experiment, "lf|pasmc"),
      date = stringr::str_extract(.data$experiment, "\\d{4}-\\d{2}-\\d{2}"),
      .before = .data$experiment
    ) |>
    dplyr::select(-.data$experiment) |>
    tidyr::separate(
      .data$id,
      c("tracer", "treatment", "time", "well"),
      sep = " "
    ) |>
    tidyr::separate(.data$ion, into = c("metabolite", "isotope")) |>
    dplyr::mutate(
      peak_area = tidyr::replace_na(peak_area, 0),
      oxygen = dplyr::case_when(
        .data$treatment %in% c("21%", "dmso", "bay") ~ "21%",
        .data$treatment == "0.5%" ~ "0.5%",
        .data$treatment == "0.2%" ~ "0.2%"
      ),
      treatment = dplyr::case_when(
        .data$treatment %in% c("21%", "0.5%", "0.2%") ~ "None",
        TRUE ~ toupper(.data$treatment)
      ),
      oxygen = factor(.data$oxygen, levels = c("21%", "0.5%", "0.2%")),
      treatment = factor(.data$treatment, levels = c("None", "DMSO", "BAY")),
      time = as.double(.data$time),
      time = dplyr::case_when(
        .data$cell_type == "lf" ~ time * 24,
        .data$cell_type == "pasmc" ~ time * 12
      )
    ) |>
    dplyr::relocate(.data$oxygen, .before = .data$treatment) |>
    dplyr::select(-filename)
}

correct_mid <- function(mid_clean){
  correction_factors <-
    dplyr::mutate(correction_factors, method = "sim")

  mid_q <-
    mid_clean |>
    dplyr::left_join(
      correction_factors,
      by = c("batch", "method", "metabolite", "isotope" = "M")
    ) |>
    dplyr::mutate(
      area_corr = dplyr::case_when(
        .data$method == "sim" ~ peak_area * cf,
        .data$method == "fs" ~ peak_area * 1
      )
    ) |>
    dplyr::select(-c(.data$peak_area, .data$cf))

  correction_matrices <-
    isotope_library |>
    dplyr::mutate(
      matrix = purrr::map2(formula, polarity, mzrtools::mz_iso_quant),
      matrix = purrr::map(matrix, purrr::pluck, "prob_matrix")
    )

  mmult <- function(m, df){
    mid <- df$mid
    if (nrow(m) > length(mid)){
      m <- m[1:length(mid), 1:length(mid)]
    }
    mid_corr <- mzrtools::mz_iso_correct(m, mid)
    dplyr::bind_cols(df, mid_corr = mid_corr)
  }

  mid_na <-
    mid_q |>
    dplyr::group_by(
      .data$method,
      .data$cell_type,
      .data$date,
      .data$tracer,
      .data$oxygen,
      .data$treatment,
      .data$time,
      .data$well,
      .data$metabolite
    ) |>
    dplyr::filter(.data$metabolite %nin% c("palmitate", "sedoheptulose")) |>
    dplyr::mutate(mid = .data$area_corr / sum(.data$area_corr)) |>
    dplyr::filter(.data$mid != "NaN") |>
    dplyr::select(-.data$area_corr) |>
    tidyr::nest() |>
    dplyr::left_join(
      dplyr::select(correction_matrices, -c(.data$formula, .data$polarity)),
      by = "metabolite"
    ) |>
    dplyr::mutate(data = purrr::map2(.data$matrix, .data$data, mmult)) |>
    tidyr::unnest(c(.data$data)) |>
    dplyr::filter(.data$isotope %in% stringr::str_c("M", 0:6)) |>
    dplyr::select(
      .data$method:.data$metabolite,
      .data$isotope,
      .data$tracer:.data$well,
      .data$mid,
      .data$mid_corr
    )
}

remove_mid_outliers <- function(mid_correct){
  mid_correct |>
    dplyr::group_by(
      dplyr::across(c(.data$method:.data$time, .data$metabolite:.data$isotope))
    ) |>
    dplyr::mutate(
      outlier = replace_outliers(.data$mid_corr)
    ) |>
    dplyr::summarise(mid = mean(mid_corr, na.rm = TRUE)) |>
    dplyr::arrange(
      .data$cell_type,
      .data$metabolite,
      .data$oxygen,
      .data$treatment,
      .data$tracer
    )
}

plot_mid_curves <- function(mids){
  mids |>
    dplyr::group_by(
      .data$method,
      .data$cell_type,
      .data$tracer,
      .data$oxygen,
      .data$treatment,
      .data$metabolite
    ) |>
    tidyr::nest() |>
    dplyr::mutate(
      title = stringr::str_c(
        .data$metabolite,
        .data$method,
        .data$tracer,
        .data$cell_type,
        stringr::str_replace(.data$oxygen, "%", ""),
        .data$treatment,
        sep = "_"
      ),
      plots = purrr::map2(
        .data$data,
        .data$title,
        ~ggplot2::ggplot(.x) +
          ggplot2::aes(
            x = .data$time,
            y = .data$mid,
            color = .data$isotope
          ) +
          ggplot2::geom_point(
            alpha = 0.3
          ) +
          ggplot2::stat_summary(
            fun = "mean",
            geom = "point",
            alpha = 0.8
          ) +
          ggplot2::stat_summary(
            fun = "mean",
            geom = "line",
            alpha = 0.8,
            show.legend = FALSE
          ) +
          ggplot2::scale_x_continuous(
            name = "Time (h)",
            breaks = seq(24, 96, by = 24)
          ) +
          ggplot2::scale_color_manual(values = viridis::viridis(7)) +
          ggplot2::labs(
            y = "Mole fraction",
            title = .y
          )
      )
    )
}
