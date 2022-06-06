# 1_viability.R

clean_viability <- function(viability_file) {
  readr::read_csv(viability_file) |>
    dplyr::mutate(
      viability = 100 * live / (dead + live),
      oxygen = factor(oxygen, levels = c("21%", "0.5%"))
    ) |>
    dplyr::group_by(time, oxygen) |>
    wmo::remove_nested_outliers(viability, remove = TRUE)
}
