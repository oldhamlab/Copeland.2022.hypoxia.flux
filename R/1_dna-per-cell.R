# 1_dna-per-cell.R

clean_dna_per_cell <- function(filename) {
  filename |>
    read_multi_excel() |>
    purrr::map(clean_technical_replicates) |>
    dplyr::bind_rows(.id = "id") |>
    tidyr::separate(id, c("cell_type", "volume"), sep = "_", convert = TRUE) |>
    dplyr::mutate(cells = 1000 * cells)
}

calculate_cells_per_dna <- function(tbl) {
  tbl %>%
    dplyr::filter(!(cell_type == "lf" & volume == "100" & cells == 300000)) %>%
    dplyr::filter(!(cell_type == "pasmc" & volume == "200" & cells == 400000)) %>%
    dplyr::group_by(cell_type, volume) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      model = map(data, ~lm(conc ~ 0 + cells, data = .x, na.action = modelr::na.warn)),
      glance = map(model, broom::tidy)
    ) %>%
    tidyr::unnest(c(glance)) %>%
    dplyr::select(cell_type, volume, slope = estimate) %>%
    dplyr::mutate(slope = 1/slope)
}
