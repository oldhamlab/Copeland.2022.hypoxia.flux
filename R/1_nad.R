# nad.R

clean_nad <- function(nad_data){
  nad_data |>
    dplyr::bind_rows(.id = "metabolite") |>
    dplyr::rename(rep = replicate) |>
    clean_technical_replicates() |>
    tidyr::separate(.data$experiment, c(NA, "date"), "_")
}

finalize_nad <- function(nad_interp, cells_per_dna){

  x <-
    dplyr::filter(cells_per_dna, cell_type == "lf" & volume == 200) |>
    dplyr::pull(slope)

  nad_interp |>
    dplyr::group_by(metabolite, date, oxygen, treatment, nucleotide) |>
    dplyr::summarise(conc = mean(conc)) |>
    dplyr::mutate(
      conc = dplyr::case_when(
        metabolite == "dna" ~ conc * x,
        metabolite == "nad" ~ conc * 720
      ),
      nucleotide = replace(nucleotide, is.na(nucleotide), "Count")
    ) |>
    tidyr::pivot_wider(-c(metabolite), names_from = nucleotide, values_from = conc) |>
    dplyr::mutate(
      oxygen = factor(oxygen, levels = c("21%", "0.5%")),
      treatment = factor(treatment, levels = c("none", "DMSO", "BAY")),
      treatment = forcats::fct_recode(treatment, "None" = "none"),
      Ratio = NADH/NAD,
      dplyr::across(tidyselect::contains("NAD"), ~. / Count * 1000)
    ) |>
    dplyr::select(-Count) |>
    tidyr::pivot_longer(c(NAD, NADH, Ratio), names_to = "measurement", values_to = "value") |>
    dplyr::arrange(measurement, oxygen, treatment)

}

annot_nad <- function(nad_final){
  nad_final |>
    dplyr::filter(treatment != "None") |>
    dplyr::group_by(measurement) |>
    tidyr::nest() |>
    dplyr::mutate(
      m = purrr::map(data, ~lmerTest::lmer(value ~ treatment * oxygen + (1 | date), data = .x)),
      s = purrr::map(
        m,
        ~emmeans::emmeans(.x, ~ treatment * oxygen) |>
          pairs(simple = "each", combine = TRUE) |>
          broom::tidy()
      )
    ) |>
    tidyr::unnest(c(s)) |>
    dplyr::select(measurement, oxygen, treatment, pval = adj.p.value) |>
    dplyr::mutate(
      lab = annot_p(pval),
      y = Inf,
      vjust = 1
    )
}
