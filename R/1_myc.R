# 1_myc.R

# 1_myc.R

combine_fluxes <- function(growth_rates, fluxes, exp) {
  growth_rates |>
    dplyr::filter(experiment == exp) |>
    dplyr::rename(flux = mu) |>
    dplyr::select(-X0) |>
    dplyr::mutate(metabolite = "growth") |>
    dplyr::bind_rows(dplyr::filter(fluxes, experiment == exp)) |>
    dplyr::filter(treatment %nin% c("siHIF1A", "siHIF2A")) |>
    dplyr::mutate(
      flux = ifelse(
        experiment == "bay-myc" & metabolite == "lactate" & date %in% c("2021-09-21", "2021-11-01"),
        flux / 3,
        flux
      )
    ) |>
    dplyr::group_by(oxygen, treatment, virus) |>
    wmo::remove_nested_outliers(flux, TRUE)
}

annot_myc_fluxes <- function(df, intervention) {
  if (intervention == "treatment") {
    lvls <- c("siCTL", "siMYC")
    x <- "oxygen"
  } else if (intervention == "virus") {
    lvls <- c("YFP", "MYC")
    x <- "virus"
  }
  fo1 <- as.formula(paste("flux ~ ", x, " * treatment + (1 | date)"))
  fo2 <- as.formula(paste("pairwise ~ ", x, "* treatment"))

  int <- rlang::sym(intervention)

  df |>
    dplyr::group_by(metabolite) |>
    tidyr::nest() |>
    dplyr::mutate(
      m = purrr::map(data, ~lmerTest::lmer(fo1, data = .x)),
      res = purrr::map(m, ~emmeans::emmeans(
        .x,
        fo2,
        simple = "each",
        adjust = "mvt",
        combine = TRUE
      )[["contrasts"]]
      ),
      out = purrr::map(res, broom::tidy)
    ) |>
    tidyr::unnest(c(out)) |>
    dplyr::filter(!!int != ".") |>
    dplyr::select(metabolite, {{ intervention }}, adj.p.value) |>
    dplyr::mutate(
      group = factor(!!int, levels = lvls),
      y_pos = Inf,
      vjust = 1.5,
      lab = annot_p(adj.p.value)
    ) |>
    dplyr::select(-!!int) |>
    dplyr::rename(!!int := group)
}

plot_myc <- function(df, annot, metab, ylab, x, fill) {
  df |>
    dplyr::filter(metabolite == metab) |>
    ggplot2::ggplot() +
    ggplot2::aes(
      x = {{x}},
      y = flux
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(
        fill = {{fill}}
      ),
      geom = "col",
      fun = "mean",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.6,
      show.legend = TRUE
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(
        group = {{fill}}
      ),
      geom = "errorbar",
      fun.data = "mean_se",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.2,
      size = 0.25,
      show.legend = FALSE,
      color = "black"
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(annot, metabolite == metab),
      ggplot2::aes(
        x = {{x}},
        y = y_pos,
        vjust = vjust,
        label = lab,
      ),
      family = "Calibri",
      color = "black",
      size = 6/ggplot2::.pt,
      show.legend = FALSE
    ) +
    ggplot2::labs(
      x = NULL,
      y = ylab,
      color = NULL,
      fill = NULL
    ) +
    ggplot2::scale_fill_manual(values = clrs, limits = force) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.05, 0.1)),
      breaks = scales::pretty_breaks(n = 6)
    ) +
    # ggplot2::coord_cartesian(ylim = c(-750, 1250)) +
    theme_plots() +
    ggplot2::theme(
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    )
}