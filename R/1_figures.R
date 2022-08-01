# 1_figures.R


# setup -------------------------------------------------------------------

clrs <- c(RColorBrewer::brewer.pal(4, "Set1")[1:4], "#08306b", RColorBrewer::brewer.pal(9, "Set1")[9:8])
names(clrs) <- c("21%", "0.5%", "DMSO", "BAY", "0.2%", "siCTL", "siMYC")

theme_plots <- function() {
  list(
    wmo::theme_wmo(
      base_family = "Calibri",
      base_size = 8
    ),
    ggplot2::theme(
      # panel.border = ggplot2::element_rect(size = 0.1),
      axis.line = ggplot2::element_line(colour = "black", size = 0.25, lineend = "square"),
      panel.border = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5, 5, 5, 5),
      plot.tag = ggplot2::element_text(face = "bold"),
      axis.title.y.left = ggplot2::element_text(margin = ggplot2::margin(r = 3))
    ),
    ggplot2::coord_cartesian(clip = "off")
  )
}

theme_patchwork <- function(design = NULL, widths = NULL, heights = NULL, tags = "A", ...) {
  list(
    patchwork::plot_layout(
      design = design,
      widths = widths,
      heights = heights,
      ...
    ),
    patchwork::plot_annotation(
      tag_levels = tags,
      theme = ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
    )
  )
}

write_figures <- function(plot, filename, path = "manuscript/figures") {
  # path <- "manuscript/figures"

  gtab <- patchwork::patchworkGrob(plot)

  overall_width <-
    grid::convertWidth(
      sum(gtab$widths),
      unitTo = "in",
      valueOnly = TRUE
    )

  overall_height <-
    grid::convertHeight(
      sum(gtab$heights),
      unitTo = "in",
      valueOnly = TRUE
    )

  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    device = ragg::agg_png,
    path = path,
    width = overall_width,
    height = overall_height,
    units = "in",
    res = 300
  )

  if (file.exists("Rplots.pdf")) unlink("Rplots.pdf")

  stringr::str_c(path, "/", filename)
}

plot_time_lines <- function(
    df,
    y,
    ylab = prot,
    clr = c("oxygen", "treatment", "group")
) {
  ggplot2::ggplot(df) +
    ggplot2::aes(
      x = time,
      y = {{y}},
      color = .data[[clr]],
      fill = .data[[clr]]
    ) +
    # ggplot2::geom_line(
    #   ggplot2::aes(group = interaction(date, group)),
    #   size = 0.25,
    #   alpha = 0.25,
    #   show.legend = FALSE
    # ) +
    # ggplot2::stat_summary(
    #   geom = "linerange",
    #   fun.data = ggplot2::mean_se,
    #   size = 0.5,
    #   show.legend = FALSE
  # ) +
  ggplot2::stat_summary(
    geom = "errorbar",
    fun.data = ggplot2::mean_se,
    color = "black",
    width = 2,
    size = 0.25,
    show.legend = FALSE
  ) +
    ggplot2::stat_summary(
      geom = "line",
      fun.data = ggplot2::mean_se,
      size = 0.5,
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      geom = "point",
      fun = "mean",
      pch = 21,
      color = "white",
      size = 1.5,
      stroke = 0.2,
      show.legend = FALSE
    ) +
    ggplot2::labs(
      x = "Time (h)",
      y = ylab
    ) +
    ggplot2::scale_x_continuous(breaks = seq(0, 72, 24)) +
    ggplot2::scale_y_continuous(limits = c(0, NA)) +
    ggplot2::scale_color_manual(values = clrs) +
    ggplot2::scale_fill_manual(values = clrs) +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    theme_plots()
}

plot_growth_curve <- function(
    df,
    cell = c("lf", "df", "pasmc"),
    exper = c("02", "05", "bay")
) {
  df |>
    dplyr::filter(
      cell_type %in% cell &
        experiment %in% exper &
        metabolite == "cells" &
        time < 96
    ) |>
    dplyr::group_by(
      date,
      group,
      time
    ) |>
    dplyr::summarize(count = mean(conc, na.rm = TRUE)) |>
    plot_time_lines(
      y = count,
      ylab = "Cell count",
      clr = "group"
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_number(scale_cut = scales::cut_short_scale())
    ) +
    ggplot2::coord_cartesian(
      ylim = c(0, NA),
      clip = "off"
    )
}

plot_growth_rates <- function(
    df,
    cell = c("lf", "pasmc"),
    exper = c("02", "05", "bay")
) {
  x <-
    df |>
    dplyr::filter(cell_type %in% cell & experiment %in% exper)

  annot <-
    lmerTest::lmer(mu ~ group + (1 | date), data = x) |>
    emmeans::emmeans(~ group) |>
    pairs(adjust = "mvt") |>
    broom::tidy() |>
    dplyr::mutate(
      group = stringr::str_extract(contrast, "(?<= - ).*"),
      x = 1.5,
      y = Inf,
      vjust = 1,
      label = annot_p(p.value)
    )

  ggplot2::ggplot(x) +
    ggplot2::aes(
      x = group,
      y = mu
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = group),
      geom = "col",
      fun = "mean",
      # width = 0.6,
      show.legend = FALSE,
      alpha = 0.5
    ) +
    ggbeeswarm::geom_beeswarm(
      ggplot2::aes(fill = group),
      pch = 21,
      size = 1,
      stroke = 0.2,
      cex = 4,
      color = "white",
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      geom = "errorbar",
      fun.data = ggplot2::mean_se,
      width = 0.2,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = annot,
      ggplot2::aes(
        x = x,
        label = label,
        y = y,
        vjust = vjust
      ),
      family = "Calibri",
      size = 8/ggplot2::.pt
    ) +
    ggplot2::labs(
      x = "Treatment",
      y = "Growth rate (/h)"
    ) +
    ggplot2::scale_fill_manual(
      values = clrs,
      limits = force
    ) +
    theme_plots() +
    ggplot2::coord_cartesian(
      # ylim = c(0, NA),
      clip = "off"
    ) +
    NULL
}

plot_expression <- function(
    df,
    exp,
    prot = c("ldha", "hif1a"),
    ylab = prot
) {

  df |>
    dplyr::filter(experiment %in% exp) |>
    dplyr::filter(protein == prot) |>
    plot_time_lines(y = fold_change, ylab = ylab, clr = "group")
}

plot_cells_per_dna <- function(dna_per_cell_clean, cell = c("lf", "pasmc")) {
  dna_per_cell_clean |>
    dplyr::filter(cell_type == cell) |>
    dplyr::filter(.data$volume == 200 & cells < 400000) |>
    ggplot2::ggplot() +
    ggplot2::facet_wrap(~cell_type, labeller = ggplot2::as_labeller(toupper)) +
    ggplot2::aes(
      x = cells,
      y = conc
    ) +
    ggplot2::geom_smooth(
      formula = y ~ 0 + x,
      method = "lm",
      color = clrs[[2]],
      size = 0.5,
      se = FALSE
    ) +
    # ggplot2::stat_summary(
    #   geom = "linerange",
    #   fun.data = "mean_se",
    #   size = 0.5,
    #   show.legend = FALSE
    # ) +
    ggplot2::stat_summary(
      geom = "errorbar",
      fun.data = ggplot2::mean_se,
      color = "black",
      width = 7500,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      geom = "point",
      fun = "mean",
      pch = 21,
      color = "white",
      fill = "black",
      size = 1.5,
      show.legend = FALSE,
      stroke = 0.2
    ) +
    ggplot2::labs(
      x = "Cell count",
      y = "DNA (ng)"
    ) +
    ggplot2::scale_x_continuous(
      labels = scales::label_number(scale_cut = scales::cut_short_scale())
    ) +
    theme_plots() +
    ggplot2::coord_cartesian(xlim = c(0, NA), clip = "off") +
    NULL
}

plot_dna_count_hypoxia <- function(dna_count_hypoxia) {
  ggplot2::ggplot(dna_count_hypoxia) +
    ggplot2::aes(
      x = count,
      y = conc,
      color = oxygen,
      fill = oxygen
    ) +
    ggplot2::geom_smooth(
      method = "lm",
      size = 0.5,
      formula = y ~ x,
      se = FALSE,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      pch = 21,
      color = "white",
      # alpha = 0.3,
      size = 1.5,
      stroke = 0.2,
      show.legend = FALSE
    ) +
    ggplot2::labs(
      x = "Cell count",
      y = "DNA (ng)"
    ) +
    ggplot2::scale_x_continuous(
      labels = scales::label_number(scale_cut = scales::cut_short_scale())
    ) +
    ggplot2::scale_color_manual(values = clrs) +
    ggplot2::scale_fill_manual(values = clrs) +
    theme_plots() +
    ggplot2::coord_cartesian(
      ylim = c(0, NA),
      xlim = c(0, 300000),
      clip = "off"
    )
}

plot_evap_data <- function(evap_clean) {
  evap_clean |>
    dplyr::filter(experiment == "05" & cell_type == "lf") |>
    dplyr::mutate(
      oxygen = factor(oxygen, levels = c("21%", "0.5%"))
    ) |>
    ggplot2::ggplot() +
    ggplot2::aes(
      x = time,
      y = volume,
      color = oxygen,
      fill = oxygen
    ) +
    ggplot2::geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      size = 0.5,
      show.legend = FALSE
    ) +
    # ggplot2::stat_summary(
    #   geom = "linerange",
    #   fun.data = "mean_se",
    #   size = 0.5,
    #   show.legend = FALSE
    # ) +
    ggplot2::stat_summary(
      geom = "errorbar",
      fun.data = ggplot2::mean_se,
      color = "black",
      width = 2,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      geom = "point",
      fun = "mean",
      pch = 21,
      color = "white",
      size = 1.5,
      stroke = 0.2,
      show.legend = FALSE
    ) +
    ggplot2::labs(
      x = "Time (h)",
      y = "Volume (mL)"
    ) +
    ggplot2::scale_x_continuous(breaks = seq(0, 96, 24)) +
    ggplot2::scale_color_manual(values = clrs) +
    ggplot2::scale_fill_manual(values = clrs) +
    theme_plots()
}

plot_k <- function(degradation_rates, k) {
  annot <-
    k |>
    dplyr::select(-k) |>
    dplyr::mutate(
      label = "*",
      ypos = Inf,
      vjust = 1
    )

  degradation_rates |>
    dplyr::mutate(
      group = dplyr::case_when(
        oxygen == "21%" & treatment == "None" ~ "21%",
        oxygen == "0.5%" ~ "0.5%",
        treatment == "DMSO" ~ "DMSO"
      ),
      group = factor(group, levels = c("21%", "0.5%", "DMSO"))
    ) |>
    dplyr::left_join(annot, by = c("metabolite", "oxygen", "treatment")) |>
    ggplot2::ggplot() +
    ggplot2::aes(
      x = reorder(toupper(abbreviation), k),
      y = k
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = group),
      geom = "col",
      fun = "mean",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.6,
      show.legend = TRUE
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "black",
      lwd = 0.25
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = group),
      geom = "errorbar",
      fun.data = ggplot2::mean_se,
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.2,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        color = group,
        label = label,
        y = ypos,
        vjust = vjust
      ),
      family = "Calibri",
      size = 6/ggplot2::.pt,
      position = ggplot2::position_dodge(width = 0.6),
      show.legend = FALSE
    ) +
    ggplot2::labs(
      x = "Metabolite",
      y = "Rate constant (/h)",
      fill = NULL,
      color = NULL
    ) +
    ggplot2::scale_color_manual(
      values = clrs,
      limits = force
    ) +
    ggplot2::scale_fill_manual(
      values = clrs,
      limits = force
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(alpha = 1))
    ) +
    # ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.2)) +
    theme_plots() +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(color = "gray90"),
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

plot_high_fluxes <- function(
    df,
    cell = c("lf", "pasmc"),
    exper = c("02", "05", "bay")
) {

  x <-
    df |>
    dplyr::filter(
      cell_type %in% cell &
        experiment %in% exper &
        metabolite %in% c("lactate", "glucose")
    )

  annot <-
    x |>
    dplyr::group_by(abbreviation) |>
    tidyr::nest() |>
    dplyr::mutate(
      model = purrr::map(
        data,
        ~lmerTest::lmer(flux ~ group + (1 | date), data = .x) |>
          emmeans::emmeans(~ group) |>
          pairs() |>
          broom::tidy()
      )
    ) |>
    tidyr::unnest(c(model)) |>
    dplyr::mutate(
      group = stringr::str_extract(contrast, "(?<= - ).*"),
      y = Inf,
      vjust = 1,
      label = annot_p(p.value)
    )

  ggplot2::ggplot(x) +
    ggplot2::aes(
      x = reorder(toupper(abbreviation), flux),
      y = flux
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = group),
      geom = "col",
      fun = "mean",
      # width = 0.6,
      position = ggplot2::position_dodge2(),
      show.legend = TRUE,
      alpha = 0.5
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "black",
      lwd = 0.25
    ) +
    ggbeeswarm::geom_beeswarm(
      ggplot2::aes(fill = group),
      dodge.width = 0.9,
      pch = 21,
      size = 1,
      stroke = 0.2,
      cex = 4,
      color = "white",
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = group),
      geom = "errorbar",
      fun.data = ggplot2::mean_se,
      position = ggplot2::position_dodge(width = 0.9),
      width = 0.2,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = annot,
      ggplot2::aes(
        x = abbreviation,
        label = label,
        y = y,
        vjust = vjust
      ),
      family = "Calibri",
      size = 8/ggplot2::.pt,
      inherit.aes = FALSE
    ) +
    ggplot2::labs(
      x = "Metabolite",
      y = "Flux (fmol/cell/h)",
      fill = NULL
    ) +
    ggplot2::scale_fill_manual(values = clrs, limits = force) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = 0.2),
      breaks = scales::extended_breaks(n = 7)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(alpha = 1))
    ) +
    theme_plots() +
    ggplot2::theme(
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

plot_low_fluxes <- function(
    df,
    cell = c("lf", "pasmc"),
    exper = c("02", "05", "bay")
) {

  x <-
    df |>
    dplyr::filter(
      cell_type %in% cell &
        experiment %in% exper &
        metabolite %nin% c("lactate", "glucose")
    ) |>
    dplyr::filter(
      !(metabolite == "glutamine" & flux > 0)
    )

  annot <-
    x |>
    dplyr::group_by(abbreviation) |>
    tidyr::nest() |>
    dplyr::mutate(
      model = purrr::map(
        data,
        ~lmerTest::lmer(flux ~ group + (1 | date), data = .x) |>
          emmeans::emmeans(~ group) |>
          pairs() |>
          broom::tidy()
      )
    ) |>
    tidyr::unnest(c(model)) |>
    dplyr::mutate(
      group = stringr::str_extract(contrast, "(?<= - ).*"),
      y = Inf,
      vjust = 1,
      label = annot_p(p.value)
    )

  ggplot2::ggplot(x) +
    ggplot2::aes(
      x = reorder(toupper(abbreviation), flux),
      y = flux
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = group),
      geom = "col",
      fun = "mean",
      position = ggplot2::position_dodge2(),
      # width = 0.9,
      alpha = 0.5,
      show.legend = TRUE
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "black",
      lwd = 0.25
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = group),
      geom = "errorbar",
      fun.data = ggplot2::mean_se,
      position = ggplot2::position_dodge(width = 0.9),
      width = 0.2,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = annot,
      ggplot2::aes(
        x = abbreviation,
        y = y,
        vjust = vjust,
        label = label,
        group = group
      ),
      family = "Calibri",
      size = 8/ggplot2::.pt
    ) +
    ggplot2::labs(
      x = "Metabolite",
      y = "Flux (fmol/cell/h)",
      fill = NULL
    ) +
    ggplot2::scale_fill_manual(values = clrs, limits = force) +
    ggplot2::scale_y_continuous(
      trans = ggallin::pseudolog10_trans,
      breaks = c(-100, -10, 0, 10),
      limits = c(-250, 50)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(alpha = 1))
    ) +
    theme_plots() +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(color = "gray80", size = 0.1),
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

annot_mids_main <- function(a, formula) {
  lmerTest::lmer(
    mid ~ group * isotope + (1 | date),
    data = a
  ) |>
    emmeans::emmeans(~ group * isotope) |>
    emmeans::mvcontrast("pairwise", mult.name = "isotope") |>
    tibble::as_tibble() |>
    dplyr::select(contrast, tidyselect::contains("p.value")) |>
    dplyr::rename_with(~ "pval", .cols = tidyselect::contains("p.value")) |>
    dplyr::mutate(
      label = dplyr::case_when(
        pval < 0.05 & contrast == "21% - 0.5%" ~ "*",
        pval < 0.05 & contrast == "21% - 0.2%" ~ "*",
        pval < 0.05 & contrast == "DMSO - BAY" ~ "†",
        pval < 0.05 & contrast == "0.5% - BAY" ~ "‡"
      ),
      order = dplyr::case_when(
        contrast == "21% - 0.5%" ~ 1,
        contrast == "21% - 0.2%" ~ 1,
        contrast == "DMSO - BAY" ~ 2,
        contrast == "0.5% - BAY" ~ 3
      )
    ) |>
    dplyr::filter(!is.na(label)) |>
    dplyr::arrange(order) |>
    dplyr::pull(label) |>
    stringr::str_c(collapse = " ")
}

plot_mids <- function(df, cell, metab, t = 72, track) {
  tracer_labels <- c(
    glc2 = expression(paste("[1,2-"^13, "C"[2], "]-GLC")),
    glc6 = bquote("[U-"^13 * "C"[6] * "]-GLC →" ~ .(metab)),
    q5 = bquote("[U-"^13 * "C"[5] * "]-GLN →" ~ .(metab)),
    lac3 = expression(paste("[U-"^13, "C"[3], "]-LAC"))
  )
  tracer_levels <- c("glc2", "glc6", "q5", "lac3")

  x <-
    df |>
    dplyr::filter(
      cell_type == cell &
        metabolite == metab &
        time == t &
        tracer == track
    ) |>
    dplyr::filter(
      !(metabolite == "CIT" & isotope == "M6")
    ) |>
    dplyr::mutate(
      metabolite = factor(metabolite, levels = metab),
      tracer = factor(tracer, levels = tracer_levels, labels = tracer_labels),
      group = dplyr::case_when(
        oxygen == "21%" & treatment == "None" ~ "21%",
        oxygen == "0.5%" & treatment == "None" ~ "0.5%",
        oxygen == "21%" & treatment == "DMSO" ~ "DMSO",
        oxygen == "21%" & treatment == "BAY" ~ "BAY"
      ),
      group = factor(group, levels = c("21%", "0.5%", "DMSO", "BAY"))
    )

  annot <-
    x |>
    dplyr::group_by(tracer, metabolite) |>
    tidyr::nest() |>
    dplyr::mutate(annot = purrr::map(data, annot_mids_main)) |>
    tidyr::unnest(c(annot))

  ggplot2::ggplot(x) +
    ggplot2::aes(
      x = isotope,
      y = mid
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = group),
      geom = "col",
      fun = "mean",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.6,
      show.legend = TRUE
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "black",
      lwd = 0.1
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = group),
      geom = "errorbar",
      fun.data = "mean_se",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.2,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = annot,
      ggplot2::aes(
        x = -Inf,
        y = Inf,
        vjust = 1.5,
        hjust = -0.3,
        label = annot
      ),
      family = "Calibri",
      size = 8/ggplot2::.pt
    ) +
    ggplot2::labs(
      x = "Isotope",
      y = "Mole fraction",
      fill = NULL,
      title = tracer_labels[track]
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, 1, 0.25),
      limits = c(0, 1.1)
    ) +
    ggplot2::theme(
      strip.placement = "outside",
      legend.title = ggplot2::element_blank()
    ) +
    ggplot2::scale_fill_manual(values = clrs, limits = force) +
    theme_plots() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "grey80", size = 0.1),
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

plot_all_mids <- function(df, cell, t = 72, track) {
  tracer_labels <- c(
    expression(paste("[1,2-"^13, "C"[2], "] glucose")),
    expression(paste("[U-"^13, "C"[6], "] glucose")),
    expression(paste("[U-"^13, "C"[5], "] glutamine")),
    expression(paste("[U-"^13, "C"[3], "] lactate"))
  )
  tracer_levels <- c("glc2", "glc6", "q5", "lac3")
  metab <- c(
    "FBP",
    "3PG",
    "PYR",
    # "ALA",
    # "SER",
    "LAC",
    "CIT",
    "AKG",
    # "GLN",
    "GLU",
    "MAL"
    # "ASP"
  )

  x <-
    df |>
    dplyr::filter(
      cell_type == cell &
        time == t &
        metabolite %in% metab
    ) |>
    dplyr::filter(
      !(metabolite == "CIT" & isotope == "M6")
    ) |>
    dplyr::mutate(
      metabolite = factor(metabolite, levels = metab),
      metabolite = forcats::fct_recode(metabolite, "`3PG`" = "3PG"),
      tracer = factor(tracer, levels = tracer_levels, labels = tracer_labels),
      isotope = stringr::str_replace(isotope, "M", ""),
      group = dplyr::case_when(
        oxygen == "21%" & treatment == "None" ~ "21%",
        oxygen == "0.5%" & treatment == "None" ~ "0.5%",
        oxygen == "21%" & treatment == "DMSO" ~ "DMSO",
        oxygen == "21%" & treatment == "BAY" ~ "BAY"
      ),
      group = factor(group, levels = c("21%", "0.5%", "DMSO", "BAY"))
    )

  annot <-
    x |>
    dplyr::group_by(tracer, metabolite) |>
    tidyr::nest() |>
    dplyr::mutate(annot = purrr::map(data, annot_mids_main)) |>
    tidyr::unnest(c(annot))

  ggplot2::ggplot(x) +
    ggplot2::aes(
      x = isotope,
      y = mid
    ) +
    ggplot2::facet_grid(
      tracer ~ metabolite,
      labeller = ggplot2::label_parsed,
      # switch = "y",
      scales = "free_x",
      space = "free_x"
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = group),
      geom = "col",
      fun = "mean",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.6,
      show.legend = TRUE
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "black",
      lwd = 0.1
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = group),
      geom = "errorbar",
      fun.data = "mean_se",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.2,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = annot,
      ggplot2::aes(
        x = -Inf,
        y = Inf,
        vjust = 1.5,
        hjust = -0.3,
        label = annot
      ),
      family = "Calibri",
      size = 8/ggplot2::.pt
    ) +
    ggplot2::labs(
      x = "Isotope",
      y = "Mole fraction",
      fill = NULL
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, 1, 0.25),
      limits = c(0, 1.1)
    ) +
    ggplot2::theme(
      strip.placement = "outside",
      legend.title = ggplot2::element_blank()
    ) +
    ggplot2::scale_fill_manual(values = clrs, limits = force) +
    theme_plots() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "grey80", size = 0.1),
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10),
      axis.line = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, size = 0.25)
    ) +
    theme_patchwork(
      widths = unit(9.5, "in"),
      heights = unit(7.5, "in"),
      tags = NULL
    )
}

plot_m5_citrate <- function(df) {
  x <-
    df |>
    dplyr::filter(
      metabolite == "CIT" &
        tracer == "q5" &
        treatment == "None" &
        isotope == "M5" &
        method == "sim"
    ) |>
    dplyr::filter(
      (cell_type == "lf" & time == 48) |
        (cell_type == "pasmc" & time == 36)
    )

  annot <-
    x |>
    dplyr::group_by(cell_type) |>
    tidyr::nest() |>
    dplyr::mutate(
      model = purrr::map(
        data,
        ~lmerTest::lmer(mid ~ oxygen + (1 | date), data = .x) |>
          emmeans::emmeans(~ oxygen) |>
          pairs() |>
          broom::tidy()
      )
    ) |>
    tidyr::unnest(c(model)) |>
    dplyr::mutate(
      x = 1.5,
      y = Inf,
      vjust = 1.5,
      label = annot_p(p.value)
    )

  ggplot2::ggplot(x) +
    ggplot2::aes(
      x = oxygen,
      y = mid
    ) +
    ggplot2::facet_wrap(
      ~cell_type,
      labeller = ggplot2::as_labeller(toupper)
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = oxygen),
      geom = "col",
      fun = "mean",
      # width = 0.6,
      show.legend = FALSE,
      alpha = 0.5
    ) +
    ggbeeswarm::geom_beeswarm(
      ggplot2::aes(fill = oxygen),
      pch = 21,
      size = 1,
      stroke = 0.2,
      cex = 4,
      color = "white",
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      geom = "errorbar",
      fun.data = ggplot2::mean_se,
      width = 0.2,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = annot,
      ggplot2::aes(
        x = x,
        label = label,
        y = y,
        vjust = vjust
      ),
      family = "Calibri",
      size = 8/ggplot2::.pt
    ) +
    ggplot2::labs(
      x = "Oxygen",
      y = "M5 Citrate fraction"
    ) +
    ggplot2::scale_fill_manual(
      values = clrs,
      limits = force
    ) +
    theme_plots() +
    ggplot2::coord_cartesian(
      # ylim = c(0, NA),
      clip = "off"
    ) +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, size = 0.25)
    ) +
    NULL
}

format_time_course_mids <- function(model_mids) {
  tracer_labels <- c(
    expression(paste("[1,2-"^13, "C"[2], "] glucose")),
    expression(paste("[U-"^13, "C"[6], "] glucose")),
    expression(paste("[U-"^13, "C"[5], "] glutamine")),
    expression(paste("[U-"^13, "C"[3], "] lactate"))
  )
  tracer_levels <- c("glc2", "glc6", "q5", "lac3")

  model_mids |>
    tidyr::unnest(c(data)) |>
    dplyr::filter(metabolite %in% c("FBP", "PYR", "CIT", "MAL")) |>
    dplyr::filter(tracer != "lac3") |>
    # dplyr::filter(isotope != "M6") |>
    dplyr::mutate(
      tracer = factor(
        tracer,
        levels = tracer_levels,
        labels = tracer_labels
      ),
      metabolite = factor(metabolite, levels = c("FBP", "PYR", "CIT", "MAL"))
    )
}

plot_mid_time_course <- function(df, cells, o2, treat, color) {
  df |>
    dplyr::filter(
      cell_type == cells &
        oxygen %in% o2 &
        treatment == treat &
        isotope == "M0"
    ) |>
    ggplot2::ggplot() +
    ggplot2::aes(
      x = time,
      y = 1 - mean,
      color = oxygen,
      fill = oxygen
    ) +
    ggplot2::facet_grid(
      metabolite ~ tracer,
      labeller = ggplot2::label_parsed
    ) +
    ggplot2::geom_linerange(
      ggplot2::aes(
        ymin = 1 - mean - se,
        ymax = 1 - mean + se
      ),
      color = "black",
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_line(
      size = 0.5,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      pch = 21,
      color = "white",
      size = 1.5,
      stroke = 0.2,
      show.legend = TRUE
    ) +
    ggplot2::labs(
      x = "Time (h)",
      y = "Labeled fraction",
      color = NULL,
      fill = NULL
    ) +
    ggplot2::scale_x_continuous(breaks = seq(0, 72, 24)) +
    ggplot2::scale_color_manual(
      values = clrs,
      limits = force,
      aesthetics = c("fill", "color")
    ) +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    theme_plots() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "grey80", size = 0.1),
      legend.key.size = ggplot2::unit(0.5, units = "lines"),
      axis.line = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, size = 0.25),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    ) +
    theme_patchwork(
      widths = unit(3, "in"),
      heights = unit(3, "in"),
      tags = NULL
    )
}

plot_exch_flux <- function(df, enzyme) {
  x <-
    df |>
    dplyr::filter(id == enzyme & type == "exch") |>
    dplyr::filter(treatment %in% c("21%", "0.5%"))

  eq <- unique(x$equation)

  lhs <- stringr::str_extract(eq, "\\w+")
  rhs <- stringr::str_extract(eq, "(?<=-> ).*")

  title <- stringr::str_c(rhs, " → ", lhs)

  ggplot2::ggplot(x) +
    ggplot2::aes(
      x = toupper(cell_type),
      y = flux,
      fill = treatment
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = lb,
        ymax = ub
      ),
      position = ggplot2::position_dodge(width = 0.5),
      width = 0.2,
      size = 0.25
    ) +
    ggplot2::geom_point(
      pch = 21,
      position = ggplot2::position_dodge(width = 0.5),
      color = "white",
      size = 1.5,
      stroke = 0.2
    ) +
    ggplot2::scale_fill_manual(values = clrs, limits = force) +
    ggplot2::labs(
      x = NULL,
      y = "Exchange flux\n(fmol/cell/h)",
      fill = NULL,
      title = title
    ) +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    theme_plots() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "grey80", size = 0.1),
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

plot_lactate_mids <- function(df, cell, t = 72) {
  tracer_labels <- c(
    expression(paste("[1,2-"^13, "C"[2], "] glucose")),
    expression(paste("[U-"^13, "C"[6], "] glucose")),
    expression(paste("[U-"^13, "C"[5], "] glutamine")),
    expression(paste("[U-"^13, "C"[3], "] lactate"))
  )
  tracer_levels <- c("glc2", "glc6", "q5", "lac3")
  metab <- c(
    "LAC",
    "FBP",
    "PYR",
    "CIT"
  )

  x <-
    df |>
    dplyr::filter(
      cell_type == cell &
        time == t &
        metabolite %in% metab &
        tracer == "lac3"
    ) |>
    dplyr::filter(
      !(metabolite == "CIT" & isotope == "M6")
    ) |>
    dplyr::mutate(
      metabolite = factor(metabolite, levels = metab),
      metabolite = forcats::fct_recode(metabolite, "`3PG`" = "3PG"),
      tracer = factor(tracer, levels = tracer_levels, labels = tracer_labels),
      isotope = stringr::str_replace(isotope, "M", ""),
      group = dplyr::case_when(
        oxygen == "21%" & treatment == "None" ~ "21%",
        oxygen == "0.5%" & treatment == "None" ~ "0.5%",
        oxygen == "21%" & treatment == "DMSO" ~ "DMSO",
        oxygen == "21%" & treatment == "BAY" ~ "BAY"
      ),
      group = factor(group, levels = c("21%", "0.5%", "DMSO", "BAY"))
    )

  annot <-
    x |>
    dplyr::group_by(tracer, metabolite) |>
    tidyr::nest() |>
    dplyr::mutate(annot = purrr::map(data, annot_mids_main)) |>
    tidyr::unnest(c(annot))

  ggplot2::ggplot(x) +
    ggplot2::aes(
      x = isotope,
      y = mid
    ) +
    ggplot2::facet_grid(
      tracer ~ metabolite,
      labeller = ggplot2::label_parsed,
      # switch = "y",
      scales = "free_x",
      space = "free_x"
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = group),
      geom = "col",
      fun = "mean",
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.6,
      show.legend = TRUE
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "black",
      lwd = 0.1
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = group),
      geom = "errorbar",
      fun.data = ggplot2::mean_se,
      position = ggplot2::position_dodge(width = 0.6),
      width = 0.2,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = annot,
      ggplot2::aes(
        x = -Inf,
        y = Inf,
        vjust = 1.5,
        hjust = -0.3,
        label = annot
      ),
      family = "Calibri",
      size = 8/ggplot2::.pt
    ) +
    ggplot2::labs(
      x = "Isotope",
      y = "Mole fraction",
      fill = NULL
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, 1, 0.25),
      limits = c(0, 1.1)
    ) +
    ggplot2::theme(
      strip.placement = "outside",
      legend.title = ggplot2::element_blank()
    ) +
    ggplot2::scale_fill_manual(values = clrs, limits = force) +
    theme_plots() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "grey80", size = 0.1),
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10),
      axis.line = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, size = 0.25)
    )
}

plot_hyp_bay_fluxes <- function(df, annot, metab, ylab) {
  df |>
    dplyr::filter(metabolite == metab) |>
    ggplot2::ggplot() +
    ggplot2::aes(
      x = oxygen,
      y = flux
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = treatment),
      geom = "col",
      fun = "mean",
      # width = 0.6,
      position = ggplot2::position_dodge2(),
      show.legend = TRUE,
      alpha = 0.5
    ) +
    ggbeeswarm::geom_beeswarm(
      ggplot2::aes(fill = treatment),
      dodge.width = 0.9,
      pch = 21,
      size = 1,
      stroke = 0.2,
      cex = 4,
      color = "white",
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = treatment),
      geom = "errorbar",
      fun.data = ggplot2::mean_se,
      position = ggplot2::position_dodge(width = 0.9),
      width = 0.2,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(annot, metabolite == metab),
      ggplot2::aes(
        x = oxygen,
        y = y_pos,
        vjust = vjust,
        label = lab,
      ),
      family = "Calibri",
      size = 8/ggplot2::.pt
    ) +
    ggplot2::labs(
      x = "Oxygen",
      y = ylab,
      fill = NULL
    ) +
    ggplot2::scale_fill_manual(
      values = clrs,
      limits = force
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.05, 0.1)),
      breaks = scales::pretty_breaks(n = 6)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(alpha = 1))
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

plot_hyp_bay_densities <- function(df, annot, prot, ylab) {
  annot1 <-
    dplyr::filter(annot, !is.na(treatment)) |>
    dplyr::mutate(
      treatment = factor(treatment, levels = c("DMSO", "BAY")),
      oxygen = factor(oxygen, levels = c("21%", "0.5%"))
    )
  annot2 <-
    dplyr::filter(annot, is.na(treatment)) |>
    dplyr::mutate(
      treatment = factor(treatment, levels = c("DMSO", "BAY")),
      oxygen = factor(oxygen, levels = c("21%", "0.5%"))
    )

  df |>
    dplyr::mutate(
      treatment = factor(treatment, levels = c("DMSO", "BAY")),
      oxygen = factor(oxygen, levels = c("21%", "0.5%"))
    ) |>
    ggplot2::ggplot() +
    ggplot2::aes(
      x = oxygen,
      y = fold_change
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = treatment),
      geom = "col",
      fun = "mean",
      # width = 0.6,
      position = ggplot2::position_dodge2(),
      show.legend = TRUE,
      alpha = 0.5
    ) +
    ggbeeswarm::geom_beeswarm(
      ggplot2::aes(fill = treatment),
      dodge.width = 0.9,
      pch = 21,
      size = 1,
      stroke = 0.2,
      cex = 4,
      color = "white",
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = treatment),
      geom = "errorbar",
      fun.data = ggplot2::mean_se,
      position = ggplot2::position_dodge(width = 0.9),
      width = 0.2,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = annot1,
      ggplot2::aes(
        color = treatment,
        y = y_pos,
        label = lab,
        vjust = vjust
      ),
      family = "Calibri",
      size = 8/ggplot2::.pt,
      position = ggplot2::position_dodge(width = 0.9),
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = annot2,
      ggplot2::aes(
        y = y_pos,
        label = lab,
        vjust = vjust
      ),
      family = "Calibri",
      size = 8/ggplot2::.pt,
      color = "black",
      show.legend = FALSE
    ) +
    ggplot2::labs(
      x = "Oxygen",
      y = ylab,
      fill = NULL,
      color = NULL
    ) +
    ggplot2::scale_fill_manual(
      values = clrs,
      limits = force,
      aesthetics = c("fill", "color")
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.05, 0.1)),
      breaks = scales::pretty_breaks(n = 6)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(alpha = 1))
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

plot_nad <- function(df, annot, metab, ylab) {
  annot <- dplyr::filter(annot, measurement == metab)
  annot1 <-
    dplyr::filter(annot, oxygen != ".") |>
    dplyr::mutate(
      treatment = "BAY",
      treatment = factor(treatment, levels = c("DMSO", "BAY")),
      oxygen = factor(oxygen, levels = c("21%", "0.5%"))
    )
  annot2 <-
    dplyr::filter(annot, treatment != ".") |>
    dplyr::mutate(
      oxygen = "21%",
      treatment = factor(treatment, levels = c("DMSO", "BAY")),
      oxygen = factor(oxygen, levels = c("21%", "0.5%"))
    )

  df |>
    dplyr::filter(treatment != "None" & measurement == metab) |>
    ggplot2::ggplot() +
    ggplot2::aes(
      x = treatment,
      y = value,
      fill = oxygen
    ) +
    ggplot2::stat_summary(
      # ggplot2::aes(fill = treatment),
      geom = "col",
      fun = "mean",
      # width = 0.6,
      position = ggplot2::position_dodge2(),
      show.legend = TRUE,
      alpha = 0.5
    ) +
    ggbeeswarm::geom_beeswarm(
      # ggplot2::aes(fill = treatment),
      dodge.width = 0.9,
      pch = 21,
      size = 1,
      stroke = 0.2,
      cex = 4,
      color = "white",
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      # ggplot2::aes(group = treatment),
      geom = "errorbar",
      fun.data = ggplot2::mean_se,
      position = ggplot2::position_dodge(width = 0.9),
      width = 0.2,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = annot1,
      ggplot2::aes(
        color = oxygen,
        y = y,
        label = lab,
        vjust = vjust
      ),
      family = "Calibri",
      size = 8/ggplot2::.pt,
      position = ggplot2::position_dodge(width = 0.6),
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = annot2,
      ggplot2::aes(
        y = y,
        label = lab,
        vjust = vjust
      ),
      family = "Calibri",
      size = 8/ggplot2::.pt,
      color = "black",
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = clrs,
      limits = force,
      aesthetics = c("fill", "color")
    ) +
    ggplot2::labs(
      x = "Treatment",
      y = ylab,
      fill = NULL,
      color = NULL
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(alpha = 1))
    ) +
    theme_plots() +
    ggplot2::theme(
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    )

}

arrange_fluxes <- function(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10) {
  layout <- "
  abc
  def
  ghi
  jjj
  "

  p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + p10 +
    theme_patchwork(
      design = layout,
      widths = unit(1, "in"),
      heights = unit(1, "in"),
      guides = "collect"
    ) &
    theme(
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

arrange_s1 <- function(p1, p2, p3, p4, p5, p6) {
  layout <- "
  abc
  de#
  fff
  "

  p1 + p2 + p3 + p4 + p5 + p6 +
    theme_patchwork(
      design = layout,
      widths = unit(1, "in"),
      heights = unit(1, "in")
    )
}

arrange_m3 <- function(p1, p2, p3, p4, p5, p6) {
  A <-
    p1 + p2 + p3 +
    patchwork::plot_layout(
      guides = "collect",
      widths = c(1, 1.5, 1.5)
    ) &
    ggplot2::theme(
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    )


  row1 <- wrap_plots(A, p4, nrow = 1, widths = c(4, 1))

  row2 <-
    p5 + p6 +
    patchwork::plot_layout(guides = "collect")

  row1 / row2 +
    theme_patchwork(
      # design = layout,
      widths = ggplot2::unit(7.5, "in"),
      heights = ggplot2::unit(c(1, 4), "in")
    )

}

arrange_s6 <- function(p1, p2) {
  layout <- "
  a
  b
  "

  p1 + p2 +
    theme_patchwork(
      design = layout,
      widths = unit(2.75, "in"),
      heights = unit(3, "in")
    )
}

arrange_s7 <- function(p1, p2, p3, p4, p5) {
  layout <- "
    abc
    def
  "

  p1 + p2 + p3 + p4 + p5 + patchwork::guide_area() +
    theme_patchwork(
      design = layout,
      widths = unit(3.5, "in"),
      heights = unit(c(3.5), "in"),
      guides = "collect"
    )
}

arrange_m4 <- function(p1, p2, p3) {
  layout <- "
  ab
  cc
  "

  (
    (p1 + p2) + patchwork::plot_layout(guides = "collect")
  ) /
    p3 +
    theme_patchwork(
      # design = layout,
      widths = unit(3, "in"),
      heights = unit(c(1, 1), "in")
    ) &
    theme(
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

arrange_m5 <- function(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10) {
  row1 <- patchwork::wrap_plots(p1, p2, p3, guides = "collect")

  row2 <- patchwork::wrap_plots(p4, p5)

  row3 <- patchwork::wrap_plots(p6, p7)

  row4 <- patchwork::wrap_plots(p8, p9, p10, guides = "collect")

  row1 / row2 / row3 / row4 +
    theme_patchwork(
      widths = unit(5, "in"),
      heights = unit(c(1, 1.5, 1.5, 1), "in")
    ) &
    theme(
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

arrange_s8 <- function(p1, p2, p3, p4, p5) {
  layout <- "
  ad
  be
  ce
  "

  p1 + p2 + p3 + p4 + p5 +
    theme_patchwork(
      design = layout,
      widths = unit(2.5, "in"),
      heights = unit(c(1.5), "in")
    )
}

arrange_s9 <- function(p1, p2, p3, p4, p5, p6, p7, p8, p9) {
  layout <- "
  aeg
  beh
  #e#
  cfi
  df#
  "

  p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 +
    theme_patchwork(
      design = layout,
      widths = unit(2.5, "in"),
      heights = unit(c(1.5, 1.5, 1.1, 1.5, 1.5), "in")
    )
}

arrange_m6 <- function(p1, p2, p3, p4, p5, p6, p7) {
  layout <- "
  acdf
  bce#
  "

  p1 + p2 + p3 + p4 + p5 + p6 +
    theme_patchwork(
      design = layout,
      widths = unit(1.5, "in"),
      heights = unit(c(1.5), "in")
    )
}

arrange_m7 <- function(p1, p2, p3, p4, p5, p6, p7, p8, p9) {
  layout <- "
  abc
  def
  ghi
  "

  p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 +
    theme_patchwork(
      design = layout,
      widths = unit(1, "in"),
      heights = unit(c(1), "in")
    )
}

create_resources <- function() {
  tibble::tribble(
    ~category, ~`REAGENT or RESOURCE`, ~SOURCE, ~IDENTIFIER,
    "Antibodies", "HIF-1α", "BD Biosciences", "610958",
    "Antibodies", "c-MYC", "Cell Signaling Technologies", "D84C12",
    "Antibodies", "LDHA", "Cell Signaling Technologies", "2012",
    "Antibodies", "HRP-α-Rabbit IgG", "Cell Signaling Technologies", "7074",
    "Antibodies", "HRP-α-Mouse IgG", "Cell Signaling Technologies", "7076",
    "Bacterial and virus strains", "c-MYC adenovirus", "Vector Biolabs", "1285",
    "Bacterial and virus strains", "YFP adenovirus", "Oldham et al., 2015", "",
    "Chemicals, peptides, and recombinant proteins", "[1,2-^13^C~1~] glucose", "Cambridge Isotope Labs", "CLM-504-PK",
    "Chemicals, peptides, and recombinant proteins", "[U-^13^C~6~] glucose", "Cambridge Isotope Labs", "CLM-1396-PK",
    "Chemicals, peptides, and recombinant proteins", "[U-^13^C~5~] glutamine", "Cambridge Isotope Labs", "CLM-1822-H-PK",
    "Chemicals, peptides, and recombinant proteins", "[U-^13^C~3~] lactate", "Sigma", "485926",
    "Chemicals, peptides, and recombinant proteins", "Molidustat (BAY-85-3934)", "Cayman", "15297",
    "Critical commercial assays", "Glucose colorimetric assay kit", "Cayman", "10009582",
    "Critical commercial assays", "ʟ-Lactate assay kit", "Cayman", "700510",
    "Critical commercial assays", "Pyruvate assay kit", "Cayman", "700470",
    "Depositied data", "Raw and analyzed data", "This paper", "https://github.com/oldhamlab/Copeland.2021.hypoxia.flux",
    "Depositied data", "RNA-seq reads", "This paper", "SRA: PRJNA721596",
    "Depositied data", "Summarized RNA-seq data", "This paper", "https://github.com/oldhamlab/rnaseq.lf.hypoxia.molidustat",
    "Experimental models: Cell lines", "Normal human lung fibroblasts", "Lonza", "CC-2512",
    "Experimental models: Cell lines", "Pulmonary artery smooth muscle cells", "Lonza", "CC-2581",
    "Oligonucleotides", "ACTB (Hs03023943_g1)", "Life Technologies", "4351370",
    "Oligonucleotides", "GLUT1 (Hs00892681_m1)", "Life Technologies", "4351370",
    "Oligonucleotides", "LDHA (Hs00855332_g1)", "Life Technologies", "4351370",
    "Oligonucleotides", "MYC ON-TARGETplus siRNA", "Dharmacon", "L-003282-02-0005",
    "Oligonucleotides", "ON-TARGETplus non-targeting control pool", "Dharmacon", "D-001810-10-05"
  ) |>
    flextable::as_grouped_data(groups = c("category")) |>
    flextable::as_flextable(hide_grouplabel = TRUE) |>
    flextable::bold(j = 1, i = ~ !is.na(category), bold = TRUE, part = "body") |>
    flextable::bold(part = "header", bold = TRUE) |>
    flextable::colformat_double(
      i = ~ is.na(category),
      j = "REAGENT or RESOURCE",
      digits = 0,
      big.mark = ""
    ) |>
    flextable::compose(
      i = 11,
      j = 1,
      part = "body",
      value = flextable::as_paragraph("[1,2-", flextable::as_sup("13"), "C", flextable::as_sub("2"), "] glucose")
    ) |>
    flextable::compose(
      i = 12,
      j = 1,
      part = "body",
      value = flextable::as_paragraph("[U-", flextable::as_sup("13"), "C", flextable::as_sub("6"), "] glucose")
    ) |>
    flextable::compose(
      i = 13,
      j = 1,
      part = "body",
      value = flextable::as_paragraph("[U-", flextable::as_sup("13"), "C", flextable::as_sub("5"), "] glutamine")
    ) |>
    flextable::compose(
      i = 14,
      j = 1,
      part = "body",
      value = flextable::as_paragraph("[U-", flextable::as_sup("13"), "C", flextable::as_sub("3"), "] lactate")
    ) |>
    flextable::font(fontname = "Calibri", part = "all") |>
    flextable::fontsize(size = 9, part = "all") |>
    flextable::set_table_properties(layout = "autofit")
}

format_flux_table <- function(
    flux_differences,
    cell = c("lf", "pasmc"),
    experiment = c("0.5%", "BAY"),
    ssr_ctl = NULL,
    ssr_exp = NULL
) {

  big_border <- flextable::fp_border_default(color = "black", width = 1)
  small_border <- flextable::fp_border_default(color = "black", width = 0.25)

  conditions <- c(unique(flux_differences$ctl), unique(flux_differences$exp))

  df <-
    flux_differences |>
    dplyr::ungroup() |>
    dplyr::filter(normalization == "none" & cell_type == cell & exp == experiment)

  conditions <- c(unique(df$ctl), unique(df$exp))

  df |>
    dplyr::select(-c(normalization, cell_type, index, ctl, exp)) |>
    dplyr::arrange(desc(type)) |>
    dplyr::mutate(
      pathway = stringr::str_to_sentence(pathway),
      type = toupper(type),
      dplyr::across(tidyselect::matches("ctl|exp"), ~scales::scientific(.x))
    ) |>
    flextable::as_grouped_data(groups = c("type", "pathway")) |>
    flextable::as_flextable(hide_grouplabel = TRUE) |>
    flextable::border_remove() |>
    flextable::bold(j = 1, i = ~ !is.na(pathway), bold = TRUE, part = "body") |>
    flextable::bold(j = 1, i = ~ !is.na(type), bold = TRUE, part = "body") |>
    flextable::add_header_row(
      values = c("", conditions, ""),
      colwidths = c(2, 3, 3, 1)
    ) |>
    flextable::compose(
      i = 1,
      j = 3:5,
      part = "header",
      value = flextable::as_paragraph(conditions[[1]], flextable::as_sup("a"))
    ) |>
    flextable::compose(
      i = 1,
      j = 6:8,
      part = "header",
      value = flextable::as_paragraph(conditions[[2]], flextable::as_sup("b"))
    ) |>
    flextable::set_header_labels(
      id = "ID",
      equation = "Reaction",
      flux_ctl = "Flux",
      lb_ctl = "LB",
      ub_ctl = "UB",
      flux_exp = "Flux",
      lb_exp = "LB",
      ub_exp = "UB",
      ratio = "Ratio"
    ) |>
    flextable::align(i = 1, part = "header", align = "center") |>
    flextable::align(i = 2, j = 3:8, part = "header", align = "center") |>
    flextable::merge_h(part = "header") |>
    flextable::bold(part = "header", bold = TRUE) |>
    flextable::colformat_double(
      digits = 2,
      big.mark = ""
    ) |>
    flextable::hline_top(part = "header", border = big_border) |>
    flextable::hline_bottom(part = "all", border = big_border) |>
    flextable::hline(i = ~ !is.na(type), border = small_border) |>
    flextable::hline(i = 1, j = c(3:5, 6:8), border = small_border, part = "header") |>
    flextable::add_footer_lines(c("a", "b")) |>
    flextable::compose(
      i = 1,
      part = "footer",
      value = flextable::as_paragraph(flextable::as_sup("a"), ssr_ctl)
    ) |>
    flextable::compose(
      i = 2,
      part = "footer",
      value = flextable::as_paragraph(flextable::as_sup("b"), ssr_exp)
    ) |>
    flextable::font(fontname = "Calibri", part = "all") |>
    flextable::fontsize(size = 8, part = "all") |>
    flextable::set_table_properties(layout = "autofit")
}
