# 1_figures.R

# plot setup
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
    ggplot2::stat_summary(
      geom = "linerange",
      fun.data = ggplot2::mean_se,
      size = 0.5,
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
    pairs() |>
    broom::tidy() |>
    dplyr::mutate(
      group = stringr::str_extract(contrast, "(?<= - ).*"),
      y = Inf,
      vjust = 1,
      label = annot_p(p.value)
    )

  x |>
    dplyr::group_by(date) |>
    dplyr::mutate(
      grand_mean = mean(x$mu),
      group_mean = mean(mu),
      adj = group_mean - grand_mean,
      mu_corr = mu + adj
    ) |>
    ggplot2::ggplot() +
    ggplot2::aes(
      x = group,
      y = mu_corr
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = group),
      geom = "col",
      fun = "mean",
      # width = 0.6,
      show.legend = FALSE,
      alpha = 0.5
    ) +
    ggplot2::stat_summary(
      geom = "errorbar",
      fun.data = ggplot2::mean_se,
      width = 0.2,
      size = 0.25,
      show.legend = FALSE
    ) +
    ggbeeswarm::geom_beeswarm(
      ggplot2::aes(fill = group),
      pch = 21,
      size = 1,
      stroke = 0.25,
      cex = 4,
      color = "white",
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = annot,
      ggplot2::aes(
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
      ylim = c(0, NA),
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
    ggplot2::stat_summary(
      geom = "linerange",
      fun.data = "mean_se",
      size = 0.5,
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      geom = "point",
      fun = "mean",
      pch = 21,
      color = "white",
      fill = "black",
      size = 1.5,
      show.legend = FALSE
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
      show.legend = FALSE) +
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
  evap_clean %>%
    dplyr::filter(experiment == "05" & cell_type == "lf") %>%
    dplyr::mutate(
      oxygen = factor(oxygen, levels = c("21%", "0.5%"))
    ) %>%
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
    ggplot2::stat_summary(
      geom = "linerange",
      fun.data = "mean_se",
      size = 0.5,
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      geom = "point",
      fun = "mean",
      pch = 21,
      color = "white",
      size = 1.5,
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
    k %>%
    dplyr::select(-k) %>%
    dplyr::mutate(
      label = "*",
      ypos = Inf,
      vjust = 1.5
    )

  degradation_rates %>%
    dplyr::mutate(
      group = dplyr::case_when(
        oxygen == "21%" & treatment == "None" ~ "21%",
        oxygen == "0.5%" ~ "0.5%",
        treatment == "DMSO" ~ "DMSO"
      ),
      group = factor(group, levels = c("21%", "0.5%", "DMSO"))
    ) %>%
    dplyr::left_join(annot, by = c("metabolite", "oxygen", "treatment")) %>%
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
      fun.data = "mean_se",
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
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.2)) +
    theme_plots() +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(color = "gray90"),
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.position = "bottom",
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

arrange_fluxes <- function(p1, p2, p3, p4, p5, p6, p7, p8, p9) {
  layout <- "
  abc
  def
  gh#
  "

  p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 +
    theme_patchwork(
      design = layout,
      widths = unit(1, "in"),
      heights = unit(1.25, "in"),
      guides = "collect"
    ) &
    theme(
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(0.5, "lines"),
      legend.key.height = ggplot2::unit(0.5, "lines"),
      legend.box.margin = ggplot2::margin(t = -10)
    )
}

arrange_s1 <- function(p1, p2, p3, p4, p5) {
  layout <- "
  abc
  de#
  "

  p1 + p2 + p3 + p4 + p5 +
    theme_patchwork(
      design = layout,
      widths = unit(1, "in"),
      heights = unit(1.25, "in")
    )
}
