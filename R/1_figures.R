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
      panel.border = ggplot2::element_rect(size = 0.25),
      plot.margin = ggplot2::margin(5, 5, 5, 5),
      plot.tag = ggplot2::element_text(face = "bold"),
      axis.title.y.left = ggplot2::element_text(margin = ggplot2::margin(r = 3))
    ),
    ggplot2::coord_cartesian(clip = "off")
  )
}
