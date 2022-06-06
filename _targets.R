# setup -------------------------------------------------------------------

devtools::load_all()
library(targets)
library(tarchetypes)

invisible(
  lapply(
    list.files(path = "R", pattern = "\\.R$", full.names = TRUE),
    source
  )
)

conflicted::conflict_prefer("filter", "dplyr")

options(
  tidyverse.quiet = TRUE,
  usethis.quiet = TRUE,
  dplyr.summarise.inform = FALSE
)

future::plan(future.callr::callr(workers = future::availableCores() - 1))

# target-specific options
tar_option_set(
  packages = c("wmo", "tidyverse", "patchwork"),
  # packages = c("tidyverse", "patchwork", "xcms"),
  # imports = c("rnaseq.lf.hypoxia.molidustat"),
  format = "qs"
)

# list of targets ---------------------------------------------------------

list(

  # dna per cell ------------------------------------------------------------

  tar_target(
    dna_per_cell_file,
    path_to_data("dna-per-cell-number.xlsx"),
    format = "file"
  ),
  tar_target(
    dna_per_cell_raw,
    clean_dna_per_cell(dna_per_cell_file)
  ),
  tar_target(
    dna_per_cell_std,
    make_std_curves(dna_per_cell_raw)
  ),
  tar_target(
    dna_per_cell_clean,
    interp_data(dna_per_cell_raw, dna_per_cell_std)
  ),
  tar_target(
    cells_per_dna,
    calculate_cells_per_dna(dna_per_cell_clean)
  ),
  tar_render(
    dna_per_cell_report,
    path = path_to_reports("dna-per-cell.Rmd"),
    output_dir = system.file("analysis/pdfs", package = "Copeland.2022.hypoxia.flux")
  ),

  # extracellular fluxes ----------------------------------------------------

  tar_target(
    fluxes_meta_files,
    path_to_data("(lf|pasmc)_.*_meta\\.csv"),
    format = "file"
  ),
  tar_target(
    fluxes_meta,
    clean_flux_meta(fluxes_meta_files)
  ),
  tar_target(
    fluxes_data_files,
    path_to_data("(lf|pasmc)_.*_[a-z]_\\d{4}-\\d{2}-\\d{2}\\.xlsx"),
    format = "file"
  ),
  tar_target(
    fluxes_data,
    assemble_flux_data(fluxes_data_files)
  ),
  tar_target(
    conc_raw,
    clean_fluxes(fluxes_data)
  ),
  tar_target(
    conc_std,
    make_std_curves(conc_raw)
  ),
  tar_target(
    conc_std_plots,
    print_plots(conc_std$plots, conc_std$title, "fluxes/01_standard_curves"),
    format = "file"
  ),
  tar_target(
    conc_std_clean_fld,
    make_std_curves(dplyr::filter(conc_raw, !(detector == "fld" & conc > 900)))
  ),
  tar_target(
    conc_std_clean,
    clean_flux_std(conc_raw)
  ),
  tar_target(
    conc_interp,
    interp_data(conc_raw, conc_std_clean)
  ),
  tar_target(
    conc_with_missing,
    fill_missing_fluxes(conc_interp, fluxes_meta)
  ),
  tar_target(
    conc_clean,
    filter_assays(conc_with_missing)
  ),
  tar_target(
    evap_raw,
    assemble_evap_data(fluxes_data)
  ),
  tar_target(
    evap_clean,
    fill_missing_evap(evap_raw, conc_clean)
  ),
  tar_target(
    flux_measurements,
    assemble_flux_measurements(conc_clean, evap_clean)
  ),
  tar_target(
    growth_curves,
    plot_growth_curves(flux_measurements)
  ),
  tar_target(
    growth_curve_plots,
    print_plots(growth_curves$plots, growth_curves$title, "fluxes/02_growth_curves"),
    format = "file"
  ),
  tar_target(
    growth_rates,
    calculate_growth_rates(growth_curves)
  ),
  tar_target(
    degradation_curves,
    plot_degradation_curves(flux_measurements)
  ),
  tar_target(
    degradation_curve_plots,
    print_plots(degradation_curves$plots, degradation_curves$title, "fluxes/03_degradation_curves"),
    format = "file"
  ),
  tar_target(
    degradation_rates,
    calculate_degradation_rates(degradation_curves)
  ),
  tar_target(
    k,
    clean_degradation_rates(degradation_rates)
  ),
  tar_target(
    mass_curves,
    plot_mass_curves(flux_measurements)
  ),
  tar_target(
    mass_curve_plots,
    print_plots(mass_curves$plots, mass_curves$title, "fluxes/04_mass_curves"),
    format = "file"
  ),
  tar_target(
    flux_curves,
    plot_flux_curves(mass_curves, k, growth_rates)
  ),
  tar_target(
    flux_curve_plots,
    print_plots(flux_curves$plots, flux_curves$title, "fluxes/05_flux_curves"),
    format = "file"
  ),
  tar_target(
    fluxes,
    calculate_fluxes(flux_curves)
  ),
  # tar_target(
  #   fluxes_pairwise_annot,
  #   annot_pairwise(fluxes)
  # ),
  tar_render(
    extracellular_fluxes_report,
    path = path_to_reports("extracellular-fluxes.Rmd"),
    output_dir = system.file("analysis/pdfs", package = "Copeland.2022.hypoxia.flux")
  ),

  # write manuscript --------------------------------------------------------

  # tar_target(
  #   template,
  #   system.file("manuscript/template.docx", package = "McGarrity.2022.hypoxia.omics"),
  #   format = "file"
  # ),
  # tar_target(
  #   pkgs,
  #   system.file("manuscript/packages.bib", package = "McGarrity.2022.hypoxia.omics"),
  #   format = "file"
  # ),
  # tar_target(
  #   bib,
  #   system.file("manuscript/library.json", package = "McGarrity.2022.hypoxia.omics"),
  #   format = "file"
  # ),
  # tar_target(
  #   csl,
  #   system.file("manuscript/cell-metabolism.csl", package = "McGarrity.2022.hypoxia.omics"),
  #   format = "file"
  # ),
  # tar_render(
  #   manuscript,
  #   path = path_to_manuscript("manuscript.Rmd"),
  #   output_dir = path_to_manuscript(""),
  #   output_format = bookdown::word_document2(
  #     reference_docx = template,
  #     df_print = "kable",
  #     fig_caption = TRUE,
  #     number_sections = FALSE,
  #     pandoc_args = c(
  #       "--lua-filter=scholarly-metadata.lua",
  #       "--lua-filter=author-info-blocks.lua",
  #       "--lua-filter=pagebreak.lua"
  #     )
  #   ),
  #   params = list(
  #     bibliography = c(bib, pkgs),
  #     csl = csl
  #   )
  # ),
  # tar_render(
  #   supplement,
  #   path = path_to_manuscript("supplement.Rmd"),
  #   output_dir = path_to_manuscript(""),
  #   output_format = bookdown::word_document2(
  #     reference_docx = template,
  #     df_print = "kable",
  #     fig_caption = TRUE,
  #     number_sections = FALSE,
  #     pandoc_args = c(
  #       "--lua-filter=scholarly-metadata.lua",
  #       "--lua-filter=author-info-blocks.lua",
  #       "--lua-filter=pagebreak.lua",
  #       "--lua-filter=multiple-bibliographies.lua"
  #     )
  #   ),
  #   params = list(
  #     bibliography = c(bib, pkgs),
  #     csl = csl
  #   )
  # ),

  NULL
)