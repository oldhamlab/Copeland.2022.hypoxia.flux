# setup -------------------------------------------------------------------

devtools::load_all()
library(targets)
library(tarchetypes)

extrafont::loadfonts(quiet = TRUE)

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
  packages = c("tidyverse", "patchwork", "scales", "grid"),
  # packages = c("tidyverse", "patchwork", "xcms"),
  # imports = c("rnaseq.lf.hypoxia.molidustat"),
  format = "qs",
  # error = "continue"
  error = NULL,
  # debug = "gsea_deg_oxphos"
)

# targets -----------------------------------------------------------------

list(

  # dna curve ---------------------------------------------------------------

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

  # fluxes ------------------------------------------------------------------

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

  # qbias -------------------------------------------------------------------

  tar_target(
    qbias_files,
    path_to_data("q-bias-correction"),
    format = "file"
  ),
  tar_target(
    qbias_ratios,
    import_qbias(qbias_files)
  ),
  tar_target(
    pred_ratios,
    calculate_predicted_ratios()
  ),
  tar_target(
    correction_factors,
    calculate_correction_factors(qbias_ratios, pred_ratios)
  ),
  tar_render(
    qbias_correction_factor_report,
    path = path_to_reports("qbias-correction-factors.Rmd"),
    output_dir = system.file("analysis/pdfs", package = "Copeland.2022.hypoxia.flux")
  ),

  # mids --------------------------------------------------------------------

  tar_target(
    mid_files,
    path_to_data("(a|b)_(fs|sim)_(lf|pasmc)_.*\\.csv"),
    format = "file"
  ),
  tar_target(
    mid_clean,
    clean_mids(mid_files)
  ),
  tar_target(
    mid_correct,
    correct_mid(mid_clean)
  ),
  tar_target(
    mids,
    remove_mid_outliers(mid_correct)
  ),
  tar_target(
    mid_curves,
    plot_mid_curves(mids)
  ),
  tar_target(
    mid_curve_plots,
    print_plots(mid_curves$plots, mid_curves$title, "mids"),
    format = "file"
  ),
  tar_render(
    mid_report,
    path = path_to_reports("mass-isotope-distributions.Rmd"),
    output_dir = system.file("analysis/pdfs", package = "Copeland.2022.hypoxia.flux")
  ),
  # tar_target(
  #   pasmc_m5,
  #   get_m5_citrate(pruned_mids)
  # ),

  # biomass -----------------------------------------------------------------

  tar_target(
    biomass_file,
    path_to_data("cell-mass.csv"),
    format = "file"
  ),
  tar_target(
    biomass_clean,
    clean_biomass(biomass_file)
  ),
  tar_target(
    biomass,
    calculate_biomass(biomass_clean)
  ),
  tar_target(
    biomass_equations,
    calculate_biomass_equations(biomass)
  ),
  tar_target(
    biomass_equations_out,
    write_matlab_input(biomass_equations, coefs, "_biomass.csv"),
    format = "file"
  ),
  tar_render(
    biomass_report,
    path = path_to_reports("biomass.Rmd"),
    output_dir = system.file("analysis/pdfs", package = "Copeland.2022.hypoxia.flux")
  ),

  # matlab ------------------------------------------------------------------

  tar_target(
    reactions_file,
    path_to_reports("modeling/matlab-input/reactions.csv"),
    format = "file"
  ),
  tar_target(
    model_reactions,
    format_reactions(reactions_file)
  ),
  tar_target(
    model_fluxes,
    format_fluxes(growth_rates, fluxes)
  ),
  tar_target(
    model_fluxes_out,
    write_matlab_input(model_fluxes, data, "_fluxes.csv"),
    format = "file"
  ),
  tar_target(
    pruned_mids,
    format_mids(mids)
  ),
  tar_target(
    model_mids,
    summarize_mids(pruned_mids)
  ),
  tar_target(
    model_mids_out,
    write_matlab_input(model_mids, data, "_mids.csv"),
    format = "file"
  ),

  # viability ---------------------------------------------------------------

  tar_target(
    viability_file,
    path_to_data("cell-viability.csv"),
    format = "file"
  ),
  tar_target(
    viability,
    clean_viability(viability_file)
  ),

  # blots -------------------------------------------------------------------

  tar_target(
    blot_files,
    path_to_data("immunoblots"),
    format = "file"
  ),
  tar_target(
    blot_raw,
    read_data(blot_files)
  ),
  tar_target(
    blot_norm,
    normalize_densities(blot_raw)
  ),

  # mrna --------------------------------------------------------------------

  tar_target(
    mrna_files,
    path_to_data("mrna"),
    format = "file"
  ),
  tar_target(
    mrna_raw,
    read_data(mrna_files)
  ),
  tar_target(
    mrna_norm,
    normalize_qpcr(mrna_raw)
  ),

  # model -------------------------------------------------------------------

  tar_target(
    map_flux_files,
    path_to_data("model\\.csv"),
    format = "file"
  ),
  tar_target(
    map_fluxes,
    clean_model_fluxes(map_flux_files, model_reactions)
  ),
  tar_target(
    map_flux_differences,
    assemble_flux_differences(map_fluxes)
  ),
  tar_render(
    map_flux_difference_report,
    path = path_to_reports("flux-differences.Rmd"),
    output_dir = system.file("analysis/pdfs", package = "Copeland.2022.hypoxia.flux")
  ),
  tar_target(
    node_file,
    path_to_data("nodes\\.csv"),
    format = "file"
  ),
  tar_target(
    nodes,
    readr::read_csv(node_file)
  ),
  tar_target(
    lf_hypoxia_graph,
    make_graph(map_flux_differences, nodes, cell = "lf", treat = "21%", normalizer = "none")
  ),
  tar_target(
    lf_hypoxia_graph_ratio_plot,
    plot_ratio_network(lf_hypoxia_graph, "Hypoxia/Normoxia")
  ),
  tar_target(
    bay_graph,
    make_graph(map_flux_differences, nodes, cell = "lf", treat = "DMSO", normalizer = "none")
  ),
  tar_target(
    bay_graph_ratio_plot,
    plot_ratio_network(bay_graph, "BAY/DMSO")
  ),
  tar_target(
    lf_normoxia_graph_plot,
    plot_normoxia_network(lf_hypoxia_graph, "LF\nNormoxia")
  ),
  tar_target(
    lf_hypoxia_growth_graph,
    make_graph(map_flux_differences, nodes, cell = "lf", treat = "0.5%", normalizer = "growth")
  ),
  tar_target(
    lf_hypoxia_growth_graph_plot,
    plot_ratio_network(lf_hypoxia_growth_graph, "Hypoxia/Normoxia\nGrowth Rate Normalized", edges = FALSE)
  ),
  tar_target(
    pasmc_hypoxia_graph,
    make_graph(map_flux_differences, nodes, cell = "pasmc", treat = "21%", normalizer = "none")
  ),
  tar_target(
    pasmc_normoxia_graph_plot,
    plot_normoxia_network(pasmc_hypoxia_graph, "PASMC\nNormoxia")
  ),
  tar_target(
    pasmc_hypoxia_graph_plot,
    plot_ratio_network(pasmc_hypoxia_graph, "PASMC\nHypoxia/Normoxia")
  ),
  tar_target(
    lf_pasmc_normoxia_ratio_fluxes,
    make_cell_ratio_graph(map_fluxes, nodes)
  ),
  tar_target(
    lf_pasmc_normoxia_ratio_plot,
    plot_ratio_network(lf_pasmc_normoxia_ratio_fluxes, "PASMC/LF")
  ),

  # nad ---------------------------------------------------------------------

  tar_target(
    nad_files,
    path_to_data("nad-assay_.*\\.xlsx"),
    format = "file"
  ),
  tar_target(
    nad_data,
    assemble_flux_data(nad_files)
  ),
  tar_target(
    nad_raw,
    clean_nad(nad_data)
  ),
  tar_target(
    nad_conc_std,
    make_std_curves(nad_raw, fo = ~MASS::rlm(value ~ poly(conc, 2, raw = TRUE), data = .x, , maxit = 1000))
  ),
  tar_target(
    nad_conc_std_plots,
    print_plots(nad_conc_std$plots, nad_conc_std$title, "nad/01_standard_curves"),
    format = "file"
  ),
  tar_target(
    nad_interp,
    interp_data(nad_raw, nad_conc_std)
  ),
  tar_target(
    nad_final,
    finalize_nad(nad_interp, cells_per_dna)
  ),
  tar_target(
    nad_annot,
    annot_nad(nad_final)
  ),
  tar_target(
    nad_plot,
    plot_nad(nad_final, nad_annot, "NAD", "NAD\n(nmol/cell)")
  ),
  tar_target(
    nadh_plot,
    plot_nad(nad_final, nad_annot, "NADH", "NADH\n(nmol/cell)")
  ),
  tar_target(
    nadh_ratio_plot,
    plot_nad(nad_final, nad_annot, "Ratio", "NADH/NAD ratio")
  ),

  # rnaseq ------------------------------------------------------------------

  tar_target(
    dds,
    count_rnaseq()
  ),
  tar_target(
    pca_data,
    vst_rnaseq(dds)
  ),
  tar_target(
    rnaseq_pca,
    plot_rnaseq_pca(pca_data)
  ),
  tar_target(
    hallmark_pathways,
    get_msigdb_pathways(category = "H")
  ),
  tar_target(
    tfea,
    run_tfea(dds)
  ),
  tar_target(
    tfea_fit,
    fit_tfea(dds, tfea)
  ),
  tar_map(
    values = list(
      names = c("hyp", "bay", "hyp_bay", "int"),
      comp = list(
        "h.dmso - n.dmso",
        "n.bay - n.dmso",
        "h.bay - n.bay",
        "(h.dmso - n.dmso) - (n.bay - n.dmso)"
      ),
      xlab = c(
        "Hypoxia/Normoxia",
        "BAY/DMSO",
        "Hypoxia/Normoxia in BAY",
        "??Hypoxia/??BAY"
      ),
      colors = list(clrs[2:1], clrs[4:3], clrs[2:1], clrs[c(2, 4)]),
      filename = stringr::str_c("gsea_", c("hyp", "bay", "hyp_bay", "int"))
    ),
    names = names,
    tar_target(
      deg,
      identify_deg(dds, rlang::expr(comp))
    ),
    tar_target(
      rnaseq_vol,
      plot_rnaseq_volcano(deg, xlab = xlab)
    ),
    tar_target(
      gsea,
      run_gsea(deg, hallmark_pathways)
    ),
    tar_target(
      gsea_table,
      plot_gsea_table(gsea, title = xlab, clr = colors)
    ),
    tar_target(
      gsea_table_file,
      write_table(gsea_table, filename),
      format = "file"
    ),
    tar_target(
      gsea_table_plot,
      patchwork::wrap_elements(full = plot_image(gsea_table_file, vjust = 0)) +
        ggplot2::coord_fixed() +
        theme_plots()
    ),
    tar_target(
      tfea_res,
      index_tfea(tfea_fit, names)
    ),
    tar_target(
      tfea_plot,
      plot_tfea_volcanoes(tfea_res, xlab = xlab, colors = colors, nudge = 10)
    ),
    NULL
  ),
  tar_target(
    rnaseq_venn,
    plot_rnaseq_venn(deg_hyp, deg_bay)
  ),
  tar_target(
    gsea_venn,
    plot_gsea_venn(gsea_hyp, gsea_bay)
  ),
  tar_target(
    tfea_venn,
    plot_tfea_venn(tfea_res_hyp, tfea_res_bay)
  ),
  tar_render(
    rnaseq_report,
    path = path_to_reports("rnaseq.Rmd"),
    output_dir = system.file("analysis/pdfs", package = "Copeland.2022.hypoxia.flux")
  ),

  # metab -------------------------------------------------------------------

  tar_target(
    metab_tar_files,
    path_to_data("lf_05-bay_metabolomics-targeted.xlsx"),
    format = "file"
  ),
  tar_target(
    metab_tar_raw,
    format_metab_tar(metab_tar_files)
  ),
  tar_target(
    metab_tar_clean,
    remove_missing_metab(metab_tar_raw) |>
      correct_drift() |>
      quality_control() |>
      impute_missing() |>
      pqn() |>
      log_transform() |>
      annot_metabs()
  ),
  tar_target(
    metab_tar_pca,
    calc_metab_pca(metab_tar_clean)
  ),
  tar_target(
    metab_tar_pca_plot,
    plot_metab_pca(metab_tar_clean, metab_tar_pca)
  ),
  tar_target(
    metab_tar_limma,
    fit_metab_limma(metab_tar_clean)
  ),
  tar_target(
    metab_pathways_file,
    path_to_data("metabolites.tab"),
    format = "file"
  ),
  tar_target(
    metab_pathways,
    read_metab_pathways(metab_pathways_file)
  ),
  tar_map(
    values = list(
      names = list("hyp", "bay", "hyp_bay", "int"),
      colors = list(clrs[2:1], clrs[4:3], clrs[2:1], clrs[c(2, 4)]),
      xlab = list("Hypoxia/Normoxia", "BAY/DMSO", "Hypoxia/Normoxia in BAY", "??Hypoxia/??BAY"),
      title = list("Hypoxia/Normoxia", "BAY/DMSO", "Hypoxia/Normoxia in BAY", "??Hypoxia/??BAY"),
      filename = stringr::str_c("msea_", list("hyp", "bay", "hyp_bay", "int"))
    ),
    names = names,
    tar_target(
      metab_tar_res,
      index_metab_limma(metab_tar_clean, metab_tar_limma, names)
    ),
    tar_target(
      metab_tar_vol,
      plot_metab_volcano(metab_tar_res, colors = colors, xlab = xlab)
    ),
    tar_target(
      metab_tar_msea,
      run_msea(metab_tar_res, metab_pathways)
    ),
    tar_target(
      metab_tar_msea_table,
      plot_msea_table(metab_tar_msea, title, colors, names)
    ),
    tar_target(
      msea_table_file,
      write_table(metab_tar_msea_table, filename),
      format = "file"
    ),
    tar_target(
      msea,
      patchwork::wrap_elements(full = plot_image(msea_table_file, vjust = 0)) +
        ggplot2::coord_fixed() +
        theme_plots()
    ),
    NULL
  ),
  tar_target(
    metab_venn,
    plot_metab_venn(metab_tar_res_hyp, metab_tar_res_bay)
  ),
  tar_render(
    metabolomics_report,
    path = path_to_reports("metabolomics-targeted.Rmd"),
    output_dir = system.file("analysis/pdfs", package = "Copeland.2022.hypoxia.flux")
  ),
  tar_target(
    tca_leading_edge,
    plot_leading_edge(
      metab_tar_res_hyp_bay,
      metab_pathways,
      "KEGG | Citrate cycle (TCA cycle) - Homo sapiens (human)"
    )
  ),

  # myc ---------------------------------------------------------------------

  tar_map(
    values = list(
      names = c("simyc", "oemyc"),
      exp = list("05-simyc", "bay-myc"),
      intervention = list("treatment", "virus"),
      x = rlang::syms(c("treatment", "virus")),
      y = rlang::syms(c("oxygen", "treatment"))
    ),
    names = names,
    tar_target(
      myc_fluxes,
      combine_fluxes(growth_rates, fluxes, exp = exp)
    ),
    tar_target(
      myc_fluxes_annot,
      annot_myc_fluxes(myc_fluxes, intervention)
    ),
    tar_target(
      myc_growth_plot,
      plot_myc(myc_fluxes, myc_fluxes_annot, "growth", "Growth rate (/h)", x = x, fill = y)
    ),
    tar_target(
      myc_lactate_plot,
      plot_myc(myc_fluxes, myc_fluxes_annot, "lactate", "Lactate\n(fmol/cell/h)", x = x, fill = y)
    ),
    NULL
  ),

  # cosmos ------------------------------------------------------------------

  tar_target(
    carnival_options,
    set_carnival_options()
  ),
  tar_target(
    cosmos_network,
    get_cosmos_network()
  ),
  tar_map(
    values = list(
      names = c("hyp", "bay", "hyp_bay", "int"),
      tfs = rlang::syms(stringr::str_c("tfea_res_", c("hyp", "bay", "hyp_bay", "int"))),
      metab = rlang::syms(stringr::str_c("metab_tar_res_", c("hyp", "bay", "hyp_bay", "int"))),
      deg = rlang::syms(stringr::str_c("deg_", c("hyp", "bay", "hyp_bay", "int")))
    ),
    names = names,
    tar_target(
      cosmos_tf,
      format_tf(tfs)
    ),
    tar_target(
      cosmos_metab,
      format_metab(metab)
    ),
    tar_target(
      cosmos_deg,
      format_deg(deg)
    ),
    tar_target(
      cosmos_prep_forward,
      preprocess("forward", cosmos_network, cosmos_tf, cosmos_metab, cosmos_deg, carnival_options)
    ),
    tar_target(
      cosmos_prep_reverse,
      preprocess("reverse", cosmos_network, cosmos_tf, cosmos_metab, cosmos_deg, carnival_options)
    ),
    tar_target(
      cosmos_forward,
      run_cosmos("forward", cosmos_prep_forward, carnival_options)
    ),
    tar_target(
      cosmos_reverse,
      run_cosmos("reverse", cosmos_prep_reverse, carnival_options)
    ),
    tar_target(
      cosmos_res,
      format_cosmos(cosmos_forward, cosmos_reverse)
    ),
    NULL
  ),

  # m1 ----------------------------------------------------------------------

  tar_target(
    lf_hyp_05_timeline_png,
    system.file("manuscript/ai/lf_hyp_05_timeline.png", package = "Copeland.2022.hypoxia.flux"),
    format = "file"
  ),
  tar_target(
    lf_hyp_05_timeline,
    plot_image(lf_hyp_05_timeline_png, scale = 1.6, hjust = 0.2, vjust = 0.1)
  ),
  tar_target(
    lf_hyp_05_growth_curve,
    plot_growth_curve(flux_measurements, cell = "lf", exp = "05")
  ),
  tar_target(
    lf_hyp_05_growth_rate,
    plot_growth_rates(growth_rates, cell = "lf", exp = "05")
  ),
  tar_target(
    lf_hyp_05_blot_png,
    path_to_manuscript("ai/lf_05_hif1a-ldha-blots.png"),
    format = "file"
  ),
  tar_target(
    lf_hyp_05_blot,
    plot_image(lf_hyp_05_blot_png, scale = 1.3, hjust = 0.2, vjust = 0)
  ),
  tar_target(
    lf_hyp_05_hif1a_prot,
    plot_expression(blot_norm, "lf_05", "hif1a", "HIF-1?? protein\n(normalized)")
  ),
  tar_target(
    lf_hyp_05_ldha_prot,
    plot_expression(blot_norm, "lf_05", "ldha", "LDHA protein\n(normalized)")
  ),
  tar_target(
    lf_hyp_05_glut1_rna,
    plot_expression(mrna_norm, "lf_05", "glut1", "GLUT1 mRNA\n(normalized)")
  ),
  tar_target(
    lf_hyp_05_ldha_rna,
    plot_expression(mrna_norm, "lf_05", "ldha", "LDHA mRNA\n(normalized)")
  ),
  tar_target(
    lf_hyp_05_high,
    plot_high_fluxes(fluxes, "lf", "05")
  ),
  tar_target(
    lf_hyp_05_low,
    plot_low_fluxes(fluxes, "lf", "05")
  ),
  tar_target(
    m1,
    arrange_fluxes(
      lf_hyp_05_timeline,
      lf_hyp_05_growth_curve,
      lf_hyp_05_growth_rate,
      lf_hyp_05_blot,
      lf_hyp_05_hif1a_prot,
      lf_hyp_05_glut1_rna,
      lf_hyp_05_ldha_rna,
      lf_hyp_05_ldha_prot,
      lf_hyp_05_high,
      lf_hyp_05_low
    )
  ),
  tar_target(
    m1_figure,
    write_figures(m1, "m1.png"),
    format = "file"
  ),

  # s1 ----------------------------------------------------------------------

  tar_target(
    viability_plot,
    plot_time_lines(viability, y = viability, ylab = "Cell viability (%)", clr = "oxygen")
  ),
  tar_target(
    lf_dna_curve,
    plot_cells_per_dna(dna_per_cell_clean, "lf")
  ),
  tar_target(
    pasmc_dna_curve,
    plot_cells_per_dna(dna_per_cell_clean, "pasmc")
  ),
  tar_target(
    dna_count_hypoxia_file,
    path_to_data("dna-count-hypoxia.csv"),
    format = "file"
  ),
  tar_target(
    dna_count_hypoxia,
    clean_dna_count_hypoxia(dna_count_hypoxia_file)
  ),
  tar_target(
    dna_count_hypoxia_plot,
    plot_dna_count_hypoxia(dna_count_hypoxia)
  ),
  tar_target(
    evap_plot,
    plot_evap_data(evap_clean)
  ),
  tar_target(
    k_plot,
    plot_k(degradation_rates, k)
  ),
  tar_target(
    s1,
    arrange_s1(
      viability_plot,
      lf_dna_curve,
      pasmc_dna_curve,
      dna_count_hypoxia_plot,
      evap_plot,
      k_plot
    )
  ),
  tar_target(
    s1_figure,
    write_figures(s1, "s1.png"),
    format = "file"
  ),

  # s2 ----------------------------------------------------------------------

  tar_target(
    lf_hyp_02_timeline_png,
    system.file("manuscript/ai/lf_hyp_02_timeline.png", package = "Copeland.2022.hypoxia.flux"),
    format = "file"
  ),
  tar_target(
    lf_hyp_02_timeline,
    plot_image(lf_hyp_02_timeline_png, scale = 1.6, hjust = 0.2, vjust = 0.1)
  ),
  tar_target(
    lf_hyp_02_growth_curve,
    plot_growth_curve(flux_measurements, cell = "lf", exp = "02")
  ),
  tar_target(
    lf_hyp_02_growth_rate,
    plot_growth_rates(growth_rates, cell = "lf", exp = "02")
  ),
  tar_target(
    lf_hyp_02_blot_png,
    path_to_manuscript("ai/lf_02_hif1a-ldha-blots.png"),
    format = "file"
  ),
  tar_target(
    lf_hyp_02_blot,
    plot_image(lf_hyp_02_blot_png, scale = 1.3, hjust = 0.2, vjust = 0)
  ),
  tar_target(
    lf_hyp_02_hif1a_prot,
    plot_expression(blot_norm, "lf_02", "hif1a", "HIF-1?? protein\n(normalized)")
  ),
  tar_target(
    lf_hyp_02_ldha_prot,
    plot_expression(blot_norm, "lf_02", "ldha", "LDHA protein\n(normalized)")
  ),
  tar_target(
    lf_hyp_02_glut1_rna,
    plot_expression(mrna_norm, "lf_02", "glut1", "GLUT1 mRNA\n(normalized)")
  ),
  tar_target(
    lf_hyp_02_ldha_rna,
    plot_expression(mrna_norm, "lf_02", "ldha", "LDHA mRNA\n(normalized)")
  ),
  tar_target(
    lf_hyp_02_high,
    plot_high_fluxes(fluxes, "lf", "02")
  ),
  tar_target(
    lf_hyp_02_low,
    plot_low_fluxes(fluxes, "lf", "02")
  ),
  tar_target(
    s2,
    arrange_fluxes(
      lf_hyp_02_timeline,
      lf_hyp_02_growth_curve,
      lf_hyp_02_growth_rate,
      lf_hyp_02_blot,
      lf_hyp_02_hif1a_prot,
      lf_hyp_02_glut1_rna,
      lf_hyp_02_ldha_rna,
      lf_hyp_02_ldha_prot,
      lf_hyp_02_high,
      lf_hyp_02_low
    )
  ),
  tar_target(
    s2_figure,
    write_figures(s2, "s2.png"),
    format = "file"
  ),

  # s3 ----------------------------------------------------------------------

  tar_target(
    pasmc_hyp_05_timeline_png,
    system.file("manuscript/ai/pasmc_hyp_05_timeline.png", package = "Copeland.2022.hypoxia.flux"),
    format = "file"
  ),
  tar_target(
    pasmc_hyp_05_timeline,
    plot_image(pasmc_hyp_05_timeline_png, scale = 1.6, hjust = 0.2, vjust = 0.1)
  ),
  tar_target(
    pasmc_hyp_05_growth_curve,
    plot_growth_curve(flux_measurements, cell = "pasmc", exp = "05")
  ),
  tar_target(
    pasmc_hyp_05_growth_rate,
    plot_growth_rates(growth_rates, cell = "pasmc", exp = "05")
  ),
  tar_target(
    pasmc_hyp_05_blot_png,
    path_to_manuscript("ai/pasmc_05_hif1a-ldha-blots.png"),
    format = "file"
  ),
  tar_target(
    pasmc_hyp_05_blot,
    plot_image(pasmc_hyp_05_blot_png, scale = 1.3, hjust = 0.2, vjust = 0)
  ),
  tar_target(
    pasmc_hyp_05_hif1a_prot,
    plot_expression(blot_norm, "pasmc_05", "hif1a", "HIF-1?? protein\n(normalized)")
  ),
  tar_target(
    pasmc_hyp_05_ldha_prot,
    plot_expression(blot_norm, "pasmc_05", "ldha", "LDHA protein\n(normalized)")
  ),
  tar_target(
    pasmc_hyp_05_glut1_rna,
    plot_expression(mrna_norm, "pasmc_05", "glut1", "GLUT1 mRNA\n(normalized)")
  ),
  tar_target(
    pasmc_hyp_05_ldha_rna,
    plot_expression(mrna_norm, "pasmc_05", "ldha", "LDHA mRNA\n(normalized)")
  ),
  tar_target(
    pasmc_hyp_05_high,
    plot_high_fluxes(fluxes, "pasmc", "05")
  ),
  tar_target(
    pasmc_hyp_05_low,
    plot_low_fluxes(fluxes, "pasmc", "05")
  ),
  tar_target(
    s3,
    arrange_fluxes(
      pasmc_hyp_05_timeline,
      pasmc_hyp_05_growth_curve,
      pasmc_hyp_05_growth_rate,
      pasmc_hyp_05_blot,
      pasmc_hyp_05_hif1a_prot,
      pasmc_hyp_05_glut1_rna,
      pasmc_hyp_05_ldha_rna,
      pasmc_hyp_05_ldha_prot,
      pasmc_hyp_05_high,
      pasmc_hyp_05_low
    )
  ),
  tar_target(
    s3_figure,
    write_figures(s3, "s3.png"),
    format = "file"
  ),

  # m2 ----------------------------------------------------------------------

  tar_target(
    lf_bay_timeline_png,
    system.file("manuscript/ai/lf_bay_timeline.png", package = "Copeland.2022.hypoxia.flux"),
    format = "file"
  ),
  tar_target(
    lf_bay_timeline,
    plot_image(lf_bay_timeline_png, scale = 1.6, hjust = 0.2, vjust = 0.1)
  ),
  tar_target(
    lf_bay_growth_curve,
    plot_growth_curve(flux_measurements, cell = "lf", exp = "bay")
  ),
  tar_target(
    lf_bay_growth_rate,
    plot_growth_rates(growth_rates, cell = "lf", exp = "bay")
  ),
  tar_target(
    lf_bay_blot_png,
    path_to_manuscript("ai/lf_bay_hif1a-ldha-blots.png"),
    format = "file"
  ),
  tar_target(
    lf_bay_blot,
    plot_image(lf_bay_blot_png, scale = 1.3, hjust = 0.2, vjust = 0)
  ),
  tar_target(
    lf_bay_hif1a_prot,
    plot_expression(blot_norm, "lf_bay", "hif1a", "HIF-1?? protein\n(normalized)")
  ),
  tar_target(
    lf_bay_ldha_prot,
    plot_expression(blot_norm, "lf_bay", "ldha", "LDHA protein\n(normalized)")
  ),
  tar_target(
    lf_bay_glut1_rna,
    plot_expression(mrna_norm, "lf_bay", "glut1", "GLUT1 mRNA\n(normalized)")
  ),
  tar_target(
    lf_bay_ldha_rna,
    plot_expression(mrna_norm, "lf_bay", "ldha", "LDHA mRNA\n(normalized)")
  ),
  tar_target(
    lf_bay_high,
    plot_high_fluxes(fluxes, "lf", "bay")
  ),
  tar_target(
    lf_bay_low,
    plot_low_fluxes(fluxes, "lf", "bay")
  ),
  tar_target(
    m2,
    arrange_fluxes(
      lf_bay_timeline,
      lf_bay_growth_curve,
      lf_bay_growth_rate,
      lf_bay_blot,
      lf_bay_hif1a_prot,
      lf_bay_glut1_rna,
      lf_bay_ldha_rna,
      lf_bay_ldha_prot,
      lf_bay_high,
      lf_bay_low
    )
  ),
  tar_target(
    m2_figure,
    write_figures(m2, "m2.png"),
    format = "file"
  ),

  # m3 ----------------------------------------------------------------------

  tar_target(
    mid_glc6_pyr,
    plot_mids(pruned_mids, "lf", "PYR", track = "glc6")
  ),
  tar_target(
    mid_glc6_cit,
    plot_mids(pruned_mids, "lf", "CIT", track = "glc6")
  ),
  tar_target(
    mid_q5_cit,
    plot_mids(pruned_mids, "lf", "CIT", track = "q5")
  ),
  tar_target(
    mid_q5_m5_cit,
    plot_m5_citrate(pruned_mids)
  ),
  tar_target(
    m3,
    arrange_m3(
      mid_glc6_pyr,
      mid_glc6_cit,
      mid_q5_cit,
      mid_q5_m5_cit,
      lf_hypoxia_graph_ratio_plot,
      bay_graph_ratio_plot
    )
  ),
  tar_target(
    m3_figure,
    write_figures(m3, "m3.png"),
    format = "file"
  ),

  # s4 ----------------------------------------------------------------------

  tar_target(
    lf_mids,
    plot_all_mids(pruned_mids, "lf")
  ),
  tar_target(
    s4_figure,
    write_figures(lf_mids, "s4.png")
  ),


  # s5 ----------------------------------------------------------------------

  tar_target(
    pasmc_mids,
    plot_all_mids(pruned_mids, "pasmc", t = 36)
  ),
  tar_target(
    s5_figure,
    write_figures(pasmc_mids, "s5.png")
  ),

  # s6 ----------------------------------------------------------------------

  tar_target(
    time_course_mids,
    format_time_course_mids(model_mids)
  ),
  tar_target(
    mid_time_course,
    plot_mid_time_course(time_course_mids, "lf", c("21%", "0.5%"), "None")
  ),
  tar_target(
    s6_figure,
    write_figures(mid_time_course, "s6.png"),
    format = "file"
  ),

  # s7 ----------------------------------------------------------------------

  tar_target(
    s7,
    arrange_s7(
      lf_normoxia_graph_plot,
      pasmc_normoxia_graph_plot,
      lf_pasmc_normoxia_ratio_plot,
      pasmc_hypoxia_graph_plot,
      lf_hypoxia_growth_graph_plot
    )
  ),
  tar_target(
    s7_figure,
    write_figures(s7, "s7.png")
  ),

  # m4 ----------------------------------------------------------------------

  tar_target(
    rc_fluxes,
    plot_exch_flux(map_fluxes, "IDH")
  ),
  tar_target(
    mct_fluxes,
    plot_exch_flux(map_fluxes, "MCT")
  ),
  tar_target(
    lactate_mids,
    plot_lactate_mids(pruned_mids, "lf")
  ),
  tar_target(
    m4,
    arrange_m4(
      rc_fluxes,
      mct_fluxes,
      lactate_mids
    )
  ),
  tar_target(
    m4_figure,
    write_figures(m4, "m4.png")
  ),

  # m5 ----------------------------------------------------------------------

  tar_target(
    hyp_bay_fluxes_stats,
    analyze_hyp_bay_fluxes(growth_rates, fluxes)
  ),
  tar_target(
    hyp_bay_fluxes_growth,
    plot_hyp_bay_fluxes(hyp_bay_fluxes_stats$data, hyp_bay_fluxes_stats$annot, "growth", "Growth Rate (/h)")
  ),
  tar_target(
    hyp_bay_fluxes_glc,
    plot_hyp_bay_fluxes(hyp_bay_fluxes_stats$data, hyp_bay_fluxes_stats$annot, "glucose", "Glucose\n(fmol/cell/h)") + ggplot2::scale_y_reverse()
  ),
  tar_target(
    hyp_bay_fluxes_lac,
    plot_hyp_bay_fluxes(hyp_bay_fluxes_stats$data, hyp_bay_fluxes_stats$annot, "lactate", "Lactate\n(fmol/cell/h)")
  ),
  tar_target(
    m5,
    arrange_m5(
      hyp_bay_fluxes_growth,
      hyp_bay_fluxes_glc,
      hyp_bay_fluxes_lac,
      metab_tar_pca_plot,
      metab_tar_vol_hyp_bay,
      msea_hyp_bay,
      tca_leading_edge,
      nad_plot,
      nadh_plot,
      nadh_ratio_plot
    )
  ),
  tar_target(
    m5_figure,
    write_figures(m5, "m5.png")
  ),

  # s8 ----------------------------------------------------------------------

  tar_target(
    s8,
    arrange_s8(
      metab_tar_vol_hyp,
      metab_tar_vol_bay,
      metab_venn,
      msea_hyp,
      msea_bay
    )
  ),
  tar_target(
    s8_figure,
    write_figures(s8, "s8.png")
  ),

  # s9 ----------------------------------------------------------------------

  tar_target(
    s9,
    arrange_s9(
      rnaseq_vol_hyp,
      rnaseq_vol_bay,
      rnaseq_venn,
      gsea_venn,
      gsea_table_plot_hyp,
      gsea_table_plot_bay,
      tfea_plot_hyp,
      tfea_plot_bay,
      tfea_venn
    )
  ),
  tar_target(
    s9_figure,
    write_figures(s9, "s9.png")
  ),

  # m6 ----------------------------------------------------------------------

  tar_map(
    values = list(
      sets = list(
        c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT"),
        c("HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2"),
        "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
      ),
      titles = list(
        "E2F Targets and\nG2/M Checkpoint",
        "MYC Targets",
        "Oxidative Phosphorylation"
      ),
      names = list(
        "e2f",
        "myc",
        "oxphos")
    ),
    names = names,
    tar_target(
      gsea_deg,
      plot_pathway_volcanoes(
        deg_hyp_bay,
        hallmark_pathways,
        sets = sets,
        title = titles
      )
    )
  ),
  tar_target(
    m6,
    arrange_m6(
      rnaseq_pca,
      rnaseq_vol_hyp_bay,
      gsea_table_plot_hyp_bay,
      gsea_deg_e2f,
      gsea_deg_myc,
      tfea_plot_hyp_bay
    )
  ),
  tar_target(
    m6_figure,
    write_figures(m6, "m6.png")
  ),

  # m7 ----------------------------------------------------------------------

  tar_target(
    myc_blot_file,
    path_to_manuscript("ai/lf_05-bay_myc-blots.png"),
    format = "file"
  ),
  tar_target(
    myc_blot,
    plot_blot(myc_blot_file)
  ),
  tar_target(
    myc_blot_quant,
    analyze_hyp_bay_densities(blot_norm, "myc")
  ),
  tar_target(
    myc_blot_plot,
    plot_hyp_bay_densities(myc_blot_quant$data, myc_blot_quant$annot, "myc",  "MYC protein\n(normalized)")
  ),
  tar_target(
    model_image_file,
    path_to_manuscript("ai/working-model.png"),
    format = "file"
  ),
  tar_target(
    model_image,
    plot_image(model_image_file, scale = 1.7, hjust = 0.25, vjust = 0.1)
  ),
  tar_target(
    simyc_image_file,
    path_to_manuscript("ai/lf_05-simyc_myc-blots.png"),
    format = "file"
  ),
  tar_target(
    simyc_image,
    plot_image(simyc_image_file)
  ),
  tar_target(
    oemyc_image_file,
    path_to_manuscript("ai/lf_bay-myc_myc-blots.png"),
    format = "file"
  ),
  tar_target(
    oemyc_image,
    plot_image(oemyc_image_file)
  ),
  tar_target(
    m7,
    arrange_m7(
      myc_blot,
      myc_blot_plot,
      model_image,
      simyc_image,
      myc_growth_plot_simyc + ggplot2::geom_hline(yintercept = 0, size = 0.25),
      myc_lactate_plot_simyc,
      oemyc_image,
      myc_growth_plot_oemyc,
      myc_lactate_plot_oemyc
    )
  ),
  tar_target(
    m7_figure,
    write_figures(m7, "m7.png")
  ),

  # tables ------------------------------------------------------------------

  tar_target(
    resources_table,
    create_resources()
  ),
  tar_target(
    lf_hypoxia_table,
    format_flux_table(map_flux_differences, "lf", "0.5%", " SSR 391.7 [311.2-416.6] (95% CI, 362 DOF)", " SSR 334.3 [311.2-416.6] (95% CI, 362 DOF)")
  ),
  tar_target(
    lf_bay_table,
    format_flux_table(map_flux_differences, "lf", "BAY", " SSR 393.5 [311.2-416.6] (95% CI, 362 DOF)", " SSR 392.4 [308.4-413.4] (95% CI, 359 DOF)")
  ),
  tar_target(
    pasmc_hypoxia_table,
    format_flux_table(map_flux_differences, "pasmc", "0.5%", " SSR 575.6 [499.1-630.6] (95% CI, 563 DOF)", " SSR 521.3 [482.2-611.6] (95% CI, 545 DOF)")
  ),

  # manuscript --------------------------------------------------------------

  tar_target(
    template,
    system.file("manuscript/template.docx", package = "Copeland.2022.hypoxia.flux"),
    format = "file"
  ),
  tar_target(
    pkg_citations,
    knitr::write_bib(c(.packages(all.available = TRUE)), path_to_manuscript("packages.bib")),
    cue = tar_cue(mode = "always")
  ),
  tar_target(
    pkgs,
    system.file("manuscript/packages.bib", package = "Copeland.2022.hypoxia.flux"),
    format = "file"
  ),
  tar_target(
    csl,
    system.file("manuscript/embo.csl", package = "Copeland.2022.hypoxia.flux"),
    format = "file"
  ),
  tar_render(
    manuscript,
    path = path_to_manuscript("manuscript.Rmd"),
    output_dir = path_to_manuscript(""),
    output_file = "copeland_2022_flux_manuscript.docx",
    output_format = bookdown::word_document2(
      reference_docx = template,
      df_print = "kable",
      fig_caption = TRUE,
      number_sections = FALSE,
      pandoc_args = c(
        "--lua-filter=scholarly-metadata.lua",
        "--lua-filter=author-info-blocks.lua",
        "--lua-filter=pagebreak.lua"
      )
    )
  ),

  # supplement --------------------------------------------------------------

  tar_render(
    supplement,
    path = path_to_manuscript("supplement.Rmd"),
    output_dir = path_to_manuscript(""),
    output_file = "copeland_2022_flux_supplement.docx",
    output_format = bookdown::word_document2(
      reference_docx = template,
      df_print = "kable",
      fig_caption = TRUE,
      number_sections = FALSE,
      pandoc_args = c(
        "--lua-filter=scholarly-metadata.lua",
        "--lua-filter=author-info-blocks.lua",
        "--lua-filter=pagebreak.lua",
        "--lua-filter=multiple-bibliographies.lua"
      )
    )
  ),
  NULL
)
