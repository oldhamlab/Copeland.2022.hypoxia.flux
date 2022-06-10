# 1_cosmos.R

set_carnival_options <- function(){
  path <- "analysis/carnival"
  carnival_options <- cosmosR::default_CARNIVAL_options(solver = "cplex")
  carnival_options$solverPath <- "/Applications/CPLEX_Studio201/cplex/bin/x86-64_osx/cplex"
  carnival_options$outputFolder <- path
  carnival_options$workdir <- path
  carnival_options$threads <- 0
  carnival_options$clonelog <- -1
  carnival_options$poolIntensity <- 2
  carnival_options$mipGap <- 0.2
  carnival_options
}

get_cosmos_network <- function(){
  data(meta_network, package = "cosmosR")
  meta_network
}

format_tf <- function(tfea) {
  tfea |>
    dplyr::mutate(activity = ifelse(adj.P.Val >= 0.05, 0, sign(t))) |>
    dplyr::select(tf, activity) |>
    tibble::deframe()
}

format_metab <- function(metab) {
  metab |>
    dplyr::select(hmdb, t, P.Value) |>
    dplyr::mutate(activity = dplyr::if_else(P.Value < 0.05, sign(t), 0)) |>
    dplyr::select(hmdb, activity) |>
    tibble::deframe() |>
    cosmosR::prepare_metab_inputs(compartment_codes = c("c", "m"))
}

format_deg <- function(deg) {
  deg |>
    dplyr::filter(!is.na(symbol) & symbol != "") |>
    dplyr::filter(!is.na(padj)) |>
    dplyr::select(symbol, stat, padj) |>
    dplyr::group_by(symbol) |>
    dplyr::arrange(abs(stat), .by_group = TRUE) |>
    dplyr::slice_tail(n = 1) |>
    dplyr::mutate(activity = dplyr::if_else(padj < 0.05, sign(stat), 0)) |>
    dplyr::select(symbol, activity) |>
    tibble::deframe()
}

preprocess <- function(
    route = c("forward", "reverse"),
    network,
    signaling,
    metabolites,
    transcripts,
    carnival_options
) {
  carnival_options$timelimit <- 3600

  network_members <- unique(c(network$source, network$target))
  signaling <- signaling[names(signaling) %in% network_members]
  metabolites <- metabolites[names(metabolites) %in% network_members]

  route <- match.arg(route)
  if (route == "forward") {
    fun <- cosmosR::preprocess_COSMOS_signaling_to_metabolism
    remove_nodes <- TRUE
  } else if (route == "reverse") {
    fun <- cosmosR::preprocess_COSMOS_metabolism_to_signaling
    remove_nodes <- TRUE
  }

  args <-
    list(
      meta_network = network,
      signaling_data = signaling,
      metabolic_data = metabolites,
      diff_expression_data = transcripts,
      diff_exp_threshold = 1,
      maximum_network_depth = 6,
      remove_unexpressed_nodes = remove_nodes,
      CARNIVAL_options = carnival_options
    )

  do.call(fun, args)
}

run_cosmos <- function(
    route = c("forward", "reverse"),
    prep,
    carnival_options
) {
  carnival_options$timelimit <- 4 * 3600

  if (route == "forward") {
    fun <- cosmosR::run_COSMOS_signaling_to_metabolism
  } else if (route == "reverse") {
    fun <- cosmosR::run_COSMOS_metabolism_to_signaling
  }

  args <- list(data = prep, CARNIVAL_options = carnival_options)

  do.call(fun, args)
}
