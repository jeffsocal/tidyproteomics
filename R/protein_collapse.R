#' Convert peptide quantitative data into protein quantitative data
#'
#' @description
#' `collapse()` produces a protein based tidyproteomics data-object from a peptide based tidyproteomics data-object.
#'
#' @param data a tidyproteomics data-object
#' @param collapse_to a character string
#' @param assign_by the method to by which to combine peptides into proteins
#' @param top_n a numeric to indicate the N number of peptides summed account for
#' the protein quantitative value, this assumes that peptides have been summed across
#' charge states
#' @param fasta_path if supplied, it will be used to fill in annotation values such as
#' description, protein_name and gene_name
#' @param split_abundance a boolean to indicate if abundances for razor peptides should
#' be split according to protein prevalence
#' @param .verbose a boolean
#' @param .function an assignable protein abundance summary function, fsum, fmean,
#' fgeomean and fmedian have constructed as NAs must be removed.
#' Example: fmedian <- function(x){stats::median(x, na.rm = TRUE)}
#' Example: fquantile <- function(x){stats::quantile(x, .75, na.rm = TRUE)}
#'
#' @return a tidyproteomics data-object
#' @importFrom rlang :=
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#' # data <- hela_peptides %>% collapse()
#' # data %>% summary("sample")
#'
collapse <- function(
    data = NULL,
    collapse_to = 'protein',
    assign_by = c("all-possible","razor-local","razor-global","non-homologous"),
    top_n = Inf,
    split_abundance = FALSE,
    fasta_path = NULL,
    .verbose = TRUE,
    .function = fsum
) {

  # visible bindings
  organism <- NULL
  description <- NULL
  peptides <- NULL
  accession <- NULL
  .data <- NULL
  protein <- NULL
  peptide <- NULL
  modifications <- NULL
  abundance_raw <- NULL
  imputed <- NULL
  n_peptides <- NULL
  num_proteins <- NULL
  abundance <- NULL
  num_psms <- NULL
  num_peptides <- NULL
  identifier <- NULL

  abundance_pro <- NULL
  abundance_shared <- NULL


  function_str <- gsub("\\.", "_", as.character(methods::functionBody(.function))[2])
  assign_by <- rlang::arg_match(assign_by)
  if(mode(top_n) != 'numeric') { cli::cli_abort(c("x" = "top_n must be of type numeric")) }
  if(top_n < 1) { cli::cli_abort(c("x" = "top_n must larger than zero, not `{top_n}`")) }
  if(is.infinite(top_n)) { top_n_str <- 'ALL' } else { top_n_str <- as.character(top_n) }
  if(data$analyte != 'peptides'){
    cli::cli_abort(c("x" = "Analyte type must be {.emph peptides}",
                     "i" = "This dataset is of type {.emph {data$analyte}}"))
  }
  if(!is.logical(split_abundance) && !split_abundance %in% c(TRUE, FALSE)){
    cli::cli_abort(c("x" = "split_abundance must be TRUE or FALSE, not `{split_abundance}`"))
  }
  if(!is.logical(.verbose) && !.verbose %in% c(TRUE, FALSE)){
    cli::cli_abort(c("x" = ".verbose must be TRUE or FALSE, not `{.verbose}`"))
  }

  check_data(data)

  collapse_possible <- data$identifier

  tb_fasta <- c()
  if(!is.null(fasta_path)) {

    count_peptides <- function(x){ x$peptides <- length(x$peptides); return(x) }

    tb_fasta <- fasta_parse(fasta_path)
    tb_fasta <- parallel::mclapply(tb_fasta, fasta_digest)
    tb_fasta <- parallel::mclapply(tb_fasta, count_peptides)
    tb_fasta <- tb_fasta %>%
      dplyr::bind_rows() %>%
      dplyr::select(!c('sequence')) %>%
      dplyr::rename(num_peptides_possible = peptides,
                    protein = accession)

    tb_names <- names(tb_fasta)
    tb_names <- tb_names[which(grepl('name|description|organism', tb_names))]

    collapse_possible <- c(collapse_possible, tb_names)
  }

  collapse_to <- rlang::arg_match(collapse_to, collapse_possible)

  if(.verbose == TRUE){
    cli::cli_alert_info("Collapsing by top {.emph {top_n_str}} peptides into proteins for {.emph {assign_by}} by {.emph {function_str}}")
    cli::cli_progress_bar("...", type = "tasks")
  }

  if(is.infinite(top_n)) { top_n <- 1e5 } # Inf doesn't work, just make the number big

  # one of identity
  experiment_ids <- c("sample_id", "import_file", "sample_file", "sample", "replicate")
  merge_by <- c(collapse_to, experiment_ids)

  if(.verbose == TRUE) {cli::cli_progress_step(' ... gathering peptides')}

  # get the melded peptide table
  tb_peptides <- data %>% meld(single_quant_source = TRUE)

  if(!is.null(fasta_path)) {
    join_fasta_rm <- setdiff(intersect(names(tb_fasta), names(tb_peptides)), data$identifier)
    if(length(join_fasta_rm) > 1) {tb_peptides <- tb_peptides %>% dplyr::select(!dplyr::all_of(join_fasta_rm))}
    join_fasta <- intersect(names(tb_fasta), names(tb_peptides))
    tb_peptides <- tb_peptides %>% dplyr::left_join(tb_fasta, by = join_fasta)

    if(!collapse_to %in% data$identifier){

      join_primary <- intersect(names(tb_fasta), data$identifier)
      join_headers <- intersect(names(tb_fasta), c(collapse_to, data$identifier))

      n_before <- length(unique(unlist(tb_peptides[,join_primary])))

      w <- which(is.na(tb_peptides[,collapse_to]))
      if(length(w) > 0) {tb_peptides <- tb_peptides[-w,]}

      n_after <- length(unique(unlist(tb_peptides[,join_primary])))

      if(n_before != n_after) {cli::cli_alert_danger("Number of proteins droped form {n_before} to {n_after}")}
    }
  }


  # test and indicate that imputation was not accounted for
  if(length(intersect(c("match_between_runs","imputed"), names(tb_peptides))) == 0) {
    tb_peptides <- tb_peptides %>% dplyr::mutate(imputed = 0)
  }

  # merge MBR and imputed into a single is.imputed value
  tb_peptides <- tb_peptides %>%
    tidyr::pivot_longer(
      dplyr::matches("match_between_runs|imputed"),
      names_to = 'method',
      values_to = 'imputed'
    )  %>%
    dplyr::select(!c('method'))

  tb_peptides <- tb_peptides %>%
    dplyr::group_by_at(setdiff(names(tb_peptides), 'imputed')) %>%
    dplyr::summarise(
      imputed = max(imputed),
      .groups = 'drop'
    )

  # add in protein count
  tb_peptides <- tb_peptides %>%
    dplyr::left_join(
      tb_peptides %>%
        dplyr::group_by_at(c(experiment_ids, 'peptide')) %>%
        dplyr::summarise(
          num_proteins = length(unique(protein)),
          .groups = 'drop'
        ),
      by = c(experiment_ids, 'peptide')
    )

  if(.verbose == TRUE) {cli::cli_progress_step(' ... accounting for homology')}

  # recombine for accurate accounting
  if(assign_by == "all-possible"){
    tb_prot_new <- tb_peptides
  } else {

    # slice out razor peptides
    razor_by <- c(setdiff(data$identifier, collapse_to), setdiff(merge_by, data$identifier))
    if(assign_by == "razor-global") {razor_by <- c(setdiff(data$identifier, collapse_to))}

    # remove quant value for shared peptides
    if(assign_by == "non-homologous") {tb_peptides <- tb_peptides %>% dplyr::mutate(abundance = ifelse(num_proteins > 1, 0, abundance))}

    tb_assg <- tb_peptides %>%
      dplyr::group_by_at(merge_by) %>%
      dplyr::summarise(
        n_peptides = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::inner_join(tb_peptides, by = merge_by) %>%
      dplyr::group_by_at(razor_by) %>%
      dplyr::slice_max(n_peptides, n = 1, with_ties = TRUE) %>%
      dplyr::ungroup() %>%
      dplyr::select(c(collapse_to, razor_by)) %>%
      unique()

    tb_prot_new <- tb_peptides %>%
      dplyr::select(!collapse_to) %>%
      unique() %>%
      dplyr::left_join(tb_assg, by = razor_by) %>%
      dplyr::arrange(protein, peptide, sample, replicate) %>%
      dplyr::relocate(protein)
  }

  if(split_abundance == TRUE){
    if(.verbose == TRUE) {cli::cli_progress_step(' ... computing shared peptide abundances')}
    # calculate the protein abundance
    tb_pro_quant <- tb_prot_new %>%
      dplyr::arrange(dplyr::desc(abundance)) %>%
      dplyr::group_by_at(merge_by) %>%
      dplyr::summarise(abundance_pro = .function(abundance[1:top_n]),
                       .groups = 'drop')

    # calculate the shared peptide abundance
    tb_pro_quant_shared <- tb_prot_new %>%
      dplyr::full_join(tb_pro_quant, by = merge_by) %>%
      dplyr::group_by_at(c(experiment_ids, 'peptide', 'modifications')) %>%
      dplyr::mutate(abundance_shared = abundance * abundance_pro / sum(abundance_pro),
                    n_shared = dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(dplyr::desc(abundance)) %>%
      dplyr::mutate(abundance = ifelse(num_proteins > 1, abundance_pro, abundance)) %>%
      dplyr::select(!dplyr::matches('abundance_pro'))

  } else {
    tb_pro_quant_shared <- tb_prot_new %>%
      dplyr::arrange(dplyr::desc(abundance)) %>%
      dplyr::mutate(abundance_shared = abundance)
  }

  if(.verbose == TRUE) {cli::cli_progress_step(' ... computing protein stats')}
  tb_pro_quant_summed <- tb_pro_quant_shared %>%
    dplyr::arrange(dplyr::desc(abundance_shared)) %>%
    dplyr::group_by_at(merge_by) %>%
    # calculate the shared protein abundance
    dplyr::summarise(abundance = .function(abundance_shared[1:top_n]),
                     peptides = paste(sort(unique(peptide)), collapse = "; "),
                     num_peptides = dplyr::n(),
                     num_unique_peptides = length(which(num_proteins == 1)),
                     imputed = sum(imputed) / num_peptides,
                     .groups = 'drop') %>%
    dplyr::rename(identifier := collapse_to) %>%
    # pull into protein groups
    dplyr::group_by_at(c(experiment_ids, 'abundance', 'peptides', 'num_peptides', 'num_unique_peptides', 'imputed')) %>%
    dplyr::summarise(
      num_identifiers = length(unique(identifier)),
      identifiers_grouped = paste(sort(identifier), collapse = "; "),
      identifier = sort(identifier)[1],
      .groups = 'drop'
    ) %>%
    dplyr::relocate(identifier) %>%
    dplyr::mutate(abundance = ifelse(abundance == 0, NA, abundance))

  # rename the abundance column
  colnames(tb_pro_quant_summed)[which(colnames(tb_pro_quant_summed) == 'abundance')] <- paste('abundance', data$quantitative_source, sep="_")

  tb_pro_quant_summed <- as.data.frame(tb_pro_quant_summed)
  tb_pro_quant_summed[,collapse_to] <- tb_pro_quant_summed$identifier
  tb_pro_quant_summed[,paste0("num_", collapse_to)] <- tb_pro_quant_summed$num_identifiers
  tb_pro_quant_summed[,paste0(collapse_to, "_group")] <- tb_pro_quant_summed$identifiers_grouped

  dat_pro <- tb_pro_quant_summed %>%
    tibble::as_tibble() %>%
    dplyr::select(!dplyr::matches('identifier')) %>%
    dplyr::relocate(dplyr::matches(collapse_to))

  if(!is.null(fasta_path)) {dat_pro <- dat_pro %>% dplyr::left_join(tb_fasta, by = collapse_to)}

  dat_pro <- dat_pro %>% codify(identifier = collapse_to, annotations = names(tb_fasta))

  # the output object
  out <- list(
    origin = data$origin,
    analyte = paste0(collapse_to, "s"),
    identifier = collapse_to,
    quantitative_source = data$quantitative_source,
    operations = append(data$operations, glue::glue("Top {top_n_str} peptides summed to proteins according to '{assign_by} by {function_str}'."))
  )

  out <- append(out, dat_pro)
  class(out) <- 'tidyproteomics'

  cli::cli_progress_done()

  return(out)
}
