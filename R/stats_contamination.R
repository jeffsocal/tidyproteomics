#' Assess the relative amount of protein contamination
#'
#' @description
#' `stats_contamination()` is an analysis function that can take a regular
#' expression as a means to assign subsets of proteins as contaminant.
#'
#' @param data tidyproteomics data object
#' @param pattern character string, regular expression
#'
#' @return a tibble
#'
stats_contamination <- function(
    data = NULL,
    pattern = "CRAP"
){

  # visible bindings
  sample <- NULL
  replicate <- NULL
  abundance <- NULL
  file_id <- NULL
  protein_origin <- NULL
  native <- NULL

  check_data(data)

  if( data$analyte != 'proteins'){
    cli::cli_abort(c("x" = "currently only proteins can be analyzed for contamination"))
  }

  # filter out the pattern contaminant proteins
  data_crap <- c()
  data_test <- data %>% extract()
  data_rest <- data %>% subset(description %!like% !!pattern, .verbose = FALSE) %>%
    extract() %>%
    dplyr::mutate(protein_origin = 'native')

  if(data_test %>% nrow() > data_rest %>% nrow()) {
    if(pattern == 'CRAP'){
      data_temp <- data %>% subset(description %like% !!pattern, .verbose = FALSE)

      if(data_temp$quantitative %>% nrow() > 0) {
        data_crap <- list(
          data_temp %>% subset(description %like% "Keratin", .verbose = FALSE) %>%
            extract() %>% dplyr::mutate(protein_origin = 'Keratin'),
          data_temp %>% subset(description %like% "Trypsin", .verbose = FALSE) %>%
            extract() %>% dplyr::mutate(protein_origin = 'Trypsin'),
          data_temp %>% subset(description %like% "BSA", .verbose = FALSE) %>%
            extract() %>% dplyr::mutate(protein_origin = 'BSA'),
          data_temp %>% subset(description %!like% "Keratin|Trypsin|BSA", .verbose = FALSE) %>%
            extract() %>% dplyr::mutate(protein_origin = 'Other')
        ) %>% dplyr::bind_rows()
      } else {
        data_crap <- data_temp
      }
    } else {
      data_crap <- data %>% subset(description %like% !!pattern,
                                      .verbose = FALSE) %>%
        extract() %>% dplyr::mutate(protein_origin = pattern)
    }
  }

  data_sum <- list(data_rest, data_crap) %>% dplyr::bind_rows()

  data_pr <- data_sum %>%
    dplyr::group_by(sample, replicate, protein_origin) %>%
    dplyr::summarise(
      abundance = sum(abundance, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    dplyr::group_by(sample, replicate) %>%
    dplyr::mutate(abundance = abundance/sum(abundance, na.rm = TRUE)) %>%
    dplyr::ungroup()

  data_pc <- data_pr %>%
    dplyr::mutate(abundance = abundance * 100,
                  abundance = paste0(signif(abundance,3), "%")) %>%
    tidyr::pivot_wider(
      names_from = 'protein_origin',
      values_from = 'abundance'
    ) %>%
    dplyr::full_join(data$experiments, by = c("sample", "replicate")) %>%
    dplyr::relocate(file_id) %>%
    dplyr::relocate(native, .after = 'replicate')

  return(data_pc)
}
