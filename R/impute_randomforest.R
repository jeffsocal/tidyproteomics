#' Imputes missing values based on the missForest function
#'
#' @param matrix a matrix with some NAs
#' @param cores the number of threads used to speed the calculation
#'
#' @return a matrix with imputed values
#' @export
#'
impute.randomforest <- function(
    matrix = NULL,
    cores = 2
){

  if(!is.matrix(matrix)) { cli::cli_abort('imput data must be a matrix')}
  if(!is.numeric(cores)) { cli::cli_abort('cores must be an integer, not {cores}')}

  cores = min(floor(cores), ncol(matrix))

  tryCatch({

    suppressWarnings({
      doParallel::registerDoParallel(cores=cores)

      this_matrix <- matrix %>%
        missForest::missForest(
          parallelize = 'variables',
          variablewise = TRUE)
    })

  }, error = function(err) {
    cli::cli_abort(c("x" = as.character(as.vector(err))))
  })

  # cli::cli_alert_info("impute.randomforest OOBEr {this_matrix$OOBerror}")

  return(this_matrix$ximp)
}
