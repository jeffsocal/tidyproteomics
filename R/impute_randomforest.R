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

  if (!requireNamespace("missForest", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg missForest} is required for random forest imputation. Please install it.")
  }
  cores = min(floor(cores), ncol(matrix))

  if (cores > 1) {
    if (!requireNamespace("doParallel", quietly = TRUE)) {
      cli::cli_abort("Package {.pkg doParallel} is required for parallel random forest imputation. Please install it.")
    }

    on.exit({
      if (requireNamespace("doParallel", quietly = TRUE)) {
        try(doParallel::stopImplicitCluster(), silent = TRUE)
      }
      if (requireNamespace("foreach", quietly = TRUE)) {
        try(foreach::registerDoSEQ(), silent = TRUE)
      }
    }, add = TRUE)
  }

  tryCatch({

    suppressWarnings({
      if (cores > 1) {
        doParallel::registerDoParallel(cores = cores)
        parallelize <- 'variables'
      } else {
        parallelize <- 'no'
      }

      this_matrix <- matrix %>%
        missForest::missForest(
          parallelize = parallelize,
          variablewise = TRUE)
    })

  }, error = function(err) {
    cli::cli_abort(c("x" = as.character(as.vector(err))))
  })

  # cli::cli_alert_info("impute.randomforest OOBEr {this_matrix$OOBerror}")

  return(this_matrix$ximp)
}
