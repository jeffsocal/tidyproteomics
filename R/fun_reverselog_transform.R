#' Reverse the plot axis for log transformation
#'
#' @param base a numeric
#'
#' @return a ggplot scale transformation
#'
reverselog_transformation <- function(base = 10) {
  scales::trans_new(
    name = paste0("reverselog-", base),
    transform = function(x) -log(x, base),
    inverse = function(x) base^(-x)
  )
}

#' Generate breaks for log transformation
#'
#' @param n a numeric for the approximate number of breaks
#'
#' @return a function that generates breaks
#'
log_breaks <- function(n = 11) {
  function(x) {
    # Filter out 0 or negative values to prevent -Inf in log10
    x <- x[x > 0]

    # Get range in log space
    rng <- log10(range(x, na.rm = TRUE))

    min_val <- floor(rng[1])
    max_val <- ceiling(rng[2])

    # Calculate optimal step size based on the desired number of breaks (n)
    # The max(1, ...) ensures we don't end up with a step size < 1 for small ranges
    step_size <- max(1, round((max_val - min_val) / n))

    # Create a sequence of integers covering that range with the new step size
    breaks <- 10^seq(min_val, max_val, by = step_size)
    return(breaks)
  }
}
