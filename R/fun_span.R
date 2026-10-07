#' Helper function for plotting a human readable axis
#'
#' @param data tidyproteomics data object
#' @param pairs the list of vector doublets
#'
#' @return list of vectors
#'
axis_breaks <- function(x) {
  # 1. Define the 'nice' base steps
  base_steps <- c(1, 2, 2.5, 5)

  # 2. Create a range of possible magnitudes (from 0.01 to 100)
  magnitudes <- 10^(-2:2)

  # 3. Generate all possible intervals (0.1, 0.2, 0.25, 0.5, 1, 2, etc.)
  all_intervals <- outer(base_steps, magnitudes, "*") |> as.vector() |> sort()

  # 4. Calculate how many elements each interval would create for [-x, x]
  # Formula: (2 * x) / interval + 1
  counts <- (2 * x) / all_intervals + 1

  # 5. Find the interval that keeps the count between 5 and 11
  # We prioritize the one closest to 7 or 9 for a 'balanced' look
  best_interval <- all_intervals[counts >= 5 & counts <= 11]

  # Fallback: if none fit perfectly, take the one closest to the range
  if (length(best_interval) == 0) {
    best_interval <- all_intervals[which.min(abs(counts - 8))]
  } else {
    best_interval <- max(best_interval) # Larger interval = fewer elements
  }

  # 6. Generate the sequence
  # We use 'round' to clean up floating point errors (e.g., 0.30000000004)
  seq_vec <- seq(from = -x, to = x, by = best_interval)

  return(seq_vec)
}
