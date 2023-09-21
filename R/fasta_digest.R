#' Proteolytic digest a parsed fasta list
#'
#' @description
#' `fasta_digest()` Generates peptide sequences based on *enzyme* and *partial* inputs.
#' Only works with the "list" output of the `parse()` function
#'
#' @param protein as character string
#' @param ... parameters for `peptides()`
#'
#' @return a list
#'
#' @examples
#' #\dontrun{
#' proteins <- fasta_parse("~/Local/data/fasta/ecoli_UniProt.fasta")
#' proteins <- fasta_digest(proteins, enzyme = "[K]", partial = 2)
#' }
#'
fasta_digest <- function(
    protein = NULL,
    ...
){
  if(mode(protein) != "list") {cli::cli_abort(c("x" = "protein is `{mode(protein)}`, should be an list"))}
  if(!"sequence" %in% names(protein)) { return("wtf") }
  protein[['peptides']] <- fasta_peptides(protein$sequence, ...)
  return(protein)
}
