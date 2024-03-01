#'
#' @param suffStat data.frame, of which the first column is the factor genotype,
#'                  and subsequent columns contain the traits. The name of the
#'                 first column should be genotype
#'
#' @param covariates data.frame containing covariates, that should always be used
#'                   in each conditional independence test. Should be either NULL (default)
#'                   or a data.frame with the same number of rows as suffStat
#'
#' @return A data-frame with the residuals
#'
#' @author Willem Kruijer and Pariya Behrouzi.
#'         Maintainers: Willem Kruijer \email{willem.kruijer@wur.nl} and
#'        Pariya Behrouzi \email{pariya.behrouzi@gmail.com}
#'
#' @references A paper on arxiv
#'
#' @export
#'
getMeans <- function(suffStat, covariates=NULL) {
  # suffStat <- d; covariates <- data.frame(c1 = rnorm(400))
  names(suffStat)[1] <- 'genotype'
  suffStat$genotype  <- factor(as.character(suffStat$genotype))

  p           <- ncol(suffStat)
  g           <- unique(suffStat$genotype)
  ng          <- length(g)
  m           <- matrix(NA, ng, p-1)

  rownames(m) <- g
  colnames(m) <- names(suffStat)[-1]


  if (!is.null(covariates)) {
    covariates <- as.data.frame(covariates)
    names(covariates) <- paste0('CoVaRiaTe_',1:ncol(covariates))
    stopifnot(nrow(covariates)==nrow(suffStat))
    suffStat <- cbind(suffStat, covariates)
    cv <- paste('~ -1 + ',paste(names(covariates), collapse ='+'), '+ genotype')
  } else {
    cv <- '~ -1 + genotype'
  }

  for (j in 2:p) {
    lm.obj      <- lm(as.formula(paste(names(suffStat)[j], cv)), data = suffStat)
    cff         <- coefficients(lm.obj)#[-(1+nc)]
    if (!is.null(covariates)){cff <- cff[-(1:ncol(covariates))]}
    names(cff)  <- substr(names(cff), start = 9, stop=1000)
    m[names(cff), j-1] <- as.numeric(cff)
  }

  m <- as.data.frame(m)
  m <- data.frame(G = rownames(m), m)
return(m)
}
