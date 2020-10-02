# TODO 
# - check trim usage
# - add documentation
# - add tests

# Element-wise functions ------------------------------------------------------

#' Odds of the probability vector x
#'
#' @export
odds <- function (x) {
    x / (1-x)
}

#' Logit (log-odds) of the probability vector x
#'
#' @export
logit <- function (x) {
    log(x) - log(1-x)
}

# Summaries of relative abundance vectors--------------------------------------

#' Geometric mean of x
#'
#' @export
gm_mean <- function(x, na.rm = FALSE) {
    exp(mean(log(x), na.rm = na.rm))
}

#' Geometric standard deviation of x
#'
#' Note, uses denominator n-1
#'
#' @export
gm_sd <- function(x, na.rm = FALSE) {
    exp(sd(log(x), na.rm = na.rm))
}

#' Geometric range of x
#'
#' @export
gm_range <- function(x) {
    max(x) / min(x)
}

#' Geometric absolute value of x
#'
#' @export
gm_abs <- function(x) {
    exp(abs(log(x)))
}

#' Geometric (multiplicative) version of the function f
#'
#' @export
gm <- function(f) {
    function (x, ...) exp(f(log(x), ...))
}

# Transformations for relative abundance vectors ------------------------------

#' Close the elements of x to proportions
#'
#' @export
close_elts <- function (x, na.rm = FALSE) {
    x / sum(x, na.rm = na.rm)
}

#' Geometrically center the elements of x
#'
#' @export
center_elts <- function (x, na.rm = FALSE) {
    exp(log(x) - mean(log(x), na.rm = na.rm))
}

#' Compute the centered log-ratio transform of x
#' 
#' @param x Vector of abundances.
#' @param base Base for logarithm
#'
#' @export
clr <- function(x, base = exp(1), na.rm = FALSE) {
    log(x, base = base) - mean(log(x, base = base), na.rm = na.rm)
}

# Distance / dissimilarity between samples ------------------------------------

#' Distance or dissimilarity between relative abundance vectors x and y
#'
#' @param method Distance/dissimilarity measure.
#' @param trim Should x and y be reduced to their common positive elements
#' before computing the Aitchison distance (otherwise, the distance will be
#' Inf).
#' 
#' method == "aitchison" -> Aitchison distance
#' method == "bray" -> Bray-Curtis dissimilarity between x and y. Note,
#' converts x and y to proportions before computing.
#'
#' @export
# Bray method is equal to `vegan::vegdist(rbind(close_elts(x), close_elts(y)))[1]`
xydist <- function(x, y, method = "aitchison", trim = FALSE) {
    if (length(x) != length(y)) {
        stop("x and y have different lengths")
    }
    if (method == "aitchison") {
        if (trim) {
            idx <- (x > 0) & (y > 0)
            x <- x[idx]
            y <- y[idx]
        }
        sqrt( sum( (clr(x) - clr(y))^2 ) )
    } else if (method == "bray") {
        x <- close_elts(x)
        y <- close_elts(y)
        Cxy <- purrr::map2_dbl(x, y, min) %>% sum
        1 - Cxy
    }
}

#' Aitchison norm of x
#'
#' @param na.rm Whether to remove NAs and NaNs before calculating.
#'
#' @export
anorm <- function (x, na.rm = FALSE) {
    x %>%
        clr(., na.rm = na.rm) %>%
        {.^2} %>%
        sum(., na.rm = na.rm) %>%
        sqrt(.)
}
