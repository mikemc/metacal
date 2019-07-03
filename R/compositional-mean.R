# Center (Compositional mean) -------------------------------------------------

# The function `center()` is the only exported function. It allows for multiple
# methods of computing the center, specified with the `method` argument, and
# for specifying the input and output scales and output format, and does some
# preformating to feed to the functions that compute by the specific methods.

# The functions implementing the various methods are `log_center_*()`; these
# take a matrix on the log scale, and provide as output a named vector on the
# log scale.

# TODO: test these functions behave as expected when there is an instance of
# O(s)_i/A(s)_i = \inf or 0; and add better options in this case (e.g., should
# not center the elements)

#' Compute the center (compositional mean) of a set of compositions
#'
#' Unobserved values should be marked as NaN. 
#'
#' @param .data - all-numeric data frame or matrix with taxa as columns
#' @param weights - sample (row) weights
#' @param method - "proj", "gm", or "rss"
#' @param in_scale - "linear" (default) or "log"
#' @param out_scale - "linear" (default) or "log"
#' @param denom - Taxa to use in the denominator; if NULL, use all taxa.
#' @param enframe - whether to return the bias estimate as a two-column tibble
#'
#' @export
center <- function(.data, weights = rep(1, nrow(.data)), method = "proj",
    in_scale = "linear", out_scale = "linear", denom = NULL, enframe = FALSE, 
    components = FALSE) {
    if (!(in_scale %in% c("linear", "log")))
        stop('`in_scale` must be "linear" or "log"')
    if (!(out_scale %in% c("linear", "log")))
        stop('`out_scale` must be "linear" or "log"')
    if (!(method %in% c("proj", "gm", "rss")))
        stop('`method` must be "proj", "gm", or "rss"')

    # The log_center_*()'s require a matrix of log compositions
    mat <- .data %>% as("matrix")
    if (in_scale == "linear")
        mat <- log(mat)

    # Currently, this function will not give meaningful results if `mat`
    # contains any +/- Inf or NA values.
    if ( ! all(is.finite(mat) | is.nan(mat)) )
        stop('Elements of log compositions must be finite or NaN')

    # Check that there is only one component in the co-occurence graph (and
    # that it includes all taxa)
    cmps <- cooccurrence_components(mat)
    if ( !components && (length(unique(cmps)) > 1) )
        # stop('Co-occurrance graph has multiple components')
        stop('The center is not fully determined')

    if (method == "proj") {
        b <- log_center_proj(mat, weights)
    } else if (method == "gm") {
        b <- log_center_gm(mat, weights)
    } else if (method == "rss") {
        b <- log_center_rss(mat, weights, bound = bound)
    }

    if (!is.null(denom))
        b <- b - mean(b[denom])

    if (out_scale == "linear")
        b <- exp(b)

    if (components) {
        b <- left_join(
            tibble::enframe(b, "Taxon", "Center"),
            tibble::enframe(cmps, "Taxon", "Component"), 
            by = "Taxon"
        )
    } else if (enframe) {
        b <- tibble::enframe(b, "Taxon", "Center")
    }

    b
}

#' Projection matrix P_M
#' 
#' Used internally for the "proj" method of computing the compositional mean.
#' See vandenBoogaart2006 and Bren2008.
#'
#' @param K number of elements (taxa)
#' @param M set of missing elements (default is none missing)
proj_mat <- function(K, M = c()) {
    mat <- diag(nrow = K) - 1/(K - length(M))
    mat[M,] <- 0
    mat[,M] <- 0
    mat
}

# Elements are centered w.r.t. the non-NaN elements.
# TODO: add a test for this ^
log_center_gm <- function(mat, weights = rep(1, nrow(mat))) {
    clr_center <- mat %>%
        apply(2, weighted.mean, weights) %>%
        {. - mean(., na.rm = TRUE)}
    clr_center
}

# Method from vandenBoogaart2006; described in vandenBoogaart2013 and at
# https://core.ac.uk/download/pdf/132548286.pdf (Bren2008)
log_center_proj <- function(mat, weights = rep(1, nrow(mat))) {

    # Proper normalization appears to handled by the ginv matrix, so is
    # unnecessary to normalize the weights.

    K <- ncol(mat)

    mat0 <- mat
    mat0[is.nan(mat0)] <- 0

    P_sum <- diag(0, nrow = K)
    v_sum <- rep(0, K)
    for (i in seq(nrow(mat))) {
        M <- which(is.nan(mat[i,]))
        v <- mat0[i,]
        P <- proj_mat(K, M) * weights[i]
        P_sum <- P_sum + P
        v_sum <- v_sum + P %*% v
    }

    clr_center <- (MASS::ginv(P_sum) %*% v_sum) %>% c
    names(clr_center) <- colnames(mat)
    clr_center
}

#' @param bound - single number giving the lower and upper bound on the alr
#'   efficiencies ("rss" only) 
log_center_rss <- function(mat, weights = rep(1, nrow(mat)), bound = 10) {

    # This function computes the squared Aitchison norm of x on the
    # subcomposition of taxa that are observed.
    # 
    # x is a vector of the log compositional error of a sample
    row_rss <- function (x) {
        x %>%
            {. - mean(., na.rm = TRUE)} %>%
            {.^2} %>%
            sum(., na.rm = TRUE)
    }

    # mat is a matrix of the compositional errors on the log scale
    rss <- function(alr_bias, mat, weights) {
        # Substract the alr_bias from each row (sample)
        swept <- sweep(mat, 2, c(alr_bias, 0), "-")
        # Compute and sum the RSS's, weighting samples by `weights`
        apply(swept, 1, row_rss) %>%
            {. * weights} %>%
            sum
    }

    K <- ncol(mat)

    # Initial alr bias values to be passed to `par` arg of `optim()`
    initial <- rep(0, K-1)

    # Compute minimizer
    alr_center <- optim(initial, fn = rss, mat = mat,
        method = "L-BFGS-B", weights = weights, lower = -bound, upper = bound)

    # Format and return results
    clr_center <- c(alr_center$par, 0) %>%
        {. - mean(.)}
    names(clr_center) <- colnames(mat)
    clr_center
}

# Boostrapping ----------------------------------------------------------------

# TODO: 
# - add tests
# - add documentation
# - Test handling of log scale and non-default denominator
# - Warn if not all samples have all taxa and dist = "multinomial"

# Choice of two weight distributions - Multinomial(N, rep(1/S, S)), as in the
# classic bootstrap, and Dirichlet(rep(1, S)), as in the Bayesian bootstrap of
# Rubin.

# N - the choice of N for multinomial samplign.

#' Generate bootstrap replicates of the sample center
#'
#' @export
bootrep_center <- function(.data, R = 4000, N = nrow(.data), method = "proj",
    dist = "dirichlet", in_scale = "linear", out_scale = "linear",
    denom = NULL) {

    if (!(dist %in% c("multinomial", "dirichlet")))
        stop('`dist` must be "multinomial" or "dirichlet"')

    if (N != nrow(.data) & dist == "dirichlet")
        stop('`N != nrow(.data)` only supported with dist = "multinomial"')

    # Test on the full dataset. Ensures that the bias estimate is fully-defined
    # if all samples get positive weight.
    bhat <- center(.data, method = method, in_scale = in_scale, 
        out_scale = out_scale, denom = denom)

    mat <- .data %>% as("matrix")
    if (in_scale == "linear")
        mat <- log(mat)

    # List of weights for each replicate
    if (dist == "multinomial") {
        # Weights ~ Multinomial(N, rep(1, N0) / N0)
        N0 <- nrow(.data)
        wmat <- rmultinom(R, N, rep(1, N0))
        wlist <- lapply(seq(R), function(i) wmat[,i])
    } else if (dist == "dirichlet") {
        # Weights ~ Dirichlet(rep(1, N))
        # Can get a single draw from a Dirichlet(rep(1, N)) by normalizing a
        # vector of K iid exponential random variables to sum to 1; however,
        # normalization is optional when weights are passed to the weighted
        # center functions.
        wlist <- rep(N, R) %>%
            map(rexp, rate = 1)
    }
    names(wlist) <- seq_along(wlist)

    # Define the function to compute Bhat
    if (method == "proj") {
        log_center <- function(w, mat) log_center_proj(mat, w)
    } else if (method == "gm") {
        log_center <- function(w, mat) log_center_gm(mat, w)
    } else if (method == "rss") {
        log_center <- function(w, mat) log_center_rss(mat, w)
    }
    # Compute Bhat for each bootstrap replicate
    reps <- map(wlist, log_center, mat = mat)
    # Set the denominator
    if (!is.null(denom))
        reps <- map(reps, ~ . - mean(.[denom]))
    # Join into a single data frame
    reps <- map_dfr(reps, enframe, "Taxon", "Center", .id = ".id")

    if (out_scale == "linear")
        reps <- mutate(reps, Center = exp(Center))

    reps
}

