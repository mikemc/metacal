# Functions for computing ratios

# TODO: ensure that comp_vars default won't pick group_vars; add ability to use
# symbols as arguments; Add ability to specify the taxon_var

#' Compute taxon ratios from a tidy data frame
#'
#' Assumes the taxon variable name is `Taxon`
#'
#' @param drop Whether to drop the individual taxon quantities
#' @export
compute_ratios <- function(.data, 
    # taxon_var = "Taxon",
    group_vars = c("Sample"),
    comp_vars = setdiff(names(dplyr::select_if(.data, is.numeric)), 
        group_vars),
    drop = TRUE) {

    .data <- .data %>%
        dplyr::select(Taxon, group_vars, comp_vars)

    # Get all pairs of taxa and join the abundances.
    tb <- tidyr::expand(.data, tidyr::nesting(!!!rlang::syms(group_vars)), 
        Taxon.x = Taxon, Taxon.y = Taxon) %>%
        # mutate(Pair = paste(Taxon.x, Taxon.y, sep = ":")) %>%
        dplyr::left_join(.data, by = c(group_vars, "Taxon.x" = "Taxon")) %>%
        dplyr::left_join(.data, by = c(group_vars, "Taxon.y" = "Taxon"))
    # Each var in `comp_vars` now appears twice, with a ".x" and a ".y" suffix. 

    # Compute the ratios for each var in `comp_vars`:
    for (v in comp_vars) {
        tb[v] <- tb[paste0(v, ".x")] / tb[paste0(v, ".y")]
    }

    if (drop == TRUE) {
        tb <- tb %>%
            dplyr::select(Taxon.x, Taxon.y, group_vars, comp_vars)
    }

    tb
}

# Experimental ----------------------------------------------------------------
# These are experimental and so not exported right now

#' Pairwise ratios of the elements of x
#'
pw_ratios <- function (x) {
    if (is.null(names(x))) 
        names(x) <- seq(x)
    tb <- tidyr::crossing(i = names(x), j = names(x)) %>%
        dplyr::mutate(Name = paste(i, j, sep = ":"))
    v <- purrr::map2_dbl(tb$i, tb$j, ~ x[.x] / x[.y])
    names(v) <- tb$Name
    v
}

# Turn a matrix of compositional data into a tidy data frame of ratio
# observations

#' Pairwise ratios of taxa from matrix of multiple samples
#'
pw_ratios_matrix <- function (x, na.rm = FALSE) {
    if (is.null(colnames(x))) 
        colnames(x) <- seq(x)
    tb <- tidyr::crossing(i = colnames(x), j = colnames(x)) %>%
        dplyr::mutate(Name = paste(i, j, sep = ":"))
    v <- purrr::map2(tb$i, tb$j, ~ x[,.x] / x[,.y])
    names(v) <- tb$Name
    # v
    # tidy df
    tb <- purrr::map(v, tibble::enframe, "Sample", "Ratio") %>%
        dplyr::bind_rows(.id = "Pair")
    if (na.rm) {
        tb <- dplyr::filter(tb, !is.na(Ratio))
    }
    tb
}

#' Ratio matrix
#'
ratio_matrix <- function (x, diag = TRUE, upper = TRUE) {
    K <- length(x)
    mat <- matrix(x, nrow = K, ncol = K) / 
        matrix(x, nrow = K, ncol = K, byrow = TRUE)
    rownames(mat) <- colnames(mat) <- names(x)
    if (!diag) 
        diag(mat) <- NA
    if (!upper) 
        upper.tri(mat) <- NA
    mat
}
