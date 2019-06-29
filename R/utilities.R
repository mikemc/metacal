#' @name corner
#' @importFrom useful corner
#' @export
useful::corner

# Helpers for working with "tidy" microbiome data -----------------------------

#' Mutate within groups
#'
#' @export
mutate_by <- function(.data, group_vars, ...) {
    gvs <- rlang::enquos(group_vars)
    .data %>%
        group_by_at(vars(!!!gvs)) %>%
        mutate(...) %>%
        ungroup
}

#' Create a matrix from columns of a tidy data frame
#'
#' Note: Setting `fill` will replace NAs and NaNs, along with elements
#' corresponding to missing rows, with the value of `fill`.
#'
#' @param .data A data frame with columns `rows`, `cols`, and `elts`.
#' @param rows Column that will become the rownames of the matrix.
#' @param cols Column that will become the rownames of the matrix.
#' @param elts Column that will become the matrix elements.
#' @param fill Value to use for missing elements.
#' @export
build_matrix <- function(.data, rows, cols, elts, fill = NULL) {
    rows <- rlang::enquo(rows)
    cols <- rlang::enquo(cols)
    elts <- rlang::enquo(elts)
    tb <- .data %>% select(!!rows, !!cols, !!elts)
    if (!is.null(fill)) {
        tb <- tb %>% 
            complete(!!rows, !!cols, fill = rlang::list2(!!elts := fill))
    }
    tb <- tb %>% spread(!!cols, !!elts)
    mat <- tb %>% select(-!!rows) %>% as("matrix")
    rownames(mat) <- tb %>% pull(!!rows)
    mat
}

#' Coerce a (wide) tibble to a matrix
#'
#' @export
as_matrix <- function(tb, rownames = NULL) {
    rownames <- rlang::enquo(rownames)
    if (rlang::quo_is_null(rownames)) {
        mat <- tb %>% as("matrix")
    } else {
        mat <- tb %>% select(-!!rownames) %>% as("matrix")
        rownames(mat) <- tb %>% dplyr::pull(!!rownames)
    }
    mat
}


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
    tb <- tidyr::expand(.data, tidyr::nesting(!!!syms(group_vars)), 
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
