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
        dplyr::group_by_at(vars(!!!gvs)) %>%
        dplyr::mutate(...) %>%
        dplyr::ungroup()
}

#' Create a matrix from columns of a tidy data frame
#'
#' Note: Setting `fill` will replace NAs and NaNs, along with elements
#' corresponding to missing rows, with the value of `fill`.
#'
#' @param .data A data frame with columns `rows`, `cols`, and `elts`.
#' @param rows Column that will become the rownames of the matrix.
#' @param cols Column that will become the colnames of the matrix.
#' @param elts Column that will become the matrix elements.
#' @param fill Value to use for missing elements.
#' @export
build_matrix <- function(.data, rows, cols, elts, fill = NULL) {
    rows <- rlang::enquo(rows)
    cols <- rlang::enquo(cols)
    elts <- rlang::enquo(elts)
    tb <- .data %>% 
      dplyr::ungroup() %>%
      dplyr::select(!!rows, !!cols, !!elts)
    if (!is.null(fill)) {
        tb <- tb %>% tidyr::complete(!!rows, !!cols,
            fill = rlang::list2(!!elts := fill))
    }
    tb <- tb %>% tidyr::pivot_wider(names_from = !!cols, values_from = !!elts)
    mat <- tb %>% dplyr::select(-!!rows) %>% as("matrix")
    rownames(mat) <- tb %>% dplyr::pull(!!rows)
    mat
}

#' Coerce a (wide) data frame to a matrix
#'
#' @export
as_matrix <- function(.data, rownames = NULL) {
    rownames <- rlang::enquo(rownames)
    if (rlang::quo_is_null(rownames)) {
        mat <- .data %>% as("matrix")
    } else {
        mat <- .data %>% dplyr::select(-!!rownames) %>% as("matrix")
        rownames(mat) <- .data %>% dplyr::pull(!!rownames)
    }
    mat
}
