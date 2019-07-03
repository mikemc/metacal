# Functions that can help understand bias estimates

# Taxon co-occurrence network -------------------------------------------------

#' Taxon co-occurrence network
#'
#' Functions for probing the taxon co-occurrence network, used for checking if
#' the bias estimated by the `center()` function is fully-determined.
#'
#' The edges are weighted by the number of samples the pair of taxa co-occur in
#' and so provide some information about the precision of the pairwise bias
#' estimate.
#'
#' @param mat A compositional error matrix
#' @param all Whether to include taxa that aren't in any multi-taxon samples
#' @param enframe Whether to "enframe" the returned vector
#'
#' @name cooccurrence
NULL

#' @describeIn cooccurrence Adjacency matrix of the co-occurrence network
#' @export
cooccurrence_matrix <- function(mat, all = TRUE) {
    ratios <- mat %>%
        as_tibble(rownames = "Sample") %>%
        gather("Taxon", "Value", -Sample) %>%
        compute_ratios()
    if (!all) {
        ratios <- filter(ratios, !is.na(Value))
    }
    adj_mat <- ratios %>%
        filter(Taxon.x != Taxon.y) %>%
        group_by(Taxon.x, Taxon.y) %>%
        summarize(n = sum(!is.na(Value))) %>%
        ungroup %>%
        build_matrix(Taxon.x, Taxon.y, n, 0)
    adj_mat
}

#' @describeIn cooccurrence Co-occurrence network as an `igraph` object
#' @export
cooccurrence_network <- function(mat, all = TRUE) {
    adj_mat <- cooccurrence_matrix(mat, all = all)
    igraph::graph.adjacency(adj_mat, mode = "undirected", weighted = TRUE)
}

#' @describeIn cooccurrence Connected components of the co-occurrence network
#' @export
cooccurrence_components <- function(mat, all = TRUE, enframe = FALSE) {
    g <- cooccurrence_network(mat, all = all)
    components <- igraph::components(g)$membership

    if (enframe)
        components <- tibble::enframe(components, "Taxon", "Component")

    components
}

# TODO: Add example showing usage. Include plotting the graph with the edge
# weights shown
# TODO: Allow a matrix of the actual composition to be passed instead of the
# error matrix 