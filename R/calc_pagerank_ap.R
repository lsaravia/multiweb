#' PageRank importance of species in a food web following Allesina & Pascual (2009)
#'
#' Computes the "food-web PageRank" of Allesina & Pascual (2009). Unlike the
#' standard PageRank applied directly to a food web (e.g. `igraph::page_rank()`
#' with links from prey to predator), which accumulates importance in top
#' predators, this version measures how much other species *depend* on each
#' species, so that importance flows from consumers to their resources.
#'
#' The algorithm:
#' 1. Adds a "root" node (external energy/nutrient source) that is the only
#'    resource of the basal species (species without prey).
#' 2. Adds a link from every species to the root, representing loss of matter
#'    (detritus/recycling). This makes the graph irreducible and primitive, so
#'    no damping factor is needed (damping = 1), as argued by the authors.
#' 3. Each consumer splits its importance equally (or proportionally to link
#'    weights, if `weights = TRUE`) among its resources; the importance vector
#'    is the dominant eigenvector of the resulting column-stochastic matrix.
#'
#' The root is removed from the output and the scores of the species are
#' rescaled to sum 1.
#'
#' The method is defined for trophic networks only: the root node represents
#' the flow of matter, which has no meaning for non-trophic layers.
#'
#' @param g igraph directed food web with links from prey to predator (the
#'   multiweb convention).
#' @param weights logical, if TRUE and `g` has a `weight` edge attribute, each
#'   consumer splits its importance proportionally to the link weights instead
#'   of equally. The links added to/from the root have weight 1 for basal
#'   species and, for the detritus links, the mean weight of the network.
#'   The original method is unweighted (default FALSE).
#' @param remove_loops logical, remove cannibalistic self-loops (default TRUE).
#'
#' @return A data frame sorted by decreasing importance with columns:
#' \item{Species}{species name (or vertex index if the graph has no names).}
#' \item{PageRank}{importance, summing 1 over species (root excluded).}
#' \item{Rank}{rank of the species (1 = most important; ties get the minimum rank).}
#' \item{Basal}{TRUE for species without prey, that are linked to the root.}
#'
#' @references Allesina, S. & Pascual, M. (2009). Googling food webs: can an
#'   eigenvector measure species' importance for coextinctions?
#'   PLoS Computational Biology, 5(9), e1000494. \doi{10.1371/journal.pcbi.1000494}
#'
#' @seealso [calc_centrality()]
#'
#' @examples
#' # Linear chain: resource -> herbivore -> predator
#' g <- igraph::graph_from_literal(Resource -+ Herbivore, Herbivore -+ Predator)
#' calc_pagerank_ap(g)
#'
#' # Compare with the standard PageRank, which ranks the predator first
#' igraph::page_rank(g)$vector
#'
#' @import igraph
#' @export
calc_pagerank_ap <- function(g, weights = FALSE, remove_loops = TRUE) {

  stopifnot(inherits(g, "igraph"))
  if (!igraph::is_directed(g))
    stop("The food web must be directed, with links from prey to predator")

  if (remove_loops)
    g <- igraph::simplify(g, remove.multiple = FALSE, remove.loops = TRUE)

  S <- igraph::vcount(g)
  sp_names <- igraph::V(g)$name
  if (is.null(sp_names)) sp_names <- as.character(seq_len(S))

  use_w <- weights && !is.null(igraph::edge_attr(g, "weight"))
  if (use_w) {
    A <- as.matrix(igraph::as_adjacency_matrix(g, attr = "weight", sparse = FALSE))
    if (any(A < 0)) stop("Link weights must be non-negative")
    w_root <- mean(A[A > 0])
  } else {
    A <- as.matrix(igraph::as_adjacency_matrix(g, sparse = FALSE))
    A[A > 0] <- 1                     # collapse multiple links
    w_root <- 1
  }
  # A[i, j] > 0 means i (prey) -> j (predator)

  basal <- colSums(A) == 0            # species without prey

  # Augmented matrix with the root as node S + 1
  B <- matrix(0, S + 1, S + 1)
  B[1:S, 1:S] <- A
  B[S + 1, which(basal)] <- 1         # root -> basal species
  B[1:S, S + 1] <- w_root             # every species -> root (loss of matter)

  # Each consumer j splits its importance among its resources i (column j)
  M <- sweep(B, 2, colSums(B), "/")   # column-stochastic

  ev <- eigen(M)
  idx <- which.max(Re(ev$values))
  v <- abs(Re(ev$vectors[, idx]))

  v <- v[1:S]                         # drop root
  v <- v / sum(v)

  out <- data.frame(Species  = sp_names,
                    PageRank = v,
                    Rank     = rank(-round(v, 12), ties.method = "min"),
                    Basal    = basal,
                    stringsAsFactors = FALSE)
  out <- out[order(out$Rank), ]
  rownames(out) <- NULL
  out
}
