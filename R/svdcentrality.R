#' Singular Value Decomposition (SVD) Analysis of Networks
#'
#' This function performs Singular Value Decomposition (SVD) on a network adjacency matrix.
#' It computes species importance based on the dominant singular values and returns entropy, rank, and key visualizations.
#'
#' @param A A square adjacency matrix where rows represent prey and columns represent predators.
#'          The values represent interaction strengths.
#'
#' @return A list containing:
#' \item{Rank}{The numerical rank of the matrix after small values are rounded for numerical stability.}
#' \item{Entropy}{The entropy of the singular value distribution.}
#' \item{Prey_Importance}{A data frame ranking prey species by their contribution to the largest singular value.}
#' \item{Predator_Importance}{A data frame ranking predator species by their contribution to the largest singular value.}
#' \item{Singular_Values}{A numeric vector of the singular values of the matrix.}
#' \item{Plot_Singular_Values}{A ggplot object showing the distribution of singular values.}
#' \item{Plot_Prey_Importance}{A ggplot object showing the top prey species contributing to the largest singular value.}
#' \item{Plot_Predator_Importance}{A ggplot object showing the top predator species contributing to the largest singular value.}
#'
#' @seealso \code{\link{calc_svd_entropy}} for a lighter version returning only entropy and rank.
#'
#' @examples
#'
#' results <- calc_svd_entropy_importance(netData[[1]])
#' print(results$Rank)
#' print(results$Entropy)
#' print(results$Plot_Singular_Values)
#' print(results$Plot_Prey_Importance)
#' print(results$Plot_Predator_Importance)
#'
#' @import ggplot2 dplyr
#' @export
calc_svd_entropy_importance <- function(A, threshold_factor = 1e-6) {

  # Detect if input is an igraph object and extract adjacency matrix
  if (inherits(A, "igraph")) {
    if(is.null(edge_attr(A, "weight"))) {
      A <- as_adjacency_matrix(A, sparse = FALSE)
    } else {
      A <- as_adjacency_matrix(A, attr="weight", sparse = FALSE)
    }
  }

  # Check if A is a valid matrix
  if (!is.matrix(A) || nrow(A) != ncol(A)) {
    stop("A must be a square adjacency matrix or an igraph object.")
  }
  # Compute SVD
  svd_result <- svd(A)

  # Extract row and column names from matrix A

  prey_names <- rownames(A)
  predator_names <- colnames(A)
  if(is.null(prey_names)) {
    prey_names <- 1:nrow(A)
  }
  if(is.null(predator_names)) {
    predator_names <- 1:nrow(A)
  }

  # Normalize singular values
  p <- svd_result$d / sum(svd_result$d)

  # Compute effective rank (rounding small values)
  threshold <- max(svd_result$d) * threshold_factor
  rank_approx <- sum(svd_result$d > threshold)

  # Compute SVD entropy
  svd_entropy <- -log(rank_approx) *sum(p[p > 0] * log(p[p > 0]))

  # Find species contributing to the largest singular value
  first_singular_vector_prey <- abs(svd_result$u[, 1])  # Prey contributions
  first_singular_vector_pred <- abs(svd_result$v[, 1])  # Predator contributions

  # Assign names
  names(first_singular_vector_prey) <- prey_names
  names(first_singular_vector_pred) <- predator_names

  # Sort species by importance
  prey_importance <- sort(first_singular_vector_prey, decreasing = TRUE)
  predator_importance <- sort(first_singular_vector_pred, decreasing = TRUE)

  # Convert to data frame for ggplot
  df_singular_values <- tibble(Index = seq_along(svd_result$d), Value = svd_result$d)
  df_prey_importance <- tibble(Prey = names(prey_importance), Importance = prey_importance)
  df_predator_importance <- tibble(Predator = names(predator_importance), Importance = predator_importance)

  # Plot Singular Values Distribution
  p1 <- ggplot(df_singular_values[1:rank_approx,], aes(x = Index, y = Value)) +
    geom_point(color = "blue") +
    geom_line(color = "blue") +
    labs(title = "Singular Value Spectrum", x = "Index", y = "Singular Value") +
    theme_bw()

  # Plot Top Contributing Prey
  p2 <- ggplot(df_prey_importance[1:10,], aes(x = reorder(Prey, Importance),
                                              y = Importance, fill=Importance)) +
    geom_col() +
    scale_fill_viridis_c(option = "D", direction = -1) +
    guides(fill = "none") +
    coord_flip() +
    labs(x = "Prey", y = "Contribution") +
    theme_bw()

  # Plot Top Contributing Predators
  p3 <- ggplot(df_predator_importance[1:10,], aes(x = reorder(Predator, Importance),
                                                  y = Importance, fill=Importance)) +
    geom_col() +
    scale_fill_viridis_c(option = "D", direction = -1) +
    guides(fill = FALSE) +
    coord_flip() +
    labs(x = "Predator", y = "Contribution") +
    #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme_bw()

  # Return results
  list(
    Rank = rank_approx,
    Entropy = svd_entropy,
    ImportancePrey = df_prey_importance,
    ImportancePredators = df_predator_importance,
    Plot_Singular_Values = p1,
    Plot_Prey_Importance = p2,
    Plot_Predator_Importance = p3
  )
}

#' Compute Centrality for a Weighted/Unweighted Network
#'
#' This function calculates the centrality of species in a network
#' based on an adjacency matrix or igraph object using a specified centrality measure.
#'
#' @param A An igraph object or an adjacency matrix.
#'          The values represent interaction strengths.
#' @param centrality_func A function to compute centrality. Default is `eigen_centrality`.
#'        Other options: `page_rank`, `betweenness`, `degree`, etc.
#'
#' @return A list containing:
#' \item{CentralityScores}{A data frame with species names and their centrality scores.}
#' \item{Plot_Centrality}{A ggplot object showing the top 10 species by centrality.}
#'
#' @examples
#' # Example usage with an adjacency matrix A
#' results <- calc_centrality(A, centrality_func = page_rank)
#' print(results$CentralityScores)
#' print(results$Plot_Centrality)
#'
#' @import igraph ggplot2 viridis
#' @export
calc_centrality <- function(A, centrality_func = eigen_centrality) {

  if (inherits(A, "igraph")) {
    graph_A <- A
  } else {
    # Convert matrix A to an igraph object
    graph_A <- graph_from_adjacency_matrix(A, mode = "directed", weighted = TRUE)
  }

  # Compute centrality using the specified function
  centrality_result <- centrality_func(graph_A, directed = TRUE)

  # Extract the vector if function returns a list
  if (is.list(centrality_result) && "vector" %in% names(centrality_result)) {
    centrality_result <- centrality_result$vector
  }

  # Create a data frame with species names and their centrality scores
  if (is.null(names(centrality_result))) {
    names(centrality_result) <- 1:vcount(graph_A)
  }

  df_centrality <- tibble(Species = names(centrality_result), Centrality = centrality_result)

  # Sort by centrality
  df_centrality <- df_centrality[order(-df_centrality$Centrality), ]

  # Plot top 10 species by centrality
  p <- ggplot(df_centrality[1:min(10, nrow(df_centrality)),],
              aes(x = Centrality, y = reorder(Species, Centrality), fill = Centrality)) +
    geom_col() +
    scale_fill_viridis_c(option = "D", direction = -1) +
    labs(x = "Centrality Score", y = "Species") +
    theme_minimal() + guides(fill = FALSE) +
    theme(axis.text.y = element_text(angle = 45, hjust = 1))

  # Return results
  list(
    CentralityScores = df_centrality,
    Plot_Centrality = p
  )
}





#' Singular Value Decomposition (SVD) Entropy and Rank Calculation
#'
#' This function calculates the entropy and rank of a network using Singular Value Decomposition (SVD).
#' It accepts either an adjacency matrix or an igraph object and processes it accordingly.
#'
#' @param A An adjacency matrix (square, weighted or unweighted) or an igraph object.
#' @param threshold_factor A small value (default = 1e-6) to determine numerical rank approximation.
#'
#' @return A list containing:
#' \item{Rank}{Effective rank of the matrix.}
#' \item{Entropy}{SVD entropy of the network.}
#'
#' @import igraph
#' @export
calc_svd_entropy <- function(A, threshold_factor = 1e-6) {
  #library(igraph)

  # Convert igraph object to adjacency matrix if necessary
  if (inherits(A, "igraph")) {
    if(is.null(edge_attr(A, "weight"))) {
       A <- as.matrix(as_adjacency_matrix(A, sparse = FALSE))
    } else {
      A <- as.matrix(as_adjacency_matrix(A, attr = "weight", sparse = FALSE))
    }
  }

  # Ensure A is a valid matrix
  if (!is.matrix(A)) stop("Input must be an adjacency matrix or an igraph object.")

  # Compute SVD
  svd_result <- svd(A)

  # Normalize singular values
  p <- svd_result$d / sum(svd_result$d)

  # Compute effective rank (rounding small values)
  threshold <- max(svd_result$d) * threshold_factor
  rank_approx <- sum(svd_result$d > threshold)

  # Compute SVD entropy
  svd_entropy <- -1/log(rank_approx) * sum(p[p > 0] * log(p[p > 0]))

  # Return results
  list(
    Rank = rank_approx,
    Entropy = svd_entropy
  )
}


