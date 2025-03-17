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
    guides(fill = FALSE) +
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

#' Compute Eigenvector Centrality for a weigthed/unweigthed Network
#'
#' This function calculates the eigenvector centrality of species in a network
#' based on an adjacency matrix and visualizes the top species contributing to network structure.
#'
#' @param A An igraph object or an adjacency matrix.
#'          The values represent interaction strengths.
#'
#' @return A list containing:
#' \item{EigenCentrality}{A data frame with species names and their centrality scores.}
#' \item{Plot_Centrality}{A ggplot object showing the top 10 species by eigenvector centrality.}
#'
#' @examples
#' # Example usage with an adjacency matrix A
#' results <- calculate_eigencentrality(A)
#' print(results$EigenCentrality)
#' print(results$Plot_Centrality)
#'
#' @import igraph ggplot2 viridis
#' @export
calc_eigencentrality <- function(A) {

  if (inherits(A, "igraph")) {
    graph_A <- A
  } else {
    # Convert matrix A to an igraph object
    graph_A <- graph_from_adjacency_matrix(A, mode = "directed", weighted = TRUE)
  }

  # Compute eigenvector centrality
  eigencentrality_result <- eigen_centrality(graph_A, directed = TRUE)$vector
  #eigencentrality_result <- page_rank(graph_A, directed = TRUE)$vector

  # Create a data frame with species names and their centrality scores
  if(is.null(names(eigencentrality_result))) {
    names(eigencentrality_result) <- 1:vcount(A)
  }

  df_centrality <- tibble(Species = names(eigencentrality_result), Centrality = eigencentrality_result)

  # Sort by centrality
  df_centrality <- df_centrality[order(-df_centrality$Centrality), ]

  # Plot top 10 species by eigencentrality
  p <- ggplot(df_centrality[1:10,], aes(x = Centrality, y = reorder(Species, Centrality), fill = Centrality)) +
    geom_col() +
    scale_fill_viridis_c(option = "D", direction = -1) +
    labs(x = "Centrality Score", y = "Species") +
    theme_minimal() +  guides(fill = FALSE) +
    theme(axis.text.y = element_text(angle = 45, hjust = 1))

  # Return results
  list(
    EigenCentrality = df_centrality,
    Plot_Centrality = p
  )
}



#' Shuffle a network while preserving node degree and calculates the svd/modularity variation
#'
#' @param A Original adjacency matrix (directed or undirected, weighted or unweighted).
#' @param delta Number of link swaps per iteration.
#' @param max_iterations Maximum number of iterations until network metrics stabilize.
#' @param tolerance Tolerance for metric stability (minimum change in modularity or SVD entropy).
#' @param directed Boolean: TRUE if the network is directed, FALSE if undirected.
#' @param modularity_func R function to calculate modularity, default is `cluster_infomap`.
#'
#' @return A list with:
#' \item{New_A}{Shuffled adjacency matrix.}
#' \item{Metrics}{Data frame tracking metric evolution across iterations, including the original network.}
#'
#' @import igraph
#' @export
shuffle_network_deg_svd <- function(input_graph, delta = 10, max_iterations = 100,
                                    tolerance = 1e-3, directed = TRUE, weighted=TRUE, modularity_func = cluster_infomap) {

  # Convert igraph object to adjacency matrix if needed
  if (inherits(input_graph, "igraph")) {
    if(weighted) {
      if(is.null(edge_attr(input_graph, "weight"))) {
        stop( "Weighted graph must have 'weight' attribute on edges.")
      }
      A <- as.matrix(as_adjacency_matrix(input_graph, attr = "weight", sparse = FALSE))
    } else {
      A <- as.matrix(as_adjacency_matrix(input_graph, sparse = FALSE))
    }
  } else {
    A <- input_graph  # Assume it's already an adjacency matrix
  }
  original_A <- A  # Save original matrix
  n <- nrow(A)  # Number of nodes

  mode <- ifelse(directed, "directed", "undirected")
  if(weighted)
    G <- graph_from_adjacency_matrix(A, mode = mode, weighted = TRUE, diag = TRUE)
  else
    G <- graph_from_adjacency_matrix(A, mode = mode, weighted = NULL, diag = TRUE)

  modularity_values <- c(modularity(modularity_func(G)))
  entropy_values <- c(calc_svd_entropy(A)$Entropy)

  for (iter in 1:max_iterations) {
    for (d in 1:delta) {

      # Select two distinct edges (i, j) and (i', j') to swap
      edges <- which(A > 0, arr.ind = TRUE)  # Get all existing edges
      if (nrow(edges) < 2) break  # Ensure at least two edges exist

      selected_edges <- sample(1:nrow(edges), 2)
      i <- edges[selected_edges[1], 1]
      j <- edges[selected_edges[1], 2]
      i_prime <- edges[selected_edges[2], 1]
      j_prime <- edges[selected_edges[2], 2]

      # Ensure valid swap: i ≠ i', j ≠ j', and no duplicate edges
      if (i != i_prime && j != j_prime && A[i, j_prime] == 0 && A[i_prime, j] == 0) {

        # Swap edges
        A[i, j_prime] <- A[i, j]
        A[i_prime, j] <- A[i_prime, j_prime]

        # Remove old edges
        A[i, j] <- 0
        A[i_prime, j_prime] <- 0
      }
    }
    if(weighted) {
      r <- dim(original_A)[1]
      c <- dim(original_A)[2]


      for(i in seq_len(c) ){
        ss <- sample(original_A[,i])
        ss <- ss[ss>0]
        k <- 1
        for( j in seq_len(r)){
          if(A[j,i]>0 ) {
            A[j,i] <- ss[k]
            k <- k+1
          }
        }
      }
    }

    # Compute new modularity and entropy
    G_new <- graph_from_adjacency_matrix(A, mode = mode, weighted = weighted, diag = TRUE)
    modularity_values <- c(modularity_values, modularity(modularity_func(G_new)))
    entropy_values <- c(entropy_values, calc_svd_entropy(A)$Entropy)

    # Check convergence
    if (iter > 1) {
      if (abs(modularity_values[iter] - modularity_values[iter - 1]) < tolerance &&
          abs(entropy_values[iter] - entropy_values[iter - 1]) < tolerance) {
        break
      }
    }
  }

  # Store results
  metrics <- data.frame(
    Iteration = 0:(length(modularity_values) - 1),
    Modularity = modularity_values,
    SVD_Entropy = entropy_values
  )

  return(list(New_A = A, Metrics = metrics))
}



#' Degree-Preserving Network Shuffling via Edge Swaps
#'
#' This function shuffles a directed network while preserving node degrees.
#' It follows a controlled randomization process by swapping edges iteratively.
#'
#' @param A A square adjacency matrix (directed, weighted or unweighted).
#' @param delta Number of edge swaps to perform.
#' @param max_attempts Number of times to attempt a valid swap before stopping.
#'
#' @return A shuffled adjacency matrix preserving in-degree and out-degree.
#'
#' @import igraph
#' @export
shuffle_network_deg <- function(input_graph, delta = 100, directed = TRUE,weighted=TRUE) {
  if(weighted) {
    w_attr <- "weight"
  } else {
    w_attr <- NULL
  }

  # check weighted graph
  if(weighted && is.null(edge_attr(input_graph, "weight"))) {
    stop( "Weighted graph must have 'weight' attribute on edges.")
  }

  # Convert igraph object to adjacency matrix if needed
  if (inherits(input_graph, "igraph")) {
    A <- as.matrix(as_adjacency_matrix(input_graph, attr = w_attr, sparse = FALSE))
  } else {
    A <- input_graph  # Assume it's already an adjacency matrix
  }

  mode <- ifelse(directed, "directed", "undirected")
  original_A <- A   # Save the original network

  for (iter in 1:delta) {

    # Select two distinct edges (i, j) and (i', j') to swap
    edges <- which(A > 0, arr.ind = TRUE)  # Get all existing edges
    if (nrow(edges) < 2) break  # Ensure at least two edges exist

    selected_edges <- sample(1:nrow(edges), 2)
    i <- edges[selected_edges[1], 1]
    j <- edges[selected_edges[1], 2]
    i_prime <- edges[selected_edges[2], 1]
    j_prime <- edges[selected_edges[2], 2]

    # Ensure valid swap: i ≠ i', j ≠ j', and no duplicate edges
    if (i != i_prime && j != j_prime && A[i, j_prime] == 0 && A[i_prime, j] == 0) {

      # Swap edges
      A[i, j_prime] <- A[i, j]
      A[i_prime, j] <- A[i_prime, j_prime]

      # Remove old edges
      A[i, j] <- 0
      A[i_prime, j_prime] <- 0
    }
  }
  if(weighted) {
    r <- dim(original_A)[1]
    c <- dim(original_A)[2]


    for(i in seq_len(c) ){
      ss <- sample(original_A[,i])
      ss <- ss[ss>0]
      k <- 1
      for( j in seq_len(r)){
        if(A[j,i]>0 ) {
          A[j,i] <- ss[k]
          k <- k+1
        }
      }
    }
  }
  return(A)
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
