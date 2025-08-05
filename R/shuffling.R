#' Generate a Sequence of Shuffled Networks and Track SVD/Modularity Evolution
#'
#' This function generates a sequence of shuffled networks from an original graph by applying
#' incremental shuffling steps. It tracks modularity and SVD entropy throughout the process,
#' and stops when the change in these metrics falls below a specified `tolerance` threshold.
#'
#' @param input_graph An igraph object or adjacency matrix (directed or undirected, weighted or unweighted).
#' @param delta Integer. Number of link swaps per iteration (default: 10).
#' @param max_iterations Integer. Maximum number of iterations until network metrics stabilize (default: 100).
#' @param tolerance Numeric. Tolerance for metric stability - minimum change in modularity or SVD entropy (default: 1e-3).
#' @param directed Logical. Whether the network is directed (default: TRUE).
#' @param weighted Logical. Whether the network is weighted (default: TRUE).
#' @param modularity_func Function. Function to calculate modularity (default: \code{cluster_infomap}).
#' @param shuffle_func Function. Function to perform network shuffling (default: \code{shuffle_network_deg}).
#'
#' @return A list containing:
#' \describe{
#'   \item{New_A}{Shuffled adjacency matrix}
#'   \item{Metrics}{Data frame tracking metric evolution across iterations, including the original network}
#' }
#'
#' @import igraph
#' @export
generate_shuffled_seq_tol <- function(input_graph, delta = 10, max_iterations = 100,
                                    tolerance = 1e-3, directed = TRUE, weighted=TRUE,
                                    modularity_func = cluster_infomap,
                                    shuffle_func = shuffle_network_deg) {

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
    A <- shuffle_func(A, delta = delta, directed = directed, weighted = weighted)
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
    SVD_Entropy = entropy_values,
    Modularity = modularity_values
  )

  return(list(New_A = A, Metrics = metrics))
}



#' Degree-Preserving Network Shuffling via Edge Swaps
#'
#' This function shuffles a directed network while preserving in-degree and out-degree.
#' It follows a controlled randomization process by swapping edges iteratively.
#' Similar to the `curve_ball()` function but more computationally demanding.
#' This function is useful for generating increasingly randomized networks,
#' to generate fully randomized networks, the `curbe_ball()` function is preferred.
#'
#' @param A A square adjacency matrix (directed, weighted or unweighted).
#' @param delta Number of edge swaps to perform.
#' @param max_attempts Number of times to attempt a valid swap before stopping.
#'
#' @return A shuffled adjacency matrix preserving in-degree and out-degree.
#'
#' @references
#'
#' Huaylla, C. A., Nacif, M. E., Coulin, C., Kuperman, M. N., & Garibaldi, L. A. (2021).
#' Decoding information in multilayer ecological networks: The keystone species case.
#' Ecological Modelling, 460, 109734. https://doi.org/10.1016/j.ecolmodel.2021.109734
#'
#' Strona, G., Nappo, D., Boccacci, F., Fattorini, S., & San-Miguel-Ayanz, J. (2014).
#' A fast and unbiased procedure to randomize ecological binary matrices with fixed row and column totals.
#' Nature Communications, 5, 4114. https://doi.org/10.1038/ncomms5114
#'
#' @import igraph
#' @export
shuffle_network_deg <- function(input_graph, delta = 100, directed = TRUE,weighted=TRUE) {
  if(weighted) {
    w_attr <- "weight"
  } else {
    w_attr <- NULL
  }

  # Check for weighted graphs
  if (weighted && inherits(input_graph, "igraph") && is.null(edge_attr(input_graph, "weight"))) {
    stop("Weighted graph must have 'weight' attribute on edges.")
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


#' Generate a Sequence of Shuffled Networks and Track Metric Evolution
#'
#' This function generates a sequence of shuffled networks from an original graph, applying
#' incremental shuffling steps while tracking modularity, SVD rank and entropy .
#'
#' @param original_graph An igraph object representing the original network.
#' @param max_delta Integer. Number of total shuffling steps to perform.
#' @param delta Integer. Number of swaps per shuffling step.
#' @param directed Logical. Whether the network is directed (default: TRUE).
#' @param weighted Logical. Whether the network is weighted (default: TRUE).
#' @param shuffle_func Function. A network shuffling function (default: `shuffle_network_ws`).
#' @param modularity_func Function. A modularity calculation function (default: `cluster_infomap`).
#'
#' @return A list with:
#' \item{Networks}{A list of igraph objects representing the shuffled networks at each step.}
#' \item{Metrics}{A tibble containing the step number, SVD entropy, SVD Rank, and modularity score.}
#'
#' @import igraph dplyr
#' @export
generate_shuffled_seq <- function(original_graph, max_delta = 10, delta = 10,
                                  directed = TRUE, weighted = TRUE,
                                  shuffle_func = shuffle_network_ws,
                                  modularity_func = cluster_infomap) {
  # Initialize
  shuffled_networks <- list()
  results <- tibble()

  A_current <- as.matrix(as_adjacency_matrix(original_graph, sparse = FALSE))  # Start from original
  shuffled_networks[[1]] <- original_graph

# Calculate metrics
  svd <- calc_svd_entropy(A_current)
  entropy <- svd$Entropy
  rank <- svd$Rank
  modularity_score <- modularity(modularity_func(original_graph))

  results <- bind_rows(results, tibble(Step = 1, SVD_Entropy = entropy, SVD_Rank = rank,
                                       Modularity = modularity_score))

  # Generate sequence of networks
  for (d in 1:max_delta) {
    A_shuffled <- shuffle_func(A_current, delta = delta, directed = directed, weighted = weighted)  # Apply shuffle
    g_shuffled <- graph_from_adjacency_matrix(A_shuffled, mode = "directed", weighted = weighted)

    # Store network
    shuffled_networks[[d + 1]] <- g_shuffled

    # Calculate metrics

    svd <- calc_svd_entropy(A_shuffled)
    entropy <- svd$Entropy
    rank <- svd$Rank

    modularity_score <- modularity(modularity_func(g_shuffled))

    results <- bind_rows(results, tibble(Step = d+1, SVD_Entropy = entropy, SVD_Rank = rank,
                                         Modularity = modularity_score))

    # Update A_current for next iteration
    A_current <- A_shuffled
  }

  return(list(Networks = shuffled_networks, Metrics = results))
}

#' Random Network Rewiring Without Preserving Degree Distribution
#'
#' This function randomly rewires a directed network while preserving
#' the total number of links but NOT the degree distribution, this procedure is adapted from
#' Huaylla et al. (2024) and is based on an approach described by Watts and Strogatz (1998) for small-world networks.
#'
#' The procedure is: we select a single pair of connected nodes and disconnect one of the nodes at one end of the link.
#' We then connect the free end to another node that is randomly selected.
#'
#' @param input_graph A square adjacency matrix (directed, weighted or unweighted) or an igraph object.
#' @param delta Number of rewiring attempts to perform.
#' @param directed Logical, whether the network is directed (default: TRUE).
#' @param weighted Logical, whether the network is weighted (default: TRUE).
#'
#' @return A shuffled adjacency matrix preserving the number of links but not degree distribution.
#'
#' @references
#' Watts, D. J., & Strogatz, S. H. (1998). Collective dynamics of 'small-world' networks.
#' Nature, 393(6684), 440-442. \doi{10.1038/30918}
#'
#' Huaylla, C. A., Kuperman, M. N., & Garibaldi, L. A. (2024). Comparison of two statistical
#' measures of complexity applied to ecological bipartite networks.
#' Physica A: Statistical Mechanics and Its Applications, 642, 129764.
#' https://doi.org/10.1016/j.physa.2024.129764
#'
#' @import igraph
#' @export
#'
shuffle_network_ws <- function(input_graph, delta = 100, directed = TRUE, weighted = TRUE) {
  if (weighted) {
    w_attr <- "weight"
  } else {
    w_attr <- NULL
  }

  # Check for weighted graphs
  if (weighted && inherits(input_graph, "igraph") && is.null(edge_attr(input_graph, "weight"))) {
    stop("Weighted graph must have 'weight' attribute on edges.")
  }

  # Convert igraph object to adjacency matrix if needed
  if (inherits(input_graph, "igraph")) {
    A_orig <- as.matrix(as_adjacency_matrix(input_graph, attr = w_attr, sparse = FALSE))
  } else {
    A_orig <- input_graph
  }

  mode <- ifelse(directed, "directed", "undirected")

  repeat {
    A <- A_orig  # Start from original each time

    # Step 1: Get positions of existing links (1s) and absent links (0s)
    ones_list <- which(A > 0, arr.ind = TRUE)
    zeros_list <- which(A == 0, arr.ind = TRUE)

    for (iter in seq_len(delta)) {
      if (nrow(ones_list) == 0 || nrow(zeros_list) == 0) break

      # Step 3: Choose a random existing edge
      selected_edge <- sample(nrow(ones_list), 1)
      i <- ones_list[selected_edge, 1]
      j <- ones_list[selected_edge, 2]

      # Step 4: Check if the node has more than one neighbor
      if (sum(A[i, ] > 0) > 1) {

        # Step 5: Choose a random absent link (to replace the existing one)
        selected_zero <- sample(nrow(zeros_list), 1)
        i_new <- zeros_list[selected_zero, 1]
        j_new <- zeros_list[selected_zero, 2]


        # Step 6: Perform the swap in the adjacency matrix
        A[i_new, j_new] <- A[i, j]  # Move weight (if weighted) or assign 1 (if unweighted)
        A[i, j] <- 0  # Remove the original link

        # Update ones_list and zeros_list accordingly
        ones_list <- which(A > 0, arr.ind = TRUE)
        zeros_list <- which(A == 0, arr.ind = TRUE)
      }
    }

    # Check connectivity
    g <- graph_from_adjacency_matrix(A, mode = mode, weighted = weighted)

    if (components(g, mode = "weak")$no == 1) {
      return(A)
    }

    # else: repeat again
  }
}
