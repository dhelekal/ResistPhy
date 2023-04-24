#' Simulate coalescent process times using Poisson Process Thinning
#' Note: Ne function is interpolated as a step function
#' @param samp_t a vector of sampling times
#' @param n_samp a vector containing the numbers of samples taken at corresponding sampling time.
#' @param Ne_vals a vector of Ne(t) values.
#' @param Ne_times a vector of time points corresponding to Ne(t) values.
#' @export
sim_coal <- function(samp_t, n_samp, Ne_vals, Ne_times) {
    n <- length(Ne_times)
    samp_total <- sum(n_samp)
    dt <- Ne_times[2]-Ne_times[1]
    stopifnot("Ne_vals and Ne_times must have same length"=length(Ne_vals)==n)
    coal_times <- sim_coal_fast(samp_t, n_samp, Ne_vals, Ne_times, samp_total, dt, n)
    return(list(samp_times=samp_t, n_samp=n_samp, coal_times=coal_times))
}

#' Generate a Newick string representation of a tree corresponding to a realisation of a coalescent process.
#' Note: the topology is RANDOM, If you wish for the topology to be consistent, use set.seed before calling this function.
#' @param samp_t sampling times.
#' @param n_samp a vector containing the numbers of samples taken at corresponding sampling time.
#' @param coal_t times of coalescent events / internal nodes.
#' @param node_name_prefix Prefix given to node names.
#' @return an ape phylogeny with appropriate times
#' @export
build_coal_tree <- function(samp_t, n_samp, coal_t, node_name_prefix="N") {
    n_leaves <- sum(n_samp)
    samp_t_expanded <- c()

    for(i in 1:length(samp_t)) {
        samp_t_expanded <- c(samp_t_expanded, rep(samp_t[i], n_samp[i]))
    }

    samp_t_expanded <- samp_t_expanded

    tree_nodes <- seq(1, n_leaves)
    tree_nodes <- sapply(tree_nodes, function (x) return (paste0("S", x)))

    extant_entries <- c(1)
    extant_times <- c(samp_t_expanded[1])
    
    coal_idx <- 1
    s_idx <- 1
    t <- samp_t_expanded[1]
    if (length(coal_t) > 0) {
        while (coal_idx <= length(coal_t)) {
            if (s_idx < length(samp_t_expanded) && (samp_t_expanded[s_idx + 1] < coal_t[coal_idx])) {
                s_idx <- s_idx + 1
                t <- samp_t_expanded[s_idx]
                extant_entries <- c(extant_entries, s_idx)
                extant_times <- c(extant_times, t)
            } else {
                t <- coal_t[coal_idx]
                coal_node1_idx <- trunc(runif(1, 1, length(extant_entries) + 1))
                
                coal_node1 <- extant_entries[coal_node1_idx]
                ct1 <- extant_times[coal_node1_idx]
                
                extant_times <- extant_times[-coal_node1_idx]
                extant_entries <- extant_entries[-coal_node1_idx]
                
                br_len_1 <- t - ct1
                entry1 <- tree_nodes[coal_node1]
                
                coal_node2_idx <- trunc(runif(1, 1, length(extant_entries) + 1))
                coal_node2 <- extant_entries[coal_node2_idx]
                ct2 <- extant_times[coal_node2_idx]
                br_len_2 <- t - ct2
                entry2 <- tree_nodes[coal_node2]
                
                tree_nodes[coal_node2] <- paste0("(", entry1,
                                                ":", br_len_1,
                                                ",",
                                                entry2,
                                                ":", br_len_2, ")",node_name_prefix,coal_idx)
                extant_times[coal_node2_idx] <- t
                coal_idx <- coal_idx + 1
            }
        }
    }
    tree_str <- paste0(tree_nodes[extant_entries[1]], ";")
    return(read.tree(text=tree_str))
} 