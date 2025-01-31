library(NMF)

plapply <- function(..., progress=FALSE, n.cores=parallel::detectCores(), mc.preschedule=FALSE, mc.allow.recursive=TRUE, fail.on.error=FALSE) {
  if (progress) {
    result <- pbmcapply::pbmclapply(..., mc.cores=n.cores, mc.preschedule=mc.preschedule, mc.allow.recursive=mc.allow.recursive)
    iter <- list(...)[[1]]
    if (!is.null(names(result)) && ("value" %in% names(result)) &&
       (!("value" %in% names(iter)) || ((length(result) != length(iter)) && !("value" %in% names(result$value))))) {
      # In some cases, pbmcapply changes the list structure, storing results in the result$value and warning in result$warning.
      # There is no good way to check if it was changed, though.
      if (!is.null(result$warning)) {
        warning(result$warning)
      }

      result <- result$value
    }
  } else if(n.cores > 1) {
    result <- parallel::mclapply(..., mc.cores=n.cores, mc.preschedule=mc.preschedule, mc.allow.recursive=mc.allow.recursive)
  } else {
    result <- lapply(...)
  }

  if(fail.on.error) {
    is.error <- (sapply(result, class) == "try-error")
    if (any(is.error)) {
      stop(paste("Errors in plapply:", result[is.error]))
    }
  }

  return(result)
}


knn_adjacency_matrix = function(df, k, weighted=FALSE) {
    nns = df[,c("x", "y", "z")] %>% FNN::get.knn(k=k)
    xn = nns$nn.index

    i = rep(1:nrow(xn), each=ncol(xn))
    j = t(xn) %>% as.vector()
    x = if(weighted) {t(nns$nn.dist) %>% as.vector()} else rep(1, length(i))
    adj_mat = sparseMatrix(i=i, j=j, x=x, dims=c(nrow(xn),nrow(xn)))
    return(adj_mat)
}

knn_adjacency_matrix_xenium = function(df, k, weighted=FALSE) {
    nns = df[,c("x_location", "y_location", "z_location")] %>% FNN::get.knn(k=k)
    xn = nns$nn.index

    i = rep(1:nrow(xn), each=ncol(xn))
    j = t(xn) %>% as.vector()
    x = if(weighted) {t(nns$nn.dist) %>% as.vector()} else rep(1, length(i))
    adj_mat = sparseMatrix(i=i, j=j, x=x, dims=c(nrow(xn),nrow(xn)))
    return(adj_mat)
}


#' This is a much faster implementation than `knn_count_mat`
knn_count_matrix = function(df, k, all.genes=NULL) {
    adj_mat = knn_adjacency_matrix(df, k)
    j = if (is.null(all.genes)) factor(df$gene) else factor(df$gene, levels=all.genes)
    ncm = adj_mat %*% sparseMatrix(
        i=1:length(j), j=as.integer(j), x=rep(1, length(j)),
        dims=c(nrow(adj_mat), length(levels(j)))
    )
    colnames(ncm) = levels(j)
    rownames(ncm) = rownames(df)

    return(ncm)
}

knn_count_matrix_xenium = function(df, k, all.genes=NULL) {
    adj_mat = knn_adjacency_matrix_xenium(df, k)
    j = if (is.null(all.genes)) factor(df$feature_name) else factor(df$feature_name, levels=all.genes)
    ncm = adj_mat %*% sparseMatrix(
        i=1:length(j), j=as.integer(j), x=rep(1, length(j)),
        dims=c(nrow(adj_mat), length(levels(j)))
    )
    colnames(ncm) = levels(j)
    rownames(ncm) = rownames(df)

    return(ncm)
}

#' Get KNN counts for each molecule
#' @param df Data frame with molecule coordinates
#' @param k Neighborhood size
#' @param ncores Number of cores to use
#' @return Data frame with KNN counts for each molecule
get_knn_counts = function(df, k, include_i = FALSE, ncores = 1) {
    k = as.integer(k)
    cells = df %>% count(cell) %>% filter(n > k) %>% pull(cell)
    genes = df$gene %>% unique() %>% sort()

    N = df %>%
        filter(cell %in% cells) %>%
        split(.$cell) %>%
        lapply(knn_count_matrix, k=k, all.genes=genes) %>%
        do.call(rbind, .)

    # N[is.na(N)] = 0

    # N = as.data.frame(N)

    # N %<>% tibble::column_to_rownames('id') %>% as.matrix()# %>% .[,2:ncol(.)]

    return(N)
}

get_knn_counts_xenium = function(df, k, include_i = FALSE, ncores = 1) {
    k = as.integer(k)
    cells = df %>% count(cell_id) %>% filter(n > k) %>% pull(cell_id)
    genes = df$feature_name %>% unique() %>% sort()

    N = df %>%
        filter(cell_id %in% cells) %>%
        split(.$cell_id) %>%
        lapply(knn_count_matrix_xenium, k=k, all.genes=genes) %>%
        do.call(rbind, .)

    # N[is.na(N)] = 0

    # N = as.data.frame(N)

    # N %<>% tibble::column_to_rownames('id') %>% as.matrix()# %>% .[,2:ncol(.)]

    return(N)
}

##### Factorization #####

run_weighted_nmf = function(X, k, n.downdonor = 4000) {
    # downsampling X for for faster runtime for the demonstration
    if (n.downdonor > 0) {
        rows_select <- donor(1:nrow(X), n.downdonor, replace=FALSE)
        X <- X[rows_select,]
    }

    X <- X[rowSums(X) > 0,] # we can't remove genes, as we need them all labeled

    # computing the weights matrix to use in weighted NMF
    Z <- matrix(rep(1/colSums(X), nrow(X)), nrow = nrow(X), byrow = TRUE)
    Z[,colSums(X) == 0] <- 1e-4

    res <- NMF::nmf(X, rank = k, method = 'ls-nmf', weight = Z, nrun = 30, .opt='vmp5', seed=0)

    return(res)
}

#### MRF segmentation #####

#' Get graph from cell
#' @param df_cell Data frame with molecule coordinates
#' @param h Neighborhood size
#' @param distinct_edges Whether to keep only distinct edges
#' @return Graph object
get_graph = function(df_cell, h = 10, distinct_edges = FALSE) {

    h = as.integer(h)

    if (h >= nrow(df_cell)) {
        h = min(h, nrow(df_cell) - 1)
        message('h is too large, setting to ', h, '...')
    }

    if (nrow(df_cell) <= h) {
        return(NULL)
    }

    xn = df_cell %>% select(x, y, z) %>% FNN::get.knn(k = h) %>% .$nn.index

    gene_dict = df_cell %>% mutate(id = 1:n()) %>% {setNames(.$gene, .$mol_id)}
    id_dict = df_cell %>% mutate(id = 1:n()) %>% {setNames(.$mol_id, .$id)}

    xn = apply(xn, 2, function(i){id_dict[i]})
    rownames(xn) = df_cell$mol_id

    edges = xn %>% as.data.frame() %>%
        tibble::rownames_to_column('x') %>%
        reshape2::melt(id.var = 'x', value.name = 'y') %>%
        select(x, y)

    G = graph_from_data_frame(edges, directed = FALSE)

    # only keep distinct edges
    if (distinct_edges) {
        G = simplify(G)
    }

    V(G)$gene = gene_dict[names(V(G))]

    V(G)$id = names(V(G))

    # V(G)$centrality = estimate_betweenness(G, directed = F, cutoff = -1)

    return(G)
}


get_graph_xenium = function(df_cell, h = 10, distinct_edges = FALSE) {

    h = as.integer(h)

    if (h >= nrow(df_cell)) {
        h = min(h, nrow(df_cell) - 1)
        message('h is too large, setting to ', h, '...')
    }

    if (nrow(df_cell) <= h) {
        return(NULL)
    }

    xn = df_cell %>% select(x_location, y_location, z_location) %>% FNN::get.knn(k = h) %>% .$nn.index

    gene_dict = df_cell %>% mutate(id = 1:n()) %>% {setNames(.$feature_name, .$transcript_id)}
    id_dict = df_cell %>% mutate(id = 1:n()) %>% {setNames(.$transcript_id, .$id)}

    xn = apply(xn, 2, function(i){id_dict[i]})
    rownames(xn) = df_cell$transcript_id

    edges = xn %>% as.data.frame() %>%
        tibble::rownames_to_column('x_location') %>%
        reshape2::melt(id.var = 'x_location', value.name = 'y_location') %>%
        select(x_location, y_location)

    G = graph_from_data_frame(edges, directed = FALSE)

    # only keep distinct edges
    if (distinct_edges) {
        G = simplify(G)
    }

    V(G)$gene = gene_dict[names(V(G))]

    V(G)$id = names(V(G))

    # V(G)$centrality = estimate_betweenness(G, directed = F, cutoff = -1)

    return(G)
}

#' Run MRF
#' @param df_cell Data frame with molecule coordinates
#' @param expr Data frame with gene expression
#' @param k Number of factors
#' @param h Neighborhood size
#' @param t Transition probability
#' @return Data frame with MRF segmentation
#' @export
run_crf = function(df_cell, expr, h = 10, t = 0.05, verbose=TRUE, ...) {
    G = df_cell %>% get_graph(h = h)

    expr <- expr[unique(df_cell$gene),]
    non_exp_mask <- (rowSums(expr) < 1e-10)
    if (any(non_exp_mask)) {
        warning('Expression matrix contains zero rows for genes: ', paste(rownames(expr)[non_exp_mask], collapse=', '))
        expr[non_exp_mask,] <- 1 / nrow(expr)
    }

    V = igraph::as_data_frame(G, 'vertices') %>% cbind(expr[.$gene,])

    adj = as_adjacency_matrix(G)

    k = ncol(expr)
    crf = CRF::make.crf(adj, k)

    crf$node.pot = V %>%
        select(colnames(expr)) %>%
        as.matrix()

    A = matrix(rep(t, k^2), nrow = k)
    diag(A) = 1-(k-1)*t

    for (i in 1:length(crf$edge.pot)) {
        crf$edge.pot[[i]] = A
    }

    out <- CRF::decode.lbp(crf, verbose=verbose, ...)

    df_out = data.frame(id = V$name, factor = out) %>%
        mutate(cell = unique(df_cell$cell))

    return(df_out)
}

run_crf_per_cell <- function(df_ct, H, h=5, ncores=1) {
    crf_res = df_ct %>%
      split(.$cell, drop=TRUE) %>%
      .[sapply(., nrow) > 1] %>%
      sccore::plapply(run_crf, H, h=h, mc.preschedule=TRUE, fail.on.error=TRUE, n.cores=ncores, progress=TRUE) %>%
    #   lapply(run_crf, H, h=h) %>%
      bind_rows()

    rownames(df_ct) <- df_ct$mol_id
    crf_res <- cbind.data.frame(df_ct[crf_res$id,],crf_res)
    crf_res <- crf_res[,1:(ncol(crf_res)-1)]
    return(crf_res)
}



run_crf_xenium = function(df_cell, expr, h = 10, t = 0.05, verbose=TRUE, ...) {
    G = df_cell %>% get_graph(h = h)

    expr <- expr[unique(df_cell$feature_name),]
    non_exp_mask <- (rowSums(expr) < 1e-10)
    if (any(non_exp_mask)) {
        warning('Expression matrix contains zero rows for genes: ', paste(rownames(expr)[non_exp_mask], collapse=', '))
        expr[non_exp_mask,] <- 1 / nrow(expr)
    }

    V = igraph::as_data_frame(G, 'vertices') %>% cbind(expr[.$feature_name,])

    adj = as_adjacency_matrix(G)

    k = ncol(expr)
    crf = CRF::make.crf(adj, k)

    crf$node.pot = V %>%
        select(colnames(expr)) %>%
        as.matrix()

    A = matrix(rep(t, k^2), nrow = k)
    diag(A) = 1-(k-1)*t

    for (i in 1:length(crf$edge.pot)) {
        crf$edge.pot[[i]] = A
    }

    out <- CRF::decode.lbp(crf, verbose=verbose, ...)

    df_out = data.frame(id = V$name, factor = out) %>%
        mutate(cell = unique(df_cell$cell_id))

    return(df_out)
}

run_crf_per_cell_xenium <- function(df_ct, H, h=5, ncores=1) {
    crf_res = df_ct %>%
      split(.$cell_id, drop=TRUE) %>%
      .[sapply(., nrow) > 1] %>%
      sccore::plapply(run_crf_xenium, H, h=h, mc.preschedule=TRUE, fail.on.error=TRUE, n.cores=ncores, progress=TRUE) %>%
    #   lapply(run_crf, H, h=h) %>%
      bind_rows()

    rownames(df_ct) <- df_ct$transcript_id
    crf_res <- cbind.data.frame(df_ct[crf_res$id,],crf_res)
    crf_res <- crf_res[,1:(ncol(crf_res)-1)]
    return(crf_res)
}