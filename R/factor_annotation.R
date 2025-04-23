## bridge test

get_adj_pairs <- function(df_sub,cell.annot,ncm.ct,ncm_cell,ct1,ct2) {
  # get molecules belonging to first cell type that are in cells_into vec
  ct1_mol_id <- df_sub$mol_id[df_sub$celltype==ct1]
  ncm_ct1 <- ncm.ct[ct1_mol_id,]
  # get only the molecules that are adjacent to the other cell type
  mol_near_ct2 <- rownames(ncm_ct1)[ncm_ct1[,ct2]!=0]
  if (length(mol_near_ct2)<5) {
    return(NULL)
  }

  ct_labs <- cell.annot[colnames(ncm_cell),'celltype']
  ncm_cell_sub <- ncm_cell[mol_near_ct2,which(ct_labs==ct2),drop=FALSE]
  ncm_cell_sub <- ncm_cell_sub>0
  if (ct1==ct2) { # need to make sure we exclude self cells from adjacent list
    # further limiting to only molecules that are linked to multiple cells
    ndx_keep <- which(Matrix::rowSums(ncm_cell_sub>0)>1)
    if (length(ndx_keep)==0) {
      return(NULL)
    }
    ncm_cell_sub <- ncm_cell_sub[ndx_keep,]
    self_cells <- df_sub[rownames(ncm_cell_sub),'cell']
    cn_match <- match(self_cells,colnames(ncm_cell_sub))
    non_zero_entries <- Matrix::summary(ncm_cell_sub)

    # Split the data by row index to find non-zero columns for each row
    non_zero_cols <- split(non_zero_entries$j, non_zero_entries$i)
    non_zero_cols <- lapply(1:length(non_zero_cols), function(r_ndx) {
      c_ndx <- non_zero_cols[[r_ndx]]
      c_ndx <- c_ndx[!(c_ndx %in% cn_match[r_ndx])]
      return(colnames(ncm_cell_sub)[c_ndx][1])
    })
    adj_cells <- unlist(non_zero_cols)
  } else {
    non_zero_entries <- Matrix::summary(ncm_cell_sub)
    non_zero_cols <- split(non_zero_entries$j, non_zero_entries$i)
    non_zero_cols <- lapply(non_zero_cols, function(cols) {
      return(colnames(ncm_cell_sub)[cols][1]) # just return the first adjacent cell if there are several
    })
    adj_cells <- unlist(non_zero_cols)
    self_cells <- df_sub[rownames(ncm_cell_sub),'cell']
  }

  # only keep the first instance of duplicated cells
  cell_pair_df <- cbind.data.frame(self_cells,adj_cells)
  ndx_keep <- which(!duplicated(cell_pair_df[,'self_cells']))
  cell_pair_df <- cell_pair_df[ndx_keep,]
  ndx_keep <- which(!duplicated(cell_pair_df[,'adj_cells']))
  cell_pair_df <- cell_pair_df[ndx_keep,]

  if (ct1==ct2) {
    sub1 <- cell_pair_df[,1]
    sub2 <- cell_pair_df[,2]

    ndx_keep <- which(!(sub2 %in% sub1))
    cell_pair_df <- cell_pair_df[ndx_keep,]
  }

  return(cell_pair_df)
}

get_mol_link <- function(df_cell1,df_cell2,knn_k=10) {
  df_tmp <- rbind(df_cell2,df_cell1)
  if (nrow(df_tmp)<=knn_k) {
    knn_k <- nrow(df_tmp)-1
  }
  # Nearest neighbors for points in cluster_b
  nn_result <- RANN::nn2(data = df_tmp, query = df_cell1, k = knn_k)

  index_switch <- nrow(df_cell2)

  # Check if any nearest neighbor belongs to cluster_a
  bridging_count <- sum(apply(nn_result$nn.idx, 1, function(indices) {
    any(indices<=index_switch)  # 1 means it belongs to cluster_a
  }))

  return(bridging_count)
}

prune_cells_test <- function(cell_pair_df,ncells_samp=200) {
  f_counts <- cell_pair_df[,3:ncol(cell_pair_df)]
  all_cells_test <- c()
  for (i in 1:ncol(f_counts)) {
    ndx_keep <- which(f_counts[,i]>=5)
    cells <- rownames(f_counts)[ndx_keep]
    already_included <- cells[cells %in% all_cells_test]
    cells <- cells[!(cells %in% already_included)]
    num_included <- length(already_included)
    if (num_included>ncells_samp) {
      cells_samp <- c()
    } else if ((length(cells)+num_included)>ncells_samp) {
      num_needed <- ncells_samp-num_included
      cells_samp <- sample(cells,num_needed)
    } else if (length(ndx_keep)<5) {
      cells_samp <- c()
    } else {
      cells_samp <- cells
    }
    all_cells_test <- c(all_cells_test,cells_samp)
  }
  all_cells_test <- unique(all_cells_test)
  cell_pair_df_pruned <- cell_pair_df[all_cells_test,]
  return(cell_pair_df_pruned)
}

move_cells_adj <- function(df_cell1,df_cell2,num_cross_from_target,
                           knn_k=10,step = 0.02,max_iter = 100) {

  c1_mean1 <- mean(df_cell1[,'x'])
  c1_mean2 <- mean(df_cell1[,'y'])
  r1 <- (max(df_cell1[,'x'])-min(df_cell1[,'x']))/2
  r2 <- (max(df_cell1[,'y'])-min(df_cell1[,'y']))/2
  max_c1_r <- max(r1,r2)

  c2_mean1 <- mean(df_cell2[,'x'])
  c2_mean2 <- mean(df_cell2[,'y'])
  r1 <- (max(df_cell2[,'x'])-min(df_cell2[,'x']))/2
  r2 <- (max(df_cell2[,'y'])-min(df_cell2[,'y']))/2
  max_c2_r <- max(r1,r2)

  r <- max_c1_r + max_c2_r  # radius of the circle

  # Generate a random angle between 0 and 2*pi
  theta <- runif(1, min = 0, max = 2 * pi)

  # Calculate the x and y coordinates for the point
  x2 <- c1_mean1 + r * cos(theta)
  y2 <- c1_mean2 + r * sin(theta)

  # Calculate the shift needed
  shift_x <- x2 - c2_mean1
  shift_y <- y2 - c2_mean2

  # Apply the shift to all points
  df_cell2$x <- df_cell2$x + shift_x
  df_cell2$y <- df_cell2$y + shift_y

  for (i in 1:max_iter) {
    # for (i in 1:10) {
    # Current bridging counts
    count13 <- get_mol_link(df_cell1, df_cell2, knn_k)

    # If Cluster 3 has reached the same bridging count as Cluster 2, stop
    if (count13 >= num_cross_from_target) break

    # Direction to move Cluster 3 closer to Cluster 1
    direction <- Matrix::colMeans(df_cell1) - Matrix::colMeans(df_cell2)
    myshift <- step * direction
    df_cell2 <- sweep(df_cell2,MARGIN = 2,STATS = myshift,FUN = '+')
  }

  if (i==max_iter) {
    return(NULL)
  } else {
    return(df_cell2)
  }
}


compute_bridge_scores <- function(df_crf_diff_fracs,factor_test,df_crf,cell_pair_df,knn_k=10,abs_diff=TRUE) {
  all_pair_scores <- lapply(1:nrow(cell_pair_df),function(pair_ndx) {
    # print(pair_ndx)
    cell1 <- cell_pair_df[pair_ndx,1]
    cell2 <- cell_pair_df[pair_ndx,2]

    # first compute difference in expression of the factor between cells
    fracs_c1 <- df_crf_diff_fracs[cell1,factor_test]
    fracs_c2 <- df_crf_diff_fracs[cell2,factor_test]
    if (abs_diff) {
      fracs_diff <- abs(fracs_c2 - fracs_c1)
    } else {
      fracs_diff <- fracs_c2 - fracs_c1
    }

    ### get frac of f mols in cell1 that cross to cell2
    df_cell1 <- df_crf[df_crf$cell %in% c(cell1),c('x','y','z','factor')]
    f_labs <- which(df_cell1$factor==factor_test)
    df_cell1$factor <- NULL
    df_cell2 <- df_crf[df_crf$cell %in% c(cell2),c('x','y','z')]
    df_tmp <- rbind(df_cell1,df_cell2)
    if (nrow(df_tmp)<=knn_k) {
      tmp_k <- nrow(df_tmp)-1
      nn_result <- RANN::nn2(data = df_tmp, query = df_cell1, k = tmp_k)
    } else {
      nn_result <- RANN::nn2(data = df_tmp, query = df_cell1, k = knn_k)
    }
    index_switch <- nrow(df_cell1)

    bridging_count <- sum(apply(nn_result$nn.idx[f_labs,], 1, function(indices) {
      any(indices>index_switch)
    }))
    frac_cross <- bridging_count / length(f_labs)

    pair_score <- fracs_diff * frac_cross
    return(pair_score)
  })
  all_pair_scores <- unlist(all_pair_scores)
  names(all_pair_scores) <- cell_pair_df[,1]

  return(all_pair_scores)
}

compute_bridge_scores_null <- function(df_crf_diff_fracs,cell_pair_df,df_crf,knn_k,ncells_samp=200) {
  f_dist_all <- as.matrix(dist(df_crf_diff_fracs[cell_pair_df[,2],]))

  cell_vec <- df_crf$cell
  df_cells <- df_crf %>%
    select(x,y,z) %>%
    split(cell_vec, drop = T)

  df_cells_full <- df_crf %>%
    split(cell_vec, drop = T)

  cells_match <- c()
  for (i in 1:nrow(cell_pair_df)) {
    cell2_orig <- cell_pair_df[i,2]
    c_dist_pairs <- f_dist_all[cell2_orig,]
    c_dist_pairs <- c_dist_pairs[order(c_dist_pairs,decreasing=FALSE)]
    c_dist_pairs <- c_dist_pairs[2:length(c_dist_pairs)]
    c_dist_pairs_top <- names(c_dist_pairs)[1:ceiling((length(c_dist_pairs)*.1)+1)]
    cell_select <- sample(c_dist_pairs_top,1)

    cells_match <- c(cells_match,cell_select)
  }
  cell_pair_samp <- cbind.data.frame(cell_pair_df[,1],cells_match)
  # print('cell match complete')

  colnames(cell_pair_samp) <- c('self_cells','adj_cells')
  perm_list1 <- list()
  perm_list2 <- list()
  pairs_keep <- c()
  for (i in 1:nrow(cell_pair_samp)) {
    cell1 <- cell_pair_samp[i,1]
    cell2 <- cell_pair_samp[i,2]

    # print(paste0(cell1,' ',cell2))

    df_cell1 <- df_cells[[cell1]]
    df_cell2 <- df_cells[[cell2]]

    # determine number of molecules in "admix into" cell cross over to it's real adj pair
    cell2_real <- cell_pair_df[i,2]
    df_cell2_real <- df_cells[[cell2_real]]
    num_cross_from_target <- get_mol_link(df_cell1,df_cell2_real,knn_k=knn_k)

    # now move cell1 and cell2 close enough together to have roughly same number of linked molecules
    df_cell2_moved <- move_cells_adj(df_cell1,df_cell2,num_cross_from_target,knn_k=knn_k)

    if (!is.null(df_cell2_moved)) {
      pairs_keep <- c(pairs_keep,i)
    }

    df_cell1$cell <- cell1
    df_cell2_moved$cell <- cell2
    df_cell1$factor <- df_cells_full[[cell1]]$factor
    df_cell2_moved$factor <- df_cells_full[[cell2]]$factor

    perm_list1[[i]] <- df_cell1
    perm_list2[[i]] <- df_cell2_moved
  }
  # print('cell randomization complete')

  perm_list1 <- perm_list1[pairs_keep]
  perm_list2 <- perm_list2[pairs_keep]
  cell_pair_df <- cell_pair_df[pairs_keep,]
  cell_pair_samp <- cell_pair_samp[pairs_keep,]

  df_crf_perm <- do.call(rbind,c(perm_list1,perm_list2))

  all_factors <- 1:ncol(df_crf_diff_fracs)
  scores_null_all <- c()
  for (factor_test in all_factors) {
    ndx_test <- which(cell_pair_df[,factor_test+2]>=5)
    cell_pair_samp_f <- cell_pair_samp[ndx_test,] # only test cells that have at lest 5 molecule of the factor being tested
    if (nrow(cell_pair_samp_f)<5) { # require that we have at least 5 cells total
      scores_null_all[[factor_test]] <- NA
      next
    } else if (nrow(cell_pair_samp_f)>ncells_samp) { # if there are many pairs, downsample for speed
      cell_pair_samp_f <- cell_pair_samp_f[sample(1:nrow(cell_pair_samp_f),ncells_samp),]
    }
    scores_out_null <- compute_bridge_scores(df_crf_diff_fracs,factor_test,df_crf_perm,cell_pair_samp_f,knn_k=knn_k,
                                             abs_diff = TRUE)
    scores_null_all[[factor_test]] <- scores_out_null
  }

  return(scores_null_all)
}

#' Run cell-bridging test to determine factors that represent admixture into various cell types
#'
#' @param df Data Frame of transcript locations, mol_id, cell, and celltype columns
#' @param crf.f.labs Data Frame of CRF results
#' @param cell.annot Data Frame of spatial metadata
#' @param ncm.ct Matrix of all KNN counts for molecules to cell types
#' @param ncm.cell Matrix of all KNN counts for molecules to cells
#' @param ncells.samp Number of cells to subsample to run the test with for speed (default=200)
#' @param knn.k Number of nearest neighbors to use for KNN (default=20)
#' @param n.cores Number of cores to run parallel with (default=1)
#'
#' @return a list of p-values and top quantiles for each pair of cell types tested
#' @export
run_bridge_test <- function(
    df, crf.f.labs, cell.annot, ncm.ct, ncm.cell, ncells.samp=200, knn.k=20, n.cores=1
  ) {
  df$cell <- as.character(df$cell)

  ct_counts <- table(cell.annot$celltype)
  all.ctypes <- names(ct_counts)[ct_counts>1]
  ct_perms <- gtools::permutations(n = length(all.ctypes), r = 2, v = all.ctypes)
  same_ct_lst <- lapply(all.ctypes, \(x) c(x,x))
  ct_perms <- do.call(rbind,c(same_ct_lst,list(ct_perms)))

  combo.res <- sccore::plapply(1:nrow(ct_perms),function(i) {
    # print(i)
    ct1 <- ct_perms[i,1]
    ct2 <- ct_perms[i,2]
    cell_pair_df <- get_adj_pairs(df,cell.annot,ncm.ct,ncm.cell,ct1=ct1,ct2=ct2)

    if (is.null(cell_pair_df)) {
      return(NA)
    }

    # if crf.f.labs has more than 1 column (if it was run per ct), need to load the correct vector
    if (ncol(crf.f.labs)>1) {
      df$factor <- as.factor(crf.f.labs[,ct1])
    } else {
      df$factor <- as.factor(crf.f.labs[,1])
    }

    ## append factor counts for the admix in cells
    cells_all <- c(unlist(cell_pair_df))
    df_crf <- df[df$cell %in% cells_all,]
    df_crf_diff_tbl <- table(df_crf[,c('cell','factor')])
    df_crf_diff_fracs <- df_crf_diff_tbl / Matrix::rowSums(df_crf_diff_tbl)
    df_crf_diff_tbl <- as.data.frame.matrix(df_crf_diff_tbl)
    cell_pair_df <- cbind.data.frame(cell_pair_df,df_crf_diff_tbl[cell_pair_df$self_cells,])
    rownames(cell_pair_df) <- cell_pair_df$self_cells

    cell_pair_df <- prune_cells_test(cell_pair_df, ncells_samp=ncells.samp)

    if (nrow(cell_pair_df)<5) {
      return(NA)
    }

    # cells_all <- c(cell_pair_df[,1],cell_pair_df[,2])
    # df_crf_diff_fracs <- df_crf_diff_fracs[cells_all,]
    # df_crf <- df_crf[df_crf$cell %in% cells_all,]

    # doing five permuted iterations per cell for more stable results
    scores_null_iter <- lapply(1:3, function(null_iter){
      compute_bridge_scores_null(
        df_crf_diff_fracs, cell_pair_df, df_crf,
        knn_k=knn.k, ncells_samp=ncells.samp
      )
    })

    scores_null_all <- lapply(1:length(scores_null_iter[[1]]),function(f){
      cells_all <- lapply(scores_null_iter, \(nr) names(nr[[f]]))
      cells_int <- cells_all[[1]]
      for (null_iter in 2:length(cells_all)) {
        cells_int <- intersect(cells_int,cells_all[[null_iter]])
      }
      nfres_all <- lapply(scores_null_iter, \(nr) nr[[f]][cells_int])
      nfres_all <- do.call(cbind,nfres_all)
      cmeans <- Matrix::rowMeans(nfres_all)
      return(cmeans)
    })

    all_factors <- 1:ncol(df_crf_diff_fracs)
    factor_stats <- list()
    for (factor_test in all_factors) {
      factor_null_vals <- scores_null_all[[factor_test]]
      if (is.na(factor_null_vals[1])) {
        factor_stats[[factor_test]] <- NA
        next
      }
      cell_pair_df_f <- cell_pair_df[names(factor_null_vals),]
      scores_real <- compute_bridge_scores(
        df_crf_diff_fracs, factor_test, df_crf, cell_pair_df_f,
        knn_k=knn.k, abs_diff=TRUE
      )
      q3_val <- quantile(scores_real,.75)
      pval <- wilcox.test(scores_real,factor_null_vals, paired = TRUE, alternative = "greater")[["p.value"]]

      factor_stats[[factor_test]] <- c(pval,q3_val)
    }

    return(factor_stats)
  },mc.preschedule=FALSE,n.cores=n.cores,progress=TRUE)
  names(combo.res) <- paste0(ct_perms[,1],'_',ct_perms[,2])

  return(combo.res)
}


collapse_row_col <- function(mats,c.type) {
  all_pval_mats_summary <- lapply(mats,function(x){
    if (c.type=='r_max') {
      n_iter <- nrow(x)
    } else if (c.type=='c_mean') {
      n_iter <- ncol(x)
    }
    all_max <- sapply(1:n_iter,function(i) {
      if (c.type=='r_max') {
        xvec <- x[i,]
      } else if (c.type=='c_mean') {
        xvec <- x[,i]
      }
      xvec <- xvec[!is.na(xvec)]

      if (c.type=='r_max') {
        if (length(xvec)>0) {
          return(max(xvec))
        } else {
          return(NA)
        }
      } else if (c.type=='c_mean') {
        if (length(xvec)>0) {
          return(mean(xvec))
        } else {
          return(NA)
        }
      }
    })
    if (c.type=='r_max') {
      names(all_max) <- rownames(x)
    } else if (c.type=='c_mean') {
      names(all_max) <- colnames(x)
    }
    return(all_max)
  })
  all_pval_mats_summary <- do.call(cbind,all_pval_mats_summary)
  return(all_pval_mats_summary)
}

select_bridge_remove <- function(all_pval_mats_rs,all_pval_mats_cs,pthresh) {
  # now getting list of factor-ct combos to remove
  f_ct_rm_bridge <- list()
  f_from_bridge <- list()
  for (f in 1:ncol(all_pval_mats_rs)) {
    f_nm <- paste0('f_',f)
    pval <- all_pval_mats_rs[,f]
    pval <- pval[!is.na(pval)]
    ct_rm <- names(pval)[which(pval>(-log10(pthresh)))]

    pval <- all_pval_mats_cs[,f]
    pval <- pval[!is.na(pval)]
    if (length(pval)==0) {
      next
    }
    ct_from <- names(pval)[which(pval==max(pval,na.rm=TRUE))]

    ct_rm <- ct_rm[!(ct_rm %in% ct_from)]
    if (length(ct_rm)>0 & !is.na(ct_rm[1])) {
      f_ct_rm_bridge[[f_nm]] <- ct_rm
    }
    if (length(ct_from)>0 & !is.na(ct_rm[1]) & length(ct_rm)>0) {
      f_from_bridge[[f_nm]] <- ct_from
    }
  }
  return(list(f_ct_rm_bridge,f_from_bridge))
}

extract_name_combos <- function(combo.res, all.ctypes) {
  ct.combos <- names(combo.res)
  escaped.uniques <- all.ctypes %>%
    stringr::str_replace_all("([\\^\\$\\.\\|\\?\\*\\+\\(\\)\\[\\{\\\\])", "\\\\\\1") %>%
    paste(collapse = "|")
  pattern <- paste0("^(", escaped.uniques, ")_(", escaped.uniques, ")$")
  matches <- stringr::str_match(ct.combos, pattern)
  ct.combos <- lapply(1:nrow(matches), \(i) c(matches[i,2],matches[i,3]))
  return(ct.combos)
}

#' Extract results from bridge test result into a list of factor-ctype pairs for removal
#'
#' @param combo.res List of results from bridge tests
#' @param all.ctypes Vector of cell types
#' @param n.factors Number of factors tested
#' @param nmf.type Character whether the cell type was run jointly with all cell types
#' together or separately per cell type ('joint' or 'ct')
#' @param p.thresh pvalue threshold below which to select factor-ct combinations for removal
#' @param adj.pvals Whether to adjust pvalues with FDR or not
#' @param get.p.mat Whether to get the p-value matrix instead of the list of pairs
#'
#' @return either a list of factor-ctype pairs or a matrix of pvalues.
#' @export
extract_bridge_res <- function(
    combo.res, all.ctypes, n.factors,
    nmf.type='joint', p.thresh=.01, adj.pvals=TRUE, get.p.mat=FALSE
  ) {
  nmf.type <- match.arg(nmf.type, c('joint', 'ct'))
  all_qsc_mats <- NULL
  all_pval_mats <- NULL

  for (factor_test in 1:n.factors) {
    ct_combos <- extract_name_combos(combo.res, all.ctypes)

    mat_qsc <- matrix(NA,ncol=length(all.ctypes),nrow=length(all.ctypes))
    rownames(mat_qsc) <- colnames(mat_qsc) <- all.ctypes
    mat_pval <- mat_qsc

    # loop through pairs and extract results for the factor
    for (pair_ndx in 1:length(combo.res)) {
      pair_res <- combo.res[[pair_ndx]]
      ct1 <- ct_combos[[pair_ndx]][1]
      ct2 <- ct_combos[[pair_ndx]][2]

      if (length(pair_res)==1) {
        next
      }

      pval <- pair_res[[factor_test]][1]
      qsc <- pair_res[[factor_test]][2]

      mat_qsc[ct1,ct2] <- qsc
      mat_pval[ct1,ct2] <- pval
    }

    all_qsc_mats[[factor_test]] <- mat_qsc

    if (adj.pvals) {
      mat_pval_fdr <- p.adjust(mat_pval,method='fdr')
      mat_pval_fdr <- matrix(mat_pval_fdr,nrow=nrow(mat_pval),ncol=ncol(mat_pval))
      rownames(mat_pval_fdr) <- rownames(mat_pval)
      colnames(mat_pval_fdr) <- colnames(mat_pval)
      mat_pval_fdr <- -log10(mat_pval_fdr)
      all_pval_mats[[factor_test]] <- mat_pval_fdr
    } else {
      mat_pval <- -log10(mat_pval)
      all_pval_mats[[factor_test]] <- mat_pval
    }
  }

  if (nmf.type=='joint') {
    all_pval_mats_rs <- collapse_row_col(all_pval_mats,c.type = 'r_max')
    all_pval_mats_cs <- collapse_row_col(all_pval_mats,c.type = 'c_mean')
    all_score_mats_rs <- collapse_row_col(all_qsc_mats,c.type = 'r_max')
    all_score_mats_cs <- collapse_row_col(all_qsc_mats,c.type = 'c_mean')
    if (get.p.mat) {
      return(all_pval_mats)
    }
    bridge_res <- select_bridge_remove(all_pval_mats_rs,all_pval_mats_cs,p.thresh)

  } else {
    bridge_res <- list()
    for (ct in rownames(all_pval_mats[[1]])) {
      ct_from <- lapply(1:length(all_pval_mats),function(mat_ndx) { # for a ct, get results from all it's factors
        x1 <- all_pval_mats[[mat_ndx]]
        admix_for_ct <- x1[ct,]
        admix_for_ct <- admix_for_ct[!is.na(admix_for_ct)]
        if (length(admix_for_ct)==0) {
          return(NA)
        }

        ct_rm <- names(admix_for_ct)[which(admix_for_ct>(-log10(p.thresh)))]

        if (length(ct_rm)>0) {
          ct_res <- names(admix_for_ct)[admix_for_ct==max(admix_for_ct)]
          return(ct_res)
        } else {
          return(NA)
        }
      })

      names(ct_from) <- paste0('f_',1:length(all_pval_mats))
      ndx_not_na <- which(!is.na(ct_from))
      ct_from <- ct_from[ndx_not_na]
      ct_to <- lapply(ct_from,function(x){
        return(ct)
      })
      bridge_res[[ct]] <- list(ct_to,ct_from)
    }
  }

  return(bridge_res)
}

## membrane test

line_coordinates_interpolated <- function(r0, c0, r1, c1) {
  n <- max(abs(r1 - r0), abs(c1 - c0)) + 1
  rows <- round(seq(r0, r1, length.out = n))
  cols <- round(seq(c0, c1, length.out = n))
  coords <- cbind(rows, cols)
  colnames(coords) <- c("row", "col")
  return(coords)
}


compute_summary <- function(p0, p1, im) {
  p0 <- unlist(p0)
  p1 <- unlist(p1)

  z0 <- p0[3]
  z1 <- p1[3]
  if (z0 != z1) stop("z indices must match")

  p0 <- floor(p0)
  p1 <- floor(p1)
  interpolated_xy <- line_coordinates_interpolated(p0[1], p0[2], p1[1], p1[2])

  vals <- diag(im[z0+1,interpolated_xy[,2]+1,interpolated_xy[,1]+1])
  return(max(vals))
}


#' Normalize images for membrane score tests
#'
#' @param im image matrix with dimensions of z by y by x
#'
#' @return the normalized image
#' @export
normalize_images <- function(im) {

  # 1. Subtract the mean for each channel
  channel_means <- apply(im, 1, mean)  # Mean for each channel
  im_centered <- sweep(im, 1, channel_means, FUN = "-")  # Subtract mean along channels

  # 2. Divide by the standard deviation for each channel
  channel_sds <- apply(im_centered, 1, sd)  # Standard deviation for each channel
  im_normalized <- sweep(im_centered, 1, channel_sds, FUN = "/")  # Normalize each channel

  # 3. Shift the global minimum to zero
  global_min <- min(im_normalized)
  im_shifted <- im_normalized - global_min

  # 4. Scale so the global maximum is 1
  global_max <- max(im_shifted)
  im_final <- im_shifted / global_max
  return(im_final)
}


compute_cell_fact_score <- function(coords_1,cur_df_coords,cell_cent,im_cell,end_dists,
                                    min_mols=5,max_cells_per_fov=10,c_dist_expand=.025) {

  im_test <- im_cell

  ## getting the combinations of molecules to test
  # need to do this within z-slices
  all_z <- as.numeric(unique(coords_1$z_index))

  cur_df_coords_z <- lapply(all_z,function(z_val){
    return(cur_df_coords[cur_df_coords$z_index==z_val,])
  })
  names(cur_df_coords_z) <- paste0('z_',all_z)


  all_z <- as.numeric(unique(cur_df_coords$z_index))

  cell_xmin <- min(cur_df_coords$x)
  cell_xmax <- max(cur_df_coords$x)
  cell_ymin <- min(cur_df_coords$y)
  cell_ymax <- max(cur_df_coords$y)
  cell_bounds <- c(cell_xmin,cell_xmax,cell_ymin,cell_ymax)

  ## testing balancing distances
  all_cent_dists <- apply(cur_df_coords, 1, function(coords) {
    x <- coords['x']
    y <- coords['y']
    sqrt((x-cell_cent[1])**2 + (y-cell_cent[2])**2)
  })
  d_thresh <- max(all_cent_dists)*c_dist_expand # threshold to pick random molecules within

  cell_cent <- cell_cent - 1 # need to subtract 1, since compute_summary adds 1
  scores <- c()
  scores_perm <- c()
  for (rndx in 1:nrow(coords_1)) {
    cell_cent_same_z <- c(cell_cent,coords_1[rndx,'z_index'])
    coords_tmp <- unlist(coords_1[rndx, ])
    coords_tmp[1:2] <- coords_tmp[1:2] - 1 # need to subtract 1, since compute_summary adds 1 and need to stay in im bounds
    if (all(coords_tmp==cell_cent_same_z)) {
      next
    } else {
      mol_test <- rownames(coords_1)[rndx]
      score <- compute_summary(coords_tmp, cell_cent_same_z, im_test)
      score <- score * end_dists[mol_test]

      # now pick a random other molecule at a matched distance
      same_z_ndx <- which(cur_df_coords[,3]==coords_tmp[3])
      d_real <- all_cent_dists[mol_test]
      d_range <- c(d_real - d_thresh,d_real + d_thresh)
      all_cent_dists_same_z <- all_cent_dists[same_z_ndx]
      all_cent_dists_same_z <- all_cent_dists_same_z[all_cent_dists_same_z!=0]
      ndx_sel <- which(all_cent_dists_same_z>d_range[1] & all_cent_dists_same_z<d_range[2])
      if (length(ndx_sel)==1) {
        next
      }
      mol_sel <- names(all_cent_dists_same_z)[ndx_sel]
      mol_sel <- mol_sel[mol_sel != mol_test]
      mol_sel <- sample(mol_sel,1)
      coords_tmp <- unlist(cur_df_coords[mol_sel,])
      coords_tmp[1:2] <- coords_tmp[1:2] - 1
      score_perm <- compute_summary(coords_tmp, cell_cent_same_z, im_test)
      score_perm <- score_perm * end_dists[mol_sel]
    }
    scores <- c(scores, score)
    scores_perm <- c(scores_perm, score_perm)
  }

  if (length(scores)==0) { # if only have 1 molecule per z
    return(NA)
  }

  score_rel <- mean(log(scores/scores_perm))

  return(c(score_rel,NA))
}

compute_membrane_sep_per_cell_factor <- function(
    df_fov_sub, im, cell_pair_df, min_mols=5, max_cells_per_fov=10, c_dist_expand=.025
  ) {
  results <- list()
  non_na_count <- 0
  for (cur_cell_ndx in 1:nrow(cell_pair_df)) {
    # print(cur_cell_ndx)
    if (non_na_count>=max_cells_per_fov) {
      results[[cur_cell_ndx]] <- NA
      next
    }

    cur_cell <- cell_pair_df[cur_cell_ndx,1]
    adj_cell <- cell_pair_df[cur_cell_ndx,2]
    cur_df <- df_fov_sub[[cur_cell]]
    rownames(cur_df) <- paste0('mol',1:nrow(cur_df))
    adj_df <- df_fov_sub[[adj_cell]]
    rownames(adj_df) <- paste0('mol',1:nrow(adj_df))

    if (sum(cur_df$is_admixture)<min_mols) {
      results[[cur_cell_ndx]] <- NA
      next
    }

    cur_cent <- c(mean(cur_df$x),mean(cur_df$y))
    adj_cent <- c(mean(adj_df$x),mean(adj_df$y))

    # for each factor molecule get a measure of whether it's between the centroids of each cell
    cur_xy <- cur_df[,c('x','y')]

    cur_ends <- apply(cur_xy,MARGIN = 1,FUN = function(coords){
      direction_vector <- cur_cent-coords
      vector_length <- sqrt(sum(direction_vector^2))
      unit_vector <- direction_vector / vector_length
      new_point <- coords + unit_vector
      return(new_point)
    })
    cur_ends <- t(cur_ends)
    adj_ends <- apply(cur_xy,MARGIN = 1,FUN = function(coords){
      direction_vector <- adj_cent-coords
      vector_length <- sqrt(sum(direction_vector^2))
      unit_vector <- direction_vector / vector_length
      new_point <- coords + unit_vector
      return(new_point)
    })
    adj_ends <- t(adj_ends)

    # now get distances between the two ends
    end_dists <- sapply(1:nrow(cur_xy),function(i){
      mydist <- sqrt(sum((cur_ends[i,]-adj_ends[i,])^2))
    })

    end_dists <- end_dists/2 # dividing by 2 since that is the max
    names(end_dists) <- rownames(cur_df)

    cur_df_coords <- floor(cur_df[,c('x','y','z_index')])

    im_cell <- im[,min(cur_df_coords$y):max(cur_df_coords$y),min(cur_df_coords$x):max(cur_df_coords$x),drop=FALSE]
    cur_df_coords$y <- cur_df_coords$y - min(cur_df_coords$y) + 1
    cur_df_coords$x <- cur_df_coords$x - min(cur_df_coords$x) + 1
    mols_admix_coords <- cur_df_coords[which(cur_df$is_admixture),]
    cell_cent <- round(c(mean(cur_df_coords$x),mean(cur_df_coords$y)))

    score_admix <- compute_cell_fact_score(mols_admix_coords, cur_df_coords, cell_cent,
                                           im_cell, end_dists, min_mols=min_mols,
                                           max_cells_per_fov=max_cells_per_fov,
                                           c_dist_expand=c_dist_expand)

    if (!is.na(score_admix[1])) {
      non_na_count <- non_na_count + 1
    }

    results[[cur_cell_ndx]] <- score_admix
  }

  names(results) <- paste0(cell_pair_df[,1],':',cell_pair_df[,2])
  results <- results[!is.na(results)]

  return(results)
}



#' Run membrane test to determine factor-ctype pairs for removal
#'
#' @param df Data Frame of molecule coordinates, with mol_id, cell, gene columns
#' @param crf.f.labs Data Frame of CRF results
#' @param im.load.fn A function to load the images for a give field of view or the whole dataset
#' @param num.factors Number of factors
#' @param min.mols Minimum number of admixture molecules per cell (default=5)
#' @param knn.k Neighborhood size for knn graph used to find adjacent cells (default=10)
#' @param max.cells.per.fov Maximum number of cells to test per field of view (default=10)
#' @param all.fov Vector of field of views to test or NULL if all should be tested (default=NULL)
#' @param use.block.knn Logical whether to run the KNN in blocks over the image.
#' This can be helpful to reduce runtime for large datasets. (default=FALSE)
#' @param c.dist.expand A fraction of the biggest molecule-centroid distance. Used to
#' select random molecules against for comparison based on current molecule
#' centroid dist +- the expand. May need to increase this for datasets with less
#' coverage (default=.025).
#' @param ncm.ct knn counts matrix from molecules to cell types. Can only input if
#' there is a single fov for the whole dataset (default=NULL)
#' @param ncm.cell knn counts matrix from molecules to cells. Can only input if
#' there is a single fov for the whole dataset (default=NULL)
#' @param n.cores Number of cores to run parallel with (default=1)
#'
#' @return A list with the membrane test scores per cell for each ctype comparison
#' @export
run_memb_test <- function(
    df, crf.f.labs, im.load.fn, num.factors, min.mols=5, max.cells.per.fov=10, knn.k=10,
    all.fov=NULL, use.block.knn=FALSE, c.dist.expand=.025, ncm.ct=NULL, ncm.cell=NULL, n.cores=1
  ) {

  all_ct <- as.character(unique(df$celltype))
  if (is.null(all.fov)) {
    all.fov <- unique(df$fov)
  }
  all.fov <- all.fov[order(all.fov,decreasing=FALSE)]

  # get cell type pairs
  ct_perms <- gtools::permutations(n = length(all_ct), r = 2, v = all_ct)
  same_ct_lst <- lapply(all_ct, \(x) c(x,x))
  ct_perms <- do.call(rbind,c(same_ct_lst,list(ct_perms)))

  fov_res <- list()
  for (fov_ndx in 1:length(all.fov)) {
    fov <- all.fov[fov_ndx]

    mol_ndx_keep <- which(df$fov==fov)
    df_fov <- df[mol_ndx_keep,]
    df_fov$mol_id <- 1:nrow(df_fov)
    rownames(df_fov) <- df_fov$mol_id
    cell.annot <- unique(df_fov[,c('cell','celltype')])
    rownames(cell.annot) <- cell.annot$cell

    # get knn counts matrices if not alreay provided and if there are multiple fovs
    if (length(all.fov)>1 | is.null(ncm.ct) | is.null(ncm.cell)) {
      if (use.block.knn) {
        prep_dat <- get_mol_knn_blocks(df,knn_k=knn.k,n.cores=n.cores)
      } else {
        prep_dat <- get_mol_knn(df_fov,knn_k=knn.k)
      }
      ncm.ct <- prep_dat[[1]]
      ncm.cell <- prep_dat[[2]]
    }

    colnames(df_fov)[3] <- 'z_index'

    # load images for the fov
    im = im.load.fn(fov=fov, normalize=TRUE)

    all_z <- unique(df_fov$z_index)
    all_z_nm <- paste0('z',as.character(all_z))

    # now loop through factors since we can only test some cells with some factors
    f_res <- list()
    all_fact <- 1:num.factors
    for (f in all_fact) {
      message(paste0('fov ',fov_ndx,', factor ',f))

      # loop through cell type pairs
      results <- sccore::plapply(1:nrow(ct_perms),function(i) {
        # results <- lapply(1:nrow(ct_perms),function(i) {
        # print(i)
        ct1 <- ct_perms[i,1]
        ct2 <- ct_perms[i,2]
        cell_pair_df <- get_adj_pairs(df_fov,cell.annot,ncm.ct,ncm.cell,ct1,ct2)
        if (is.null(cell_pair_df)) {
          return(list())
        }

        if (nrow(cell_pair_df)==0) {
          return(list())
        }

        f_labs_fov <- crf.f.labs[mol_ndx_keep,,drop=FALSE]
        if (ncol(crf.f.labs)>1) {
          df_fov$factor <- f_labs_fov[,ct1]
        } else {
          df_fov$factor <- f_labs_fov[,1]
        }

        df_fov$is_admixture <- df_fov$factor == f

        df_fov_sub <- df_fov[df_fov$cell %in% unlist(c(cell_pair_df)),]
        df_fov_sub <- df_fov_sub %>%
          split(.$cell, drop = T)

        pres <- compute_membrane_sep_per_cell_factor(
          df_fov_sub, im, cell_pair_df, min_mols=min.mols,
          max_cells_per_fov=max.cells.per.fov, c_dist_expand=c.dist.expand
        )
        return(pres)
      },mc.preschedule=TRUE,n.cores=n.cores,progress=TRUE)
      # })
      names(results) <- paste0(ct_perms[,1],'_',ct_perms[,2])
      f_res[[f]] <- results
    }
    fov_res[[fov_ndx]] <- f_res
  }

  # now stitching results together
  ct_ct_all <- paste0(ct_perms[,1],'_',ct_perms[,2])
  f_ct_ct_combos <- expand.grid(all_fact,ct_ct_all,stringsAsFactors = FALSE)

  scores_final <- lapply(1:nrow(f_ct_ct_combos),function(combo_ndx){
    f <- f_ct_ct_combos[combo_ndx,1]
    ct_ct <- f_ct_ct_combos[combo_ndx,2]
    f_ct_all_scores <- lapply(1:length(all.fov),function(fov_ndx){
      scores <- fov_res[[fov_ndx]][[f]][[ct_ct]]
      return(scores)
    })
    f_ct_all_scores <- unlist(f_ct_all_scores,recursive = FALSE)
    return(f_ct_all_scores)
  })
  names(scores_final) <- paste0(f_ct_ct_combos[,1],'_',f_ct_ct_combos[,2])
  return(scores_final)
}

select_memb_remove <- function(all_pval_mats_rs_m,all_pval_mats_cs_m,p.thresh) {
  f_ct_rm_mem <- list()
  f_from_mem <- list()
  for (f in 1:ncol(all_pval_mats_rs_m)) {
    f_nm <- paste0('f_',f)
    pval <- all_pval_mats_rs_m[,f]
    ct_rm <- names(pval)[which(pval>(-log10(p.thresh)))]

    pval <- all_pval_mats_cs_m[,f]
    pval <- pval[!is.na(pval)]
    if (length(pval)==0) {
      next
    }
    ct_from <- names(pval)[which(pval==max(pval,na.rm=TRUE))]

    ct_rm <- ct_rm[!(ct_rm %in% ct_from)]
    if (length(ct_rm)>0 & !is.na(ct_rm[1])) {
      f_ct_rm_mem[[f_nm]] <- ct_rm
    }
    if (length(ct_from)>0 & !is.na(ct_rm[1]) & length(ct_rm)>0) {
      f_from_mem[[f_nm]] <- ct_from
    }
  }
  return(list(f_ct_rm_mem,f_from_mem))
}

#' Extract the membrane results into a list of factor-ctype combinations for removal
#'
#' @param combo.res List of results from bridge tests
#' @param all.ctypes Vector of cell types
#' @param nmf.type Character whether the cell type was run jointly with all cell types
#' together or separately per cell type ('joint' or 'ct')
#' @param p.thresh pvalue threshold below which to select factor-ct combinations for removal
#' @param adj.pvals Whether to adjust pvalues with FDR or not
#' @param get.p.mat Whether to get the p-value matrix instead of the list of pairs
#'
#' @return either a list of factor-ctype pairs or a matrix of pvalues.
#' @export
extract_memb_res <- function(
    combo.res, all.ctypes, nmf.type, p.thresh=.01, adj.pvals=FALSE, get.p.mat=FALSE
  ) {
  if (is.null(combo.res)) {
    return(NULL)
  }
  fmats <- list()
  pmats <- list()
  scores_f_all <- sapply(names(combo.res),function(x){
    strsplit(x,split='_')[[1]][[1]]
  })
  f_all <- unique(scores_f_all)
  for (f in f_all) {
    scores_sub <- combo.res[scores_f_all==f]

    ## new names parsing accounting for names with possible underscores in them
    combos_no_num <- sub("^\\d+_", "", names(scores_sub))
    all.ctypes <- all.ctypes[order(nchar(all.ctypes), decreasing = TRUE)]
    escaped_uniques <- stringr::str_replace_all(all.ctypes, "([\\^\\$\\.\\|\\?\\*\\+\\(\\)\\[\\{\\\\])", "\\\\\\1")
    pattern <- paste0("^(", paste(escaped_uniques, collapse = "|"), ")_(",
                      paste(escaped_uniques, collapse = "|"), ")$")
    matches <- stringr::str_match(combos_no_num, pattern)
    ct_from_all <- matches[, 3]
    ct_to_all <- matches[, 2]

    ct_from <- unique(ct_from_all)
    ct_to <- unique(ct_to_all)
    fmat <- matrix(NA,ncol=length(ct_from),nrow=length(ct_to))
    colnames(fmat) <- rownames(fmat) <- ct_from
    pmat <- fmat
    for (i in 1:length(scores_sub)) {
      scores_only <- sapply(scores_sub[[i]],function(x){
        return(x[1])
      })
      checks_only <- sapply(scores_sub[[i]],function(x){
        return(x[2])
      })

      # need to check if there are any -Inf vals caused by 0 stain vals
      if (is.list(scores_only)) {
        next
      }
      scores_only <- scores_only[is.finite(scores_only)]
      if (length(scores_only)<=1) {
        next
      }

      ndx_keep <- which(!is.nan(scores_only))
      scores_only <- scores_only[ndx_keep]
      checks_only <- checks_only[ndx_keep]

      pval <- t.test(scores_only,alternative = 'greater')$p.value
      myscore <- sum(scores_only>0)/length(scores_only)
      fmat[ct_to_all[i],ct_from_all[i]] <- myscore
      pmat[ct_to_all[i],ct_from_all[i]] <- pval
    }
    if (adj.pvals) {
      pmat_fdr <- p.adjust(pmat,method='fdr')
      pmat_fdr <- matrix(pmat_fdr,ncol=ncol(pmat),nrow=nrow(pmat))
      rownames(pmat_fdr) <- rownames(pmat)
      colnames(pmat_fdr) <- colnames(pmat)
      pmat_fdr <- -log10(pmat_fdr)
      pmats[[f]] <- pmat_fdr
    } else {
      pmat <- -log10(pmat)
      pmats[[f]] <- pmat
    }

    fmats[[f]] <- fmat
  }

  if (nmf.type=='joint') {
    all_pval_mats_rs_m <- collapse_row_col(pmats,c.type = 'r_max')
    all_pval_mats_cs_m <- collapse_row_col(pmats,c.type = 'c_mean')

    if (get.p.mat) {
      return(pmats)
    }

    memb_res <- select_memb_remove(all_pval_mats_rs_m,all_pval_mats_cs_m,p.thresh)

  } else {
    memb_res <- list()
    for (ct in rownames(pmats[[1]])) {
      ct_from <- lapply(1:length(pmats),function(mat_ndx) { # for a ct, get results from all it's factors
        x1 <- pmats[[mat_ndx]]
        admix_for_ct <- x1[ct,]
        admix_for_ct <- admix_for_ct[!is.na(admix_for_ct)]

        if (length(admix_for_ct)==0) {
          return(NA)
        }

        ct_rm <- names(admix_for_ct)[which(admix_for_ct>(-log10(p.thresh)))]

        if (length(ct_rm)>0) {
          ct_res <- names(admix_for_ct)[admix_for_ct==max(admix_for_ct)]
          return(ct_res)
        } else {
          return(NA)
        }

      })
      names(ct_from) <- paste0('f_',1:length(pmats))
      ndx_not_na <- which(!is.na(ct_from))
      ct_from <- ct_from[ndx_not_na]
      ct_to <- lapply(ct_from,function(x){
        return(ct)
      })
      memb_res[[ct]] <- list(ct_to,ct_from)
    }
  }

  return(memb_res)
}

### enrichment test

#' Run the enrichment test for each cell type's markers in each factor
#'
#' @param res NMF results object
#' @param markers Data Frame of markers with columns 'gene' and 'marker_of'
#' @param n.perm Number of permutations to run (default=10000)
#' @param p.thresh Threshold below which to select significant results for removal
#' @param adj.pvals Whether to adjust pvalues using FDR (default=TRUE)
#' @param crf.h The Neighborhood size used in CRF (default=NULL)
#' @param df Data Frame with molecule cordinates, mol_id, gene, cell, celltype columns (default=NULL)
#' @param n.cores Number of cores to run parallel with (default=1)
#'
#' @return A list of factor cell type combinations for removal
#' @export
get_enr <- function(
    res, markers, n.perm=10000, p.thresh=.1, adj.pvals=TRUE, crf.h=NULL,
    df=NULL, n.cores=1
  ) {

  # summarise gene loadings
  if (res$gene_sub) {
    if (is.null(df)) {
      stop('Need to provide df since genes were subsetted before NMF')
    }
    if (is.null(crf.h)) {
      stop('Need to provide crf.h since genes were subsetted before NMF')
    }
    H <- project_nmf(df, res, nmf.h=crf.h, n.cores=n.cores)
  } else {
    H = res@fit@H
  }

  k = nrow(H)

  H_all = t(H) %>% as.data.frame() %>%
    setNames(1:k) %>%
    tibble::rownames_to_column('gene') %>%
    reshape2::melt(id.var = 'gene', variable.name = 'factor', value.name = 'loading') %>%
    group_by(factor) %>%
    mutate(frac = loading/sum(loading)) %>%
    ungroup() %>%
    mutate(k = k)

  V = H_all %>%
    reshape2::dcast(gene ~ factor, value.var = 'frac', fill = 0) %>%
    tibble::column_to_rownames('gene')

  # need to subset markers to only the genes that were tested in the NMF!
  markers <- markers[markers$gene %in% rownames(V),]

  # filtering out cell types with very few marker genes
  marker_counts <- table(markers$marker_of)
  ct_keep <- names(marker_counts)[marker_counts>2]
  markers <- markers[markers$marker_of %in% ct_keep,]

  # marker permutation test
  marker_sim = markers %>%
    dplyr::count(marker_of) %>%
    split(.$marker_of) %>%
    parallel::mclapply(
      mc.cores = n.cores,
      function(marker_n) {

        n = marker_n$n
        marker_of = marker_n$marker_of

        lapply(
          1:n.perm,
          function(i) {
            set.seed(i)
            H_all %>% group_by(factor) %>%
              sample_n(n) %>%
              summarise(frac = sum(frac), .groups = 'drop') %>%
              mutate(i = i)
          }
        ) %>%
          bind_rows() %>%
          mutate(marker_of = marker_of)

      }
    ) %>%
    bind_rows()

  # computing pvals
  dat_marker = H_all %>%
    inner_join(markers, by = 'gene') %>%
    group_by(factor, marker_of) %>%
    summarise(frac = sum(frac), .groups = 'drop') %>%
    arrange(factor, -frac) %>%
    mutate(factor_label = paste0(factor)) %>%
    tidyr::complete(factor_label, marker_of, fill = list('frac' = 0))

  dat_marker$factor <- as.numeric(dat_marker$factor)
  marker_sim$factor <- as.numeric(marker_sim$factor)

  all_pv <- c()
  for (i in 1:nrow(dat_marker)) {
    ct <- dat_marker$marker_of[i]
    f <- dat_marker$factor[i]
    frac <- dat_marker$frac[i]
    marker_sim_sub <- marker_sim[marker_sim$factor==f,]
    marker_sim_sub <- marker_sim_sub[marker_sim_sub$marker_of==ct,]
    pval <- sum(marker_sim_sub$frac > frac) / nrow(marker_sim_sub)
    # to avoid having 0 pvals, set 0 to lowest non-zero perm
    if (pval==0) {
      pval <- 1/nrow(marker_sim_sub)
    }
    all_pv <- c(all_pv,pval)
  }
  dat_marker$pval <- all_pv

  # fdr adjusting all pvalues
  if (adj.pvals) {
    # correcting per factor, like with the other methods
    all_dat_marker_subs <- list()
    all_f <- unique(dat_marker$factor)
    for (f_ndx in 1:length(all_f)) {
      f <- all_f[f_ndx]
      dat_marker_sub <- dat_marker[dat_marker$factor==f,]
      dat_marker_sub$fdr <- p.adjust(dat_marker_sub$pval,method='fdr')
      all_dat_marker_subs[[f_ndx]] <- dat_marker_sub
    }
    dat_marker <- do.call(rbind.data.frame,all_dat_marker_subs)
    # dat_marker$fdr <- p.adjust(dat_marker$pval,method='fdr')
    dat_marker_sig <- dat_marker[dat_marker$fdr<p.thresh,]
  } else {
    dat_marker_sig <- dat_marker[dat_marker$pval<p.thresh,]
  }

  all_ct <- unique(dat_marker$marker_of)
  f_test <- unique(dat_marker_sig$factor)
  f.ct.rm <- list()
  for (f in f_test) {
    f_nm <- paste0('f_',f)
    # get cell types that were not enriched
    dat_marker_sig_sub <- dat_marker_sig[dat_marker_sig$factor==f,]
    ct_rm <- all_ct[!(all_ct %in% (dat_marker_sig_sub$marker_of))]
    f.ct.rm[[f_nm]] <- ct_rm
  }

  f.from <- list()
  for (f in f_test) {
    f_nm <- paste0('f_',f)
    dat_marker_sig_sub <- dat_marker_sig[dat_marker_sig$factor==f,]
    admix_from <- unique(dat_marker_sig_sub$marker_of)
    f.from[[f_nm]] <- admix_from
  }

  return(list(f.ct.rm,f.from))
}

### false positive checking

get_distant_cell_fracs <- function(cell_annot, cf_fracs_all, f, all_other_ct, admix_from, ncm_ct_full) {
  # first stratify the admix_to ct to those not near any admix_from
  admix_to_fracs <- list()
  for (ct in all_other_ct) {
    cells_test <- names(cell_annot)[cell_annot %in% ct]
    ncm_ct_sub <- ncm_ct_full[cells_test,admix_from,drop=FALSE]
    cells_keep_ndx <- which(Matrix::rowSums(ncm_ct_sub)==0) # picking cells not near any admix from cells
    cells_keep <- rownames(ncm_ct_sub)[cells_keep_ndx]
    if (length(cells_keep)<=1) {
      admix_to_fracs[[ct]] <- NA
      next
    }
    factor_fracs <- cf_fracs_all[cells_keep,]
    factor_fracs <- factor_fracs[,f]

    cells_bg <- cells_test[!(cells_test %in% cells_keep)]
    if (length(cells_bg)<=1) {
      admix_to_fracs[[ct]] <- NA
      next
    }
    factor_fracs_bg <- cf_fracs_all[cells_bg,]
    factor_fracs_bg <- factor_fracs_bg[,f]

    admix_to_fracs[[ct]] <- list(factor_fracs,factor_fracs_bg)
  }
  return(admix_to_fracs)
}


outlier_test <- function(fracs_all,admix_to,all_other_ct) {
  all_frac_expr <- c()
  for (ct in all_other_ct) {
    factor_fracs <- fracs_all[[ct]][[1]]
    if (is.na(factor_fracs[[1]][[1]])) {
      frac_cells_expr <- 0
    } else {
      frac_cells_expr <- sum(factor_fracs>.05) / length(factor_fracs)
    }
    all_frac_expr <- c(all_frac_expr,frac_cells_expr)
  }
  names(all_frac_expr) <- all_other_ct
  ndx_to <- which(all_other_ct==admix_to)
  iqrval <- IQR(all_frac_expr)
  q3 <- quantile(all_frac_expr,.75)
  upper_thresh <- q3 + (1.5*iqrval)
  is_outlier <- all_frac_expr[ndx_to] > upper_thresh

  is_outlier1 <- is_outlier & (all_frac_expr[ndx_to]>.1) # requiring at least 10% of cells to have frac expr over 5%
  is_outlier <- is_outlier1

  return(is_outlier)
}


#' Check for false positive removals for NMFs done per cell type
#'
#' @param mb.res List of suggested removals per cell type
#' @param ct.fracs.per.ct Factor fractions per cell per cell type
#' @inheritParams check_f_rm
#'
#' @return A vector of factor_ctypes suggested for removal and pruned of potential
#' false positives that are likely to be natively expressed factors.
#' @export
check_f_rm_per_ct <- function(
    cell.annot, mb.res, ct.fracs.per.ct, ncm.ct.full, do.clean=TRUE, median.thresh=.1
  ) {
  if (is.null(mb.res)) {
    return(NULL)
  }

  ct_res_fin <- list()
  for (ct in names(mb.res)) {
    f.ct.rm <- mb.res[[ct]][[1]]
    f.from <- mb.res[[ct]][[2]]
    cf.fracs.all <- ct.fracs.per.ct[[ct]]
    check_res <- check_f_rm(
      cell.annot, f.ct.rm, f.from, cf.fracs.all, ncm.ct.full, do.clean=do.clean,
      median.thresh=median.thresh
    )
    ct_res_fin[[ct]] <- check_res
  }
  return(ct_res_fin)
}

#' Check if any factors for removal are likely native to certain cell types and shouldn't
#' be removed.
#'
#' @param cell.annot character vector of cell types named by cell ids
#' @param f.ct.rm list of factors and the corresponding cell types for removal
#' @param f.from list of factors and the cell types they are likely sources of admixture from
#' @param cf.fracs.all Factor fractions per cell per
#' @param ncm.ct.full Nearest neighbor counts from molecules to cell types
#' @param do.clean Whether to remove detected false positives or just reformat
#' the results (default=TRUE)
#' @param median.thresh Median expression threshold above which false positives are flagged
#'
#' @return A vector of factor_ctypes suggested for removal and pruned of potential
#' false positives that are likely to be natively expressed factors.
#' @export
check_f_rm <- function(
    cell.annot, f.ct.rm, f.from, cf.fracs.all, ncm.ct.full,
    do.clean=TRUE,median.thresh=.1
  ) {
  if (is.null(f.ct.rm)) {
    return(NULL)
  }

  if (length(f.ct.rm)==0) {
    return(NA)
  }

  all.ct <- unique(cell.annot)

  ndx_keep <- sapply(1:length(f.ct.rm),function(i){
    if (!is.null(f.ct.rm[[i]])) {
      return(i)
    } else {
      return(NA)
    }
  })
  ndx_keep <- ndx_keep[!is.na(ndx_keep)]
  f.ct.rm <- f.ct.rm[ndx_keep]
  f.from <- f.from[ndx_keep]

  if (length(f.ct.rm)==0) {
    return(NA)
  }

  # remove any cases where sc had a ct that was enriched but it's not in spatial data
  ndx_keep <- sapply(1:length(f.from),function(i){
    if (f.from[[i]][1] %in% all.ct) {
      return(i)
    } else {
      return(NA)
    }
  })
  ndx_keep <- ndx_keep[!is.na(ndx_keep)]
  f.ct.rm <- f.ct.rm[ndx_keep]
  f.from <- f.from[ndx_keep]

  if (length(f.ct.rm)==0) {
    return(NA)
  }

  f_ct_rm_clean <- f.ct.rm
  for (f_nm in names(f.ct.rm)) {
    # print(f_nm)
    f <- strsplit(f_nm,split='_')[[1]][[2]]
    admix_from <- f.from[[f_nm]][1] # only use top admix from ct to test near/far analysis with
    all_admix_to <- f.ct.rm[[f_nm]]
    all_admix_to <- all_admix_to[!(all_admix_to %in% admix_from)]
    all_other_ct <- all.ct[!(all.ct %in% admix_from)]

    f_ct_rm_clean[[f_nm]] <- all_admix_to

    fracs_all <- get_distant_cell_fracs(cell.annot, cf.fracs.all, f, all_other_ct, admix_from, ncm.ct.full)

    for (admix_to in all_admix_to) {
      if (!(admix_to %in% all.ct)) { # if the ct isn't in spatial, then take it off the removal list
        ct_rm <- f_ct_rm_clean[[f_nm]]
        ct_rm <- ct_rm[ct_rm != admix_to]
        f_ct_rm_clean[[f_nm]] <- ct_rm
        next
      }

      if (is.na(fracs_all[[admix_to]][[1]][[1]])) { # if there aren't any cells far away, need to take F-CT off removal list to be safe
        ct_rm <- f_ct_rm_clean[[f_nm]]
        ct_rm <- ct_rm[ct_rm != admix_to]
        f_ct_rm_clean[[f_nm]] <- ct_rm
        next
      }
      factor_fracs <- fracs_all[[admix_to]][[1]]
      factor_fracs_bg <- fracs_all[[admix_to]][[2]]

      if (all(factor_fracs==0) & all(factor_fracs_bg==0)) {
        ct_rm <- f_ct_rm_clean[[f_nm]]
        ct_rm <- ct_rm[ct_rm != admix_to]
        f_ct_rm_clean[[f_nm]] <- ct_rm
        next
      }

      is_outlier <- outlier_test(fracs_all,admix_to,all_other_ct)

      # also checking if the far cells just express the factor to some high baseline level
      med_expr <- median(factor_fracs)

      expr_diff <- mean(factor_fracs_bg) - mean(factor_fracs)

      # remove the ct from the removal list if there is high med_expr or no DE or there is a less than 1% difference in mean expression fraction
      if (do.clean) {
        if (med_expr>median.thresh | is_outlier | expr_diff<0) {
          ct_rm <- f_ct_rm_clean[[f_nm]]
          ct_rm <- ct_rm[ct_rm != admix_to]
          f_ct_rm_clean[[f_nm]] <- ct_rm
        }
      }
    }
  }

  f_ct_rm_clean2 <- lapply(1:length(f_ct_rm_clean),function(i) {
    ct_rm <- f_ct_rm_clean[[i]]
    if (length(ct_rm)==0) {
      return(c())
    }
    f_nm <- names(f_ct_rm_clean)[i]
    f <- strsplit(f_nm,split='_')[[1]][[2]]
    return(paste0(f,'_',ct_rm))
  })
  f_ct_rm_clean2 <- unlist(f_ct_rm_clean2)
  if (is.null(f_ct_rm_clean2)) {
    return(NA)
  }

  return(f_ct_rm_clean2)
}


## other utilities

get_mol_knn <- function(df,knn_k=20) {
  # get molecule KNN matrix
  rownames(df) <- df$mol_id
  knn_mat <- knn_adjacency_matrix(df, k=knn_k)
  colnames(knn_mat) <- df$mol_id
  rownames(knn_mat) <- df$mol_id

  # now convert this matrix into a matrix of cell type that the molecule belongs to
  j = factor(df$celltype)
  ncm_ct = knn_mat %*% sparseMatrix(i=1:length(j), j=as.integer(j), x=rep(1, length(j)))
  colnames(ncm_ct) = levels(j)

  j = factor(df$cell)
  ncm_cell = knn_mat %*% sparseMatrix(i=1:length(j), j=as.integer(j), x=rep(1, length(j)))
  colnames(ncm_cell) = levels(j)

  return(list(ncm_ct,ncm_cell))
}

get_mol_knn_blocks <- function(df,knn_k=20,n.cores=1) {
  vals_test_x <- seq(from = min(df$x), to = max(df$x), length.out = 10)
  vals_test_y <- seq(from = min(df$y), to = max(df$y), length.out = 10)
  vals_comb <- expand.grid(vals_test_x, vals_test_y)
  colnames(vals_comb) <- c('x','y')

  # Define grid parameters
  num_x_bins <- 5  # Number of x divisions
  num_y_bins <- 5  # Number of y divisions

  # Compute bin breaks
  x_breaks <- seq(min(df$x), max(df$x), length.out = num_x_bins + 1)
  y_breaks <- seq(min(df$y), max(df$y), length.out = num_y_bins + 1)

  # Assign each point to a grid cell
  df$x_bin <- cut(df$x, breaks = x_breaks, include.lowest = TRUE, labels = FALSE)
  df$y_bin <- cut(df$y, breaks = y_breaks, include.lowest = TRUE, labels = FALSE)

  # ## just plotting to check that the grid works as intended
  # df_sub <- df[sample(1:nrow(df),20000),]
  # df_sub$bin_comb <- paste0(df_sub$x_bin,'_',df_sub$y_bin)
  # df_sub$bin_comb <- factor(df_sub$bin_comb,levels=sample(unique(df_sub$bin_comb)))
  # ggplot(df_sub,aes(x=x,y=y,color=bin_comb)) +
  #   geom_point()
  # ##

  # Split data into a list of smaller data frames by grid cell
  df_list <- split(df, list(df$x_bin, df$y_bin), drop = TRUE)
  all_ct <- unique(df$celltype)
  grid_res <- sccore::plapply(df_list,function(df_sub) {
    if (nrow(df_sub)==0) {
      return(NULL)
    }
    # get molecule KNN matrix
    rownames(df_sub) <- df_sub$mol_id
    knn_mat <- knn_adjacency_matrix(df_sub, k=knn_k)
    colnames(knn_mat) <- df_sub$mol_id
    rownames(knn_mat) <- df_sub$mol_id

    # now convert this matrix into a matrix of cell type that the molecule belongs to
    j = factor(df_sub$celltype)
    ncm_ct = knn_mat %*% sparseMatrix(i=1:length(j), j=as.integer(j), x=rep(1, length(j)))
    colnames(ncm_ct) = levels(j)
    missing_cols <- setdiff(all_ct, colnames(ncm_ct))
    if (length(missing_cols)>0) {
      tmp_mat <- matrix(0,ncol=length(missing_cols),nrow=nrow(ncm_ct))
      colnames(tmp_mat) <- missing_cols
      ncm_ct <- cbind(ncm_ct,tmp_mat)
    }
    ncm_ct <- ncm_ct[,all_ct]

    j = factor(df_sub$cell)
    ncm_cell = knn_mat %*% sparseMatrix(i=1:length(j), j=as.integer(j), x=rep(1, length(j)))
    colnames(ncm_cell) = levels(j)

    return(list(ncm_ct,ncm_cell))
  },mc.preschedule=TRUE,n.cores=n.cores,progress=TRUE)

  ncm_ct_all <- lapply(grid_res,function(x){
    return(x[[1]])
  })
  ncm_cell_all <- lapply(grid_res,function(x){
    return(x[[2]])
  })

  ncm_ct <- do.call(rbind,ncm_ct_all)

  # combine ncm_cell matrices into one
  all_row_names <- unique(unlist(lapply(ncm_cell_all, rownames)))
  all_col_names <- unique(unlist(lapply(ncm_cell_all, colnames)))
  row_map <- setNames(seq_along(all_row_names), all_row_names)
  col_map <- setNames(seq_along(all_col_names), all_col_names)
  triplet_list <- lapply(ncm_cell_all, function(m) {
    m_t <- as(m, "TsparseMatrix")  # Convert to triplet format (dgTMatrix)
    row_inds <- row_map[rownames(m)][m_t@i + 1]
    col_inds <- col_map[colnames(m)][m_t@j + 1]
    x_vals <- m_t@x  # Values remain the same
    data.frame(i = row_inds, j = col_inds, x = x_vals)
  })
  all_triplets <- do.call(rbind, triplet_list)
  ncm_cell <- sparseMatrix(
    i = all_triplets$i,
    j = all_triplets$j,
    x = all_triplets$x,
    dims = c(length(all_row_names), length(all_col_names)),
    dimnames = list(all_row_names, all_col_names)
  )

  # order rows the same way as df
  ncm_ct <- ncm_ct[rownames(df),]
  ncm_cell <- ncm_cell[rownames(df),]

  return(list(ncm_ct,ncm_cell))
}
