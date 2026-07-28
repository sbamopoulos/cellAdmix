#' @export
run_pagoda_de <- function(cm, groups) {
  p2 <- pagoda2::Pagoda2$new(cm, n.cores=1, min.cells.per.gene=0, verbose=FALSE)
  de <- p2$getDifferentialGenes(
    groups=groups, z.threshold=0,
    upregulated.only=FALSE, append.specificity.metrics=FALSE
  )

  return(de)
}

#' @export
de_volcano_plot <- function(de.df, tumor.marker.genes, genes.to.label) {
  ## making volcano plot for original DE
  de.df$padj <- 2*pnorm(abs(de.df$Z), mean=0, sd=1, lower.tail=FALSE) # convert z back to padj
  de.df$logpv <- -log10(de.df$padj)
  ndx_lab <- which(de.df$Gene %in% genes.to.label)
  de.df$top_mark_other <- NA
  de.df[ndx_lab,'top_mark_other'] <- de.df$Gene[ndx_lab]
  de.df$logpv[de.df$logpv > 30] <- 30
  ymax <- max(de.df$logpv, na.rm=TRUE) + 1
  de.df$mark_other <- de.df$Gene %in% tumor.marker.genes
  de.df$mark_other <- factor(de.df$mark_other, levels=c(TRUE, FALSE))
  levels(de.df$mark_other) <- c('Malignant marker', 'other')

  myColors <- c('red4', 'grey70')

  # Helper variables
  limits <- range(de.df$logpv, na.rm=TRUE, finite=TRUE)
  step <- diff(limits) * 0.1
  size <- 0.45 * step

  de.df$mark_other <- factor(de.df$mark_other, levels=c('Malignant marker', 'other'))
  bottom_lim <- min(limits[1] - as.numeric(de.df$mark_other) * step - size) - .5

  p <- ggplot(de.df, aes(x=M, y=logpv, color=mark_other, label=top_mark_other)) +
    geom_segment(aes(y=ymax, yend=0, x=log2(1.5), xend=log2(1.5)), color='red', linetype='dashed') +
    geom_segment(aes(y=ymax, yend=0, x=-log2(1.5), xend=-log2(1.5)), color='red', linetype='dashed') +
    geom_hline(yintercept=0) +
    geom_hline(yintercept=-log10(.05), color='red', linetype='dashed') +
    geom_segment(
      aes(
        color=mark_other,
        xend=M,
        y=limits[1] - as.numeric(mark_other) * step + size,
        yend=limits[1] - as.numeric(mark_other) * step - size
      )
    ) +
    ggrepel::geom_text_repel(size=4, show.legend=FALSE) +
    geom_vline(xintercept=0) +
    geom_point(data=de.df[de.df$mark_other == 'other',], alpha=.5, size=1.5) +
    geom_point(data=de.df[de.df$mark_other == 'Malignant marker',], alpha=.5, size=1.5) +
    scale_y_continuous(limits=c(bottom_lim, ymax), expand=c(0, 0)) +
    xlab('logFC') +
    ylab('-log10(Padj)') +
    scale_colour_manual(breaks=c('Malignant marker', 'other'), values=myColors) +
    theme_classic(base_line_size=.5) +
    theme(plot.title=element_text(hjust=0.5),
          legend.title=element_blank())

  return(p)
}

#' @export
plot_expression_comparison <- function(
    cm.orig, cm.clean, markers.native, markers.admix, admixture.type, native.type,
    combine=TRUE
) {
  dat_orig_norm <- sweep(cm.orig, 2, colSums(cm.orig), '/')
  dat_cln_norm <- sweep(cm.clean, 2, colSums(cm.clean), '/')

  admix.av.expr.before <- Matrix::rowMeans(dat_orig_norm[markers.admix,])
  native.av.expr.before <- Matrix::rowMeans(dat_orig_norm[markers.native,])

  admix.av.expr.after <- Matrix::rowMeans(dat_cln_norm[markers.admix,])
  native.av.expr.after <- Matrix::rowMeans(dat_cln_norm[markers.native,])

  tmp1a <- cbind.data.frame(markers.admix, admix.av.expr.before, 'original')
  colnames(tmp1a) <- c('gene', 'av_expr', 'type')
  tmp1b <- cbind.data.frame(markers.admix, admix.av.expr.after, 'cleaned')
  colnames(tmp1b) <- c('gene', 'av_expr', 'type')
  tmp1 <- rbind.data.frame(tmp1a, tmp1b)
  tmp1$type <- factor(tmp1$type, levels=c('original', 'cleaned'))

  p1 <- ggplot2::ggplot(tmp1, aes(x=gene, y=av_expr, fill=type)) +
    geom_bar(stat="identity", position=position_dodge()) +
    xlab(paste(admixture.type, 'marker')) +
    ylab('Average expression') +
    ggtitle(paste(admixture.type, 'marker expression\n in', native.type)) +
    theme(plot.title=element_text(hjust=0.5))

  tmp1a <- cbind.data.frame(markers.native, native.av.expr.before, 'original')
  colnames(tmp1a) <- c('gene', 'av_expr', 'type')
  tmp1b <- cbind.data.frame(markers.native, native.av.expr.after, 'cleaned')
  colnames(tmp1b) <- c('gene', 'av_expr', 'type')
  tmp1 <- rbind.data.frame(tmp1a, tmp1b)
  tmp1$type <- factor(tmp1$type, levels=c('original', 'cleaned'))

  p2 <- ggplot2::ggplot(tmp1, aes(x=gene, y=av_expr, fill=type)) +
    geom_bar(stat="identity", position=position_dodge()) +
    xlab(paste(native.type, 'marker')) +
    ylab('Average expression') +
    ggtitle(paste(native.type, 'marker expression\n in', native.type)) +
    theme(plot.title=element_text(hjust=0.5))

  if (!combine) {
    return(list(p1, p2))
  }

  cowplot::plot_grid(p1, p2, nrow=1)
}

plot_de_comparison <- function(de.orig, de.clean, genes.to.label, tumor.marker.genes) {
  z_thresh <- qnorm(.01/2, lower.tail=FALSE)
  tmp <- de.orig %>%
    {data.frame(orig=.[,'Z'], clean=de.clean[rownames(.),'Z'], row.names=rownames(.))}

  tmp %<>% as.matrix() %>% pmin(10) %>% pmax(-10) %>% as.data.frame()

  ## trying to color points by tumor marker set belonging
  tmp$markers <- ifelse(rownames(tmp) %in% tumor.marker.genes, 'Malignant marker', 'other') %>%
    factor(levels=c('other', 'Malignant marker'))

  ndx_lab <- match(genes.to.label, rownames(tmp))
  genes.to.label <- genes.to.label[!is.na(ndx_lab)]
  ndx_lab <- ndx_lab[!is.na(ndx_lab)]
  tmp$glab <- NA
  tmp[ndx_lab,'glab'] <- genes.to.label

  p <- ggplot(tmp, aes(x=orig, y=clean, color=markers, label=glab)) +
    geom_point(data=tmp[tmp$markers == 'other',], alpha=.5, size=1.5) +
    geom_point(data=tmp[tmp$markers == 'Malignant marker',], alpha=.5, size=1.5) +
    ggrepel::geom_text_repel(size=4, show.legend=FALSE) +
    xlab('Z-score (original)') +
    ylab('Z-score (cleaned)') +
    geom_hline(yintercept=z_thresh, color='red', linetype='dashed') +
    geom_hline(yintercept=-1*z_thresh, color='red', linetype='dashed') +
    geom_hline(yintercept=0) +
    geom_abline(slope=1, intercept=0) +
    scale_colour_manual(breaks=c('Malignant marker', 'other'), values=c('red4', 'grey70')) +
    theme_classic() +
    theme(plot.title=element_text(hjust=0.5), legend.title=element_blank())

  return(p)
}


check_fp <- function(
    df, cell_annot, crf_res, bridge_res, do_clean=TRUE, knn_k=100, median_thresh=.1
) {
  ##### first get cell neighbor knn matrix
  meta_sub <- cell_annot[,c('x','y','z')]
  knn_mat <- knn_adjacency_matrix(meta_sub, k=knn_k)
  colnames(knn_mat) <- rownames(meta_sub)
  rownames(knn_mat) <- rownames(meta_sub)
  
  # now convert this matrix into a matrix of cell type that the molecule belongs to
  j = factor(cell_annot$celltype)
  ncm_ct_full = knn_mat %*% sparseMatrix(i=1:length(j), j=as.integer(j), x=rep(1, length(j)))
  colnames(ncm_ct_full) = levels(j)
  
  celltypes <- cell_annot %$% setNames(celltype, rownames(.))

  df$factor <- crf_res[,1]
  factor_counts <- table(df[,c('cell','factor')])
  cf_fracs <- factor_counts / Matrix::rowSums(factor_counts)
  ct_fracs_per_ct <- cf_fracs

  check_res <- check_f_rm(celltypes, bridge_res[[1]], bridge_res[[2]], ct_fracs_per_ct, ncm_ct_full,
                          do.clean=do_clean, median.thresh=median_thresh
  )
  
  return(check_res)
}
