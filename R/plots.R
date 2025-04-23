#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @importFrom scales pretty_breaks
NULL

#' Get scaled average expression fractions per cell type from scRNA data
#'
#' @param sc.obj sc.obj is seurat object with a cell_type meta.data column
#'
#' @return matrix of z-transformed mean expression fractions per cell type
#' @export
prep_sc_for_plt <- function(sc.obj, cell.type.col=NULL) {
  if (is.null(cell.type.col)) {
    if ('celltype' %in% colnames(sc.obj@meta.data)) {
      cell.type.col <- 'celltype'
      cell.types <- sc.obj$celltype
    } else if ('cell_type' %in% colnames(sc.obj@meta.data)) {
      cell.type.col <- 'cell_type'
      cell.types <- sc.obj$cell_type
    } else {
      cell.type.col <- 'Idents'
      cell.types <- Seurat::Idents(sc.obj)
    }
    message('Using ', cell.type.col, ' column for cell type labels.')
  }

  if (is.null(sc.obj[['RNA']]$counts)) {
    stop("Need 'counts' slot in Seurat object.")
  }

  cell.types <- as.character(cell.types)

  ### creating a version of the fractional counts with z-transform across cell types
  lib.sizes <- Matrix::colSums(sc.obj[['RNA']]$counts)
  frac.counts <- sweep(sc.obj[['RNA']]$counts, 2, lib.sizes, '/')

  all.ct <- unique(cell.types)
  av.frac.counts <- lapply(all.ct, \(ct) {
    cells.keep <- colnames(sc.obj)[cell.types == ct]
    frac.counts.sub <- frac.counts[, cells.keep]
    Matrix::rowMeans(frac.counts.sub) # gene average fractions across cells of a given ct
  })
  av.frac.counts.all.ct <- do.call(cbind, av.frac.counts)
  av.frac.counts.scaled <- t(scale(t(av.frac.counts.all.ct)))
  colnames(av.frac.counts.scaled) <- all.ct

  # needs columns of 'gene', 'celltype', 'z'
  g <- rep(rownames(av.frac.counts.scaled), ncol(av.frac.counts.scaled))
  ct <- sapply(colnames(av.frac.counts.scaled), \(ct) {
    rep(ct, nrow(av.frac.counts.scaled))
  })
  ct <- c(ct)
  z <- c(av.frac.counts.scaled)
  ct.exp <- cbind.data.frame(g, ct, z)
  colnames(ct.exp) <- c('gene', 'celltype', 'z')
  ct.exp$celltype <- factor(ct.exp$celltype, levels=sort(all.ct))
  return(ct.exp)
}

#' Plot top NMF loadings for each factor, optionally also showing clean scRNA expression
#'
#' @param nmf.res NMF results object
#' @param ct.exp Matrix of scaled scRNA expression fractions (default=NULL)
#' @param f.show Vector of factors to plot. If NULL, will plot all (default=NULL)
#' @param nplt.row Number of rows to place all the heatmaps on (default=2)
#'
#' @return a list of ComplexHeatmap objects
#' @export
plot_nmf_loadings <- function(
    nmf.res, ct.exp=NULL, f.show=NULL, nplt.row=2, add.theme=theme_get()
  ) {
  k <- nrow(nmf.res@fit@H)
  H <- as.data.frame(nmf.res@fit@H) %>% t %>%
    as.data.frame() %>%
    setNames(1:k) %>%
    tibble::rownames_to_column('gene') %>%
    reshape2::melt(id.var='gene', variable.name='factor', value.name='loading') %>%
    group_by(factor) %>%
    mutate(frac=loading/sum(loading)) %>%
    ungroup()

  if (is.null(f.show)) {
    f.show <- c(1:k)
  }

  ggs <- lapply(f.show, \(i) {
    Hi <- H %>% filter(factor == i)

    gene_order <- Hi %>% arrange(-frac) %>% pull(gene) %>% head(10) %>% rev

    gg.bar <- Hi %>%
      filter(gene %in% gene_order) %>%
      mutate(gene=factor(gene, gene_order)) %>%
      ggplot(
        aes(x=frac, y=gene)
      ) +
      xlab('Loadings fraction') +
      ylab(paste0('Factor ', as.character(i))) +
      geom_col() +
      scale_x_continuous(breaks=round(seq(0, max(Hi$frac), max(Hi$frac)/2), digits=2), position="top") +
      theme(
        axis.title.x=element_text(size=6),
        axis.title.y=element_text(size=8, face="bold")
      )

    if (is.null(ct.exp)) {
      return(gg.bar)
    }

    gg.ct <- ct.exp %>%
      filter(gene %in% gene_order) %>%
      mutate(gene=factor(gene, gene_order)) %>%
      ggplot(
        aes(x=celltype, y=gene, fill=z)
      ) +
      geom_tile(show.legend=FALSE) +
      scale_fill_gradient2(
        low='blue', mid='white',
        high='red',
        oob=scales::oob_squish,
        limits=c(-2, 2),
        na.value='whitesmoke'
      ) +
      scale_x_discrete(drop=FALSE, expand=expansion(add=0), position='bottom') +
      scale_y_discrete(drop=FALSE, expand=expansion(add=0), position='right') +
      add.theme +
      theme(
        axis.text.x=element_text(size=6, angle=60, vjust=1, hjust=1),
        axis.text.y=element_text(size=6, hjust=1, vjust=0.5),
        axis.title=element_blank(),
        panel.border=element_rect(fill=NA, color='gray'),
        axis.line=element_blank()
      ) +
      xlab('') + ylab('')

    return(plot_grid(gg.bar, gg.ct, nrow=1, align='h', rel_widths=c(.3, .5)))

  })

  fig <- plot_grid(plotlist=ggs, nrow=nplt.row, align='hv')

  return(fig)
}

#' Plot admixture score ratios (P_admix clean)/(P_admix unclean) per cell.
#' Only implemented for joint NMF/CRF results currently.
#'
#' @param scores.orig Admixture scores data.frame for uncleaned cells
#' @param scores.cln Admixture scores data.frame for cleaned cells
#' @param f.rm Character vector of factor-ctype combinations to remove
#' @param df Data Frame of molecule coordinates, gene, cell, and mol_id columns
#' @param crf.res Data Frame of CRF molecule assignments (default=.9)
#' @param upper.lim Upper limit of score ratios to show (default=1.5)
#' @param min.mean.thresh Lower threshold, below which cell types are excluded
#' @param add.theme ggplot2 theme to add to the plots (default: current theme)
#'
#' @return a list of 4 plots including: admixture scores, score ratios, % cleaned,
#' and number of molecules before/after cleaning
#' @export
plot_cell_score_ratios <- function(
    scores.orig, scores.cln, f.rm, df, crf.res,
    upper.lim=1.5, min.mean.thresh=1e-3, add.theme=theme_get()
  ) {
  if (!('cell_type' %in% colnames(df))) {
    if ('celltype' %in% colnames(df)) {
      df$cell_type <- df$celltype
    } else {
      stop('Need a column in df called cell_type or celltype.')
    }
  }
  all.ct <- unique(scores.orig$cell_type)

  f.ct <- sapply(f.rm, \(x) {
    tmp <- strsplit(x, split='_')[[1]]
    paste0(tmp[2:length(tmp)], collapse='_')
  })
  f.f <- sapply(f.rm, \(x) strsplit(x, split='_')[[1]][[1]])

  df$factor <- crf.res[, 1]

  score.ratios <- list()
  for (ct in all.ct) {
    if (!(ct %in% f.ct)) {
      next
    }
    ndx <- which(f.ct == ct)
    f.ord <- as.numeric(f.f[ndx])
    f.ord <- as.character(f.ord[order(f.ord)])

    df.sub <- df[df$cell_type == ct, ]
    f.ord <- intersect(f.ord, df.sub$factor)
    if (length(f.ord) == 0) {
      next
    }
    t2 <- table(df.sub[, c('cell', 'factor')])
    cell.ndx.test <- which(rowSums(t2[, f.ord, drop=FALSE]) > 0)
    cells.w.f.rm <- rownames(t2)[cell.ndx.test]

    scores.orig.ct <- scores.orig[scores.orig$cell_type == ct, ]
    scores.cln.ct <- scores.cln[scores.cln$cell_type == ct, ]
    cells.both <- intersect(scores.orig.ct$cell, scores.cln.ct$cell)
    scores.orig.ct <- scores.orig.ct[scores.orig.ct$cell %in% cells.both, ]
    scores.cln.ct <- scores.cln.ct[scores.cln.ct$cell %in% cells.both, ]
    rownames(scores.orig.ct) <- scores.orig.ct$cell
    rownames(scores.cln.ct) <- scores.cln.ct$cell
    scores.cln.ct <- scores.cln.ct[rownames(scores.orig.ct), ]

    cells.fin <- intersect(rownames(scores.cln.ct), cells.w.f.rm)
    frac.cleaned <- length(cells.w.f.rm) / nrow(t2) # to show this fraction for all cells
    scores.orig.ct <- scores.orig.ct[cells.fin, ]
    scores.cln.ct <- scores.cln.ct[cells.fin, ]
    tmp <- do.call(cbind.data.frame, list(scores.orig.ct$scores,
                                       scores.cln.ct$scores,
                                       scores.cln.ct$scores/scores.orig.ct$scores))
    colnames(tmp) <- c('orig', 'clean', 'scores')
    tmp$cell_type <- ct
    tmp$cell <- scores.cln.ct$cell
    tmp$frac.cleaned <- frac.cleaned

    # calculating molecules per cell and molecules removed
    mol.counts <- rowSums(t2[cells.fin, ])
    mol.counts.rm <- rowSums(t2[cells.fin, f.ord, drop=FALSE])
    mol.counts.after <- mol.counts - mol.counts.rm
    tmp$mol.counts <- mol.counts
    tmp$mol.counts.after <- mol.counts.after
    rownames(tmp) <- tmp$cell

    if (mean(tmp$orig, trim=.1) > min.mean.thresh) {
      score.ratios[[ct]] <- tmp
    }
  }
  scores.ratios.all.ct <- do.call(rbind.data.frame, score.ratios)
  scores.ratios.all.ct <- scores.ratios.all.ct[!is.nan(scores.ratios.all.ct$scores), ]
  scores.ratios.all.ct$scores[scores.ratios.all.ct$scores > upper.lim] <- upper.lim
  scale.fun <- \(x) {sprintf("%.2f", x)}

  tmp1 <- scores.ratios.all.ct[, c('orig', 'cell_type')]
  tmp2 <- scores.ratios.all.ct[, c('clean', 'cell_type')]
  tmp1$cln_status <- 'unclean'
  tmp2$cln_status <- 'clean'
  colnames(tmp1)[1] <- 'scores'
  colnames(tmp2)[1] <- 'scores'
  tmp <- rbind.data.frame(tmp1, tmp2)

  msc <- tmp1 %>%
    group_by(cell_type) %>%
    summarise(mean_score=mean(scores, trim=.1))

  ct.ord <- msc$cell_type[order(msc$mean_score, decreasing=FALSE)]
  tmp$cell_type <- factor(tmp$cell_type, levels=ct.ord)

  p.theme <- theme(
    axis.text.y=element_text(size=7),
    axis.title.y=element_text(size=8),
    plot.title=element_text(hjust=0.5, size=8.5),
    axis.title.x=element_text(size=8),
    axis.text.x=element_text(size=7),
    legend.title=element_blank(),
    legend.text=element_text(size=7),
    legend.key.height=unit(.5, 'cm'),
    legend.key.width=unit(.5, 'cm'),
    axis.line=element_line(size=0.3),
    axis.ticks=element_line(size=0.3),
    panel.grid.major.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.y=element_blank(),
    panel.border=element_blank()
  )

  p1 <- ggplot(tmp, aes(y=cell_type, x=scores, fill=cln_status)) +
    geom_boxplot(outliers=FALSE, size=.3) +
    geom_hline(yintercept=0:length(unique(tmp$cell_type)) + .5, color='gray90') +
    ggtitle(paste0('Admixture scores')) +
    xlab('P_admix') +
    ylab('Cell type') +
    scale_x_continuous(labels=scale.fun, breaks=pretty_breaks(n=3)) +
    guides(fill=guide_legend(reverse=TRUE)) +
    add.theme +
    p.theme

  scores.ratios.all.ct$cell_type <- factor(scores.ratios.all.ct$cell_type, levels=ct.ord)
  p2 <- ggplot(scores.ratios.all.ct, aes(x=scores, y=cell_type, fill=factor(after_stat(quantile), levels=c('4', '3', '2', '1')))) +
    geom_hline(yintercept=0:length(unique(scores.ratios.all.ct$cell_type)) + .5, color='gray90') +
    ggridges::stat_density_ridges(
      scale=1.25, linewidth=.3, geom="density_ridges_gradient",
      calc_ecdf=TRUE, quantiles=4, quantile_lines=TRUE
    ) +
    geom_vline(xintercept=1, color='red', linetype='dashed') +
    scale_fill_viridis_d(name="Quartiles") +
    xlab('(P_admix after cleaning) / (P_admix before cleaning)') +
    ylab('Cell type') +
    ggtitle('Admixture score ratios per cell') +
    scale_x_continuous(
      limits=c(0, upper.lim),
      expand=expansion(mult=c(0.025, 0.025)),
      breaks=c(0, .5, 1, upper.lim),
      labels=c("0", ".5", "1", paste0('>=', upper.lim))
    ) +
    guides(fill=guide_legend(reverse=TRUE)) +
    add.theme +
    p.theme +
    theme(
      legend.title=element_text(size=8),
      legend.key.height=unit(.3, 'cm'),
      legend.key.width=unit(.3, 'cm'),
    )

  f.c.unq <- unique(scores.ratios.all.ct[, c('cell_type', 'frac.cleaned')])
  f.c.unq$cell_type <- factor(f.c.unq$cell_type, levels=ct.ord)
  f.c.unq$frac.cleaned <- f.c.unq$frac.cleaned * 100

  p3 <- ggplot(f.c.unq, aes(y=cell_type, x=frac.cleaned)) +
    geom_bar(stat='identity', fill='darkred', alpha=.8, width=.7) +
    geom_hline(yintercept=0:length(unique(f.c.unq$cell_type)) + .5, color='gray90') +
    scale_x_continuous(labels=scale.fun, breaks=pretty_breaks(n=4)) +
    ggtitle('\nPercent cells modified') +
    xlab('%') +
    ylab('Cell type') +
    add.theme +
    theme(
      plot.title=element_text(hjust=0.5, size=8),
      axis.title.x=element_text(size=7),
      axis.text.x=element_text(size=6)
    )

  tmp1 <- scores.ratios.all.ct[, c('mol.counts', 'cell_type')]
  tmp2 <- scores.ratios.all.ct[, c('mol.counts.after', 'cell_type')]
  tmp1$mol_type <- 'unclean'
  tmp2$mol_type <- 'clean'
  colnames(tmp1)[1] <- 'counts'
  colnames(tmp2)[1] <- 'counts'
  tmp <- rbind.data.frame(tmp1, tmp2)

  tmp$mol_type <- factor(tmp$mol_type, levels=c('clean', 'unclean'))
  p4 <- ggplot(tmp, aes(y=cell_type, x=counts, fill=mol_type)) +
    geom_boxplot(outliers=FALSE, size=.3) +
    geom_hline(yintercept=0:length(unique(tmp$cell_type)) + .5, color='gray90') +
    ggtitle('Molecule counts in modified cells') +
    xlab('Number of molecules per cell') +
    scale_x_continuous(labels=scale.fun, breaks=pretty_breaks(n=3)) +
    guides(fill=guide_legend(reverse=TRUE)) +
    ylab('Cell type') +
    add.theme +
    p.theme

  return(list(p1, p2, p3, p4))
}

#' Plot p-value results from bridge- or membrane-based factor annotation
#'
#' @param rs The p-value rowmax matrix for admixture into each cell types.
#' Dims are cell types x factors
#' @param ad.sources A list of reported "admixture-from" cell types for each factor
#' @param final.rm The vector of factor cell type combinations for removal, structured
#' as f_ctype.
#' @param annot.type Either 'membrane' or 'bridge' depending on the type of annotation
#' used to generate the other inputs.
#' @param p.bridge The p-value used to extract bridge results
#' @param p.memb The p-value used to extract membrane results

#' @return a ComplexHeatmap object showing p-values, as well as which cell types are
#' recommended for removal (labeled *) and which source cell types are behind each
#'factor (labeled S)
#' @export
plot_annot_hmap <- function(rs, ad.sources, final.rm, annot.type, p.bridge=.1, p.memb=.01) {
  f.final <- unique(sapply(final.rm, \(x) strsplit(x, split='_')[[1]][[1]]))
  if (annot.type == 'membrane') {
    col_fun <- circlize::colorRamp2(c(0, -log10(p.memb), 10), c("white", "white", "red"))
  } else if (annot.type == 'bridge') {
    col_fun <- circlize::colorRamp2(c(0, -log10(p.bridge), 10), c("white", "white", "red"))
  }

  colnames(rs) <- paste0('F', 1:ncol(rs))
  hmap <- ComplexHeatmap::Heatmap(
    rs, name="-log_pmax", cluster_rows=FALSE, cluster_columns=FALSE, show_row_dend=FALSE,
    show_column_dend=FALSE, col=col_fun, border=TRUE, row_names_side='left',
    column_title='Factor', row_title='Cell type', column_title_side='bottom',
    column_names_side='bottom', row_title_side='left', na_col='white',
    cell_fun=function(j, i, x, y, width, height, fill) {
      f.nm <- paste0('f_', j)
      ct.test <- rownames(rs)[i]
      f.ct <- paste0(j, '_', ct.test)
      fontsize <- min(grid::convertWidth(width, "mm", valueOnly=TRUE),
                    grid::convertHeight(height, "mm", valueOnly=TRUE)) * 2
      if (f.ct %in% final.rm) {
        # ★
        grid::grid.text(
          sprintf("\u2605"), x, y,
          gp=gpar(fontsize=fontsize, fontface="bold", fontfamily="DejaVu Sans"),
          just="center"
        )
        # grid::grid.text(latex2exp::TeX(r'$\star$'), x, y,
        #                 gp = gpar(fontsize = fontsize,markdown = TRUE),
        #                 just = "center")
      } else if (f.nm %in% names(ad.sources) & as.character(j) %in% f.final) {
        from.ct <- ad.sources[[f.nm]][1]
        ndx.star <- which(rownames(rs) == from.ct)
        if (ndx.star == i) {
          grid::grid.text(sprintf("S"), x, y, gp=gpar(fontsize=fontsize, fontface="bold"))
        }
      }
    })
  return(hmap)
}

#' Plot the preservation of gene correlation
#'
#' @param cor.cors A vector of correlation values
#' @param bins The number of bins to use in the histogram
#'
#' @return A ggplot object
#' @export
plot_correlation_preservation <- function(cor.cors, bins=30) {
  ggplot(data.frame(x=cor.cors), aes(x=x)) +
    geom_histogram(bins=bins) +
    ggtitle("Preservation of gene correlation") +
    xlab("Correlation of gene correlations") + ylab("Num. genes") +
    xlim(min(cor.cors), 1.0) +
    scale_y_continuous(expand=c(0, 0, 0.05, 0))
}