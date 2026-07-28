### make df into a counts matrix and a metadata matrix per cell

#' Transform transcript molecule df into standard seurat object with counts matrix
#'
#' @param tx.dat dataframe of transcripts with at least gene, cell, and celltype columns
#'
#' @return Seurat object with the data in counts matrix form
#' @export
get_counts_meta_seurat <- function(tx.dat, normalize=FALSE) {
  gene.factors <- factor(tx.dat$gene)
  cell.factors <- factor(tx.dat$cell)

  cm <- Matrix::sparseMatrix(
    i = as.integer(gene.factors),  # Now rows are genes
    j = as.integer(cell.factors),  # And columns are cells
    x = rep(1, nrow(tx.dat)),      # Each row in `tx.dat` contributes 1 count
    dims = c(length(levels(gene.factors)), length(levels(cell.factors))),
    dimnames = list(
      levels(gene.factors),        # Row names (gene names)
      levels(cell.factors)         # Column names (cell names)
    )
  )

  meta <- tx.dat[,c("cell", "celltype")] %>%
    data.table::as.data.table() %>%
    unique() %>%
    as.data.frame()

  rownames(meta) <- meta$cell
  meta <- meta[colnames(cm),]

  s.obj <- Seurat::CreateSeuratObject(cm, meta.data=meta)

  if (normalize) {
    s.obj <- Seurat::NormalizeData(s.obj)
  }
  return(s.obj)
}


