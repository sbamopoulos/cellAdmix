#' Example per cell meta-data for testing and demonstrating cellAdmix code.
#' It is downsampled from an NSCLC CosMx dataset, which can be downloaded at
#' https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/
#'
#' @format A data frame with 412 rows and 6 variables:
#' \describe{
#'   \item{x}{cell centroid x location}
#'   \item{y}{cell centroid y location}
#'   \item{cell}{cell ID}
#'   \item{celltype}{cell type}
#'   \item{celltype_coarse}{optional secondary cell type annotation}
#'   \item{niche}{optional tissue region annotation}
#' }
#' @source Generated using 'data-raw/data_prep_vignette.R'
"cell_meta"

#' Example spatial transcript data for testing and demonstrating cellAdmix code.
#' It is downsampled from an NSCLC CosMx dataset, which can be downloaded at
#' https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/
#'
#' @format A data frame with 100 rows and 3 variables:
#' \describe{
#'   \item{x}{transcript x location}
#'   \item{y}{transcript y location}
#'   \item{z}{transcript z location}
#'   \item{gene}{gene label}
#'   \item{cell}{cell ID}
#'   \item{celltype}{cell type}
#'   \item{mol_id}{unique molecule ID}
#' }
#' @source Generated using 'data-raw/data_prep_vignette.R'
"df"