
library(giotto)
library(dplyr)
library(magrittr)
library(Seurat)
library(readr)
library(usethis)

# load the giotto object called gem
load('./data/spatial_data/dialogue_lung/giotto_dat/SMI_Giotto_Object.RData')

# extract relevant metadata from giotto object
cell_meta <- gem@cell_metadata$rna
cell_locs <- gem@spatial_locs$raw
cell_meta <- cbind.data.frame(cell_locs[,c(1:2)],cell_meta)

# subset to a single donor/replicate
cell_meta <- cell_meta[cell_meta$Run_Tissue_name=='Lung5_Rep1',]
colnames(cell_meta)[3] <- 'cell'
colnames(cell_meta)[34] <- 'celltype'

# rename tumor cells
cell_meta$celltype %<>% {ifelse(startsWith(., 'tumor'), 'malignant', .)}

# group all immune cells under one annotation for the visualization
immune_cell_types <- c(
  'B-cell', 'NK', 'T CD4 memory', 'T CD4 naive', 'T CD8 memory', 'T CD8 naive', 'Treg', 
  'plasmablast', 'mast', 'mDC', 'monocyte', 'pDC', 'neutrophil'
)

cell_meta$cell_type_coarse <- cell_meta$celltype %>% {ifelse(. %in% immune_cell_types, 'immune other', .)}
rownames(cell_meta) <- cell_meta$cell

# load the molecule-level spatial data 
tx_dat <- as.data.frame(read_csv('./data/spatial_data/dialogue_lung/smi_raw/Lung5_Rep1/Lung5_Rep1-Flat_files_and_images/Lung5_Rep1_tx_file.csv',
                                 show_col_types = FALSE))
tx_dat$cell <- paste0('c_1_',tx_dat$fov,'_',tx_dat$cell_ID)

# subsetting data to same cells we have annotations for
df <- tx_dat[tx_dat$cell %in% cell_meta$cell,]

# append cell type annotations to molecule-level data
match_ndx <- match(df$cell,cell_meta$cell)
df$celltype <- cell_meta$celltype[match_ndx]

# Change x and y coordinate column names
colnames(df)[3:4] <- c('x','y')

# Change gene column name
colnames(df)[8] <- c('gene')

# adding a column for molecule ID
df$mol_id <- as.character(1:nrow(df))

# converting coordinate units from pixels to microns based on information from their publication
df$x <- (df$x * 180) / 1000
df$y <- (df$y * 180) / 1000

# assigning z-coordinates to an expected height in microns
df$z <- (df$z * 800) / 1000


# downsampling the data, so it takes less memory
# selecting more fibroblasts than other cells
cells_keep1 <- sample(rownames(cell_meta)[cell_meta$celltype=='fibroblast'],300)
cells_keep2 <- sample(rownames(cell_meta),200)
cells_keep <- unique(c(cells_keep1,cells_keep2))
cell_meta <- cell_meta[cells_keep,c('sdimx','sdimy','cell','celltype','cell_type_coarse','niche')]
colnames(cell_meta)[1:2] <- c('x','y')
df <- df[df$cell %in% cells_keep,c('x','y','z','gene','cell','celltype','mol_id')]

# downsampling molecules in df
df <- df[sample(1:nrow(df),25000),]
# removing cells that have less than 20 molecules
mcounts <- table(df$cell)
cells_keep <- names(mcounts)[mcounts>=20]
df <- df[df$cell %in% cells_keep,]
cell_meta <- cell_meta[cells_keep,]

setwd('./cellAdmix/')
usethis::use_data(cell_meta, compress = "xz", overwrite = TRUE)
usethis::use_data(df, compress = "xz", overwrite = TRUE)











