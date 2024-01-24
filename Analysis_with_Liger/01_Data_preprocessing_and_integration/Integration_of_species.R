# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
projects_dir <- "/home/ytamal2/Documents/2023/PhD_projects_Yasir/scExplorer/Functions/"

# Functions - Markers identification
list_1 <- list.files(paste0(projects_dir, "Library_handler"), pattern = "*.R$", full.names = TRUE) 
sapply(list_1, source, .GlobalEnv)

# Functions - Data manipulation
list_2 <- list.files(paste0(projects_dir, "Functions_matrix_manipulation"), pattern = "*.R$", full.names = TRUE) 
sapply(list_2, source, .GlobalEnv)

# Functions - Library and packages handler
list_3 <- list.files(paste0(projects_dir, "Functions_marker_identification"), pattern = "*.R$", full.names = TRUE) 
sapply(list_3, source, .GlobalEnv)


###
# Load data tables
###

###
# Protoplasting induced genes
###

# Load the table containing the list of protoplasting-induced genes.
PP_genes_table = read.csv("/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/Protoplasting_genes/Ox_Co0_leaf_protoplast_v12_final_August_2022_ortho.csv")

# Gene IDs - protoplasting-induced genes
PP_genes = PP_genes_table$GeneID

###
# Orthologues table
###

ortho_table = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Input_files/Additional_inputs/Orthologues_n_correspondence/Orthos_table.csv")

###
# WT C. hirsuta
###

# Load data - WT OX SAM (identifier - experiment #11th) - 2490 cells, 27773 genes
OX_data_1E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Ox_RNA_11th_SAM_Newest/filtered_feature_bc_matrix/")

# Load data - WT OX SAM (identifier - experiment #12th_A) - 2200 cells, 27773 genes
OX_data_2E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Ox_RNA_12th_SAM_A_Newest/filtered_feature_bc_matrix/")

# Load data - WT OX SAM (identifier - experiment #12th_B) - 2980 cells, 27773 genes
OX_data_3E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Ox_RNA_12th_SAM_B_Newest/filtered_feature_bc_matrix/")


# Convert the gene ids in the data table to ortho gene ids
OX_DF_1 = prepare_ortho_data(input_data = OX_data_1E, ortho_data = ortho_table, ortho_column_name_of_gene_ids = "C.hirsutaOX", ortho_column_name_to_assign = "A.thaliana.TAIR10")

OX_DF_2 = prepare_ortho_data(input_data = OX_data_2E, ortho_data = ortho_table, ortho_column_name_of_gene_ids = "C.hirsutaOX", ortho_column_name_to_assign = "A.thaliana.TAIR10")

OX_DF_3 = prepare_ortho_data(input_data = OX_data_3E, ortho_data = ortho_table, ortho_column_name_of_gene_ids = "C.hirsutaOX", ortho_column_name_to_assign = "A.thaliana.TAIR10")


###
# A. thaliana - SAM (Apex)
###

# Load data - WT COL0 SAM (identifier - experiment #12th) - 3050 cells, 27629 genes
COL_SAM <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Col0_RNA_12th_SAM_New/filtered_feature_bc_matrix/", gene.column = 1)


# All gene IDs - Arabidopsis Thaliana
thaliana_genes = rownames(COL_SAM)

# extracting the Cardamine IDs that are present in orthologues table 
thaliana_ortho_genes = as.character(ortho_table$A.thaliana.TAIR10)

# not all the thaliana ids are present in the ortho data - 51 thaliana genes are missing in the thaliana data
thaliana_ortho_genes = intersect(thaliana_genes, thaliana_ortho_genes)

# SAM data
COL_SAM <- COL_SAM[thaliana_ortho_genes, ]

# remove the missing genes from the data
OX_DF_1 <- OX_DF_1[thaliana_ortho_genes, ]
OX_DF_2 <- OX_DF_2[thaliana_ortho_genes, ]
OX_DF_3 <- OX_DF_3[thaliana_ortho_genes, ]

# All gene IDs - all datasets
ortho_genes = rownames(COL_SAM)

# Remove protoplasting-induced genes from the total set of hirsuta genes
genes_to_keep = setdiff(ortho_genes, PP_genes)

##### Remove the protoplasting induced genes
OX_DF_1 <- OX_DF_1[genes_to_keep, ]
OX_DF_2 <- OX_DF_2[genes_to_keep, ]
OX_DF_3 <- OX_DF_3[genes_to_keep, ]

# SAM data
COL_SAM <- COL_SAM[genes_to_keep, ]


###
# OX - 1 E
###

# First replicate - OX 1E - total cells 6640; filter out genes that are not detected in at least 13 cells
OX_1E <- CreateSeuratObject(counts = OX_DF_1, project = "OX_1E", min.features = 200)

# Add metadata information to the seurat object
OX_1E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Hirsuta", "OX-1", "WT", "SAM")

# Remove cells with a total count more than 110000
OX_1E <- subset(OX_1E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
OX_1E[["percent.mt"]] <- PercentageFeatureSet(OX_1E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
OX_1E[["percent.pt"]] <- PercentageFeatureSet(OX_1E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
OX_1E <- subset(OX_1E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
OX_1E <- NormalizeData(OX_1E, verbose = FALSE)

# Find a set of highly avariable genes - 3000 HVGs
OX_1E <- FindVariableFeatures(OX_1E, selection.method = "vst", nfeatures = 3000)

# Extract the count matrix
OX_W1 <- GetAssayData(OX_1E, assay = "RNA", slot = "counts")

# Add a string (identifier) with the barcodes
colnames(OX_W1) <- paste("O1", colnames(OX_W1), sep = "_")


###
# OX - 2 E
###

# First replicate - OX 2E - total cells 10760; filter out genes that are not detected in at least 21 cells
OX_2E <- CreateSeuratObject(counts = OX_DF_2, project = "OX_2E", min.features = 200)

# Add metadata information to the seurat object
OX_2E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Hirsuta", "OX-2", "WT", "SAM")

# Remove cells with a total count more than 110000
OX_2E <- subset(OX_2E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
OX_2E[["percent.mt"]] <- PercentageFeatureSet(OX_2E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
OX_2E[["percent.pt"]] <- PercentageFeatureSet(OX_2E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
OX_2E <- subset(OX_2E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
OX_2E <- NormalizeData(OX_2E, verbose = FALSE)

# Find a set of highly avariable genes - 3000 HVGs
OX_2E <- FindVariableFeatures(OX_2E, selection.method = "vst", nfeatures = 3000)

# Extract the count matrix
OX_W2 <- GetAssayData(OX_2E, assay = "RNA", slot = "counts")

# Add a string (identifier) with the barcodes
colnames(OX_W2) <- paste("O2", colnames(OX_W2), sep = "_")


###
# OX - 3 E
###

# First replicate - OX 3E - total cells 4100; filter out genes that are not detected in at least 8 cells
OX_3E <- CreateSeuratObject(counts = OX_DF_3, project = "OX_3E", min.features = 200)

# Add metadata information to the seurat object
OX_3E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Hirsuta", "OX-3", "WT", "SAM")

# Remove cells with a total count more than 110000
OX_3E <- subset(OX_3E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
OX_3E[["percent.mt"]] <- PercentageFeatureSet(OX_3E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
OX_3E[["percent.pt"]] <- PercentageFeatureSet(OX_3E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
OX_3E <- subset(OX_3E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
OX_3E <- NormalizeData(OX_3E, verbose = FALSE)

# Find a set of highly avariable genes - 3000 HVGs
OX_3E <- FindVariableFeatures(OX_3E, selection.method = "vst", nfeatures = 3000)

# Extract the count matrix
OX_W3 <- GetAssayData(OX_3E, assay = "RNA", slot = "counts")

# Add a string (identifier) with the barcodes
colnames(OX_W3) <- paste("O3", colnames(OX_W3), sep = "_")


###
# Intersection of highly variable genes - Cardamine hirsuta
###

# Lets find common variable features for the replicates of C. hirsuta
ox_hvgs = Reduce(intersect, list(OX_1E@assays$RNA@var.features, OX_2E@assays$RNA@var.features, OX_3E@assays$RNA@var.features))

fileGenerator(ox_hvgs, fileName = "Shared_HVG_between_replicates_CH.txt")

#### SAM
# First replicate - SAM 1E - total cells 4850; filter out genes that are not detected in at least 13 cells
COL_1E <- CreateSeuratObject(counts = COL_SAM, project = "COL_1E", min.features = 200)

# Add metadata information to the seurat object
COL_1E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Thaliana", "WT-COL-1", "WT", "SAM")

# Remove cells with a total count more than 110000
COL_1E <- subset(COL_1E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
COL_1E[["percent.mt"]] <- PercentageFeatureSet(COL_1E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
COL_1E[["percent.pt"]] <- PercentageFeatureSet(COL_1E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
COL_1E <- subset(COL_1E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
COL_1E <- NormalizeData(COL_1E, verbose = FALSE)

# Find a set of highly avariable genes - 2000 HVGs
COL_1E <- FindVariableFeatures(COL_1E, selection.method = "vst", nfeatures = 2000)

# Select HVGs for COL_1E
col_hvgs = VariableFeatures(COL_1E)

fileGenerator(col_hvgs, "Col0_HVG_single_replicate.txt")

# Extract the count table from the seurat object
COL_W1 <- GetAssayData(COL_1E, assay = "RNA", slot = "counts")

# Prepend a distinctive string to the cell labels 
colnames(COL_W1) <- paste("S1", colnames(COL_W1), sep = "_")


###
# Combine the two sets of highly variable genes - Arabidopsis thaliana and Cardamine hirsuta
###

HVGs_combined = union(ox_hvgs, col_hvgs)

fileGenerator(HVGs_combined, "Combined_HVG_intersect_union_approach.txt")

# Lets create the liger object
WT_Species <- createLiger(list(WO1 = OX_W1, WO2 = OX_W2, WO3 = OX_W3, WC1 = COL_W1), remove.missing = F)

# Add metadata information to the liger object for the replicates in the same way as it was added in the seurat object

# Add replicate information
WT_Species@cell.data$Replicates <- WT_Species@cell.data$dataset
WT_Species@cell.data$Replicates <- factor(WT_Species@cell.data$Replicates, levels = c("WO1", "WO2", "WO3", "WC1"), labels = c("WT-OX-1", "WT-OX-2", "WT-OX-3", "WT-COL-1"))

# Add species information
WT_Species@cell.data$Species <- str_sub(WT_Species@cell.data$dataset, 1, nchar(WT_Species@cell.data$dataset) - 1)
WT_Species@cell.data$Species <- factor(WT_Species@cell.data$Species, levels = c("WO", "WC"), labels = c("Hirsuta", "Thaliana"))

# Add genotype information
WT_Species@cell.data$Genotype <- "WT"
WT_Species@cell.data$Genotype <- factor(WT_Species@cell.data$Genotype)

# Normalization of the data
WT_Species@norm.data <- list(WO1 = OX_1E@assays$RNA@data,
                             WO2 = OX_2E@assays$RNA@data, 
                             WO3 = OX_3E@assays$RNA@data, 
                             WC1 = COL_1E@assays$RNA@data)

# Selecting a set of highly variable genes
# WT_Species <- selectGenes(WT_Species, num.genes = 2000, do.plot = FALSE)

WT_Species@var.genes <- HVGs_combined

WT_Species <- scaleNotCenter(WT_Species)

# Check which datasets are we integrating
table(WT_Species@cell.data$dataset)

# Run liger integration - factorization of the matrices
WT_Species <- optimizeALS(WT_Species, k = 50, nrep = 10, lambda = 5)

# Quantile normalization of the data - integration in the shared space
WT_Species <- quantile_norm(WT_Species)

# Run liger implemented UMAP
WT_Species <- runUMAP(WT_Species)

Liger_object_K_50 <- WT_Species

#
save(Liger_object_K_50, file = "integrated_species_SAM_liger.RData")

writeLines(capture.output(sessionInfo()), "Session_info_species_SAM_liger.txt")
