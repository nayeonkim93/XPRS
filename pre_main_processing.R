start_time_T <- Sys.time()
# Load libraries and source files
library(dplyr)
library(readr)
library(data.table)
library(Rcpp)
library(optparse)
library(parallel)

option_list = list(
  make_option(c("-t", "--test.file"), type="character", default=NULL, 
              help="full path to population.genomic.file", metavar="file"),
  make_option(c("-b", "--beta.file"), type="character", default=NULL, 
              help="full path to PRS.scroing.file", metavar="file"),
  make_option(c("-g", "--GWAS.file"), type="character", default=NULL, 
              help="full path to sumstat.file", metavar="file"),
  make_option(c("-a", "--annotation.file"), type="character", default="./data/ucsc_hg19.txt", 
               help="full path to annotation.file", metavar="file"),
  make_option(c("-c", "--cS2G.file"), type="character", default="./data/cS2G_annotation.txt", 
               help="full path to cS2G.file", metavar="file"),
  make_option(c("-i", "--important_genes.file"), type="character", default= NULL, 
               help="full path to GWAS.association.file", metavar="file"),      
  make_option(c("-w", "--window.size"), type="integer", default=200000, 
               help="specify the window size for gene mapping", metavar="number"),
  make_option(c("-n", "--number.cpu"), type="integer", default=8, 
               help="specify the CPU cores you are going to use", metavar="number"),
  make_option(c("-k", "--h2.cut"), type="double", default= 0.5, 
               help="specify the snp heritability cut off value", metavar="number")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$test.file)){
  print_help(opt_parser)
  stop("An arugment must be supplied for population.genomic.file", call.=FALSE)
}

if (is.null(opt$beta.file)){
  print_help(opt_parser)
  stop("An arugment must be supplied for PRS.scoring.file", call.=FALSE)
}

if (is.null(opt$important_genes.file)){
  print_help(opt_parser)
  stop("An arugment must be supplied for GWAS.association.file", call.=FALSE)
}


source("./R/h2_cut.R")
source("./R/map_snps.R")
source("./R/map_snps_no_gwas.R")
sourceCpp("./C++/BedFileReader.cpp")
setClass("BedFileReader", representation( pointer = "externalptr" ) )
BedFileReader_method <- function(name) {paste( "BedFileReader", name, sep = "__" ) }
setMethod( "$", "BedFileReader", function(x, name ) {function(...) .Call(BedFileReader_method(name),x@pointer, ... )})
setMethod("initialize","BedFileReader", function(.Object, ...) {
  .Object@pointer <-.Call(BedFileReader_method("new"), ... )
  .Object})



################################# INPUT FILE #################################

##### IMPORTANT: snpid has to be rsid, if not, you should convert cS2G file to your formatted snpid #####

##### You can change these files #####
# 1. population.genomic.file (plink file)
test.file <- opt$test.file

# 2. PRS.scoring.file : beta obatained from other PRS model such as lassosum/PRS-CS
# header should be chr snpid pos A1 A2 beta
beta.file <- opt$beta.file

# 3. GWAS summary stat :gwas summary statistics that is used to construct PRS.scoring.file
# header should be CHR SNPID POS A1 A2 P.VALUE
gwas.file <- opt$GWAS.file

# 4. GWAS.association.file : gwas association file that is downloaded from GWAS catalog associated with phenotype 
# header should be MAPPED GENES P_VALUE
important.genes.file <- opt$important_genes.file

##### Given (but if you want you can change as format given) #####
# annotation file: positioning gene mapping from UCSC hg38 reference genome
annotation.file <- opt$annotation.file

# cS2G file: combine Snp to Gene annotation file from the paper
cS2G.file <- opt$cS2G.file

##### Parameters you can change #####
# Window size: default window size is 500kb
window_size <- opt$window.size

# number of CPU cores : default CPU cores using are 8
numCores <- opt$number.cpu

# snp heritability cut off : default cut off value using is 1.0e-15
cut_off <- opt$h2.cut

# Before doing so we need to cut by snp heritability
# this function include snp match which match a1, a2 allele in test file and beta file
snp <- h2_cal_cut(test.file, beta.file, cut_off)
# snp$beta <- (snp$beta - mean(snp$beta)) / sd(snp$beta)
test.file <- "./data/new_test"

# C++ constructor
BedFileReader <- new( "BedFileReader", paste0(test.file,".fam"), paste0(test.file,".bim"), paste0(test.file,".bed"))

# input
print(paste0("test.file: ", test.file))
print(paste0("beta.file: ", beta.file))
print(paste0("GWAS.file: ", gwas.file))
print(paste0("important.genes.file: ", important.genes.file))
print(paste0("annotation.file: ", annotation.file))
print(paste0("cS2G.file: ", cS2G.file))

# parameters 
print(paste0("window size: ", window_size))
print(paste0("number of CPU: ", numCores))
print(paste0("snp heritability cut off value: ", cut_off))

################################# pre-processing #################################
start_time <- Sys.time()
# snp mapping (1. positional mapping, 2. cS2G mapping, 3. GWAS p-value mapping/ snp-heritability mapping)
if (is.null(gwas.file) == TRUE){
    print("we DON'T have gwas.file so we are going to map based on snp heritability for leftover snps")
   annotated_info <- map_snps_no_gwas(snp, beta.file, annotation.file, cS2G.file, test.file, window_size, numCores)

}
if (is.null(gwas.file) == FALSE){
    print("we have gwas.file so we are going to map based on gwas p-value for leftover snps")
    annotated_info <- map_snps(snp, beta.file, annotation.file, cS2G.file, gwas.file, window_size, numCores)
}

end_time <- Sys.time()
time_taken <- end_time - start_time
time_taken_minutes <- as.numeric(time_taken, units = "mins")
print(paste0("Time it takes to pre-processing for population is ", time_taken_minutes, " minutes"))

# head(annotated_info)

################################# main-processing : PRS calculation  #################################
# load fam file 
fam <- fread(paste0(test.file,".fam"))

# calculate PRS to full_PRS
start_time <- Sys.time()
full_PRS <- BedFileReader$calculateWholePRS(snp$snpid,snp$beta, snp$flip)


snp$beta <- snp$beta * sd(full_PRS)
full_PRS <- (full_PRS - mean(full_PRS)) / sd(full_PRS)  # standardized PRS
end_time <- Sys.time()
print(paste0("mean of PRS is ", mean(full_PRS)))
print(paste0("sd of PRS is ", sd(full_PRS)))
print(paste0("Time it takes to calculate PRS for population is ",end_time - start_time, " seconds"))

# calculate Gene Contribution Score to PRS_info

snp_map <- data.table(
  snpid = snp$snpid,
  beta = snp$beta,
  flip = snp$flip,  # Ensure this column exists
  AF = snp$AF
)
setkey(snp_map, snpid)

# Ensure the sequence is an integer sequence
num_rows <- nrow(annotated_info)
sequence <- seq(1, num_rows)

start_time <- Sys.time()
PRS_info <- mclapply(sequence, function(i) {{
    snp_ids <- annotated_info$snps_included[[i]]
    # betas <- snp_map[J(snp_ids), beta, nomatch = 0]
    # Retrieve beta and flip values for these SNPs from snp_map
    snp_subset <- snp_map[J(snp_ids), .(beta, flip, AF), nomatch = 0]
    # Handle cases where no SNPs are found
    if (nrow(snp_subset) == 0) {
        return(0)  # Or another appropriate default value
    }

    # Extract beta and flip lists
    betas <- snp_subset$beta
    flips <- snp_subset$flip
    AFs <- snp_subset$AF


    # Call the C++ function to calculate PRS with flipping
    PRS <- BedFileReader$calculatePRS(snp_ids, betas, flips, AFs)
    
    return(PRS)
}}, mc.cores = numCores, mc.set.seed = TRUE)
end_time <- Sys.time()
print(paste0("Time it takes to calculate Gene Contribution Score for individual is ", end_time - start_time, " seconds"))


PRS_info_matrix <- do.call(cbind, PRS_info)
# Calculate the mean PRS for each gene (column means)
col_means <- colMeans(PRS_info_matrix, na.rm = TRUE)
# Calculate the variance for each gene (optional)
col_variance <- apply(PRS_info_matrix, 2, var, na.rm = TRUE)
PRS_info_matrix<- sweep(PRS_info_matrix, 2, col_means, "-")
PRS_info <- as.data.frame(t(PRS_info_matrix))
colnames(PRS_info) <- fam$V2
rownames(PRS_info) <- NULL
PRS_info$Gene <- rep(annotated_info$symbol, length.out = nrow(PRS_info))


annotated_info$mean <- col_means
annotated_info$variance <- col_variance

# calculate the rank with in region by variance 
annotated_info <- annotated_info %>%
  group_by(region) %>%
  mutate(within_region_rank_by_variance = order(order(variance, decreasing=TRUE)))

# organize FULL_PRS and save important_genes
FULL_PRS <- data.frame(IID = fam$V2, PRS = full_PRS)
important_genes <- fread(important.genes.file)

# save all information to data.RDS
list_of_objects <- list(annotated_info = annotated_info, full_prs = FULL_PRS, PRS_info = PRS_info, snp_info = snp, important_genes = important_genes)

################################# post-processing : population-based  #################################
source("./R/which_symbol_important_population.R")

# 1. Within population: which region/gene is important?
start_time <- Sys.time()
fig <- which_gene_is_important_population_plot(annotated_info, important_genes)
which_gene_is_important_population_table(annotated_info, important_genes)

list_of_objects$plot <- fig
image_binary <- readBin("./static/data/population_gene.png", "raw", file.info("./static/data/population_gene.png")$size)
list_of_objects$table <- image_binary
list_of_objects$window_size <- window_size
saveRDS(list_of_objects, "./output/data.RDS")
end_time <- Sys.time()
print(paste0("Time it takes for post-processing for population is ",end_time - start_time, " seconds"))
end_time_T <- Sys.time()
time_taken_T <- end_time_T - start_time_T
time_taken_minutes_T <- as.numeric(time_taken_T, units = "mins")
print(paste0("Total time it takes is ", time_taken_minutes_T, " minutes"))
