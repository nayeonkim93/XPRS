# Load libraries and source files
# library(dplyr)
# library(readr)
# library(data.table)
# library(parallel)

map_snps <- function(snp, beta.file, annotation.file, cS2G.file, gwas.file, window_size, numCores){
  
  ####### NEED TO ADD ERROR MESSEGE #######
  
  
  # negate function
  '%notin%' <- Negate('%in%')

  
  ########## 1. POSITIONAL MAPPING: snps mapping to gene with position and given window size ##########
  
  # load annotation file
  annotation <- fread(annotation.file)
  system.time(annotation$snps_included <- mclapply(1:nrow(annotation), 
                                                   function(i) if (length(list(snp[(snp$chr == annotation$chr[i]) & (snp$pos > annotation$start[i] - window_size) & (snp$pos < annotation$end[i] + window_size),]$snpid)[[1]]!= 0))
                                                   {snp[(snp$chr == annotation$chr[i]) & (snp$pos > annotation$start[i] - window_size) & (snp$pos < annotation$end[i] + window_size),]$snpid} 
                                                   else {c()}, mc.cores = numCores))
  
  
  print("Done with 1. POSITIONAL MAPPING")
  
  ########## 2. cS2G MAPPING: snps mapping to gene with cS2G (eQTL/Hi-C etc) ##########
  # load cS2G file
  cS2G <- fread(cS2G.file)
  cS2G <- data.frame(cS2G)
  cS2G<- cS2G[cS2G$SNP %in% snp$snpid,]
 
  # remove genes that are not in annotation file
  cS2G <- cS2G[cS2G$GENE %in% annotation$symbol, ] 

  # function to annotate cS2G snps that are not mapped to genes in 1. positional mapping

   cS2G_function <- function(i){

    if (length(annotation[which((annotation$symbol == cS2G$GENE[i]) & (snp[cS2G$SNP[i] == snp$snpid,]$chr == annotation$chr)),]$snps_included[[1]]) == 0) {
      return(list(c(cS2G$SNP[i]), which((annotation$symbol == cS2G$GENE[i]) & (snp[cS2G$SNP[i] == snp$snpid,]$chr == annotation$chr))))
    }
    
    else if ((cS2G$SNP[i] %notin% annotation[which((annotation$symbol == cS2G$GENE[i]) & (snp[cS2G$SNP[i] == snp$snpid,]$chr == annotation$chr)),]$snps_included[[1]]) == TRUE){
      ls <- append(annotation[which((annotation$symbol == cS2G$GENE[i]) & (snp[cS2G$SNP[i] == snp$snpid,]$chr == annotation$chr)),]$snps_included[[1]], cS2G$SNP[i])
      return(list(sort(ls), which((annotation$symbol == cS2G$GENE[i]) & (snp[cS2G$SNP[i] == snp$snpid,]$chr == annotation$chr))))
    }

  }

  
 
  all_list  <- mclapply(1:nrow(cS2G), function (i) cS2G_function(i), mc.cores = numCores)
  # all_list  <- lapply(1:nrow(cS2G), function (i) cS2G_function(i))
  all_list = all_list[-which(sapply(all_list, is.null))]
  snp_list <- lapply(1:length(all_list), function (i) all_list[[i]][[1]])
  index_list <- unlist(lapply(1:length(all_list), function (i) all_list[[i]][[2]]))
  annotation[index_list,]$snps_included <- snp_list
  
  
  ########## post processing ##########
  
  # function to remove the genes that none of the snp is assigned 
  remove_zero <- function(i){
    if (is.null(annotation$snps_included[[i]]) == FALSE){
      return(i)
    }
  }
  
  # remove the genes that none of the snp is assigned
  ls <- unlist(lapply(1:nrow(annotation), function(i) remove_zero(i)))
  annotation <- annotation[ls,]
  
  print("Done with 2. cS2G MAPPING")
  
  # annotation$snps_included to character 
  annotation$snps_included <- as.character(annotation$snps_included)
  
  
  ########## case1: "fully overlap" grouping ##########
  # leave only 1 fully overlap gene to avoid duplicated calculation of gene PRS
  
  # one
  annotation_unique <- annotation[!duplicated(annotation$snps_included),] 
  # rest
  annotation_not_unique <- annotation[annotation$symbol %notin% annotation_unique$symbol,]
  
  # rest are assigned to annotation_unique$snps_sharing_genes 
  annotation_unique$snps_sharing_genes <- NA
  annotation_unique$snps_sharing_genes <- mclapply(1:nrow(annotation_unique), 
                                                   function (i) annotation_not_unique[annotation_not_unique$snps_included == annotation_unique$snps_included[i],]$symbol,
                                                   mc.cores = numCores)
  
  remove_character0 <- function(i){
    if ( identical(annotation_unique$snps_sharing_genes[[i]],character(0)) == TRUE){
      return(i)
    }
  }
  
  ls <- unlist(lapply(1:nrow(annotation_unique), function(i) remove_character0(i)))
  annotation_unique[ls,]$snps_sharing_genes <- NA
  
  print("Done with 'case 1: fully overlapped genes' grouping")
  

  ########## case2: "partially overlap" regionizing ##########
  annotation_unique$snps_included  <- strsplit(gsub('[c()\n "]', '', annotation_unique$snps_included),",")
  annotation_unique <- data.frame(annotation_unique)
  annotation_unique$region <- NA
  annotation_unique$snps_count <- sapply(1:nrow(annotation_unique), function(i) length(unlist(annotation_unique$snps_included[i])) ) 
  
  # Finding out the index of partially overlap
  final <- data.frame(matrix(ncol = 8, nrow = 0))
  colnames(final) <- names(annotation_unique)
  intersect_func <- function(x, y) length(intersect(unlist(x), unlist(y)))
  count <- 1
  system.time(for (c in 1:22){
    temp <- annotation_unique[annotation_unique$chr == c,]
    print(paste("Regionizing chromosome", c))
    m <- matrix(0, nrow= nrow(temp), ncol = nrow(temp))
    for (i in c(1:nrow(temp))){
      for (j in c(i:nrow(temp))){
        m[i,j] <- intersect_func(temp$snps_included[i], temp$snps_included[j])
      }
    }
    
    m[lower.tri(m)]<- t(m)[lower.tri(m)]
    
    ls_gene <- temp$symbol
    
    while (length(ls_gene) != 0){
      # maximum number of snps assigned to a single gene
      max_num <- max(temp[temp$symbol %in% ls_gene,]$snps_count)
      # gene index where maximum number of snps 
      num <- which(temp$snps_count == max_num)
      # since it can have several gene has maximum number of snps go through loop
      if (length(num >= 1)){
        for (i in c(1:length(num))){
          if (nrow(temp[which(m[num[i],] >= floor(max(max_num) *2/3)),][is.na(temp[which(m[num[i],] >= floor(max(max_num) *2/3)),]$region)== T,]) >= 1){
            temp[which(m[num[i],] >= floor(max(max_num) *2/3)),][is.na(temp[which(m[num[i],] >= floor(max(max_num) *2/3)),]$region)== T,]$region <- count
            ls_gene <-  ls_gene[ls_gene %notin% temp[which(m[num[i],] >= floor(max(max_num) *2/3)),]$symbol]
            count <- count + 1
            #print(count)
            #print(length(ls_gene))
          }
        }
      }
    }
    
    final <- rbind(final, temp)
  } )
  
  print("Done with 'case 2: partially overlapped genes' regionizing")
  
  ########## 3. GWAS p-value based mapping OR Snp heritability-based mapping  ##########


  ########## 3-1. GWAS p-value based mapping: snps mapping to gene with position and GWAS p-value ##########
  
  # load GWAS file 
  gwas <- fread(gwas.file)
  
  # snp that are left and not annotated to genes
  snp_not <- snp[snp$snpid %notin% unique(unlist(final$snps_included))]
  snp_not <- left_join(snp_not, gwas, by = c("snpid" = "SNPID"))
  snp_not <- snp_not[snp_not$chr == snp_not$CHR, ]
  snp_not <- snp_not[snp_not$pos == snp_not$POS, ]
  
  final <- as.data.frame(final)
  
  i = nrow(final) + 1
  
  count = max(final$region) + 1
  j = 1
  while (nrow(snp_not) != 0) {
    
    ls_2 <- c() # empty list
    min_snp <- snp_not[snp_not$P_VALUE == min(snp_not$P_VALUE),]
    snp_not_new <- snp_not[(snp_not$chr == min_snp$chr) & (snp_not$pos > min_snp$pos - window_size) & (snp_not$pos < min_snp$pos + window_size),]
    ls_2 <- append(ls_2, snp_not_new$snpid)
    final[i,] <- NA
    final$chr[i] <- min_snp$chr
    final$start[i] <- min_snp$pos - window_size
    final$end[i] <- min_snp$pos + window_size
    final$snps_included[i] <- list(snp_not_new$snpid)
    final$region[i] <- count
    final$symbol[i] <- paste("NA_", j, sep = "")  # Add this line
    print(final$symbol[i])

    snp_not <- snp_not[snp_not$snpid %notin% ls_2,]
    #print(paste("left over snps",nrow(snp_not)))
    #print(i)
    count = count + 1
    i = i + 1 
    j = j + 1
  }
  
  print("Done with 3-1. GWAS mapping")

  
  final$snps_count <- sapply(1:nrow(final), function(i) length(unlist(final$snps_included[i])) ) 
  
  final <- final[,c("chr", "symbol", "snps_sharing_genes", "start", "end", "snps_included", "snps_count", "region")]
  
  return(final)
  
}

