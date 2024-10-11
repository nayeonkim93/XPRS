source("./R/gene_manhattan_plot.R")
source("./R/individual_allele_info_generation.R")
source("./R/snpid_generation.R")
source("./R/jsonfile_generation.R")

which_symbol_is_important_individual <- function(person, test.file, snp, full_PRS, PRS_info, annotated_info, important_genes,output.file) {
  gene_level <- gene_manhattan_plot(person, full_PRS, PRS_info, annotated_info, important_genes)
  temp <- gene_level[[2]]
  temp$Gene <- as.character(temp$Gene)

  # Prepare data for SNP-level manhattan plot
  temp_0 <- temp[, c("Gene")]
  temp_0$Info <- paste0(temp$Chr, ":", temp$Start, ":", temp$End)

  top10 <- temp_0[1:10, ]
  bottom10 <- temp_0[11:20, ]

  write.table(top10, file = './static/data/top10.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(bottom10, file = './static/data/bottom10.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(person, file = './static/data/iid.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)

  allele_info <- individual_allele_info_generation(person, test.file, snp)

  final <- data.table(chr = character(), snpid = character(), pos = numeric(), A1 = character(), A2 = character(), beta = numeric(), PRS = numeric(), Gene = character())

  for (i in 1:nrow(temp)) {
    gene <- temp$Gene[i]
    snps_list <- temp$snps_included[[i]]
    df <- snp[snp$snpid %in% snps_list]
    
    # Ensure beta is numeric
    df$beta <- as.numeric(df$beta)
    
    # Compute PRS using a vectorized approach
    prs_values <- sapply(df$snpid, function(id) {
      if (!is.null(allele_info[[id]]) && is.numeric(allele_info[[id]])) {
        return(df$beta[df$snpid == id] * allele_info[[id]])
      } else {
        return(NA)
      }
    })
    
    df$PRS <- prs_values
    df$Gene <- gene
    final <- rbind(final, df, fill = TRUE)
  }


    final <- final[c(-1),]

  # Check SNPID format for locuszoom
  if (nrow(final[!grepl("[0-9]+:[0-9]+_[A-Z]+/[A-Z]+", final$snpid), ]) > nrow(final) / 2) {
    final <- snpid_generation(final)
  }

  jsonfile_generation(final, output.file)
  return(gene_level[[1]])
}
