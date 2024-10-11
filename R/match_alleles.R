match_alleles <- function(test.file, snp){
    # read plink bim file and beta file
    bim <- fread(paste0(test.file,".bim"))

    # Check where beta is negative
    negative_beta <- snp$beta < 0

    # Flip A1 and A2 for rows where beta is negative
    snp[negative_beta, c("A1", "A2")] <- snp[negative_beta, c("A2", "A1")]
    snp[negative_beta, "beta"] <- snp[negative_beta, beta ] * (-1)

    # merge the bim file and beta file then match the allele A1/A2
    merge <- left_join(bim, snp,  by = c("V2" = "snpid"))
    merge <- merge[is.na(merge$chr) == F,]
    merge <- merge[(merge$V1 == merge$chr),]
    merge <- merge[(merge$V4 == merge$pos),]

    # if bim A1 == PRS beta A2, bim A2 == PRS beta A1, flip = 1, if same 0, if doesn't match NA
    merge$flip <- ifelse(merge$V5 == merge$A2 & merge$V6 == merge$A1, 1, 
                            ifelse(merge$V5 == merge$A1 & merge$V6 == merge$A2, 0, NA))
    # remove the snps with flip == NA
    merge <- merge[!is.na(merge$flip), ]
    # if flip == 1, take A2 frequency. else A1
    merge[, AF := ifelse(flip == 1, A2_freq, A1_freq)]

    # save beta file and return it 
    merge <- merge[,c("V1", "V2", "V4", "A1", "A2", "beta", "AF", "flip", "h2")]
    names(merge) <- c("chr", "snpid", "pos", "A1", "A2", "beta", "AF","flip", "h2")

    print("finished matching alleles")

    return(merge)
}