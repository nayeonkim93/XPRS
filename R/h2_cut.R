source("./R/match_alleles.R")

h2_cal_cut <- function(test.file, beta.file, cut_off){

    snp <- fread(beta.file)
    write.table(snp$snpid, file = "./data/snp.txt", row.names = F, col.names=F, quote = F)
    cmd = paste0("./plink --bfile ",test.file, " --extract ./data/snp.txt --make-bed --out ./data/temp_test")
    system(cmd)


    system("cut -f2 ./data/temp_test.bim | sort | uniq -d > ./data/duplicates.txt")
    cmd <- paste0("./plink --bfile ./data/temp_test --exclude ./data/duplicates.txt --make-bed --out ./data/temp_test_no_dup")
    system(cmd)

    #########
    cmd=paste0("./plink --bfile ./data/temp_test_no_dup --freqx --out ./output/study_allele_freqx")
    system(cmd)
    freq <- fread("./output/study_allele_freqx.frqx")
    names(freq) <- c("CHR", "SNP", "A1", "A2", "Hom_A1", "Het", "Hom_A2", "Hap_A1", "Hap_A2", "Missing")

    # Calculate Total Alleles
    freq[, total_alleles := 2 * (Hom_A1 + Het + Hom_A2)]

    # Calculate Allele Frequencies
    freq[, A1_freq := (2 * Hom_A1 + Het) / total_alleles]
    freq[, A2_freq := (2 * Hom_A2 + Het) / total_alleles]

    # Calculate Minor Allele Frequency (MAF)
    freq[, MAF := pmin(A1_freq, A2_freq)]

    ######

    # Rename 'snpid' to 'SNP' in snp
    setnames(snp, old = "snpid", new = "SNP")

    # Optionally, ensure that chromosome columns have consistent naming
    # Rename 'chr' to 'CHR' in snp to match freq
    setnames(snp, old = "chr", new = "CHR")
    
    # Merge freq and snp on 'SNP'
    merged_dt <- merge(freq, snp, by = "SNP", suffixes = c("_freqdt", "_snp"), all = FALSE) 

    merged_dt[, Allele_Match := (A1_freqdt == A1_snp & A2_freqdt == A2_snp) |
                              (A1_freqdt == A2_snp & A2_freqdt == A1_snp)]
    final_dt <- merged_dt[Allele_Match == TRUE]
    final_dt
    
    final_dt <- final_dt[, .(
        CHR_snp,
        SNP,
        pos,
        A1_snp,
        A2_snp,
        A1_freq,
        A2_freq,
        MAF,
        beta
    )]

    names(final_dt) <- c("chr", "snpid", "pos", "A1", "A2", "A1_freq", "A2_freq","MAF", "beta")


    final_dt$h2 <- NA
    final_dt$h2 <- (final_dt$beta)^2*(final_dt$MAF)*(1-final_dt$MAF)

    if (nrow(final_dt) > 100000){
        final_dt <- final_dt[order(-snp$h2),]
        number <- floor(cut_off* nrow(snp))
        number <- max(number, 100000) 
        final_dt <- final_dt[c(1:number),]
        final_dt <- final_dt[order(chr, snpid),]
    }


    final <- match_alleles("./data/temp_test_no_dup", final_dt)
    # after cut_off make new_test file
    write.table(final$snpid, file = "./data/snp_cut_off.txt", row.names = F, col.names=F, quote = F)
    cmd = paste0("./plink --bfile ./data/temp_test_no_dup --extract ./data/snp_cut_off.txt --make-bed --out ./data/new_test")
    system(cmd)

    system("rm -rf ./data/snp.txt")
    system("rm -rf ./data/snp_cut_off.txt")
    system("rm -rf ./data/temp_test.*")
    system("rm -rf ./data/temp_test_no_dup*")
    system("rm -rf ./output/study_allele_freqx.*")
    system("rm -rf ./data/duplicates.txt")
    

    return(final)
}