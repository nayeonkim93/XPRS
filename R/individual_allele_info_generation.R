library(hash)


sourceCpp("./C++/BedFileReader_indi.cpp")
setClass("BedFileReader_indi", representation(pointer = "externalptr"))
BedFileReader_indi_method <- function(name) {paste("BedFileReader_indi", name, sep = "_")}
setMethod("$", "BedFileReader_indi", function(x, name) {function(...) .Call(BedFileReader_indi_method(name), x@pointer, ...)})
setMethod("initialize", "BedFileReader_indi", function(.Object, ...) {
    .Object@pointer <- .Call(BedFileReader_indi_method("new"), ...)
    .Object
})

individual_allele_info_generation <- function(person, test.file, snp) {

  BedFileReader <- new("BedFileReader_indi", paste0(test.file, ".fam"), paste0(test.file, ".bim"), paste0(test.file, ".bed"))
  bim <- fread(paste0(test.file, ".bim"))
  allele <- BedFileReader$readOneIndi(person)

  # Selective SNP index
  selective_snp_index <- match(snp$snpid, bim$V2)
  allele <- allele[selective_snp_index]

  # Ensure that 'allele' and 'snp' have the same length
  if(length(allele) != nrow(snp)) {
    stop("Length of 'allele' vector and number of rows in 'snp' data.table must match.")
  }

  # Ensure that 'flip' column is integer
  snp[, flip := as.integer(flip)]

  # Create a copy of the allele vector to preserve the original (optional)
  allele_flipped <- allele

  # Flip alleles where flip == 1
  allele_flipped[snp$flip == 1] <- 2 - allele_flipped[snp$flip == 1]

  #adjust the allele_flipped by substracting allele frequency
  allele_flipped <- allele_flipped - snp$AF

  # Save it to a data.table for efficient storage and retrieval
  allele_info <- data.table(snpid = snp$snpid, allele = allele_flipped)

  # Convert to hash for compatibility 
  allele_hash <- hash()
  for (i in seq_len(nrow(allele_info))) {
    allele_hash[[allele_info$snpid[i]]] <- allele_info$allele[i]
  }

  return(allele_hash)
}
