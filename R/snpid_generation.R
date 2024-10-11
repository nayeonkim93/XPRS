# to make the snpid as locuszoom style : {chrom}:{pos}_{ref}/{alt}

# A1 in bim is alternative
# A2 in bim is reference (not a risk allele)

snpid_generation <- function(df){
  
  df$snpid <- paste0(df$chr,":",df$pos,"_",df$A2,"/",df$A1)
  
  return(df)
  
}

