# library(dplyr)
# library(data.table)
# library(rjson)
# library(jsonlite)

jsonfile_generation <- function(df, output.file){

  df <- df[-1]
  
  assoc_temp <- list(data= NULL, lastPage = NULL)
  assoc_temp$data$variant <- df$snpid
  assoc_temp$data$position <- df$pos
  assoc_temp$data$ref_allele <- df$A2
  # assoc_temp$data$log_pvalue <- df$PRS

  # Flatten the PRS values to ensure they are a single vector
  assoc_temp$data$log_pvalue <- as.vector(df$PRS)
  
  # Convert the list to JSON
  assoc <- toJSON(assoc_temp, auto_unbox = TRUE)

  # assoc <- toJSON(assoc_temp)
  write(assoc, file=output.file)

  
}

