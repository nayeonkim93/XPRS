library(ggplot2)
library(plotly)
library(htmltools)
library(colorspace)
library(gridExtra)


which_region_is_important_population <- function(annotated_info){
  df <- annotated_info %>% 
    group_by(region) %>% 
    filter(within_region_rank_by_variance == 1) 
  
  df <- df[order(df$variance, decreasing = TRUE),][c(1:10),]
  
  df$neighboring_symbol <- NA
  df$neighboring_symbol <- sapply(1:10, function(i) setdiff(annotated_info[annotated_info$region == df$region[i],]$symbol, df$symbol[i]))
  
  df$region  <- with(df, reorder(region, variance))
  df$region <- as.factor(df$region)
  names(df)[2] <- "Gene"
  
  p <- ggplot(df, aes_string(x="region", y="variance", fill = "Gene", text= "neighboring_symbol")) +   
    geom_bar(stat = "identity") + coord_flip() + 
    labs(x = "Region",
         y = "Max(variance) in region") + 
    theme_classic() +
    scale_fill_discrete_qualitative(name = "Gene with max(variance)",'Dark 3')

  
  fig <- ggplotly(p)
  # fig
  
  df$neighboring_symbol <- sapply(1:10, function(i) 
    if (identical(df$neighboring_symbol[i][[1]], character(0))){NA} else {df$neighboring_symbol[i]})
  
  df <- df[,c("region","Gene", "chr", "start", "end", "snps_count", "variance", "neighboring_symbol", 'snps_sharing_genes')]
  names(df) <- c("Region", "Gene", "Chr", "Start", "End", "#_Snps", "Gene_Variance", "Neighboring_Gene", "Snps_Sharing_Gene")
  # View(df)

  htmltools::save_html(html = fig, file = "./static/population_region.html")
  
  return(df)
  # return(fig)
}


which_gene_is_important_population_table <- function(annotated_info, important_genes){
  df <- annotated_info %>% 
    group_by(region) %>% 
    filter(within_region_rank_by_variance == 1) 
  
  df <- df[order(df$variance, decreasing = TRUE),][c(1:10),]
  
  df$neighboring_symbol <- NA
  df$neighboring_symbol <- sapply(1:10, function(i) setdiff(annotated_info[annotated_info$region == df$region[i],]$symbol, df$symbol[i]))
  
  df$symbol <- with(df, reorder(symbol, variance))
  df$region <- as.factor(df$region)
  names(df)[2] <- "Gene"
  

  
  df$neighboring_symbol <- sapply(1:10, function(i) 
    if (identical(df$neighboring_symbol[i][[1]], character(0))){NA} else {df$neighboring_symbol[i]})
  
  df <- df[,c("region", "Gene", "chr", "start", "end", "snps_count", "variance", "neighboring_symbol", 'snps_sharing_genes')]
  names(df) <- c("Region", "Gene", "Chr", "Start", "End", "#_Snps", "Gene_Variance", "Neighboring_Gene", "Snps_Sharing_Gene")
  # View(df)

  df$Gene <- as.character(df$Gene)

  for(i in c(1:nrow(df))) {
    need <- append(unlist(df$Neighboring_Gene[i]), unlist(df$Snps_Sharing_Gene[i][[1]]))
    need <- c(df$Gene[i], unlist(need))
    df$Neighboring_Gene[i] <- list(need[!is.na(need)])
  }  

  df <- df[,c("Chr", "Neighboring_Gene",  "Start", "End", "#_Snps", "Gene_Variance")]
  names(df) <- c("Chr", "Genes", "Start", "End", "#_Snps", "Gene_Variance")


for (i in 1:nrow(df)) {
  # Filter genes from df that are in important_genes
  matched_genes <- df$Genes[[i]][df$Genes[[i]] %in% important_genes$MAPPED_GENE]

  if (length(df$Genes[[i]]) > 3) {
    
    if (length(matched_genes) > 0) {
      # Get the p-values for the matched genes and sort by P_VALUE
      gene_pvals <- important_genes[important_genes$MAPPED_GENE %in% matched_genes, ]
      sorted_genes <- gene_pvals[order(gene_pvals$P_VALUE), "MAPPED_GENE"]
      
      # Convert sorted_genes to a vector of gene names
      sorted_genes <- sorted_genes$MAPPED_GENE

      # Keep only the top 3 genes with the lowest p-values
      df$Genes[[i]] <- sorted_genes[1:min(3, length(sorted_genes))]
      
    } else {
      # If no genes are in important_genes, keep the first 3 genes from df
      df$Genes[[i]] <- df$Genes[[i]][1:3]
    }
    
  }
}


  df$Genes <- as.character(df$Genes)
  for(i in c(1:nrow(df))) {
    df$Genes[i] <- gsub('[c()/"]','',df$Genes[i])  
  }


  png("./static/data/population_gene.png", height = 25*nrow(df), width = 150*ncol(df))
  grid.table(df)
  dev.off()
}

which_gene_is_important_population_plot <- function(annotated_info, important_genes){
df <- annotated_info[order(annotated_info$variance, decreasing = T),]

df <-
  df %>% 
  group_by(region) %>% 
  filter(within_region_rank_by_variance == 1)

first_10 <-  df[c(1:10),]

df$neighboring_symbol <- NA
df$neighboring_symbol <- sapply(1:nrow(df), function(i) setdiff(annotated_info[annotated_info$region == df$region[i],]$symbol, df$symbol[i]))


df_manhattan <- df %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(start)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, start) %>%
  mutate( BPcum=start+tot) %>%
  
  # Add highlight and annotation information
  
  mutate(is_highlight = if ((symbol %in% first_10$symbol) & (chr %in% first_10$chr)){'yes_top'} 
         else {'no'}) %>%
  mutate( is_annotate=ifelse(symbol %in% first_10$symbol, "yes", "no")) 


df_manhattan[is.na(df_manhattan$symbol),]$is_highlight = "no"
# print(df_manhattan[is.na(df_manhattan$symbol),]$is_highlight)

# Prepare X axis
axisdf <- df_manhattan %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

for(i in c(1:nrow(df_manhattan))) {
  need <- append(unlist(df_manhattan$neighboring_symbol[i]), unlist(df_manhattan$snps_sharing_genes[i][[1]]))
  need <- c(df_manhattan$symbol[i], unlist(need))
  df_manhattan$symbol[i] <- list(need[!is.na(need)])
}  


for (i in 1:nrow(df_manhattan)) {
  # Filter genes from df that are in important_genes
  matched_genes <- df_manhattan$symbol[[i]][df_manhattan$symbol[[i]] %in% important_genes$MAPPED_GENE]

  if (length(df_manhattan$symbol[[i]]) > 3) {
    
    if (length(matched_genes) > 0) {
      # Get the p-values for the matched genes and sort by P_VALUE
      gene_pvals <- important_genes[important_genes$MAPPED_GENE %in% matched_genes, ]
      sorted_genes <- gene_pvals[order(gene_pvals$P_VALUE), "MAPPED_GENE"]
      
      # Convert sorted_genes to a vector of gene names
      sorted_genes <- sorted_genes$MAPPED_GENE

      # Keep only the top 3 genes with the lowest p-values
      df_manhattan$symbol[[i]] <- sorted_genes[1:min(3, length(sorted_genes))]
      
    } else {
      # If no genes are in important_genes, keep the first 3 genes from df
      df_manhattan$symbol[[i]] <- df_manhattan$symbol[[i]][1:3]
    }
    
  }
}


df_manhattan$symbol <- as.character(df_manhattan$symbol)
for(i in c(1:nrow(df_manhattan))) {
  df_manhattan$symbol[i] <- gsub('c|\\(|\\)|"|list', '',df_manhattan$symbol[i])  
}


# Prepare the text description for each point, including gene, position, and chromosome
df_manhattan$text <- paste("gene: ", df_manhattan$symbol, 
                           "\nposition: ", df_manhattan$start, 
                           "\nchromosome: ", df_manhattan$chr, 
                           "\nvariance: ", df_manhattan$variance)

# Create the plot
p <- ggplot() +
  # Plot all non-highlighted points (grey and blue)
  geom_point(data = df_manhattan[df_manhattan$is_highlight == "no",], 
             aes(x = BPcum, y = variance, color = as.factor(chr), 
                 text = paste("Gene:", symbol, 
                              "\nPosition:", start, 
                              "\nChromosome:", chr, 
                              "\nVariance:", variance)),
             alpha = 0.8, size = 1.3) +
  # Plot the highlighted points ("yes_top") on top, with their own text
  geom_point(data = df_manhattan[df_manhattan$is_highlight == "yes_top",], 
             aes(x = BPcum, y = variance, 
                 text = paste("Gene:", symbol, 
                              "\nPosition:", start, 
                              "\nChromosome:", chr, 
                              "\nVariance:", variance)), 
             color = "#CC3333", size = 2) +
  # Set color palette for chromosomes
  scale_color_manual(values = rep(c("grey", "skyblue"), 22)) +
  # Custom X axis labels
  scale_x_continuous(label = axisdf$chr, breaks = axisdf$center) +
  # Adjust Y axis limits
  scale_y_continuous(expand = c(0, 0)) +
  ylim(min(df_manhattan$variance) - 0.00001, max(df_manhattan$variance) + 0.0001) +
  # Apply theme settings
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  xlab("Chr") + ylab("Variance of Gene Contribution Score")

# Convert to plotly with custom tooltip for each point
fig <- ggplotly(p, tooltip = "text")

# Save the interactive plot
htmltools::save_html(html = fig, file = "./static/data/population_gene.html")

return(fig)


}

