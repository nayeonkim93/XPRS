gene_manhattan_plot <- function(person, full_PRS, PRS_info, annotated_info) {
  '%notin%' <- Negate('%in%')
  
 
  # Convert to data.table if not already
  if (!is.data.table(PRS_info)) setDT(PRS_info)
  if (!is.data.table(annotated_info)) setDT(annotated_info)
  
  print(person)
  PRS_info$Gene <- annotated_info$symbol
  annotated_info$attributed_value <- PRS_info[[person]]

  df <- annotated_info[order(annotated_info$attributed_value, decreasing = T),]
  # Step 1: Clean the symbol in df by removing trailing * but leave the original symbol intact
  # df$clean_symbol <- gsub("\\*+$", "", df$symbol)  # Remove trailing * from df$symbol
  # # Step 2: Merge df with important_genes on the cleaned df$symbol and important_genes$MAPPED_GENE
  # df <- merge(df, important_genes[, c("MAPPED_GENE")], 
  #             by.x = "clean_symbol", by.y = "MAPPED_GENE", all.x = TRUE)
  # # Step 3: Sort first by attributed_value and then by P_VALUE
  df <- df[order(-df$attributed_value), ]
  # df_filtered <- df[!is.na(df$P_VALUE),]
  df_filtered <- df[df$within_region_rank_by_variance == 1, ]
  df_filtered <- df_filtered %>%
  distinct()

  df <- df_filtered


  
  # Select the top and bottom 10 genes
  first_10 <- df[1:10]
  last_10 <- tail(df, 10)
  ls <- c(first_10$symbol, last_10$symbol)
  
  # Prepare data for plotting
  setDT(df)
  df_manhattan <- df[, .(chr, start, attributed_value, symbol)][order(chr, start)]
  df_manhattan[, BPcum := cumsum(as.numeric(start))]
  df_manhattan[, is_highlight := fifelse(symbol %in% first_10$symbol, 'yes_top', fifelse(symbol %in% last_10$symbol, 'yes_bottom', 'no'))]
  df_manhattan[, is_annotate := fifelse(symbol %in% ls, "yes", "no")]
  df_manhattan[, text := paste("gene:", symbol, "\nposition:", start, "\nchromosome:", chr)]

  # Prepare X axis
  axisdf <- df_manhattan[, .(center = mean(BPcum)), by = chr]
  
  # Create the plot
  p <- ggplot() +
    # Plot all points (grey and blue)
    geom_point(data = df_manhattan[is_highlight == "no"], 
              aes(x = BPcum, y = attributed_value, color = as.factor(chr), 
                  text = paste("Gene:", symbol, "\nPosition:", start, "\nChromosome:", chr, "\nAttributed Value:", attributed_value)),
              alpha = 0.8, size = 1.3) +
    # Plot the highlighted points on top, with their own text
    geom_point(data = df_manhattan[is_highlight == "yes_top"], 
              aes(x = BPcum, y = attributed_value, 
                  text = paste("Gene:", symbol, "\nPosition:", start, "\nChromosome:", chr, "\nAttributed Value:", attributed_value)), 
              color = "#CC3333", size = 2) +
    geom_point(data = df_manhattan[is_highlight == "yes_bottom"], 
              aes(x = BPcum, y = attributed_value, 
                  text = paste("Gene:", symbol, "\nPosition:", start, "\nChromosome:", chr, "\nAttributed Value:", attributed_value)), 
              color = "#3366CC", size = 2) +
    scale_color_manual(values = rep(c("grey", "skyblue"), 22)) +
    scale_x_continuous(label = axisdf$chr, breaks = axisdf$center) +
    scale_y_continuous(expand = c(0, 0)) +
    ylim(min(df_manhattan$attributed_value) - 0.0003, max(df_manhattan$attributed_value) +0.0003) +
    # ylim(min(df_manhattan$attributed_value) - 0.0001, max(df_manhattan$attributed_value) + 0.0001)


    theme_bw() +
    theme(
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    xlab("Chromosome") + ylab("Attributed Value of Gene")


  # Convert to plotly with custom tooltip for each type
  fig <- ggplotly(p, tooltip = "text")

  # Save the interactive plot
  htmltools::save_html(fig, file = "./static/data/individual_gene.html")

  
  # Save PRS data
  full_PRS <- as.data.table(full_PRS)
  person_prs <- full_PRS[IID == person, PRS]
  write.table(person_prs, file = "./static/data/prs.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # Calculate the cumulative distribution function (CDF) for the PRS values
  ecdf_PRS <- ecdf(full_PRS$PRS)

  # Compute the percentile of the person_prs
  percentile <- ecdf_PRS(person_prs) * 100

  # Determine whether it's in the top or bottom percentage
  if (percentile > 50) {
    position_text <- paste("Top", round(100 - percentile, 2), "%")
  } else {
    position_text <- paste("Bottom", round(percentile, 2), "%")
  }

  # Get the y-value of the density plot at the position of person_prs
  density_data <- density(full_PRS$PRS)
  density_at_person_prs <- approx(density_data$x, density_data$y, xout = person_prs)$y

  # Define the y-position for the annotation text, increased for larger text size
  text_y_position <- max(density_data$y) + 0.2  # Increase the value to accommodate larger text size

  # Create the density plot with an arrow marking 'person_prs'
  p_density <- ggplot(full_PRS, aes(x = PRS)) +
    geom_density(fill = "white", colour = "black", adjust = 1) +  # Draw the density plot
    theme_bw(base_size = 20) +
    annotate(
      "segment",
      x = person_prs, xend = person_prs, 
      y = 0, yend = text_y_position,  # Extend the arrow up to the new text position
      arrow = arrow(type = "closed", length = unit(0.15, "inches")),
      color = "red"
    ) +
    annotate(
      "text",
      x = person_prs, 
      y = text_y_position,  # Position the text above the plot
      label = position_text,
      color = "red",
      size = 8,  # Increase text size to 8
      vjust = -0.5,  # Adjust the text position slightly above the end of the arrow
      hjust = 0.5  # Center the text horizontally above the arrow
    ) +
    labs(x = "PRS", y = "Density") +  # Add labels for clarity
    ylim(0, max(density_data$y) + 0.7)+
    xlim(min(density_data$x) - 1, max(density_data$x) + 0.3)  # Increase x-axis limits on both sides  # Adjust y-axis limits to provide more space at the top

  # Save the modified plot
  ggsave("./static/data/dist_indi.jpg", plot = p_density)

  
  # Combine first_10 and last_10 into df_summary
  df_summary <- rbindlist(list(first_10, last_10), fill = TRUE)
  
  # Rename the second column to "gene"
  setnames(df_summary, "symbol", "Gene")
  
  # Calculate neighboring_genes
  df_summary[, neighboring_genes := sapply(1:.N, function(i) {
    genes <- setdiff(annotated_info[region == df_summary$region[i], symbol], df_summary$Gene[i])
    if (length(genes) == 0) NA_character_ else paste(genes, collapse = ";")
  })]
  
  # Reorder gene based on attributed_value and add risk column
  df_summary[, gene := factor(reorder(Gene, attributed_value))]
  df_summary[, risk := factor(fifelse(attributed_value > 0, "Increasing a risk", "Decreasing a risk"))]
  
  # Prepare text for tooltips
  df_summary[, text := paste("chr:", chr, "\nneiboring genes:", neighboring_genes, "\nsnp sharing genes:", snps_sharing_genes)]
  
  # Select and rename columns
  df_summary <- df_summary[, .(Region = region, Gene = gene, Chr = chr, Start = start, End = end, 
                               snps_included = snps_included, `#_Snps` = snps_count, Attributed_Value = attributed_value, 
                               Neighboring_Genes = neighboring_genes, Snps_Sharing_Genes = snps_sharing_genes)]
  
  first_10_df <- first_10[,c("chr", "symbol", "start", "end", "snps_count","attributed_value")]
  names(first_10_df) <- c("Chr", "Gene", "Start", "End", "#_Snps", "Attributed_Value")


  # Save the population gene table as an image
  png("./static/data/individual_gene_top.png", height = 25 * nrow(first_10_df), width = 150 * ncol(first_10_df))
  grid.table(first_10_df)
  dev.off()

  last_10_df <- last_10[,c("chr", "symbol", "start", "end", "snps_count","attributed_value")]
  names(last_10_df) <- c("Chr", "Gene", "Start", "End", "#_Snps", "Attributed_Value")


  png("./static/data/individual_gene_bottom.png", height = 25 * nrow(last_10_df), width = 150 * ncol(last_10_df))
  grid.table(last_10_df)
  dev.off()

  # Return both the Manhattan plot and the summary table
  return(list(fig, df_summary))


}
