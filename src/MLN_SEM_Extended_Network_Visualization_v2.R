################################################################################
# MLN-SEM Extended Network Visualization v2.0
#
# Streamlined visualization focusing on:
# 1. Complete validated pathways network
# 2. Correlation vs Causal network comparison with fixed layout
#
# Key improvement: Overlay causal paths on correlation network background
################################################################################

library(ggplot2)
library(ggraph)
library(igraph)
library(tidyverse)
library(RColorBrewer)
library(patchwork)

################################################################################
# HELPER FUNCTIONS
################################################################################

#' Extract all unique nodes from MLN-SEM results
extract_all_nodes <- function(results) {
  # Microbiome modules
  mz_nodes <- colnames(results$microbiome_modules$MEs)
  
  # Metabolome modules  
  mx_nodes <- colnames(results$metabolome_modules$MEs)
  
  # Phenotypes from passed pathways
  if(nrow(results$passed_pathways) > 0) {
    pheno_nodes <- unique(results$passed_pathways$outcome)
  } else {
    pheno_nodes <- character(0)
  }
  
  return(list(
    mz = mz_nodes,
    mx = mx_nodes,
    pheno = pheno_nodes
  ))
}

#' Create node data frame with attributes
create_node_df <- function(nodes_list) {
  node_df <- data.frame(
    name = c(nodes_list$mz, nodes_list$mx, nodes_list$pheno),
    type = c(
      rep("Microbiome", length(nodes_list$mz)),
      rep("Metabolome", length(nodes_list$mx)),
      rep("Phenotype", length(nodes_list$pheno))
    ),
    layer = c(
      rep(1, length(nodes_list$mz)),
      rep(2, length(nodes_list$mx)),
      rep(3, length(nodes_list$pheno))
    ),
    stringsAsFactors = FALSE
  )
  
  return(node_df)
}

#' Calculate correlation matrix for all nodes
calculate_correlation_matrix <- function(results, pheno_data, nodes_list, method = "spearman") {
  # Get module eigengenes
  ME_micro <- results$microbiome_modules$MEs
  ME_metab <- results$metabolome_modules$MEs
  
  # Get phenotype data (matching samples)
  pheno_selected <- pheno_data[, nodes_list$pheno, drop = FALSE]
  
  # Combine all data
  all_data <- cbind(ME_micro, ME_metab, pheno_selected)
  
  # Calculate correlation using specified method (default: Spearman)
  cor_mat <- cor(all_data, use = "pairwise.complete.obs", method = method)
  
  return(cor_mat)
}

################################################################################
# 1. VALIDATED PATHWAYS NETWORK (ALL PASSED PATHWAYS)
################################################################################

#' Create detailed network for ALL validated pathways
#' 
#' @param results MLN-SEM results object
#' @param pheno_data Phenotype data frame
#' @param layout Graph layout algorithm
#' @param show_labels Whether to show node labels
#' @param output_file PDF output file path
#' 
#' @return ggplot object
plot_validated_pathways_network <- function(results,
                                           pheno_data,
                                           layout = "fr",
                                           show_labels = TRUE,
                                           arrow_size = 0.25,
                                           output_file = NULL) {
  
  cat("\n════════════════════════════════════════════════════════════\n")
  cat("Creating network for ALL validated pathways...\n")
  cat("════════════════════════════════════════════════════════════\n\n")
  
  if(nrow(results$passed_pathways) == 0) {
    cat("⚠ No significant pathways found\n")
    return(NULL)
  }
  
  # Get all pathways
  all_paths <- results$passed_pathways
  
  cat(sprintf("Total validated pathways: %d\n", nrow(all_paths)))
  cat(sprintf("  Weak-F range: %.1f - %.1f\n", 
              min(all_paths$weak_F), max(all_paths$weak_F)))
  cat(sprintf("  P-value range: %.4f - %.4f\n\n", 
              min(all_paths$p_value), max(all_paths$p_value)))
  
  # Extract relevant nodes
  mx_nodes <- unique(all_paths$module)
  pheno_nodes <- unique(all_paths$outcome)
  
  # Extract relevant MZ nodes
  mz_nodes <- unique(unlist(strsplit(as.character(all_paths$instruments), ", ")))
  
  # Create edges
  # MZ → MX
  mz_mx_edges <- all_paths %>%
    mutate(instruments_list = strsplit(as.character(instruments), ", ")) %>%
    unnest(instruments_list) %>%
    mutate(
      from = instruments_list,
      to = module,
      weight = weak_F / max(weak_F),  # Normalized
      edge_type = "Instrumental",
      f_stat = weak_F
    ) %>%
    select(from, to, weight, edge_type, f_stat) %>%
    distinct()
  
  # MX → Phenotype
  mx_pheno_edges <- all_paths %>%
    mutate(
      from = module,
      to = outcome,
      weight = abs(beta),
      edge_type = "Causal",
      f_stat = weak_F
    ) %>%
    select(from, to, weight, edge_type, f_stat)
  
  all_edges <- bind_rows(mz_mx_edges, mx_pheno_edges)
  
  # Create node dataframe
  node_df <- data.frame(
    name = c(mz_nodes, mx_nodes, pheno_nodes),
    type = c(
      rep("Microbiome", length(mz_nodes)),
      rep("Metabolome", length(mx_nodes)),
      rep("Phenotype", length(pheno_nodes))
    ),
    stringsAsFactors = FALSE
  )
  
  cat(sprintf("Network composition:\n"))
  cat(sprintf("  MZ nodes: %d\n", length(mz_nodes)))
  cat(sprintf("  MX nodes: %d\n", length(mx_nodes)))
  cat(sprintf("  Phenotype nodes: %d\n", length(pheno_nodes)))
  cat(sprintf("  MZ→MX edges: %d\n", nrow(mz_mx_edges)))
  cat(sprintf("  MX→Phenotype edges: %d\n\n", nrow(mx_pheno_edges)))
  
  # Create graph
  g <- graph_from_data_frame(all_edges, directed = TRUE, vertices = node_df)
  
  # Node attributes
  V(g)$degree <- degree(g, mode = "all")
  
  # Colors
  node_colors <- c(
    "Microbiome" = "#377EB8",
    "Metabolome" = "#E41A1C",
    "Phenotype" = "#4DAF4A"
  )
  
  edge_colors <- c(
    "Instrumental" = "#984EA3",
    "Causal" = "#FF7F00"
  )
  
  # Plot
  p <- ggraph(g, layout = layout) +
    geom_edge_link(
      aes(
        width = weight,
        color = edge_type,
        alpha = weight
      ),
      arrow = arrow(length = unit(arrow_size, "cm"), type = "closed"),
      end_cap = circle(3, "mm")
    ) +
    geom_node_point(
      aes(color = type, size = degree),
      alpha = 0.9
    ) +
    scale_color_manual(values = node_colors, name = "Node Type") +
    scale_edge_color_manual(values = edge_colors, name = "Edge Type") +
    scale_size_continuous(
      range = c(5, 15), 
      name = "Node Degree"  # ADD LEGEND
    ) +
    scale_edge_width_continuous(range = c(0.5, 3), guide = "none") +
    scale_edge_alpha_continuous(range = c(0.4, 1), guide = "none") +
    theme_graph(base_family = "sans") +
    labs(
      title = sprintf("Validated Causal Pathways (n=%d)", nrow(all_paths)),
      subtitle = "Microbiome-instrumented pathways: MZ → MX → Phenotype"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      legend.position = "right"
    )
  
  # Add labels if requested
  if(show_labels) {
    if(requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        aes(x = x, y = y, label = name),
        size = 3.5,
        max.overlaps = 30,
        box.padding = 0.5,
        point.padding = 0.3,
        segment.size = 0.2,
        segment.alpha = 0.5
      )
    } else {
      p <- p + geom_node_text(
        aes(label = name),
        size = 3,
        repel = TRUE,
        check_overlap = TRUE
      )
    }
  }
  
  if(!is.null(output_file)) {
    ggsave(output_file, p, width = 14, height = 10, dpi = 300, device = cairo_pdf)
    cat(sprintf("✓ Validated pathways network saved: %s\n\n", output_file))
  }
  
  return(p)
}

################################################################################
# 2. NETWORK COMPARISON: CORRELATION (LEFT) VS CAUSAL OVERLAY (RIGHT)
################################################################################

#' Create comparison plot with correlation network and causal overlay
#' 
#' Strategy:
#' 1. Build correlation-based network (all nodes involved in passed pathways)
#' 2. Use fixed layout from correlation network
#' 3. Left panel: Standard correlation network
#' 4. Right panel: Grey correlation edges + colored causal paths
#' 
#' @param results MLN-SEM results object
#' @param pheno_data Phenotype data frame
#' @param cor_threshold Correlation threshold for edge inclusion
#' @param layout Layout algorithm (same for both panels)
#' @param show_labels Whether to show node labels
#' @param output_file PDF output file path
#' 
#' @return Combined plot object
plot_network_comparison <- function(results,
                                   pheno_data,
                                   cor_threshold = 0.3,
                                   layout = "fr",
                                   fixed_layout = NULL,
                                   correlation_method = "spearman",
                                   show_labels = TRUE,  # Changed to TRUE
                                   arrow_size = 0.2,
                                   output_file = NULL) {
  
  cat("\n════════════════════════════════════════════════════════════\n")
  cat("Creating Network Comparison: Correlation vs Causal\n")
  cat("════════════════════════════════════════════════════════════\n\n")
  
  if(nrow(results$passed_pathways) == 0) {
    cat("⚠ No significant pathways found\n")
    return(NULL)
  }
  
  # ========================================================================
  # STEP 1: Extract nodes involved in validated pathways
  # ========================================================================
  
  all_paths <- results$passed_pathways
  
  # Get all nodes
  mx_nodes <- unique(all_paths$module)
  pheno_nodes <- unique(all_paths$outcome)
  mz_nodes <- unique(unlist(strsplit(as.character(all_paths$instruments), ", ")))
  
  all_nodes <- c(mz_nodes, mx_nodes, pheno_nodes)
  
  cat(sprintf("Nodes in network: %d (MZ=%d, MX=%d, Pheno=%d)\n",
              length(all_nodes), length(mz_nodes), length(mx_nodes), length(pheno_nodes)))
  
  # ========================================================================
  # STEP 2: Calculate correlation matrix for these nodes
  # ========================================================================
  
  # Get module eigengenes
  ME_micro <- results$microbiome_modules$MEs[, mz_nodes, drop = FALSE]
  ME_metab <- results$metabolome_modules$MEs[, mx_nodes, drop = FALSE]
  pheno_selected <- pheno_data[, pheno_nodes, drop = FALSE]
  
  # Combine all data
  all_data <- cbind(ME_micro, ME_metab, pheno_selected)
  
  # Calculate correlation using specified method (Spearman by default)
  cor_mat <- cor(all_data, use = "pairwise.complete.obs", method = correlation_method)
  
  cat(sprintf("Using %s correlation\n", toupper(correlation_method)))
  
  # Create correlation edge list
  edges_idx <- which(abs(cor_mat) > cor_threshold & abs(cor_mat) < 0.999, arr.ind = TRUE)
  
  if(nrow(edges_idx) == 0) {
    cat(sprintf("⚠ No correlations found above threshold (%.2f)\n", cor_threshold))
    return(NULL)
  }
  
  cor_edge_df <- data.frame(
    from = rownames(cor_mat)[edges_idx[, 1]],
    to = colnames(cor_mat)[edges_idx[, 2]],
    correlation = cor_mat[edges_idx],
    weight = abs(cor_mat[edges_idx]),
    sign = sign(cor_mat[edges_idx]),
    stringsAsFactors = FALSE
  )
  
  # Remove duplicate edges (undirected)
  cor_edge_df <- cor_edge_df[!duplicated(t(apply(cor_edge_df[,1:2], 1, sort))), ]
  
  cat(sprintf("Correlation edges: %d (|r| > %.2f)\n\n", nrow(cor_edge_df), cor_threshold))
  
  # ========================================================================
  # STEP 3: Build causal edges
  # ========================================================================
  
  # MZ → MX edges
  mz_mx_edges <- all_paths %>%
    mutate(instruments_list = strsplit(as.character(instruments), ", ")) %>%
    unnest(instruments_list) %>%
    mutate(
      from = instruments_list,
      to = module,
      edge_type = "Instrumental",
      weight = weak_F / max(weak_F),
      f_stat = weak_F
    ) %>%
    select(from, to, edge_type, weight, f_stat) %>%
    distinct()
  
  # MX → Phenotype edges
  mx_pheno_edges <- all_paths %>%
    mutate(
      from = module,
      to = outcome,
      edge_type = "Causal",
      weight = abs(beta),
      f_stat = weak_F
    ) %>%
    select(from, to, edge_type, weight, f_stat)
  
  causal_edges <- bind_rows(mz_mx_edges, mx_pheno_edges)
  
  cat(sprintf("Causal edges: %d (MZ→MX=%d, MX→Pheno=%d)\n\n",
              nrow(causal_edges), nrow(mz_mx_edges), nrow(mx_pheno_edges)))
  
  # ========================================================================
  # STEP 4: Create node dataframe
  # ========================================================================
  
  node_df <- data.frame(
    name = all_nodes,
    type = c(
      rep("Microbiome", length(mz_nodes)),
      rep("Metabolome", length(mx_nodes)),
      rep("Phenotype", length(pheno_nodes))
    ),
    stringsAsFactors = FALSE
  )
  
  # Colors
  node_colors <- c(
    "Microbiome" = "#377EB8",
    "Metabolome" = "#E41A1C",
    "Phenotype" = "#4DAF4A"
  )
  
  edge_colors_causal <- c(
    "Instrumental" = "#984EA3",
    "Causal" = "#FF7F00"
  )
  
  # ========================================================================
  # STEP 5: Create graphs and get/use fixed layout
  # ========================================================================
  
  # Correlation graph (for layout)
  g_cor <- graph_from_data_frame(cor_edge_df, directed = FALSE, vertices = node_df)
  V(g_cor)$degree <- degree(g_cor)
  
  # Use provided layout or create new one
  if(!is.null(fixed_layout)) {
    cat("Using provided fixed layout\n")
    layout_df <- fixed_layout
  } else {
    cat(sprintf("Creating new layout using '%s' algorithm\n", layout))
    # Get layout
    set.seed(42)
    if(layout == "fr") {
      layout_coords <- layout_with_fr(g_cor)
    } else if(layout == "kk") {
      layout_coords <- layout_with_kk(g_cor)
    } else if(layout == "drl") {
      layout_coords <- layout_with_drl(g_cor)
    } else {
      layout_coords <- layout_with_fr(g_cor)
    }
    
    layout_df <- data.frame(
      name = V(g_cor)$name,
      x = layout_coords[, 1],
      y = layout_coords[, 2]
    )
  }
  
  # ========================================================================
  # STEP 6: LEFT PANEL - Correlation Network
  # ========================================================================
  
  p_left <- ggraph(g_cor, layout = "manual", x = layout_df$x, y = layout_df$y) +
    geom_edge_link(
      aes(
        width = weight,
        alpha = weight,
        color = as.factor(sign)
      )
    ) +
    geom_node_point(
      aes(color = type, size = degree),
      alpha = 0.9
    ) +
    scale_color_manual(values = node_colors, name = "Node Type") +
    scale_edge_color_manual(
      values = c("-1" = "#0571B0", "1" = "#CA0020"),
      name = "Correlation",
      labels = c("-1" = "Negative", "1" = "Positive")
    ) +
    scale_size_continuous(
      range = c(4, 12), 
      name = "Node Degree"  # ADD LEGEND
    ) +
    scale_edge_width_continuous(range = c(0.3, 2), guide = "none") +
    scale_edge_alpha_continuous(range = c(0.2, 0.8), guide = "none") +
    theme_graph(base_family = "sans") +
    labs(
      title = "Correlation Network",
      subtitle = sprintf("%s correlations (|r| > %.2f)", 
                        tools::toTitleCase(correlation_method), cor_threshold)
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "none"
    )
  
  if(show_labels) {
    if(requireNamespace("ggrepel", quietly = TRUE)) {
      p_left <- p_left + ggrepel::geom_text_repel(
        aes(x = x, y = y, label = name),
        size = 2.5,
        max.overlaps = 20
      )
    }
  }
  
  # ========================================================================
  # STEP 7: RIGHT PANEL - Causal Network with Correlation Background
  # ========================================================================
  
  # Create correlation graph for background
  g_cor_right <- graph_from_data_frame(cor_edge_df, directed = FALSE, vertices = node_df)
  
  # Identify nodes in causal pathways
  causal_nodes <- unique(c(causal_edges$from, causal_edges$to))
  
  # Add attributes
  V(g_cor_right)$degree <- degree(g_cor)[V(g_cor_right)$name]
  V(g_cor_right)$in_causal <- V(g_cor_right)$name %in% causal_nodes
  
  # Prepare causal edges with coordinates
  causal_edge_coords <- causal_edges %>%
    left_join(layout_df %>% rename(from = name, x_from = x, y_from = y), by = "from") %>%
    left_join(layout_df %>% rename(to = name, x_to = x, y_to = y), by = "to")
  
  # Create base plot with correlation edges
  p_right <- ggraph(g_cor_right, layout = "manual", x = layout_df$x, y = layout_df$y) +
    # Background: grey correlation edges
    geom_edge_link(
      aes(width = weight, alpha = weight),
      color = "grey80",
      show.legend = FALSE
    ) +
    # Foreground: colored causal edges using geom_segment with ARROWS
    geom_segment(
      data = causal_edge_coords,
      aes(x = x_from, y = y_from, xend = x_to, yend = y_to, 
          color = edge_type, linewidth = weight),
      arrow = arrow(length = unit(arrow_size * 2, "cm"), type = "closed", angle = 20),
      alpha = 0.9,
      inherit.aes = FALSE
    ) +
    # Nodes: grey for non-causal, colored for causal
    geom_node_point(
      aes(
        size = degree,
        color = ifelse(in_causal, type, "grey"),
        alpha = ifelse(in_causal, 0.9, 0.3)
      )
    ) +
    scale_color_manual(
      values = c(node_colors, edge_colors_causal, "grey" = "grey60"),
      name = "Type",
      breaks = c("Microbiome", "Metabolome", "Phenotype", "Instrumental", "Causal")
    ) +
    scale_size_continuous(
      range = c(4, 12), 
      name = "Node Degree"  # ADD LEGEND
    ) +
    scale_edge_width_continuous(range = c(0.2, 2), guide = "none") +
    scale_linewidth_continuous(range = c(0.5, 2.5), guide = "none") +
    scale_alpha_identity() +
    theme_graph(base_family = "sans") +
    labs(
      title = "Causal Pathways (on Correlation Background)",
      subtitle = sprintf("Validated causal paths (n=%d) overlaid", nrow(all_paths))
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "right"
    )
  
  if(show_labels) {
    if(requireNamespace("ggrepel", quietly = TRUE)) {
      p_right <- p_right + ggrepel::geom_text_repel(
        aes(x = x, y = y, label = name),
        size = 2.5,
        max.overlaps = 20
      )
    }
  }
  
  # ========================================================================
  # STEP 8: Combine plots
  # ========================================================================
  
  combined <- p_left + p_right +
    plot_layout(ncol = 2, widths = c(1, 1)) +
    plot_annotation(
      title = "Network Comparison: Causality vs Correlation",
      subtitle = "Left: Simple correlations | Right: Causal pathways highlighted on correlation background",
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        text = element_text(family = "sans")
      )
    )
  
  # Save if requested
  if(!is.null(output_file)) {
    ggsave(output_file, combined, width = 20, height = 10, dpi = 300, device = cairo_pdf)
    cat(sprintf("✓ Network comparison saved: %s\n\n", output_file))
  }
  
  cat("════════════════════════════════════════════════════════════\n")
  cat("✓ Network comparison complete!\n")
  cat("════════════════════════════════════════════════════════════\n\n")
  
  # Return both plot and layout for reuse
  return(list(
    plot = combined,
    layout = layout_df
  ))
}

################################################################################
# 3. MASTER VISUALIZATION FUNCTION
################################################################################

#' Generate all extended network visualizations
#' 
#' @param results MLN-SEM results object
#' @param pheno_data Phenotype data frame
#' @param output_dir Output directory for figures
#' @param prefix Filename prefix
#' 
#' @return List of plot objects
visualize_extended_networks <- function(results,
                                       pheno_data,
                                       output_dir = "figures/",
                                       prefix = "mlnsem") {
  
  cat("\n╔════════════════════════════════════════════════════════════════╗\n")
  cat("║  EXTENDED NETWORK VISUALIZATION SUITE v2.0                     ║\n")
  cat("╚════════════════════════════════════════════════════════════════╝\n\n")
  
  # Create output directory
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat(sprintf("✓ Created directory: %s\n\n", output_dir))
  }
  
  plots <- list()
  
  # 1. Validated pathways network (all) - WITH labels by default
  tryCatch({
    plots$validated <- plot_validated_pathways_network(
      results = results,
      pheno_data = pheno_data,
      show_labels = TRUE,  # Show labels by default
      output_file = file.path(output_dir, paste0(prefix, "_validated_pathways.pdf"))
    )
  }, error = function(e) {
    cat(sprintf("✗ Validated pathways network failed: %s\n", e$message))
  })
  
  # 2. Network comparison
  tryCatch({
    result <- plot_network_comparison(
      results = results,
      pheno_data = pheno_data,
      correlation_method = "spearman",
      show_labels = TRUE,  # Show labels
      output_file = file.path(output_dir, paste0(prefix, "_network_comparison.pdf"))
    )
    plots$comparison <- result$plot
    plots$layout <- result$layout  # Save layout for potential reuse
  }, error = function(e) {
    cat(sprintf("✗ Network comparison failed: %s\n", e$message))
  })
  
  cat("\n╔════════════════════════════════════════════════════════════════╗\n")
  cat("║  EXTENDED NETWORK VISUALIZATION COMPLETE!                      ║\n")
  cat("╚════════════════════════════════════════════════════════════════╝\n\n")
  
  cat("Generated visualizations:\n")
  for(name in names(plots)) {
    if(!is.null(plots[[name]])) {
      cat(sprintf("  ✓ %s\n", name))
    }
  }
  cat(sprintf("\n✓ All files saved to: %s\n\n", output_dir))
  
  return(plots)
}

################################################################################
# FOOTER
################################################################################

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║  MLN-SEM Extended Network Visualization v2.0 Loaded!           ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

cat("Available functions:\n")
cat("  • plot_validated_pathways_network() - All validated pathways (with labels)\n")
cat("  • plot_network_comparison()         - Correlation vs Causal overlay\n")
cat("  • visualize_extended_networks()     - Generate all plots\n\n")

cat("Key improvements:\n")
cat("  ✓ Streamlined to 2 core visualizations\n")
cat("  ✓ Causal paths overlaid on correlation background\n")
cat("  ✓ Fixed layout for direct comparison\n")
cat("  ✓ Grey background with colored causal highlights\n")
cat("  ✓ Proper directed arrows for causal edges\n")
cat("  ✓ Node degree legend included\n\n")

cat("Example usage:\n")
cat("  # After running MLN-SEM analysis:\n")
cat("  visualize_extended_networks(results, pheno_liver,\n")
cat("                              output_dir='figures/',\n")
cat("                              prefix='liver_study')\n\n")
