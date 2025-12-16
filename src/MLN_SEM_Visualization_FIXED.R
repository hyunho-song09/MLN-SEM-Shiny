################################################################################
# MLN-SEM Visualization Functions - FIXED VERSION
#
# Publication-ready plots for multi-omics pathway analysis
# FIXES:
#   - Font handling for cross-platform compatibility
#   - Proper ggrepel dependency check
#   - Absolute values for network edge weights
################################################################################

library(ggplot2)
library(ggraph)
library(igraph)
library(gridExtra)
library(RColorBrewer)
library(patchwork)
library(tidyverse)

# Check for optional packages
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

################################################################################
# 1. NETWORK PLOT - FIXED
################################################################################

plot_network <- function(results,
                         layout = "fr",
                         show_labels = TRUE,
                         node_size_range = c(3, 10),
                         edge_width_range = c(0.3, 2),
                         title = "Multi-Layer Network",
                         output_file = NULL) {
  
  cat("Creating network plot...\n")
  
  # Extract modules
  ME_metab <- results$metabolome_modules$MEs
  ME_micro <- results$microbiome_modules$MEs
  
  # Calculate correlations
  cor_mat <- cor(cbind(ME_metab, ME_micro), use = "pairwise.complete.obs")
  
  # Create edge list (only strong correlations)
  edges <- which(abs(cor_mat) > 0.3 & abs(cor_mat) < 1, arr.ind = TRUE)
  
  if(nrow(edges) == 0) {
    cat("No edges found with correlation > 0.3\n")
    return(NULL)
  }
  
  edge_df <- data.frame(
    from = rownames(cor_mat)[edges[, 1]],
    to = colnames(cor_mat)[edges[, 2]],
    weight = abs(cor_mat[edges]),  # CRITICAL: Use absolute values for layout
    sign = sign(cor_mat[edges]),   # Keep sign for reference
    stringsAsFactors = FALSE
  )
  
  # Remove duplicate edges (undirected graph)
  edge_df <- edge_df[!duplicated(t(apply(edge_df[,1:2], 1, sort))), ]
  
  # Create graph
  g <- graph_from_data_frame(edge_df, directed = FALSE)
  
  # Node attributes
  V(g)$type <- ifelse(grepl("^MX_", V(g)$name), "Metabolome", "Microbiome")
  V(g)$size <- degree(g)
  
  # Color palette
  colors <- c("Metabolome" = "#E41A1C", "Microbiome" = "#377EB8")
  
  # Plot
  p <- ggraph(g, layout = layout) +
    geom_edge_link(aes(edge_width = weight, edge_alpha = weight),
                   color = "grey60") +
    geom_node_point(aes(color = type, size = size)) +
    scale_color_manual(values = colors, name = "Layer") +
    scale_size_continuous(range = node_size_range, guide = "none") +
    scale_edge_width_continuous(range = edge_width_range, guide = "none") +
    scale_edge_alpha_continuous(range = c(0.2, 0.8), guide = "none") +
    theme_graph() +
    labs(title = title) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "bottom",
      text = element_text(family = "sans")  # Use safe default font
    )
  
  # Add labels conditionally
  if(show_labels) {
    if(has_ggrepel) {
      p <- p + ggrepel::geom_text_repel(
        aes(x = x, y = y, label = name),
        size = 3,
        max.overlaps = 20
      )
    } else {
      p <- p + geom_node_text(aes(label = name), size = 3, check_overlap = TRUE)
    }
  }
  
  if(!is.null(output_file)) {
    # Use Cairo device for better font handling
    tryCatch({
      cairo_pdf(output_file, width = 10, height = 8)
      print(p)
      dev.off()
      cat(sprintf("Network plot saved: %s\n", output_file))
    }, error = function(e) {
      # Fallback to regular PDF
      ggsave(output_file, p, width = 10, height = 8, dpi = 300, device = "pdf")
      cat(sprintf("Network plot saved (fallback): %s\n", output_file))
    })
  }
  
  return(p)
}

################################################################################
# 2. PATHWAY DIAGRAM
################################################################################

plot_pathways <- function(results,
                          highlight_threshold = 0.01,
                          show_coefficients = TRUE,
                          show_pvalues = TRUE,
                          title = "Significant Pathways",
                          output_file = NULL) {
  
  cat("Creating pathway diagram...\n")
  
  if(nrow(results$passed_pathways) == 0) {
    cat("No significant pathways to plot\n")
    return(NULL)
  }
  
  pathways <- results$passed_pathways
  
  # Prepare data for plotting
  plot_data <- pathways %>%
    mutate(
      pathway = paste(module, "→", outcome),
      significant = p_value < highlight_threshold,
      label = if(show_coefficients) {
        sprintf("β=%.2f\np=%.3f\nF=%.1f", beta, p_value, weak_F)
      } else {
        sprintf("p=%.3f", p_value)
      }
    )
  
  # Main plot
  p <- ggplot(plot_data, aes(x = reorder(pathway, -weak_F), y = weak_F)) +
    geom_col(aes(fill = significant), alpha = 0.8) +
    geom_hline(yintercept = 10, linetype = "dashed", color = "red", size = 0.8) +
    geom_text(aes(label = label), vjust = -0.5, size = 3) +
    scale_fill_manual(
      values = c("FALSE" = "grey60", "TRUE" = "#E41A1C"),
      name = "Significant",
      labels = c("FALSE" = "No", "TRUE" = "Yes")
    ) +
    labs(
      title = title,
      x = "Pathway",
      y = "Weak-F Statistic"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      legend.position = "bottom",
      panel.grid.major.x = element_blank(),
      text = element_text(family = "sans")
    )
  
  if(!is.null(output_file)) {
    ggsave(output_file, p, width = 10, height = 6, dpi = 300, device = cairo_pdf)
    cat(sprintf("Pathway plot saved: %s\n", output_file))
  }
  
  return(p)
}

################################################################################
# 3. COMPREHENSIVE SUMMARY FIGURE
################################################################################

plot_summary <- function(results,
                        output_file = "results_summary.pdf",
                        width = 14, height = 10) {
  
  cat("Creating comprehensive summary figure...\n")
  
  # Panel A: IV analysis overview
  iv_data <- results$iv_table %>%
    filter(!is.na(weak_F)) %>%
    mutate(
      passed = weak_F > 10 & p_value < 0.05 & sargan_p > 0.05,
      category = case_when(
        weak_F > 10 ~ "Strong",
        weak_F > 5 ~ "Moderate",
        TRUE ~ "Weak"
      )
    )
  
  p1 <- ggplot(iv_data, aes(x = weak_F, y = -log10(p_value))) +
    geom_point(aes(color = category, size = abs(beta)), alpha = 0.7) +
    geom_vline(xintercept = 10, linetype = "dashed", color = "red") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    scale_color_manual(
      values = c("Strong" = "#4DAF4A", "Moderate" = "#FF7F00", "Weak" = "#E41A1C"),
      name = "Instrument\nStrength"
    ) +
    scale_size_continuous(range = c(1, 5), name = "|β|") +
    labs(
      title = "A. IV Analysis Overview",
      x = "Weak-F Statistic",
      y = "-log10(P-value)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      legend.position = "right",
      text = element_text(family = "sans")
    )
  
  # Panel B: Pathway coefficients
  if(nrow(results$passed_pathways) > 0) {
    passed_data <- results$passed_pathways %>%
      mutate(pathway = paste(module, "→", outcome))
    
    p2 <- ggplot(passed_data, aes(x = reorder(pathway, beta), y = beta)) +
      geom_col(aes(fill = beta > 0), alpha = 0.8) +
      geom_errorbar(aes(ymin = beta - 1.96*se, ymax = beta + 1.96*se),
                    width = 0.2) +
      coord_flip() +
      scale_fill_manual(
        values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"),
        guide = "none"
      ) +
      labs(
        title = "B. Significant Pathway Coefficients",
        x = "Pathway",
        y = "Effect Size (β ± 95% CI)"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 12),
        axis.text.y = element_text(size = 10),
        text = element_text(family = "sans")
      )
  } else {
    p2 <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No significant pathways",
               size = 6) +
      theme_void()
  }
  
  # Panel C: SEM fit measures
  if(length(results$sem_results) > 0) {
    sem_data <- map_dfr(names(results$sem_results), function(name) {
      sem <- results$sem_results[[name]]
      if(sem$success) {
        fm <- sem$fit_measures
        tibble(
          pathway = name,
          CFI = fm["cfi"],
          RMSEA = fm["rmsea"],
          SRMR = fm["srmr"]
        )
      }
    })
    
    if(nrow(sem_data) > 0) {
      sem_long <- sem_data %>%
        pivot_longer(cols = c(CFI, RMSEA, SRMR),
                     names_to = "measure",
                     values_to = "value")
      
      p3 <- ggplot(sem_long, aes(x = pathway, y = value, fill = measure)) +
        geom_col(position = "dodge") +
        coord_flip() +
        scale_fill_brewer(palette = "Set2", name = "Fit Measure") +
        labs(
          title = "C. SEM Fit Measures",
          x = "Pathway",
          y = "Value"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(face = "bold", size = 12),
          axis.text.y = element_text(size = 10),
          text = element_text(family = "sans")
        )
    } else {
      p3 <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No SEM validation",
                 size = 6) +
        theme_void()
    }
  } else {
    p3 <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No SEM validation",
               size = 6) +
      theme_void()
  }
  
  # Panel D: Module sizes
  module_sizes <- data.frame(
    Layer = c("Metabolome", "Microbiome"),
    Modules = c(ncol(results$metabolome_modules$MEs),
                ncol(results$microbiome_modules$MEs))
  )
  
  p4 <- ggplot(module_sizes, aes(x = Layer, y = Modules, fill = Layer)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = Modules), vjust = -0.5, size = 5) +
    scale_fill_manual(values = c("Metabolome" = "#E41A1C", "Microbiome" = "#377EB8")) +
    labs(
      title = "D. Module Summary",
      x = "",
      y = "Number of Modules"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      legend.position = "none",
      text = element_text(family = "sans")
    )
  
  # Combine panels
  combined <- (p1 + p2) / (p3 + p4) +
    plot_annotation(
      title = "MLN-SEM Analysis Summary",
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        text = element_text(family = "sans")
      )
    )
  
  if(!is.null(output_file)) {
    ggsave(output_file, combined, width = width, height = height, dpi = 300, device = cairo_pdf)
    cat(sprintf("Summary figure saved: %s\n", output_file))
  }
  
  return(combined)
}

################################################################################
# 4. MODULE HEATMAP
################################################################################

plot_module_heatmap <- function(results,
                                layer = "metabolome",
                                top_features = 20,
                                output_file = NULL) {
  
  cat(sprintf("Creating %s module heatmap...\n", layer))
  
  if(layer == "metabolome") {
    membership <- results$metabolome_modules$membership
    prefix <- "MX_"
  } else {
    membership <- results$microbiome_modules$membership
    prefix <- "MZ_"
  }
  
  # Get top features per module
  top_feat <- membership %>%
    enframe(name = "feature", value = "module") %>%
    group_by(module) %>%
    slice_head(n = top_features) %>%
    ungroup() %>%
    mutate(module = paste0(prefix, module))
  
  # Create matrix data
  feat_mod_mat <- table(top_feat$feature, top_feat$module)
  
  # Convert to long format for ggplot
  plot_data <- as.data.frame(feat_mod_mat) %>%
    rename(feature = Var1, module = Var2, value = Freq) %>%
    filter(value > 0)
  
  p <- ggplot(plot_data, aes(x = module, y = feature)) +
    geom_tile(fill = "#377EB8", alpha = 0.8) +
    scale_x_discrete(position = "top") +
    labs(
      title = paste(tools::toTitleCase(layer), "Module Membership"),
      x = "Module",
      y = "Feature"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 0, size = 10),
      axis.text.y = element_text(size = 6),
      panel.grid = element_blank(),
      text = element_text(family = "sans")
    )
  
  if(!is.null(output_file)) {
    ggsave(output_file, p, width = 10, height = 12, dpi = 300, device = cairo_pdf)
    cat(sprintf("Heatmap saved: %s\n", output_file))
  }
  
  return(p)
}

################################################################################
# 5. EFFECT SIZE FOREST PLOT
################################################################################

plot_forest <- function(results,
                       include_all = FALSE,
                       title = "Effect Sizes of Pathways",
                       output_file = NULL) {
  
  cat("Creating forest plot...\n")
  
  if(include_all && nrow(results$iv_table) > 0) {
    plot_data <- results$iv_table %>%
      filter(!is.na(beta), !is.na(se)) %>%
      mutate(
        ci_lower = beta - 1.96 * se,
        ci_upper = beta + 1.96 * se,
        pathway = paste(module, "→", outcome),
        significant = p_value < 0.05
      )
  } else if(nrow(results$passed_pathways) > 0) {
    plot_data <- results$passed_pathways %>%
      mutate(
        ci_lower = beta - 1.96 * se,
        ci_upper = beta + 1.96 * se,
        pathway = paste(module, "→", outcome),
        significant = TRUE
      )
  } else {
    cat("No data to plot\n")
    return(NULL)
  }
  
  p <- ggplot(plot_data, aes(x = reorder(pathway, beta), y = beta)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper, color = significant),
                  width = 0.3, size = 1) +
    geom_point(aes(color = significant, size = weak_F)) +
    coord_flip() +
    scale_color_manual(
      values = c("TRUE" = "#E41A1C", "FALSE" = "grey60"),
      name = "Significant"
    ) +
    scale_size_continuous(range = c(2, 6), name = "Weak-F") +
    labs(
      title = title,
      x = "Pathway",
      y = "Effect Size (β) with 95% CI"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.y = element_text(size = 10),
      legend.position = "bottom",
      text = element_text(family = "sans")
    )
  
  if(!is.null(output_file)) {
    ggsave(output_file, p, width = 10, height = max(6, nrow(plot_data) * 0.3), dpi = 300, device = cairo_pdf)
    cat(sprintf("Forest plot saved: %s\n", output_file))
  }
  
  return(p)
}

################################################################################
# 6. DIAGNOSTIC PLOTS
################################################################################

plot_diagnostics <- function(results,
                            output_file = "diagnostics.pdf") {
  
  cat("Creating diagnostic plots...\n")
  
  iv_data <- results$iv_table %>% filter(!is.na(weak_F))
  
  # Weak-F distribution
  p1 <- ggplot(iv_data, aes(x = weak_F)) +
    geom_histogram(bins = 20, fill = "#377EB8", alpha = 0.7) +
    geom_vline(xintercept = 10, linetype = "dashed", color = "red") +
    labs(title = "A. Weak-F Distribution", x = "Weak-F Statistic", y = "Count") +
    theme_minimal() +
    theme(text = element_text(family = "sans"))
  
  # Sargan p-value distribution
  p2 <- ggplot(iv_data %>% filter(!is.na(sargan_p)),
               aes(x = sargan_p)) +
    geom_histogram(bins = 20, fill = "#4DAF4A", alpha = 0.7) +
    geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
    labs(title = "B. Sargan Test P-values", x = "P-value", y = "Count") +
    theme_minimal() +
    theme(text = element_text(family = "sans"))
  
  # Effect size vs. weak-F
  p3 <- ggplot(iv_data, aes(x = weak_F, y = abs(beta))) +
    geom_point(aes(color = p_value < 0.05), alpha = 0.7, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    scale_color_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "grey60"),
                       name = "Significant") +
    labs(title = "C. Effect Size vs. Instrument Strength",
         x = "Weak-F Statistic", y = "|Effect Size|") +
    theme_minimal() +
    theme(text = element_text(family = "sans"))
  
  # P-value distribution
  p4 <- ggplot(iv_data, aes(x = p_value)) +
    geom_histogram(bins = 20, fill = "#FF7F00", alpha = 0.7) +
    geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
    labs(title = "D. P-value Distribution", x = "P-value", y = "Count") +
    theme_minimal() +
    theme(text = element_text(family = "sans"))
  
  combined <- (p1 + p2) / (p3 + p4) +
    plot_annotation(
      title = "Diagnostic Plots",
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        text = element_text(family = "sans")
      )
    )
  
  if(!is.null(output_file)) {
    ggsave(output_file, combined, width = 12, height = 10, dpi = 300, device = cairo_pdf)
    cat(sprintf("Diagnostic plots saved: %s\n", output_file))
  }
  
  return(combined)
}

################################################################################
# 7. QUICK VISUALIZATION FUNCTION
################################################################################

visualize_all <- function(results,
                         output_dir = "figures/",
                         prefix = "mlnsem") {
  
  cat("\n╔════════════════════════════════════════════════════════════════╗\n")
  cat("║  CREATING ALL VISUALIZATIONS                                   ║\n")
  cat("╚════════════════════════════════════════════════════════════════╝\n\n")
  
  # Create output directory
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat(sprintf("Created directory: %s\n", output_dir))
  }
  
  # Track success/failure
  results_list <- list()
  
  # Generate all plots with error handling
  tryCatch({
    plot_network(results, 
                 output_file = file.path(output_dir, paste0(prefix, "_network.pdf")))
    results_list$network <- "✓"
  }, error = function(e) {
    cat(sprintf("✗ Network plot failed: %s\n", e$message))
    results_list$network <- "✗"
  })
  
  tryCatch({
    plot_pathways(results, 
                  output_file = file.path(output_dir, paste0(prefix, "_pathways.pdf")))
    results_list$pathways <- "✓"
  }, error = function(e) {
    cat(sprintf("✗ Pathway plot failed: %s\n", e$message))
    results_list$pathways <- "✗"
  })
  
  tryCatch({
    plot_summary(results, 
                 output_file = file.path(output_dir, paste0(prefix, "_summary.pdf")))
    results_list$summary <- "✓"
  }, error = function(e) {
    cat(sprintf("✗ Summary plot failed: %s\n", e$message))
    results_list$summary <- "✗"
  })
  
  tryCatch({
    plot_forest(results,
                output_file = file.path(output_dir, paste0(prefix, "_forest.pdf")))
    results_list$forest <- "✓"
  }, error = function(e) {
    cat(sprintf("✗ Forest plot failed: %s\n", e$message))
    results_list$forest <- "✗"
  })
  
  tryCatch({
    plot_diagnostics(results,
                     output_file = file.path(output_dir, paste0(prefix, "_diagnostics.pdf")))
    results_list$diagnostics <- "✓"
  }, error = function(e) {
    cat(sprintf("✗ Diagnostic plots failed: %s\n", e$message))
    results_list$diagnostics <- "✗"
  })
  
  tryCatch({
    plot_module_heatmap(results, layer = "metabolome",
                        output_file = file.path(output_dir, paste0(prefix, "_heatmap_metab.pdf")))
    results_list$heatmap_metab <- "✓"
  }, error = function(e) {
    cat(sprintf("✗ Metabolome heatmap failed: %s\n", e$message))
    results_list$heatmap_metab <- "✗"
  })
  
  tryCatch({
    plot_module_heatmap(results, layer = "microbiome",
                        output_file = file.path(output_dir, paste0(prefix, "_heatmap_micro.pdf")))
    results_list$heatmap_micro <- "✓"
  }, error = function(e) {
    cat(sprintf("✗ Microbiome heatmap failed: %s\n", e$message))
    results_list$heatmap_micro <- "✗"
  })
  
  cat("\n✓ Visualization process complete!\n")
  cat(sprintf("  Output directory: %s\n\n", output_dir))
  
  # Summary
  cat("Status summary:\n")
  for(name in names(results_list)) {
    cat(sprintf("  %s: %s\n", name, results_list[[name]]))
  }
  
  invisible(results_list)
}

cat("MLN-SEM Visualization Functions Loaded! (FIXED VERSION)\n")
cat("Available functions:\n")
cat("  - plot_network(): Network structure\n")
cat("  - plot_pathways(): Significant pathways\n")
cat("  - plot_summary(): Comprehensive results\n")
cat("  - plot_forest(): Effect size forest plot\n")
cat("  - plot_diagnostics(): Diagnostic plots\n")
cat("  - plot_module_heatmap(): Module membership\n")
cat("  - visualize_all(): Generate all plots\n\n")
cat("Fixes applied:\n")
cat("  ✓ Font handling for cross-platform compatibility\n")
cat("  ✓ Cairo device for PDF output\n")
cat("  ✓ Optional ggrepel dependency\n")
cat("  ✓ Error handling in visualize_all()\n\n")
