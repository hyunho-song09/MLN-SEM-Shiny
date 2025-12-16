################################################################################
# MLN-SEM Simple Module Assignment Extractor
#
# Extract feature-to-module assignments (simplified version)
# Uses membership information from MLN-SEM results
################################################################################

library(openxlsx)
library(dplyr)

################################################################################
# MAIN EXTRACTION FUNCTIONS
################################################################################

#' Extract microbiome module assignments (simple)
#' 
#' @param results MLN-SEM results object
#' @param microbiome Original microbiome data
#' 
#' @return Data frame with feature and module assignment
extract_microbiome_modules_simple <- function(results, microbiome) {
  
  cat("\n════════════════════════════════════════════════════════════\n")
  cat("Extracting Microbiome Module Assignments\n")
  cat("════════════════════════════════════════════════════════════\n\n")
  
  # Get microbiome features (columns starting with "s_")
  micro_cols <- grep("^s_", colnames(microbiome), value = TRUE)
  
  # Get module membership
  if(!is.null(results$microbiome_modules$membership)) {
    membership <- results$microbiome_modules$membership
    cat("✓ Found membership information\n")
  } else {
    stop("Cannot find membership in results$microbiome_modules")
  }
  
  # Create simple assignment table
  result_df <- data.frame(
    Feature = micro_cols,
    Module = membership,
    stringsAsFactors = FALSE
  )
  
  # Sort by module
  result_df <- result_df %>%
    arrange(Module, Feature)
  
  # Summary statistics
  cat("\n--- Summary ---\n")
  module_summary <- result_df %>%
    group_by(Module) %>%
    summarise(N_Features = n()) %>%
    arrange(desc(N_Features))
  
  print(module_summary)
  
  cat(sprintf("\nTotal features: %d\n", nrow(result_df)))
  cat(sprintf("Total modules: %d\n", length(unique(membership))))
  cat("\n✓ Microbiome module assignments extracted\n")
  
  return(result_df)
}

#' Extract metabolome module assignments (simple)
#' 
#' @param results MLN-SEM results object
#' @param metabolome Original metabolome data
#' 
#' @return Data frame with feature and module assignment
extract_metabolome_modules_simple <- function(results, metabolome) {
  
  cat("\n════════════════════════════════════════════════════════════\n")
  cat("Extracting Metabolome Module Assignments\n")
  cat("════════════════════════════════════════════════════════════\n\n")
  
  # Get metabolome features (columns starting with "M_")
  metab_cols <- grep("^M_", colnames(metabolome), value = TRUE)
  
  # Get module membership
  if(!is.null(results$metabolome_modules$membership)) {
    membership <- results$metabolome_modules$membership
    cat("✓ Found membership information\n")
  } else {
    stop("Cannot find membership in results$metabolome_modules")
  }
  
  # Create simple assignment table
  result_df <- data.frame(
    Feature = metab_cols,
    Module = membership,
    stringsAsFactors = FALSE
  )
  
  # Sort by module
  result_df <- result_df %>%
    arrange(Module, Feature)
  
  # Summary statistics
  cat("\n--- Summary ---\n")
  module_summary <- result_df %>%
    group_by(Module) %>%
    summarise(N_Features = n()) %>%
    arrange(desc(N_Features))
  
  print(module_summary)
  
  cat(sprintf("\nTotal features: %d\n", nrow(result_df)))
  cat(sprintf("Total modules: %d\n", length(unique(membership))))
  cat("\n✓ Metabolome module assignments extracted\n")
  
  return(result_df)
}

################################################################################
# EXPORT FUNCTIONS
################################################################################

#' Save simple module assignments to Excel
#' 
#' @param micro_assignments Microbiome assignments
#' @param metab_assignments Metabolome assignments
#' @param output_dir Output directory
#' @param prefix File prefix
#' 
#' @return NULL
save_simple_module_assignments <- function(micro_assignments,
                                          metab_assignments,
                                          output_dir = "output/",
                                          prefix = "mlnsem") {
  
  cat("\n════════════════════════════════════════════════════════════\n")
  cat("Saving Module Assignments to Excel\n")
  cat("════════════════════════════════════════════════════════════\n\n")
  
  # Create output directory if needed
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # --- Microbiome ---
  micro_file <- file.path(output_dir, paste0(prefix, "_microbiome_module_assignments.xlsx"))
  
  cat("Creating microbiome workbook...\n")
  wb_micro <- createWorkbook()
  
  # Main sheet with assignments
  addWorksheet(wb_micro, "Module_Assignments")
  writeData(wb_micro, "Module_Assignments", micro_assignments)
  
  # Format header
  headerStyle <- createStyle(
    fontSize = 12,
    fontColour = "#FFFFFF",
    halign = "center",
    fgFill = "#4F81BD",
    border = "TopBottomLeftRight",
    borderColour = "#4F81BD",
    textDecoration = "bold"
  )
  addStyle(wb_micro, "Module_Assignments", headerStyle, rows = 1, cols = 1:2, gridExpand = TRUE)
  
  # Set column widths
  setColWidths(wb_micro, "Module_Assignments", cols = 1:2, widths = c(40, 15))
  
  # Add summary sheet
  addWorksheet(wb_micro, "Summary")
  summary_micro <- micro_assignments %>%
    group_by(Module) %>%
    summarise(N_Features = n()) %>%
    arrange(desc(N_Features))
  
  writeData(wb_micro, "Summary", summary_micro)
  addStyle(wb_micro, "Summary", headerStyle, rows = 1, cols = 1:2, gridExpand = TRUE)
  setColWidths(wb_micro, "Summary", cols = 1:2, widths = "auto")
  
  # Save
  saveWorkbook(wb_micro, micro_file, overwrite = TRUE)
  cat(sprintf("✓ Microbiome assignments saved: %s\n", micro_file))
  
  # --- Metabolome ---
  metab_file <- file.path(output_dir, paste0(prefix, "_metabolome_module_assignments.xlsx"))
  
  cat("Creating metabolome workbook...\n")
  wb_metab <- createWorkbook()
  
  # Main sheet
  addWorksheet(wb_metab, "Module_Assignments")
  writeData(wb_metab, "Module_Assignments", metab_assignments)
  addStyle(wb_metab, "Module_Assignments", headerStyle, rows = 1, cols = 1:2, gridExpand = TRUE)
  setColWidths(wb_metab, "Module_Assignments", cols = 1:2, widths = c(40, 15))
  
  # Summary sheet
  addWorksheet(wb_metab, "Summary")
  summary_metab <- metab_assignments %>%
    group_by(Module) %>%
    summarise(N_Features = n()) %>%
    arrange(desc(N_Features))
  
  writeData(wb_metab, "Summary", summary_metab)
  addStyle(wb_metab, "Summary", headerStyle, rows = 1, cols = 1:2, gridExpand = TRUE)
  setColWidths(wb_metab, "Summary", cols = 1:2, widths = "auto")
  
  # Save
  saveWorkbook(wb_metab, metab_file, overwrite = TRUE)
  cat(sprintf("✓ Metabolome assignments saved: %s\n", metab_file))
  
  cat("\n════════════════════════════════════════════════════════════\n")
  cat("✓ All module assignments saved to Excel!\n")
  cat("════════════════════════════════════════════════════════════\n\n")
}

################################################################################
# MASTER FUNCTION
################################################################################

#' Extract and save simple module assignments
#' 
#' @param results MLN-SEM results object
#' @param microbiome Original microbiome data
#' @param metabolome Original metabolome data
#' @param output_dir Output directory
#' @param prefix File prefix
#' 
#' @return List with microbiome and metabolome assignments
extract_and_save_simple_modules <- function(results,
                                           microbiome,
                                           metabolome,
                                           output_dir = "output/",
                                           prefix = "mlnsem") {
  
  cat("\n")
  cat("╔════════════════════════════════════════════════════════════════╗\n")
  cat("║  SIMPLE MODULE ASSIGNMENT EXTRACTION                           ║\n")
  cat("╚════════════════════════════════════════════════════════════════╝\n")
  
  # Extract assignments with error handling
  micro_assignments <- NULL
  metab_assignments <- NULL
  
  tryCatch({
    micro_assignments <- extract_microbiome_modules_simple(results, microbiome)
  }, error = function(e) {
    cat("\n✗ Error extracting microbiome assignments:\n")
    cat(sprintf("  %s\n\n", e$message))
    stop("Microbiome module assignment extraction failed.")
  })
  
  tryCatch({
    metab_assignments <- extract_metabolome_modules_simple(results, metabolome)
  }, error = function(e) {
    cat("\n✗ Error extracting metabolome assignments:\n")
    cat(sprintf("  %s\n\n", e$message))
    stop("Metabolome module assignment extraction failed.")
  })
  
  # Save to Excel
  save_simple_module_assignments(
    micro_assignments = micro_assignments,
    metab_assignments = metab_assignments,
    output_dir = output_dir,
    prefix = prefix
  )
  
  # Return results
  return(list(
    microbiome = micro_assignments,
    metabolome = metab_assignments
  ))
}

################################################################################
# FOOTER
################################################################################

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║  MLN-SEM Simple Module Assignment Extractor Loaded!            ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

cat("Available functions:\n")
cat("  • extract_microbiome_modules_simple()   - Extract microbiome modules\n")
cat("  • extract_metabolome_modules_simple()   - Extract metabolome modules\n")
cat("  • save_simple_module_assignments()      - Save to Excel\n")
cat("  • extract_and_save_simple_modules()     - Complete workflow\n\n")

cat("Output includes:\n")
cat("  ✓ Feature name\n")
cat("  ✓ Module assignment\n")
cat("  ✓ Summary statistics (feature count per module)\n\n")

cat("Example usage:\n")
cat("  module_assignments <- extract_and_save_simple_modules(\n")
cat("    results = results,\n")
cat("    microbiome = micro_liver,\n")
cat("    metabolome = metab_liver,\n")
cat("    output_dir = 'output/',\n")
cat("    prefix = 'liver_study'\n")
cat("  )\n\n")

cat("✓ Simple version: No kME calculation, just feature-to-module mapping!\n\n")
