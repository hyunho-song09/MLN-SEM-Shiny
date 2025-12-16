################################################################################
# MLN-SEM Shiny App - Package Installer
# 
# This script installs all required packages for the MLN-SEM Shiny application
################################################################################

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║  MLN-SEM SHINY APP - PACKAGE INSTALLATION                      ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org/"))

cat("Checking R version...\n")
r_version <- R.version$major
r_minor <- R.version$minor

cat(sprintf("  R version: %s.%s\n", r_version, r_minor))

if (as.numeric(r_version) < 4) {
  stop("ERROR: R version 4.0.0 or higher is required!")
}

cat("  ✓ R version check passed\n\n")

# Define required packages
cat("Defining required packages...\n\n")

packages <- list(
  # Shiny & UI
  shiny = list(
    package = "shiny",
    description = "Web Application Framework",
    category = "UI"
  ),
  shinydashboard = list(
    package = "shinydashboard",
    description = "Dashboard Layout for Shiny",
    category = "UI"
  ),
  shinyWidgets = list(
    package = "shinyWidgets",
    description = "Custom Input Widgets",
    category = "UI"
  ),
  DT = list(
    package = "DT",
    description = "Interactive Data Tables",
    category = "UI"
  ),
  plotly = list(
    package = "plotly",
    description = "Interactive Plots",
    category = "UI"
  ),
  
  # Tidyverse ecosystem
  tidyverse = list(
    package = "tidyverse",
    description = "Data Science Packages",
    category = "Data"
  ),
  dplyr = list(
    package = "dplyr",
    description = "Data Manipulation",
    category = "Data"
  ),
  tidyr = list(
    package = "tidyr",
    description = "Data Tidying",
    category = "Data"
  ),
  purrr = list(
    package = "purrr",
    description = "Functional Programming",
    category = "Data"
  ),
  tibble = list(
    package = "tibble",
    description = "Modern Data Frames",
    category = "Data"
  ),
  stringr = list(
    package = "stringr",
    description = "String Operations",
    category = "Data"
  ),
  
  # Statistical methods
  AER = list(
    package = "AER",
    description = "Applied Econometrics with R",
    category = "Statistics"
  ),
  ivreg = list(
    package = "ivreg",
    description = "Instrumental Variables Regression",
    category = "Statistics"
  ),
  glmnet = list(
    package = "glmnet",
    description = "LASSO and Ridge Regression",
    category = "Statistics"
  ),
  sandwich = list(
    package = "sandwich",
    description = "Robust Covariance Matrix",
    category = "Statistics"
  ),
  lmtest = list(
    package = "lmtest",
    description = "Linear Model Testing",
    category = "Statistics"
  ),
  lavaan = list(
    package = "lavaan",
    description = "Structural Equation Modeling",
    category = "Statistics"
  ),
  Matrix = list(
    package = "Matrix",
    description = "Sparse and Dense Matrices",
    category = "Statistics"
  ),
  
  # Visualization
  ggplot2 = list(
    package = "ggplot2",
    description = "Grammar of Graphics",
    category = "Visualization"
  ),
  igraph = list(
    package = "igraph",
    description = "Network Analysis and Visualization",
    category = "Visualization"
  ),
  pheatmap = list(
    package = "pheatmap",
    description = "Pretty Heatmaps",
    category = "Visualization"
  ),
  
  # File I/O
  writexl = list(
    package = "writexl",
    description = "Export to Excel",
    category = "IO"
  ),
  
  # Deployment (optional)
  rsconnect = list(
    package = "rsconnect",
    description = "Deployment to ShinyApps.io",
    category = "Deployment",
    optional = TRUE
  )
)

# Function to install a package
install_package <- function(pkg_info) {
  pkg_name <- pkg_info$package
  
  if (requireNamespace(pkg_name, quietly = TRUE)) {
    cat(sprintf("  ✓ %s already installed\n", pkg_name))
    return(TRUE)
  } else {
    cat(sprintf("  → Installing %s (%s)...\n", pkg_name, pkg_info$description))
    tryCatch({
      install.packages(pkg_name, dependencies = TRUE, quiet = FALSE)
      if (requireNamespace(pkg_name, quietly = TRUE)) {
        cat(sprintf("    ✓ Successfully installed %s\n", pkg_name))
        return(TRUE)
      } else {
        cat(sprintf("    ✗ Failed to verify %s installation\n", pkg_name))
        return(FALSE)
      }
    }, error = function(e) {
      cat(sprintf("    ✗ Error installing %s: %s\n", pkg_name, e$message))
      return(FALSE)
    })
  }
}

# Install packages by category
categories <- unique(sapply(packages, function(x) x$category))

for (cat_name in categories) {
  cat("\n")
  cat(sprintf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"))
  cat(sprintf("Installing %s Packages\n", cat_name))
  cat(sprintf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n"))
  
  cat_packages <- packages[sapply(packages, function(x) x$category == cat_name)]
  
  for (pkg_info in cat_packages) {
    # Skip optional packages if user doesn't want them
    if (!is.null(pkg_info$optional) && pkg_info$optional) {
      cat(sprintf("\n%s is optional. Install? (y/n): ", pkg_info$package))
      if (interactive()) {
        response <- readline()
        if (tolower(response) != "y") {
          cat(sprintf("  ⊘ Skipped %s\n", pkg_info$package))
          next
        }
      }
    }
    
    install_package(pkg_info)
  }
}

# Final verification
cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║  VERIFYING INSTALLATION                                        ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

required_packages <- sapply(packages[sapply(packages, function(x) 
  is.null(x$optional) || !x$optional)], function(x) x$package)

missing_packages <- c()

for (pkg in required_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  ✓ %s: OK\n", pkg))
  } else {
    cat(sprintf("  ✗ %s: MISSING\n", pkg))
    missing_packages <- c(missing_packages, pkg)
  }
}

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║  INSTALLATION SUMMARY                                          ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

if (length(missing_packages) == 0) {
  cat("  ✓✓✓ All required packages installed successfully! ✓✓✓\n\n")
  cat("You can now run the Shiny app:\n")
  cat("  library(shiny)\n")
  cat("  runApp('app.R')\n\n")
  
  return(invisible(TRUE))
} else {
  cat("  ✗✗✗ Some packages failed to install ✗✗✗\n\n")
  cat("Missing packages:\n")
  for (pkg in missing_packages) {
    cat(sprintf("  - %s\n", pkg))
  }
  cat("\nPlease try installing them manually:\n")
  cat(sprintf("  install.packages(c(%s))\n", 
             paste0("'", missing_packages, "'", collapse = ", ")))
  cat("\nOr check for system dependencies (Linux):\n")
  cat("  sudo apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev\n\n")
  
  return(invisible(FALSE))
}

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("For deployment to ShinyApps.io, also run:\n")
cat("  source('deploy_shinyapps.R')\n")
cat("═══════════════════════════════════════════════════════════════\n\n")
