################################################################################
# ShinyApps.io Deployment Script for MLN-SEM Application
################################################################################

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║  MLN-SEM SHINY APP - SHINYAPPS.IO DEPLOYMENT                   ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# Check if rsconnect is installed
if (!requireNamespace("rsconnect", quietly = TRUE)) {
  cat("Installing rsconnect package...\n")
  install.packages("rsconnect")
}

library(rsconnect)

# Configuration
APP_NAME <- "mln-sem-analysis"
APP_TITLE <- "MLN-SEM Analysis Platform"

cat("Configuration:\n")
cat(sprintf("  App Name:  %s\n", APP_NAME))
cat(sprintf("  App Title: %s\n\n", APP_TITLE))

# Check account configuration
accounts <- rsconnect::accounts()

if (nrow(accounts) == 0) {
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat("NO ACCOUNT CONFIGURED\n")
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")
  
  cat("To configure your ShinyApps.io account:\n\n")
  
  cat("1. Go to https://www.shinyapps.io/ and create an account\n")
  cat("2. Click Account → Tokens\n")
  cat("3. Click 'Show' then 'Show Secret'\n")
  cat("4. Run this command in R with your credentials:\n\n")
  
  cat("   rsconnect::setAccountInfo(\n")
  cat("     name = 'your-account-name',\n")
  cat("     token = 'your-token',\n")
  cat("     secret = 'your-secret'\n")
  cat("   )\n\n")
  
  cat("5. Then run this script again.\n\n")
  
  stop("Account configuration required. See instructions above.")
}

# Display configured accounts
cat("Configured accounts:\n")
for (i in 1:nrow(accounts)) {
  cat(sprintf("  %d. %s (%s)\n", i, accounts$name[i], accounts$server[i]))
}
cat("\n")

# Select account
if (nrow(accounts) == 1) {
  ACCOUNT_NAME <- accounts$name[1]
  cat(sprintf("Using account: %s\n\n", ACCOUNT_NAME))
} else {
  cat("Which account do you want to use? ")
  if (interactive()) {
    account_idx <- as.integer(readline())
    ACCOUNT_NAME <- accounts$name[account_idx]
  } else {
    ACCOUNT_NAME <- accounts$name[1]
  }
  cat(sprintf("Using account: %s\n\n", ACCOUNT_NAME))
}

# Check for existing deployments
cat("Checking for existing deployments...\n")
existing <- rsconnect::deployments(appName = APP_NAME)

if (nrow(existing) > 0) {
  cat(sprintf("  Found existing deployment: %s\n", APP_NAME))
  cat("  This will UPDATE the existing app.\n\n")
} else {
  cat("  No existing deployment found.\n")
  cat("  This will CREATE a new app.\n\n")
}

# Confirm deployment
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("READY TO DEPLOY\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

cat("This will deploy the application to ShinyApps.io.\n")
cat("The deployment process may take several minutes.\n\n")

if (interactive()) {
  cat("Continue? (y/n): ")
  response <- readline()
  if (tolower(response) != "y") {
    cat("\nDeployment cancelled.\n")
    quit(save = "no")
  }
}

cat("\n")
cat("Starting deployment...\n\n")

# Deploy
tryCatch({
  
  rsconnect::deployApp(
    appDir = "D:/programming/Github/MLN-SEM",
    appName = APP_NAME,
    appTitle = APP_TITLE,
    account = ACCOUNT_NAME,
    server = "shinyapps.io",
    forceUpdate = TRUE,
    launch.browser = FALSE,
    logLevel = "verbose",
    lint = FALSE  # Skip linting for faster deployment
  )
  
  cat("\n")
  cat("╔════════════════════════════════════════════════════════════════╗\n")
  cat("║  DEPLOYMENT SUCCESSFUL!                                        ║\n")
  cat("╚════════════════════════════════════════════════════════════════╝\n\n")
  
  app_url <- sprintf("https://%s.shinyapps.io/%s/", ACCOUNT_NAME, APP_NAME)
  cat(sprintf("Your app is now live at:\n  %s\n\n", app_url))
  
  cat("Next steps:\n")
  cat("  1. Visit the URL above to access your app\n")
  cat("  2. Go to https://www.shinyapps.io/admin/#/applications\n")
  cat("  3. Configure instance size (recommended: Large or XLarge)\n")
  cat("  4. Set timeout limits (recommended: 3600 seconds)\n")
  cat("  5. Enable logging for debugging\n\n")
  
  if (interactive()) {
    cat("Open app in browser now? (y/n): ")
    response <- readline()
    if (tolower(response) == "y") {
      browseURL(app_url)
    }
  }
  
}, error = function(e) {
  cat("\n")
  cat("╔════════════════════════════════════════════════════════════════╗\n")
  cat("║  DEPLOYMENT FAILED                                             ║\n")
  cat("╚════════════════════════════════════════════════════════════════╝\n\n")
  
  cat(sprintf("Error: %s\n\n", e$message))
  
  cat("Common issues:\n")
  cat("  • Missing packages: Run source('install_packages.R')\n")
  cat("  • Account token expired: Reconfigure with setAccountInfo()\n")
  cat("  • File size too large: Check .Rdata files\n")
  cat("  • Timeout: Try deploying again\n\n")
  
  cat("For detailed logs, check:\n")
  cat("  ~/.rsconnect/\n\n")
})

cat("═══════════════════════════════════════════════════════════════\n\n")
