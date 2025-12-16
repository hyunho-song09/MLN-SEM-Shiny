# MLN-SEM Shiny Application

[![Shiny App](https://img.shields.io/badge/Shiny-Live%20Demo-blue?style=for-the-badge&logo=r)](https://zufqvh-0-0.shinyapps.io/mln-sem-analysis/)
[![License](https://img.shields.io/badge/License-Academic-green?style=for-the-badge)](LICENSE)
[![R Version](https://img.shields.io/badge/R-%E2%89%A5%204.0.0-blue?style=for-the-badge&logo=r)](https://www.r-project.org/)

> **Multi-Layer Network Structural Equation Modeling (MLN-SEM)** web application for integrative microbiome-metabolome-clinical data analysis with causal inference.

## 🌐 Live Demo

**Try it now:** [https://zufqvh-0-0.shinyapps.io/mln-sem-analysis/](https://zufqvh-0-0.shinyapps.io/mln-sem-analysis/)

No installation required - runs directly in your browser!

---

## 📋 Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)

---

## 🔬 Overview

MLN-SEM is a comprehensive Shiny web application that enables researchers to:
- Perform **causal inference** on multi-omics data
- Integrate **microbiome, metabolome, and clinical** datasets
- Validate pathways using **instrumental variable (IV) regression**
- Generate **publication-ready visualizations**
- All through an intuitive **web interface** - no coding required!

### Key Methodology

- **Module Construction**: clustering and dimension reduction
- **Instrumental Variables**: Robust IV-SEM framework
- **Pathway Validation**: Sargan tests and weak-F statistics
- **Network Visualization**: Correlation vs. causal relationship comparison

---

## ✨ Features

### 🎯 Complete Analysis Workflow

```
Data Upload → Parameter Configuration → Analysis Execution → Results Visualization → Download
```

### 📊 Analysis Components

- ✅ **Automated Module Construction**
  - clustering and dimension reduction
  - Eigengene calculation
  - Optimal parameter suggestions

- ✅ **IV-SEM Analysis**
  - Multi-instrumental variable regression
  - Weak-F statistics for instrument validation
  - Sargan tests for overidentification

- ✅ **Advanced Visualizations**
  - Validated pathways network
  - Correlation vs. causal comparison
  - Module heatmaps
  - Effect size plots
  - Weak-F distributions

### 🚀 User-Friendly Interface

- 📁 **Drag-and-drop data upload**
- 💡 **Automatic parameter recommendations**
- 📊 **Real-time progress tracking**
- 📥 **One-click result downloads**
- 📦 **Built-in example dataset** (120 samples)

---

## 💻 Installation

### Prerequisites

- **R** ≥ 4.0.0 ([Download R](https://cran.r-project.org/))
- **RStudio** (recommended) ([Download RStudio](https://posit.co/download/rstudio-desktop/))
