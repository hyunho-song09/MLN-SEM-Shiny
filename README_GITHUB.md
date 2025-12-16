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
- [Quick Start](#quick-start)
- [Usage](#usage)
- [Documentation](#documentation)
- [Example Data](#example-data)
- [Citation](#citation)
- [Authors](#authors)
- [License](#license)

---

## 🔬 Overview

MLN-SEM is a comprehensive Shiny web application that enables researchers to:
- Perform **causal inference** on multi-omics data
- Integrate **microbiome, metabolome, and clinical** datasets
- Validate pathways using **instrumental variable (IV) regression**
- Generate **publication-ready visualizations**
- All through an intuitive **web interface** - no coding required!

### Key Methodology

- **Module Construction**: WGCNA-based eigengene analysis
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
  - WGCNA-based clustering
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

### Option 1: Clone Repository

```bash
git clone https://github.com/YOUR-USERNAME/MLN-SEM-Shiny.git
cd MLN-SEM-Shiny
```

### Option 2: Download ZIP

1. Click the green **"Code"** button above
2. Select **"Download ZIP"**
3. Extract to your desired location

### Install Required Packages

```r
# Automatic installation
source("install_packages.R")

# Or manually
install.packages(c(
  "shiny", "shinydashboard", "shinyWidgets", "DT", "plotly",
  "tidyverse", "dplyr", "tidyr", "purrr", "tibble",
  "AER", "ivreg", "glmnet", "sandwich", "lmtest", "lavaan",
  "ggplot2", "igraph", "pheatmap", "writexl"
))
```

---

## 🚀 Quick Start

### 3-Step Launch

```r
# 1. Set working directory
setwd("/path/to/MLN-SEM-Shiny")

# 2. Load Shiny
library(shiny)

# 3. Run app
runApp("app.R")
```

The app will automatically open in your default web browser!

### Using Example Data

1. Click **"Load Example Data"** button in the app
2. Wait 3 seconds for automatic loading
3. Proceed to **"Parameters"** tab
4. Click **"Get Suggested Parameters"**
5. Run analysis!

**Total time:** ~5 minutes to complete workflow 🎉

---

## 📖 Usage

### Step-by-Step Workflow

#### 1. **Data Upload**

Upload three CSV files:
- **Clinical data** (`01_cp.csv`): Patient phenotypes and clinical parameters
- **Metabolome data** (`02_met.csv`): Metabolite abundance (prefix columns with `M_`)
- **Microbiome data** (`03_micro.csv`): Microbial abundance (prefix columns with `s_`)

Or use the built-in example dataset (120 samples, liver disease study).

#### 2. **Configure Parameters**

- Select **outcome variables** to predict
- Choose **covariates** to adjust for (Age, Gender, BMI)
- Set **module counts** (k) for dimensionality reduction
- Adjust **quality thresholds** (Weak-F, p-value)

💡 **Tip:** Click "Get Suggested Parameters" for automatic optimization!

#### 3. **Run Analysis**

- Click **"▶ Run Complete Analysis"**
- Monitor progress in real-time
- Typical runtime: 3-10 minutes

The app performs:
- Data preprocessing & CLR transformation
- Module construction
- IV regression analysis
- Pathway validation
- Network visualization

#### 4. **Explore Results**

**Results Tab:**
- Module construction summary
- IV analysis statistics
- Top validated pathways
- Complete results table

**Visualizations Tab:**
- Validated pathways network
- Correlation vs. causal comparison
- Module heatmaps
- Effect size plots
- Weak-F statistics

#### 5. **Download Results**

Available formats:
- 📊 **RDS**: Complete analysis object
- 📄 **CSV**: Pathway tables and IV results
- 📑 **XLSX**: Module assignments
- 📈 **PDF**: Individual visualizations
- 📦 **ZIP**: All plots bundled
- 📰 **HTML**: Automated report

---

## 📚 Documentation

| Document | Description |
|----------|-------------|
| [한글 사용가이드](한글_사용가이드.md) | Complete Korean user manual |
| [QUICKSTART.md](QUICKSTART.md) | 5-minute quick start guide |
| [DEPLOYMENT_GUIDE.md](DEPLOYMENT_GUIDE.md) | Deployment instructions |
| [FILE_STRUCTURE.md](FILE_STRUCTURE.md) | File organization guide |
| [EXAMPLE_DATA_GUIDE.md](EXAMPLE_DATA_GUIDE.md) | Example data documentation |
| [CHANGELOG.md](CHANGELOG.md) | Version history |

---

## 📦 Example Data

Built-in dataset included (no download required):
- **Samples**: 120 patients
- **Study**: Liver disease (MASLD/ALD)
- **Clinical**: 18 variables (AST, ALT, BMI, etc.)
- **Metabolome**: ~200 metabolites
- **Microbiome**: ~150 microbial species

### Data Format Requirements

**Clinical Data:**
```csv
Index,Sample,Group,Gender,Age,BMI,AST,ALT,...
1,patient_1,NC,2,69,19.13,22,14,...
```

**Metabolome Data:**
```csv
Index,M_Glucose,M_Lactate,M_Pyruvate,...
1,125.3,15.2,8.7,...
```

**Microbiome Data:**
```csv
Index,s_Lactobacillus,s_Bacteroides,...
1,0.15,0.23,...
```

---

## 📊 Screenshots

### Main Interface
![Data Upload](screenshots/data_upload.png)
*Intuitive data upload with example data loading*

### Analysis Results
![Results](screenshots/results.png)
*Comprehensive statistical results and validated pathways*

### Network Visualization
![Network](screenshots/network.png)
*Publication-ready network visualizations*

---

## 🔬 Citation

If you use this application in your research, please cite:

```bibtex
@software{mln_sem_shiny,
  author = {HHSONG and Lee, Do Yup},
  title = {MLN-SEM: Multi-Layer Network Structural Equation Modeling Shiny Application},
  year = {2025},
  url = {https://github.com/YOUR-USERNAME/MLN-SEM-Shiny},
  note = {Metabolomics Lab, Seoul National University}
}
```

---

## 👥 Authors

**Developer & Researcher:**
- HHSONG (PhD Candidate)
- Metabolomics Lab, Seoul National University

**Principal Investigator:**
- Prof. Do Yup Lee
- Seoul National University

---

## 📄 License

This software is provided for **academic and research purposes**.

For commercial use, please contact the authors.

---

## 🤝 Contributing

Contributions are welcome! Please feel free to:
- 🐛 Report bugs via [Issues](https://github.com/YOUR-USERNAME/MLN-SEM-Shiny/issues)
- 💡 Suggest new features
- 🔧 Submit pull requests

---

## 📞 Support

- **Documentation**: See [한글_사용가이드.md](한글_사용가이드.md) for detailed Korean guide
- **Issues**: [GitHub Issues](https://github.com/YOUR-USERNAME/MLN-SEM-Shiny/issues)
- **Email**: [Your Email]
- **Lab**: Metabolomics Lab, Seoul National University

---

## 🙏 Acknowledgments

- **WGCNA**: Module construction methodology
- **AER Package**: Instrumental variable regression
- **Shiny Team**: Web application framework
- **Community**: Testing and feedback

---

## 📈 Project Status

![Status](https://img.shields.io/badge/Status-Active-brightgreen)
![Version](https://img.shields.io/badge/Version-1.2.0-blue)
![Maintained](https://img.shields.io/badge/Maintained-Yes-green)

**Latest Version:** v1.2.0 (2025-01-16)
- ✅ Example data embedded
- ✅ Simplified deployment
- ✅ Enhanced error handling
- ✅ Complete documentation

---

## 🗺️ Roadmap

### Upcoming Features
- [ ] Multiple example datasets
- [ ] Batch analysis mode
- [ ] Custom visualization themes
- [ ] R Markdown report generation
- [ ] API endpoint access

---

<div align="center">

**⭐ If you find this useful, please star the repository! ⭐**

Made with ❤️ at Seoul National University

[🌐 Live Demo](https://zufqvh-0-0.shinyapps.io/mln-sem-analysis/) | [📖 Documentation](한글_사용가이드.md) | [🐛 Report Issue](https://github.com/YOUR-USERNAME/MLN-SEM-Shiny/issues)

</div>
