# Botany_Upatras_2026
A 2-day standalong workshop developed for the postgraduate students enrolled in the Applied Ecology and Environmental Management postgraduate studies programme offered by the Department of Biology, Faculty of Sciences, University of Patras

# Spatial Data Retrieval, Biodiversity and Biogeographical Patterns in R

[![DOI](https://zenodo.org/badge/1248075132.svg)](https://doi.org/10.5281/zenodo.20369680)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R version](https://img.shields.io/badge/R-%E2%89%A54.5.3-276DC3.svg)](https://www.r-project.org/)
[![Publish](https://github.com/kkougiou/Botany_Upatras_2026/actions/workflows/publish.yml/badge.svg)](https://github.com/kkougiou/Botany_Upatras_2026/actions/workflows/publish.yml)
[![Code style: tidyverse](https://img.shields.io/badge/code%20style-tidyverse-1F77B4.svg)](https://style.tidyverse.org/)

> A two-day reproducible methodological seminar delivered at the **Department of Biology, University of Patras (2026)**, covering spatial data retrieval, cleaning, biodiversity analysis and species distribution modelling in R, with a Mediterranean focus.

🌐 **Live slides:**

- Day 1: <https://kkougiou.github.io/Botany_Upatras_2026/Biogeography.html>
- Day 2: <https://kkougiou.github.io/Botany_Upatras_2026/Spatial-data.html>

📄 **DOI:** <https://doi.org/10.5281/zenodo.20369680>  
✉️ **Contact:** kkougiou@aua.gr  

---

## 📚 What's inside

| Day | Topic | File |
|:---:|---|---|
| **1** | Biodiversity metrics, β-diversity | `Biogeography.Rmd` | 
| **2** | Spatial data retrieval, cleaning, harmonisation (SDMs) | `Spatial data.Rmd` |



---

## 🚀 Quick start

### Option A — Clone

```bash
# 1. Clone
git clone [https://github.com/kkougiou/Botany_Upatras_2026.git](https://github.com/kkougiou/Botany_Upatras_2026.git)
cd Botany_Upatras_2026
```

### Option B — Browse only
Just visit the [landing page](https://kkougiou.github.io/Botany_Upatras_2026).

### Option C — Run in the cloud
[Open in GitHub Codespaces](https://codespaces.new/kkougiou/Botany_Upatras_2026)

---

## 🧱 System requirements

* **R** ≥ 4.5.3
* A [GBIF account](https://www.gbif.org/) (free) for the credentialed downloads

---

## 📂 Repository layout

```text
seminar-patras/
├── .github/workflows/      # CI: render + deploy to GitHub Pages
├── R scripts/              # Helper functions sourced by targets
├── RDS/                    # Several rds files needed to run time-consuming parts of the analyses
├── Figures/                # Static assets
├── Excel/                  # Excel files
├── Rasters/                # Tif files
├── Shapefiles/             # Shp files
├── .zenodo.json            # Zenodo metadata
└── CITATION.cff            # Citation metadata
```

---

## 🛠️ Pre-Workshop Setup Instructions

Welcome to the Spatial Data and Biodiversity Patterns Workshop! To ensure we can dive straight into the analytics and avoid technical delays, please complete this setup **before** arriving at the workshop.

We are using `renv` for this project. This tool ensures that everyone is running the exact same package versions, guaranteeing that the code works flawlessly on every machine.

### Step 1: System Prerequisites
Before downloading the workshop materials, please ensure your system is up to date:

* **R (version 4.5.3):** Please update to this exact version. [Download R here](https://cran.r-project.org/).
* **RStudio:** Ensure you have a recent version installed.
* **⚠️ Windows Users ONLY - Rtools45:** You **must** install Rtools to compile certain spatial packages from source. 
  1. Download the Rtools45 installer from the CRAN website.
  2. Run the installer and **leave all settings on their defaults** (just keep clicking Next).
  3. Restart your computer or RStudio after installation.
  *(Mac and Linux users can skip this step).*

### Step 2: Download the Project
Clone this repository using Git, or click the green **"<> Code"** button at the top of this page and select **"Download ZIP"**. Extract the folder somewhere easily accessible on your computer.

### Step 3: Restore the Environment
Now, let's let `renv` do the heavy lifting and install the required packages.

1. Open the project folder and double-click the **`.Rproj`** file to open RStudio.
2. `renv` should automatically activate (you will see a message in the console).
3. In the R console, type the following command and press Enter:
   ```R
   renv::restore()
   ```
4. Type `y` (yes) when it asks you to proceed.

**☕ Note:** This process will automatically download and install several heavy spatial analysis packages (like `sf`, `terra`, and `phyloregion`). It may take 10–15 minutes depending on your internet connection and machine. Let it run completely until the console returns to the `>` prompt. 

Once it is done, you are 100% ready for the workshop!

---

## 🎓 Learning outcomes

By the end of the seminar, participants will be able to:

* Programmatically retrieve occurrence and environmental data from various databases.
* Clean and harmonise geospatial data (CRS, taxonomy, coordinate quality).
* Compute and visualise α- and β-diversity across spatial grids.
* Fit and spatially cross-validate species distribution models.


---

## 📖 How to cite

Kougioumoutzis, K. (2026). *Spatial Data Retrieval and Biodiversity Patterns: A Two-Day Seminar in R* (Version 1.0.0) [Educational material]. Zenodo. https://doi.org/10.5281/zenodo.20369680

**BibTeX:**
```bibtex
@misc{kougioumoutzis2026seminar,
  author       = {Kougioumoutzis, Konstantinos},
  title        = {Spatial Data Retrieval and Biodiversity Patterns: A Two-Day Seminar in R},
  year         = 2026,
  publisher    = {Zenodo},
  version      = {1.0.0},
  doi          = {10.5281/zenodo.20369680},
  url          = {[https://doi.org/10.5281/zenodo.20369680](https://doi.org/10.5281/zenodo.20369680)}
}
```

---

## 🤝 Contributing

Pull requests are welcome — especially additions for other Mediterranean or European case studies, fixes to broken data endpoints, or translations of the slides. Please open an issue first if your change is substantial.

By participating, you agree to abide by the [Contributor Covenant Code of Conduct](https://www.contributor-covenant.org/).

---

## 📜 License

* **Code** — MIT
* **Slides, exercises, written content** — CC-BY-4.0

---

## 🙏 Acknowledgements

Materials build upon the work of the rOpenSci community, the GBIF Secretariat and the developers of `sf`, `terra`, `rgbif` and many other R packages.
