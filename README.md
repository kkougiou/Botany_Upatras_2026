# Biodiversity_Biogeography_Upatras
A 2-day standalong workshop developed for the postgraduate students enrolled in the Applied Ecology and Environmental Management postgraduate studies programme offered by the Department of Biology, Faculty of Sciences, University of Patras

# Spatial Data Retrieval, Biodiversity and Biogeographical Patterns in R

[![DOI](https://zenodo.org/badge/1248075132.svg)](https://doi.org/10.5281/zenodo.20363115)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Quarto](https://img.shields.io/badge/Made%20with-Quarto-39729E.svg)](https://quarto.org/)
[![R version](https://img.shields.io/badge/R-%E2%89%A54.5.3-276DC3.svg)](https://www.r-project.org/)
[![renv](https://img.shields.io/badge/dependencies-renv-2C7FB8.svg)](https://rstudio.github.io/renv/)
[![targets](https://img.shields.io/badge/pipeline-targets-FF7043.svg)](https://docs.ropensci.org/targets/)
[![Publish](https://github.com/kkougiou/Biodiversity_Biogeography_Upatras/actions/workflows/publish.yml/badge.svg)](https://github.com/kkougiou/Biodiversity_Biogeography_Upatras/actions/workflows/publish.yml)
[![Code style: tidyverse](https://img.shields.io/badge/code%20style-tidyverse-1F77B4.svg)](https://style.tidyverse.org/)

> A two-day reproducible methodological seminar delivered at the **Department of Biology, University of Patras (2026)**, covering spatial data retrieval, cleaning, biodiversity analysis and species distribution modelling in R, with a Mediterranean focus.

🌐 **Live slides:** <https://kkougiou.github.io/Biodiversity_Biogeography_Upatras/>  
📄 **DOI:** <https://doi.org/10.5281/zenodo.20363115>  
✉️ **Contact:** kkougiou@aua.gr  

---

## 📚 What's inside

| Day | Topic | File |
|:---:|---|---|
| **1** | Spatial data retrieval, cleaning, harmonisation | `day1-spatial-data.qmd` |
| **2** | Biodiversity metrics, β-diversity, SDMs | `day2-biodiversity.qmd` |



---

## 🚀 Quick start

### Option A — Clone

```bash
# 1. Clone
git clone [https://github.com/kkougiou/Biodiversity_Biogeography_Upatras.git](https://github.com/kkougiou/Biodiversity_Biogeography_Upatras.git)
cd Biodiversity_Biogeography_Upatras
```

### Option B — Browse only
Just visit the [live slides](https://kkougiou.github.io/Biodiversity_Biogeography_Upatras/).

### Option C — Run in the cloud
[Open in GitHub Codespaces](https://codespaces.new/kkougiou/Biodiversity_Biogeography_Upatras)

---

## 🧱 System requirements

* **R** ≥ 4.5.3
* A [GBIF account](https://www.gbif.org/) (free) for the credentialed downloads

---

## 📂 Repository layout

```text
seminar-patras/
├── .github/workflows/      # CI: render + deploy to GitHub Pages
├── R/                      # Helper functions sourced by targets
├── data/{raw,processed}/   # gitignored — populated by the pipeline
├── images/                 # Logos, static assets
├── .zenodo.json            # Zenodo metadata
└── CITATION.cff            # Citation metadata
```

---

## 🎓 Learning outcomes

By the end of the seminar, participants will be able to:

* Programmatically retrieve occurrence and environmental data from various databases.
* Clean and harmonise geospatial data (CRS, taxonomy, coordinate quality).
* Compute and visualise α- and β-diversity across spatial grids.
* Fit and spatially cross-validate species distribution models.


---

## 📖 How to cite

Kougioumoutzis, K. (2026). *Spatial Data Retrieval and Biodiversity Patterns: A Two-Day Seminar in R* (Version 1.0.0) [Educational material]. Zenodo. https://doi.org/10.5281/zenodo.XXXXXXX

**BibTeX:**
```bibtex
@misc{kougioumoutzis2026seminar,
  author       = {Kougioumoutzis, Konstantinos},
  title        = {Spatial Data Retrieval and Biodiversity Patterns: A Two-Day Seminar in R},
  year         = 2026,
  publisher    = {Zenodo},
  version      = {1.0.0},
  doi          = {10.5281/zenodo.20363115},
  url          = {[https://doi.org/10.5281/zenodo.20363115](https://doi.org/10.5281/zenodo.20363115)}
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

Materials build upon the work of the rOpenSci community, the GBIF Secretariat, the Copernicus programme, and the developers of `sf`, `terra`, `rgbif`, `targets`, `renv` and Quarto. Logo and theme adapted from the official University of Patras visual identity.
