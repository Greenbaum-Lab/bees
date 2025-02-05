# Analysis of Social Complexity Dataset

This repository contains R scripts and data used for the analysis of a social complexity dataset as described in Peled et al. (2025). If you have any questions or need assistance, please contact:

**Ohad Peled**  
[ohad.peled@mail.huji.ac.il](mailto:ohad.peled@mail.huji.ac.il)

---

## Repository Structure

The repository is organized into the following directories and files:

- **data/**  
  The original data collected from an extensive literature survey.

- **imputed_data/**  
  Data after imputation of missing values using an iterative PCA method.

- **tree/**  
  A phylogenetic tree containing 80 bee species.

- **pca_scores/**  
  Results of the phylogenetically corrected PCA used in downstream analyses.

- **dimensionality_reduction/**  
  Scripts for:
  - Data preprocessing
  - Imputation of missing data
  - Analysis of PCA and UMAP
  - Visualization of dimensionality reduction results

- **phylosig_cor/**  
  Scripts for:
  - Calculating correlations between traits
  - Hierarchical clustering
  - Estimating phylogenetic signal using Pagel's lambda and Blomberg's K

- **evolutionary_routes/**  
  Scripts for:
  - Phylogenetic reconstruction
  - Analysis of adaptive regime shifts

- **directionality/**  
  Scripts for:
  - Analyzing directional evolution based on PCA values
  - Calculating phenotypic diversification rates
  - Evaluating turn angles in evolutionary trajectories

---

## Getting Started

1. **Data Preparation:**  
   - The original dataset is located in the `data/` folder.
   - For analyses that require imputed data, refer to the `imputed_data/` folder.

2. **Phylogenetic and PCA Analyses:**  
   - The phylogenetic tree is stored in the `tree/` folder.
   - The PCA results used for downstream analyses are in the `pca_scores/` folder.

3. **Analysis Scripts:**  
   - Use the scripts in the `dimensionality_reduction/`, `phylosig_cor/`, `evolutionary_routes/`, and `directionality/` folders to reproduce the analyses presented in the manuscript.

For detailed instructions on running each analysis, please refer to the header comments within the individual R scripts.

---

## Citation

If you use this data or any of the code in your research, please cite the following:

