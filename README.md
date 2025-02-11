# Analysis of Social Complexity Dataset

This repository contains R scripts and data used for the analysis of a social complexity dataset as described in Peled et al. (2025). 

The full article can be found in:
https://doi.org/10.1016/j.cub.2025.01.009

If you have any questions or need assistance, please contact:

**Ohad Peled**  
[ohad.peled@mail.huji.ac.il](mailto:ohad.peled@mail.huji.ac.il)

---

## Repository Structure

The repository is organized into the following directories and files:

- **data**  
  The original data collected from an extensive literature survey.

- **imputed_data**  
  Data after imputation of missing values using an iterative PCA method.

- **tree**  
  A phylogenetic tree containing 80 bee species.

- **pca_scores**  
  Results of the phylogenetically corrected PCA used in downstream analyses.
  
- **reconstruced_data**  
  Data of ancetral state reconstrion including internal nodes. used in directionality and phenotypic diversifcation analyses.

- **dimensionality_reduction**  
  Scripts for:
  - Data preprocessing
  - Imputation of missing data
  - Analysis of PCA and UMAP
  - Visualization of dimensionality reduction results

- **phylosig_cor**  
  Scripts for:
  - Calculating correlations between traits
  - Hierarchical clustering
  - Estimating phylogenetic signal using Pagel's lambda and Blomberg's K

- **evolutionary_routes**  
  Scripts for:
  - Phylogenetic reconstruction
  - Analysis of adaptive regime shifts

- **directionality**  
  Scripts for:
  - Analyzing directional evolution based on PCA values
  - Calculating phenotypic diversification rates
  - Evaluating turn angles in evolutionary trajectories

- **supplementary_figures**  
  Scripts for:
  - Principal Component Analysis (PCA) of four PCs.
  - UMAP dimensionality reduction and plotting for multiple parameters.
  - Ancestral state reconstruction of four PCs.

- **sensitivity_analyses**  
  Scripts for:
  - Confidence ellipses using leave-one-out PCA.
  - Testing for taxonomic sampling bias.
  - Phylogenetic tree sensitivity by sampling alternative trees.

---

## Getting Started
   - Use the scripts to reproduce the analyses presented in the manuscript.

For detailed instructions on running each analysis, please refer to the header comments within the individual R scripts.

---

## Citation

If you use this data or any of the code in your research, please cite the following:

