# Pepper *et al.* 2026

[![DOI:10.3389/fimmu.2026.1855464](https://zenodo.org/badge/DOI/10.3389/fimmu.2026.1855464.svg)](https://doi.org/10.3389/fimmu.2026.1855464)
[![Citation Badge](https://api.juleskreuer.eu/citation-badge.php?doi=10.3389/fimmu.2026.1855464)](https://juleskreuer.eu/citation-badge/)
[![PyPI license](https://img.shields.io/pypi/l/ansicolortags.svg)](https://pypi.python.org/pypi/ansicolortags/)

This repository contains the code and processed data supporting the paper:

> Pepper JS, Granland CM, Clark SL, Thornton RB, Bayliss J, Foo MZE, Billingham W, Scott N, Fulurija A, Strickland DH, Wiertsema SP, Richmond PC, Seppanen EJ, Kirkham L-AS and Tjiam MC (2026) Mucosal immune cell priming by intranasally delivered Haemophilus haemolyticus is associated with heterologous protection against influenza and nontypeable Haemophilus influenzae. Front. Immunol. 17:1855464. doi: 10.3389/fimmu.2026.1855464

## Overview

This repository provides the analysis workflow, figure-generation scripts, and final processed outputs associated with the study. The goal is to make the analysis and derived figures inspectable and reproducible from the available project materials.

## Repository structure

- `Analysis Pipeline.Rmd` — main analysis workflow notebook/report file.
- `ClusterID_App.R` — interactive or supporting clustering application script.
- `data/` — project metadata and raw data inputs.
- `Figures/` — figure-specific analysis and plotting scripts organized by figure type.
- `output/` — processed summary tables and intermediate outputs used for downstream plotting and comparisons.

## Data availability

Processed data required for inspection and graph generation has been included in this repository under the `output/` directory.

The raw data supporting the conclusions of this article are not stored directly in this repository. Access to raw FCS and associated source data may be requested from the corresponding author or first author, in accordance with the publication statement.

## Reproducibility notes

To reproduce the analysis or regenerate figures:

1. Open the project in RStudio using the included `Pepper-2026.Rproj` file.
2. Ensure the required R packages are installed.
3. Run the analysis scripts in the `Figures/` subfolders and the primary analysis R Markdown file as appropriate.
4. Use the processed files in `output/` for inspection, comparisons, and figure assembly.

## Citation

Please cite the paper using the citation information above if you reuse the analysis, figures, or processed data in this repository.

## License

This repository is distributed under the MIT license found in the repository root.
