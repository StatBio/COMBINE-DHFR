# Allostery without Large Motions: Molecular Dissection of a Minimal-Shift MWC Allosteric Regulation

**Authors:** Paul Guénon, Marion Chauveau, Karolina Filipowska, Clément Nizak, Guillaume Stirnemann, Kimberly Reynolds, Olivier Rivoire, Damien Laage

---

## Overview

This repository contains the code used for the prediction of the compensatory mutation within the allosteric network of *E. coli* DHFR. The method is named COMBINE for COmpensatory Mutations via Boltzmann machine INference of Epistasis.

---

## Contents

- **`Data/`** — Contains the MSA and the models used for prediction. Note that the MSA is taken from Thompson et al. 2020 (https://doi.org/10.7554/eLife.53476).
- **`utils_COMBINE.py`** — Includes the main functions used to predict compensatory mutations.  
- **`COMBINE_DHFR.ipynb`** — Jupyter notebook that reproduces the analyses and generates the figures shown in the paper.

---

## Model Inference

The models provided in this repository were inferred using the code available in a separate repository: https://github.com/marionchv/SBM.  
Only the resulting models are included here. 

## Citation

If you use this code or data, please cite the associated publication:

> Guénon, P., Chauveau, M., Filipowska, K., Nizak, C., Stirnemann, G., Reynolds, K., Rivoire, O., & Laage, D.  
> *Allostery without Large Motions: Molecular Dissection of a Minimal-Shift MWC Allosteric Regulation.* 