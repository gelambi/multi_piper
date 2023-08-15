[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8250520.svg)](https://doi.org/10.5281/zenodo.8250520)

# Data from: "Multiscale Variability in Nutrients and Secondary Metabolites in a Bat-Dispersed Neotropical Fruit"

Data from: "Multiscale Variability in Nutrients and Secondary Metabolites in a Bat-Dispersed Neotropical Fruit" by Gelambi, M. & Whitehead, S. R. Published in ___Ecology and Evolution___, 2023.
This repository contains scripts and data files related to the analysis of chemical traits, including nutrients and defensive metabolites, in ripe *P. sancti-felicis* fruits expressed as mg/100mg dry weight. The study focuses on exploring multiscale variability at the intraspecific and intraindividual levels and investigates the association between nutrients and secondary metabolites using various statistical approaches.

## Scripts

### 1. `metaMS_alkenylphenols.R`

This script processes the alkenylphenol raw chromatogram using the R package metaMS. It provides the necessary steps to prepare the data for subsequent analysis.

### 2. `Objective1.R` and `Objective2.R`

These two scripts conduct comprehensive analysis of the different chemical traits. They include data exploration, VCA (Variance Components Analysis), several Generalized Linear Mixed Models (GLMMs), NMDS (Non-metric Multidimensional Scaling), and correlation matrix computation. The scripts delve into the association between nutrients and secondary metabolites at both the intraspecific and intraindividual levels.

## Data files

### 1. `data_FINAL.csv`

This data file contains the quantified chemical traits, including nutrients and defensive metabolites, in ripe *P. sancti-felicis* fruits, expressed as mg/100mg dry weight. The columns in the dataset are as follows:

- `sampleID`: A unique label per fruit in the format "Tree Number, Fruit Number."
- `Tree`: Tree ID associated with each fruit.
- `dglucose`: The content of glucose in the fruit, expressed in mg/100mg dry weight.
- `totalphenolics`: The content of total phenolics in the fruit, expressed in mg/100mg dry weight.
- `fructose`: The content of fructose in the fruit, expressed in mg/100mg dry weight.
- `proteins`: The content of total proteins in the fruit, expressed in mg/100mg dry weight.
- `F`: Content of alkenylphenol "F" in the fruit, expressed in mg/100mg dry weight.
- `B`: Content of alkenylphenol "B" in the fruit, expressed in mg/100mg dry weight.
- `D`: Content of alkenylphenol "D" in the fruit, expressed in mg/100mg dry weight.
- `E`: Content of alkenylphenol "E" in the fruit, expressed in mg/100mg dry weight.
- `A`: Content of alkenylphenol "A" in the fruit, expressed in mg/100mg dry weight.
- `G`: Content of alkenylphenol "G" in the fruit, expressed in mg/100mg dry weight.
- `H`: Content of alkenylphenol "H" in the fruit, expressed in mg/100mg dry weight.
- `I`: Content of alkenylphenol "I" in the fruit, expressed in mg/100mg dry weight.
- `J`: Content of alkenylphenol "J" in the fruit, expressed in mg/100mg dry weight.
- `C`: Content of alkenylphenol "C" in the fruit, expressed in mg/100mg dry weight.
- `alkenylphenols_total`: The total content of alkenylphenols (sum of individual alkenylphenols), expressed in mg/100mg dry weight.

### 2. `data_average.csv`

This data file contains the average content of nutrients and defensive metabolites in ripe *P. sancti-felicis* fruits, expressed as mg/100mg dry weight. 

## Figures folder

The 'Figures' folder contains various output files generated from the analyses conducted in the main scripts. These figures visually represent the results and insights obtained from the data exploration and statistical modeling.