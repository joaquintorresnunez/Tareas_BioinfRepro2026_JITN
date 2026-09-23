Tarea 1.3 BioinfRepro2026 JITN

Bioinformatics and Reproducible Genomic Analysis

Author: Joaquín Ignacio Torres Núñez

Course: Bioinformatics and Reproducible Genomic Analysis

Professor: Ricardo Verdugo

Unit 1 – Session 3

General Description

This repository contains the assignments and solutions for practical R workflows in the Bioinformatics and Reproducible Genomic Analysis course. The main goal of this unit is to master R fundamentals, implement control structures (for loops, if/next), construct custom functions, execute landscape genetics analyses using simple and Partial Mantel tests for Isolation by Resistance (IBR) in Berberis alpina, and manipulate genomic metadata using dplyr and tidyverse.

Repository Contents
Exercise 1–3: Basic R operations calculations, vector operations, and logical indexing on numerical ranges.
Exercise 4: Data loading importing maize and teosinte metadata (maizteocintle_SNP50k_meta_extended.txt) into R data frames.
Exercise 5: Loops and control structures writing for loops with conditional filtering (next) and dynamic data frame storage (rbind).
Exercise 6: Script interpretation analysis of the 1.IBR_testing.r script, identifying required R packages (ade4, ggplot2, sp, vegan) and input files.
Exercise 7: Custom function creation defining calc.tetha(Ne, u) to calculate population diversity parameter $\theta = 4N_e\mu$.
Exercise 8: Partial Mantel Test (IBR in B. alpina) extending 1.IBR_testing.r using vegan::mantel.partial to test $F_{ST}$ against present and LGM climate matrices while controlling for topographic flat distance (B.flat), marked with # Elefante blanco.
Exercise 9: Maize metadata exploration (ExplorandoMaiz.R) R script for data exploration, summary statistics, filtering, and CSV export using dplyr and readr.ToolsTutorials:Bioinformatics and Reproducible Genomic Analysis course repository: https://github.com/u-genoma/BioinfinvReproR Project: https://www.r-project.org/Vegan R Package: https://cran.r-project.org/web/packages/vegan/index.htmlDplyr / Tidyverse: https://dplyr.tidyverse.org/Usage instructionsTo clone this repository:

To clone this repository:
```
git clone https://github.com/joaquintorresnunez/Tareas_BioinfRepro2026_JITN.git
cd Tareas_BioinfRepro2026_JITN
```
