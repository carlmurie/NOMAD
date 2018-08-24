# NOMAD  

## Introduction
The NOMAD (NOrmalization for MAss spectrometry Data) R package implements an ANOVA normalization method designed for iTRAQ mass spectrometry data in a computationally efficient manner. 
 
The NOMAD package provides two main functions. The first (nomadNormalization) applies an ANOVA model to remove the bias of multiple factors and produces normalized peptide abundances. The second (nomadProteinAsssembly) combines the normalized peptide abundances into summary protein abundances and the user is given multiple options as to how the proteins are assembled. 

The structure of the factorial design ANOVA model used in NOMAD allows for a simple algebraic solution identical to the more computationally demanding matrix solution of a linear model. NOMAD scales exceptionally well and provides ridiculously faster computation times than the standard linear models provided in R such as lm.

## Installation

```{r, eval = FALSE}
# install.packages("devtools")
devtools::install_github("carlmurie/NOMAD")
```
## Citation

If you find NOMAD useful in your research please consider citing the following paper.  
  
[Normalization of mass spectrometry data (NOMAD)](https://www.ncbi.nlm.nih.gov/pubmed/29174395) Murie C, Sandri B, Sandberg AS, Griffin TJ, Lehtiö J, Wendt C, Larsson O. Adv Biol Regul. 2018 Jan;67:128-133
