# Joint Symmetric NonNegative Matrix Factorization with Laplacian Regularization
## Abstract
This repository implements the integrative mathematical framework presented in the paper "Clustering and Integrating of Heterogeneous Microbiome Data by Joint Symmetric Nonnegative Matrix Factorization with Laplacian Regularization" (Ma Y. et al, 2020).

This framework uses joint/multi-view NonNegative Matrix factorization on symmetrical similarity matrices of the data points, with an additional laplacian regularization. 

Evaluation metrics (accuracy, NMI, adjusted rand index) provided, when 'ground truth labels' are used.

Example: lung cancer dataset from database "TCGA" (https://portal.gdc.cancer.gov/), the two 'ground truth' groups are Lung Adenocarcinoma and Lung Squamous Carcinoma. Omic 'views' used are: gene expresion, protein expression, DNA methylation, miRNA expression. All necessary initial matrices provided.

Recommended hyperparameter values:
beta= 0.005, gamma=0.001, delta=0.8

