# CellMem-Power-Seek
## Power-Seek 
An algorithm to discover genes whose expression levels are inheritable across multiple cell divisions, from a single sample scRNA-seq dataset. 
The theory behind this algorithm is described in detail here: https://doi.org/10.1101/2025.01.15.633183  


## Brief summary of Power-Seek implementation on an scRNA-seq dataset 
1. Pre-process the data with PowerSeekDataPrepV1.RMD   
2. Further pre-process with PowerSeekDataPrepV2.ipynb (NOTE: log transforms should *NOT* be performed during normalization!)
3. Find the number of eigenvalues in the power law regime and then generate a ranked list of genes using Power-Seek with PowerSeekAlgorithmOnData.ipynb


 

## Detailed example of the workflow
Here we illustrate how we can process scRNA-seq count matrix data (obtained from 10x Genomics Cell Ranger pipeline or similar methods) to apply Power-Seek and generate a ranked list of genes based on inheritance time-scale. The dataset used in this example is the WM989 melanoma cell line from Harmange et al., Nat Comm 2023.



### Libraries and software needed
R library: Seurat v5

Python library: Numpy (version 2.2), pandas (2.3), scipy (1.15), matplotlib (3.10), Python os module (Python's standard built-in module). Details are in the [requirements.txt](https://github.com/RajuLab/CellMem-Power-Seek/blob/main/requirements.txt).

Software: To run.RMD notebook, RStudio Desktop, and for .ipynb notebook, VS Code have been used.

Installation: A guide to install R and RStudio can be found [here](https://rstudio-education.github.io/hopr/starting.html). VS Code installation guide can be found [here](https://code.visualstudio.com/download). For easy installation of necessary Python libraries, Anaconda distribution can be installed (installation guide [here](https://www.anaconda.com/docs/getting-started/anaconda/install)).



### Step 01 - Initial preprocessing (R-based)

In the notebook '[PowerSeekDataPrepV1.RMD](https://github.com/RajuLab/CellMem-Power-Seek/blob/main/code/PowerSeekDataPrepV1.RMD)', we use data from the WM989 melanoma cell line from Harmange et al., Nat Comm 2023. For the 
convenience of readers, the dataset is provided in this GitHub repo as '[filtered_feature_bc_matrix_WM983.zip](https://github.com/RajuLab/CellMem-Power-Seek/blob/main/filtered_feature_bc_matrix_WM983.zip)'. Please extract. 

Here, using Seurat v5, we perform basic filtering based on the number of mitochondrial and ribosomal genes, the number of unique genes (features) per cell, and the total number of raw UMI counts per cell. 

This step produces a cell-by-gene (features) dataframe (CSV file) for downstream processing.



### Step 02 - Preprocessing (Python-based)

In the notebook '[PowerSeekDataPrepV2.ipynb](https://github.com/RajuLab/CellMem-Power-Seek/blob/main/code/PowerSeekDataPrepV2.ipynb)', we use processed data from step 01 in the previous step (provided in this GitHub repo as '[full_count_mat_batch_1.csv](https://github.com/RajuLab/CellMem-Power-Seek/blob/main/full_count_mat_batch_1.zip)').

Here, we remove genes which have been expressed in cells less than a threshold value (sparsity-based filtering), normalize the expression of individual genes based on the total number of raw counts in each cell (library size normalization). Note that this step should *not* involve a log transform, since this transform can destroy lineage correlations (none of the following steps require assumptions of homoskedasticity or normality of the expression distributions). We also systematically remove genes that show very high expression in a very small number of cells. Then, we select the highly variable genes by comparing the coefficient of variation of each gene. 



### Step 03 - Power-Seek algorithm (Python-based)

In the notebook '[PowerSeekAlgorithmOnData.ipynb](https://github.com/RajuLab/CellMem-Power-Seek/blob/main/code/PowerSeekAlgorithmOnData.ipynb)', we take the output of Step 02 as input. For this particular example, we have used the top 3000 highly variable genes from the previous step (provided in this GitHub repo as ' HighVariable3000CountMatrix_batch_1.csv').

Importantly, in our example, the output of the previous step is in gene-by-cell data matrix form. We convert it to a cell-by-gene data matrix, where all gene expressions are centralized around their mean and normalized with the standard deviation across cells. Then we calculate the eigenspectra of the cell covariance matrix, and plot the log(eigenvalues) versus log(rank of eigenvalue). 

Thereafter, we determine the number of top eigenvalues that follow a power-law by fitting a straight line to the log(eigenvalues) versus log(rank of eigenvalue) plot and checking the value of goodness of fit, R^2. As the threshold increases, the goodness of fit goes up initially within the power-law regime and then decreases. However, beyond a certain point after the power-law regime ends, goodness of fit increases again because of the increased number of eigenvalues in the noisy region. Therefore, we choose the number of eigenvales around the first peak as the **threshold** (n) here. (In our experience, some fluctuation in threshold value does not affect memory gene detection much.) 

Thereafter, we use this threshold (n) to implement the Power-Seek algorithm. This step generates a metric $\delta$ for each gene. For $\delta > 0$, a larger $\delta$ corresponds to a longer memory timescale and therefore, it allows identification of top rank memory genes.
 
## Brief description of other coding files

**PowerSeekDataPrepV1_BreastTissue.Rmd**: Similar as PowerSeekDataPrepV1_WM989.RMD file where we processed breast cancer tissue data from A. Janesick et al., 2023. 

**SimulationDatamat.ipynb**: Generates a cell-by-gene count matrix with lineage correlation.

**PowerSeekSimulation.ipynb**: Applying Power-Seek algorithm to simulated cell-by-gene count matrix.

**functionsV2.py**: Contains function written by us for simulating cell-by-gene count matrix. 






