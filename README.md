DNA Methylation Age in Normal vs Adjacent Normal vs Tumor Breast Tissue
====================

This repository contains code that trains a model that predicts the age of tissue samples based on patterns in their respective DNA Methylation data. This approach follows that of Steve Horvath in [1]. A grand majority of the code in this repository is purposed to clean/format DNA Methylation data obtained from the Illumina TruSeq Methyl Capture EPIC Library Prep Kit [2] to a format that can be readily used to train/use the model. The following outlines the pipeline for this analysis and notes what should be edited before running each step.

### 1. data/makeMeth.R 
This code reads in individual sample files that contain columns for chromosome, locus position (bp), number of sequenced Cytosine reads, and the total number of reads for the locus and merges each into a single file. To ensure a high-precision methylation beta value, CpGs are dropped if the number of reads is less than 10. In addition, the code reads in an annotation file and creates a covariate file that contains such items as sample age, batch, and various enviornmental factors. The following variables should be inspected for accuracy before running this code:
* setwd -> This should point to the EpigeneticAge/data directory
* missingness.fraction -> Set the maximum fraction of samples missing a particular CpG. If this threshold is exceeded, the CpG is dropped from the data set
* tissue.type -> Set as "K" for normal, "N" for adjacent normal, or "T" for tumor tissue
* DATADIR -> This should point to the directory storing the individual sample files

### 2. data/keep_cpgs_in_KNT.py
This code reads in merged methylation files for K, N, and T tissue types and removed CpGs that are not present in each individual file. The purpose of this is to guarantee that the "Clock CpGs" selected by the elastic net regression model will be present in each data set. This will make life easier when comparing prediction results across different tissue types. The following variables should be inspected for accuracy before running this code:
* meth_file_K -> Name of the merged methylation data set for normal tissue
* meth_file_N -> Name of the merged methylation data set for adjacent normal tissue
* meth_file_T -> Name of the merged methylation data set for tumor tissue
* out_file_K -> Name of the output file for the reduced methylation data set for normal tissue
* out_file_N -> Name of the output file for the reduced methylation data set for adjacent normal tissue
* out_file_T -> Name of the output file for the reduced methylation data set for tumor tissue

### 3. data/imputeMissingData.R
This code utilizes the K-Nearest Neighbor [3,4] technique to impute missing methylation data. One would only need to run this code if the missingness.fraction variable in makeMeth.R was greater than 0. If you are unsure if you should run this, run the following command on your methylation data set and if it returns a number greater than 0, then impute.
```
grep -o "NA" <methylation file> | wc -l 
```
The following variables should be inspected for accuracy before running this code:
* meth.file.in  -> The name of the input methylation data
* meth.file.out -> The name of the output file containing the methylation data plus imputed missing values

### 4. data/split_train_vali.py
This code reads in the methylation data set and splits it into two separate data sets for training and testing the model, respectively. The following variables should be inspected for accuracy before running this code:
* random.seed()
* train_fraction -> The fraction of samples that will go into the training data set
* cov_filename  -> The name of the covariate file for normal tissue samples
* meth_filename -> The name of the input methylation dataset for normal tissue samples

### 5. trainModel.R
This code reads in the training data set and utilizes the glmnet elastic net linear regression [5-7] algorithm to train the model and select clock CpGs. The code saves .RData files that contain the lambda selected from 10-fold cross validation and the final model. The following variables should be inspected for accuracy before running this code:
* setwd -> This should point to the EpigeneticAge directory
* cov.train -> Name of the input covariate file
* meth.train -> Name of the input methylation training data set

### 6. validateModel.R
This code reads in the output trainModel.R as well as the normal tissue test data and outputs some quick sanity check plots to allow one to quickly assess the quality of the trained model. In addition, the trained model coefficients are stored in a .csv file for convenient access/later use. The following variables should be inspected for accuracy before running this code:
* cov.vali  -> The name of the covariate file for the normal tissue test data set
* meth.vali -> The name of the methylation file for the normal tissue test data set 
* model.dir -> The directory in which the output .RData from the previous step are stored

### 7. data/get_clock_cpgs.sh
This code reads the output .csv from the previous step and generates a text file with the names of all the Clock CpGs selected by the elastic net algorithm. The following variables should be inspected for accuracy before running this code:
* model_output -> Point to the output .csv from the previous step
* out_file -> Desired name of the output text file

### 8. data/keep_clock_cpgs.sh
This code reads in a text file containing the names of the Clock CpGs as well as a methylation data set and outputs a reduced methylation data set that only contains the Clock CpGs. This will greatly increase the speed and reduce resource consumption of subsequent analyses. The following variables should be inspected for accuracy before running this code:
* cpg_file -> Point to the text file containing the names of the clock cpgs
* meth_file -> Point to the target methylation data set to be reduced
* out_file -> Desired name of the output reduced methylation data set

### 9. DNAmage_KvsNvsT.R
This code reads in reduced data sets for normal, adjacent normal, and tumor tissue samples (and respective covariate files) and assesses the relationship of DNA methylation age and sample age. The following variables should be inspected for accuracy before running this code:
* setwd()
* model.dir -> Point to the directory where the model .RData is stored (end string with a "/")
* meth.file.K -> Point to the reduced methylation file for normal tissue
* meth.file.N -> Point to the reduced methylation file for adjacent normal tissue
* meth.file.T -> Point to the reduced methylation file for tumor tissue
* cov.file.K -> Point to the covariate file for normal tissue
* cov.file.N -> Point to the covariate file for adjacent normal tissue
* cov.file.T -> Point to the covariate file for tumor tissue

### 10. DNAmage_CancerSubtype.R
This code reads in reduced data sets for tumor tissue samples (and respective covariate file) and assesses the relationship of DNA methylation age and sample age for different cancer subtypes. The following variables should be inspected for accuracy before running this code:
* setwd()
* model.dir -> Point to the directory where the model .RData is stored (end string with a "/")
* meth.file.T -> Point to the reduced methylation file for tumor tissue
* cov.file.T -> Point to the covariate file for tumor tissue

References
---------------------
1. Horvath S. DNA methylation age of human tissues and cell types. Genome Biology. 2013;14(10):R115. doi:10.1186/gb-2013-14-10-r115.
2. Illumina, TruSeq methyl capture EPIC library prep kit. 2016: p. 1-8.
3. Hastie, T., Tibshirani, R., Sherlock, G., Eisen, M., Brown, P. and Botstein, D., Imputing Missing Data for Gene Expression Arrays, Stanford University Statistics Department Technical report (1999), http://www-stat.stanford.edu/~hastie/Papers/missing.pdf Olga Troyanskaya, Michael Cantor, Gavin Sherlock, Pat Brown, Trevor Hastie, Robert Tibshirani, David Botstein and Russ B. Altman, Missing value estimation methods for DNA microarrays BIOINFORMATICS Vol. 17 no. 6, 2001 Pages 520-525 
4. https://www.rdocumentation.org/packages/impute/versions/1.46.0/topics/impute.knn
5. Friedman, J., Hastie, T. and Tibshirani, R. (2008) Regularization Paths for Generalized Linear Models via Coordinate Descent, https://web.stanford.edu/~hastie/Papers/glmnet.pdf Journal of Statistical Software, Vol. 33(1), 1-22 Feb 2010 http://www.jstatsoft.org/v33/i01/ Simon, N., Friedman, J., Hastie, T., Tibshirani, R. (2011) Regularization Paths for Cox's Proportional Hazards Model via Coordinate Descent, Journal of Statistical Software, Vol. 39(5) 1-13 http://www.jstatsoft.org/v39/i05/ Tibshirani, Robert., Bien, J.
6. Friedman, J.,Hastie, T.,Simon, N.,Taylor, J. and Tibshirani, Ryan. (2012) Strong Rules for Discarding Predictors in Lasso-type Problems, JRSSB vol 74, http://statweb.stanford.edu/~tibs/ftp/strong.pdf Stanford Statistics Technical Report Glmnet Vignette https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html
7. https://www.rdocumentation.org/packages/glmnet/versions/2.0-16/topics/glmnet
