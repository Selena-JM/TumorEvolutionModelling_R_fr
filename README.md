# Description

**This project was done as part of a course on modeling biological systems in fifth year of engineering school.**

The goal of this project was to reproduce and analyze a model in a paper. The paper this project is based on is : [Can the Kuznetsov Model Replicate and Predict Cancer
Growth in Humans?](https://link.springer.com/article/10.1007/s11538-022-01075-7)

Everything is explained in detail in our report (in french) saved in this repository, but here is a summary and a simple code explanation.

## Goal of the study

In this paper, the authors use a slightly simplified Kuznetsov model to try to predict the tumor evolution for patients in this database : [Classical mathematical models for prediction of response to chemotherapy and immunotherapy](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009822#pcbi.1009822.s006). This database contains the evolution of tumor long diameter of patients from five different studies designed to assess the efficiency of Atezolizumab, an immunotherapy, on Non Small Cell Lung Cancer and bladder cancer. 

## Method

The article can be separated in 4 main parts : 

- OP1 : find the optimal parameters of the model to fit the measure points of each patient. This problem is solved for each patient separately.
- OP2 : identifiability analysis -> finding the maximal and minimal parameters so that the resulting curve is close to the optimal curve (close by 20%). This problem is solved 6 times (1 time for each parameter) for each patient.
- predictions : solve OP1 without taking into consideration the last two points of each patient. The curve at the time of the last two points is thus a prediction.
- OP3 : this problem aims to find the worse predictions the model could do. This problem gives two curves that are probable on the n-2 first points (not farther that 10% of the predictive curve) but that are as far away as possible from each other at the last point. Those two curves give a sort of confidence interval of the prediction.

We reproduced all the optimization problems and analyses in this paper and we conducted additionnal analyses, like an uncertainty and a sensitivity analysis as well as an extension of the model with clustering methods to study the possibility of finding global parameters. 

## Running the code 

If you have the data frames : set the parameter "create_data_frames" to FALSE at the very beginning of the "Main.R" file, and you can run the whole code. The data frames are loaded and every function can be used. To visualize results just set the patient number of the patient you want to see in the parameter pat_nb just before each plotting function, in the same code section. There are also more general results that can be printed (like the goodness of fit analysis for example), eveything is already written so you just need to run the wanted section of the code. 

If you don't have the data frames : set the parameter "create_data_frames" to TRUE at the very beginning of the "Main.R" file, and every data frame will be created. This process takes a lot of time because there are a lot of optimization problems to solve. Then every function can be used. 

## Code explanation

The main script is the file "Main.R", the other files contain functions that are called in "Main.R". The script is divided in 5 parts : 

- Data preprocessing and processing : calling the functions to create the data frames used in the rest of the study : Data and Data_converted. Data contains the imported Data with our preprocessing and Data_converted contains the same patients but the long diameter of the tumor is converted to a number of tumor cell. Then the functions to solve each of the problems on each of the patients and create the right data frames are called. The following data frames are created: Data_OP1, Data_OP2, Data_Pred and Data_OP3 for the Data, and Fit and Fit_OP3 for the goodness of fit study.

- Plots : calling the adequate plotting functions. Every plot showed in the report was done thanks to one of those functions.

- Uncertainty analysis : we used a uniforme distribution between -2mm and 2mm to model the measure uncertainties, and we then propagated those uncertainties to analyse their effect on the predictions and their confidence interval.

- Sensitivity analysis : local unidimensional and global multidimensional sensitivity analysis (using Sobol indices for the latter)

- Global optimization : trying clustering methods to see if patients can be clustered thanks to their parameters and their acceptable intervals.

For the last 3 parts, the plots are done in each part separately.
