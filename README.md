## Fitness_Assay_Package

### To use the Lang Lab Fitness Assay Package: 

from the R console 
`library(devtools)`  

`install_github("LangYeastEvoLab/FitnessAssayAnalysis")`

- devtools is installed on the lab computer, but may need to be installed on your personal computer

### Package dependencies
- reshape2
- ggplot2

### Input data format

This script is based on inputs generated by FlowJo and a user supplied well key that identifies the competitions. 
*Note - when reading in data - stringsAsFactors=FALSE!!*<br/>

A Well_key and FC_data template can be found as txt files in this repository. Well_key column names *must* match the provided template. 


### Functions

### Spot_check_data

`output<- Spot_check_data(your well_key, your_flow_data)`  
- goes through every sample (every well) and determines standard error. 
- creates two plots in working directory: 
      - std_error_hisotgram.pdf *histogram of sample standard errors*
      - raw_data_plots.pdf *plots time course ln(exp/ref) for each well*
      
      `Eliminate poor samples? (TRUE/FALSE)` 
      
- If `TRUE` the user will be asked to specifiy a standard error cutoff above which samples (wells) will be removed from the dataset. The output return will contain only those wells that pass the cutoff. The final column in the returned dataframe will be standard error. *Note- if samples are removed from the dataset, the well key must be updated to reflect those by running the `Adjust_well_key(your_well_key, your_reduced_flow_data)` with the output of `Spot_check_data()` as the second argument*

- If `FALSE` function returns the flow_data supplied as an argument unaltered. The final column in the dataframe will be standard error. 

### Analyze_Fitness_Data
`output <- Analyze_Fitness_Data(your well_key, your_flow_data)` 
- cacluates the coefficient of selection of a query strain using time-course compettive growth assay flow cytometry data
- output is a table of selection coefficient, standard error, and 95% CI for each competition or replicate group

      `Group replicates? (TRUE/FALSE):` 

- The well key supplied will be used to look for any duplicate competitions. If `TRUE`, grouping replicates will treat all competitions with identical strain identifiers as biological replicates (columns 2-4 in the well  key are identical). The coefficient of selection for each replicate group will be found by fitting a linear model based on all data points within a replicate. 
- If `FALSE`, each sample will be treated individually. Coefficient of selection will be found by fitting a linear model for each individual sample. 
*Replicates MUST be named identically in the well key.*

`Calculate error? (TRUE/FALSE):` 

- If you have multiple replicates (`TRUE` for `Group replicates?`), the script will find the standard error of regression and 95% CI by fitting a linear model based on all data points within a replicate. 
*Replicates MUST be named identically in the well key.*
- If you do not have replicates of individual competitions (`FALSE` for `Group replicates?`), this method will assign each fitness measurement an error based on the error of regression for each individual sample. 

      `Plot results? (TRUE/FALSE):` 

- For no replicates, each sample will be plotted individually with error bars corresponding to 95% CI. 
- If you have replicates, each replicate will be plotted as one point with error bars corresponding to 95% CI

 ### Fitness_ANCOVA
 `output<-Fitness_ANCOVA(your well_key, your_flow_data)` 
- Uses an ANCOVA (analysis of covariance) to look for statistical differences in slope between competitions.
- User options remain the same as `Analyze_Fitness_Data`. 
- Output is a summary of the ANCOVA on the linear models of all competitions. If there are significant differences amongst competitions, the console will print out a message indicating that. 

### Fitness_ANCOVA_post_hoc
`output<- Fitness_ANCOVA_post_hoc(output_of_Fitness_ANCOVA)`
- output of `Fitness_ANCOVA` can be fed into `Fitness_ANCOVA_post_hoc` *if* the former finds a significant interaction of competition and generation on fitness. 
- will perform a TUKEY HSD post hoc test to find pairwise differences in your data 
- output generated will be a table of all pairwise comparisons and their individual p values 
- significant results will also print to the console. 



