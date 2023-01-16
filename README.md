# Repository for containerization of Tumor-TIL alignment and downstream analysis

## This code is affiliated with the following publications
- [Saltz et al. 2018](https://www.cell.com/cell-reports/fulltext/S2211-1247(18)30447-9?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2211124718304479%3Fshowall%3Dtrue)
- [Le et al. 2020](https://ajp.amjpathol.org/article/S0002-9440(20)30188-7/fulltext)
- [Fassler et al. 2022](https://pubmed.ncbi.nlm.nih.gov/35565277/)

## Requirements
- Docker (to build and run pipeline)
- Output from SBU BMI Cancer and Lymphocyte prediction pipelines

## Building and running the pipeline
### Clone the repo

Once you clone/fork/download this repository, set it as your working directory and build the container as below

```
cd REPO_NAME
docker build -t til_analyses .
```

You can then call the alignment and analytics functions. 
# Available functions and their outputs

The functionality for this repository is split into two functions, `callAlign.sh` and `callAnalytics.sh`. They should be run in this order, as `callAnalytics.sh` requires output from `callALign.sh`. Please see details below. You can also call their BaSH level help functions using the -h flag after building the container.
## `callAlign.sh`

This script calls the alignment portion of the pipeline. Assuming appropriate input it will:
- Read and parse cancer and lymphocyte prediction probabilities
- Threshold data into yes/no calls based on algorithm specific thresholds
- Overlay lymphocyte and cancer maps
- Return invasion metrics

### Input
- You will need a base directory with two subfolders, each containing the outputs of our prediction pipelines (either txt or json). **Note: If the subdirs are not named exactly like below, you will need to specify the folder name in the align call (see help for how to update path)**
  * tilPreds 
  * cancPreds
- If the file names within the directories are __exactly__ the same, you don't need more information. If they vary, you will need
  * A csv in in basedir, specified with a -s flag in the align call, with
    - Col 1 is sample name
    - Col 2 is lymph prediction file name for that sample
    - Col 3 is cancer prediction file name for the same sample
  * For a sample basedir, look at "data_for_sample_run/"
### Output (will be written to baseDirectory/outputs)
- Histogram of Percent Invasion
- A csv with 7 columns
  * Sample Name: Col 1 in sampFile or file name
  * n_Canc_patch: total n of cancer patches in WSI
  * n_TIL_patch: total n of lymph patches in WSI
  * n_TIL_patch_overlap: total n of patches containing both Canc and TIL
  * percent_pos: n_TIL_patch_overlap / n_Canc_patch
  * scaled_PP: percent_pos / standard deviation of percent_pos
  * TIL_Class: Binned low vs high around mean invasion value

### Command line options and defaults 

   * -t -- tilPreds: "tilPreds" ## A directory of lymph predictions (see below)
   * -T -- tilThresh: .5 ## The probability threshold that separates a Lymph call from a no-lymph call (Inception, 0.1; VGG, 0.42; ResNet, 0.5, overwritten by -a)
   * -c -- cancPreds: "cancPreds" ## A directory of cancer predictions (see below)
   * -C -- cancThresh: .5 ## The probability threshold that separates a cancer call from a no-lymph call (0.5)
   * -s -- sampFile: "" ## Optional file with structure broken down below (to be used if file names differ between canc and lymph)
   * -o -- outputFile: "Percent_Invasion.csv" (file name for csv)
   * -O -- outputDir: "outputs" ## folder within mounted volume to save all outputs
   * -w -- writePNG ##If passed as a flag, will make a folder PNGs/ within outputDir and write thresholded pngs for all images
   * -h -- help ## Print help

### Sample call

To see help, run

```
docker run til_analyses callAlign.sh -h
```
 To run with default parameters (only works if all prediction files have exact same names), run
```
docker run -v /PATH/TO/BASEDIR:/data til_analyses callAlign.sh
```
If you know the algorithm your lymphocyte predictions were made using, run:
```
docker run -v /PATH/TO/BASEDIR:/data til_analyses callAlign.sh -a [First letter of algorithm (i,v, or r)]
```
If you need to specify sample pairs between Lymph and Canc predictions AND want to write the overlaid maps, run 
```
docker run -v /PATH/TO/BASEDIR:/data til_analyses callAlign.sh -s /path/to/csv (can be relative) -w
```
## `callanalytics.sh`

This script will take the invasion metrics calculated by `callAlign.sh` and calculate:
- Descriptive statistics:
  * Overall Invasion distribution
  * Continuous invasion distribution faceted by variables of interest
  * Invasion class calls faceted by variables of interest
- Survival correlations:
  * Univariable Kaplan-Meier and Cox regression of class calls
  * Univariable Cox regression of scaled invasion (Percent Invasion scaled by standard deviation)
  * Bivariable KM of class calls and variables of interest
  * Bivariable Cox regressions of scaled invasion and variables of interest
  * Fully modified Cox regression using scaled invasion and all variables of interest

### Input
You will need a csv with the following columns
- scaled_PP (output from callAlign)
- TIL_Class (output from callAlign)
- variables of interest (can be any valid column name, code will grab all)
  * Note: this code will use ALL columns (excludes survival columns from descriptive stats), so please trim any excess columns before running
- Column of 0s and 1s indicating outcome censor status names **survCensor**
  * This code expects __right__ censored data, with 0's indicating censor and 1's indicating event
- Column of time to event (numeric), whatever your endpoint of choice may be, named **survTime**
  * You may use different colnames for censor and time, but they must be specific in the call (see samples)
  
  ### Command line options and defaults 

   * -p -- csvPath: "" ## Path to csv to analyze (will assume it is in basedir [/data] of mounted volume)
   * -s -- excludeSurv ## If passed, code will skip survival analyses portion
   * -t -- survTime: "survivalA" ## Column name in csv with time to event
   * -c -- survCensor: .5 ## Column name in csv indicating if patient is censored (0) or has an event (1)
   * -h -- help ## Print help
   
### Output
- A pdf or html of descriptive statistics and survival correlations (if requested)
- A csv of stratification metrics
  * Do we want to write test outcomes for each plot?
  * Do we want to write Hazard ratios and confidence intervals for each group?

### Sample calls

To run using all defaults, run
```
docker run -v /PATH/TO/BASEDIR:/data til_analyses callAnalytics.sh -p CSVFILENAME 
```

To run with different survival column names
```
docker run -v /PATH/TO/BASEDIR:/data til_analyses callAnalytics.sh -p CSVFILENAME -c censorColname -t timeColName
```

To run without survival information
```
docker run -v /PATH/TO/BASEDIR:/data til_analyses callAnalytics.sh  -p CSVFILENAME -s FALSE 
```


Happy analyzing! Please report any issues you may find through the issues tab.
