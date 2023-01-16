#!/usr/bin/env Rscript --vanilla

args = commandArgs(trailingOnly=TRUE)
params = list(csvPath = args[1],
              excludeSurv = args[2],
              survTime = args[3],
              survCensor = args[4],
              outputFormat = args[5])
print(params)

## Read in dataset -- check approx num of vars to set figure height in output
setwd("/data")
df =  as.data.frame(suppressMessages(readr::read_csv(file = params$csvPath)))
params$nVars = sum(which(!colnames(df) %in% c("scaled_PP","TIL_Class",params$survCensor, params$survTime)))

## Call render
if(params$excludeSurv == "true"){
   params$excludeSurv = TRUE
   rmarkdown::render("/code/Descriptive_Statistics.rmd", params = list(csvPath = params$csvPath,
                                                                 nVars = params$nVars,
                                                                 survTime = params$survTime,
                                                                 survCensor = params$survCensor),
                     output_format = params$outputFormat,
                     output_file = "/data/Descriptive_Statistics", knit_root_dir = "/data")
} else {
   params$excludeSurv = FALSE
   rmarkdown::render("/code/Descriptive_Statistics.rmd", params = list(csvPath = params$csvPath,
                                                                 nVars = params$nVars),
                     output_format = params$outputFormat,
                     output_file = "/data/Descriptive_Statistics", knit_root_dir = "/data")
}

