#!/usr/bin/env -S Rscript --vanilla

args = commandArgs(trailingOnly=TRUE)
params = list(csvPath = args[1],
              survTime = args[2],
              survCensor = args[3],
              outputFormat = args[4])
print(params)
setwd("/data")

## Read in dataset info -- check approx num of vars to set figure height in output (this is output from callAlign.sh)
df =  as.data.frame(suppressMessages(readr::read_csv(file = params$csvPath)))

## ==== Drop TIL analytics that are not used ====
# -- n_Canc_patch
# -- n_TIL_patch
# -- n_TIL_patch_overlap
# -- patch_ratio
# -- percent_pos
droppable = c(grep("*_Canc_patch", colnames(df)), ## *_ regex used to be flexible with subtype specific patch metrics
              grep("*_TIL_patch*", colnames(df)),
              grep("patchRatio", colnames(df)),
              grep("*_percent_pos", colnames(df)))

df = df[,-droppable]
params$nVars = sum(which(!colnames(df) %in% c("slideID","scaled_PP","TIL_Class",params$survCensor, params$survTime)))

## Call render
rmarkdown::render("/tmp/rmarkdowndir/Descriptive_Statistics.rmd",
                  params = list(csvPath = params$csvPath,
                                nVars = params$nVars,
                                survTime = params$survTime,
                                survCensor = params$survCensor),
                  # Specify the intermediates dir because otherwise it will use /code.
                  # Using /code will fail in singularity because images are immutable.
                  # It will also fail in non-sudo docker because /code is owned by root.
                  intermediates_dir = tempdir(),
                  output_format = params$outputFormat,
                  output_file = "/data/Descriptive_Statistics", knit_root_dir = "/data")
