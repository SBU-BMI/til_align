#!/usr/bin/env Rscript --vanilla

library(abind,quietly = T)
library(plyr,quietly = T)
library(dplyr,quietly = T)
library(ggplot2,quietly = T)
library(readr,quietly = T)
library(raster,quietly = T)


setwd("/data")
args = commandArgs(trailingOnly=TRUE)
params = list(algorithm = args[1],
              tilDir = args[2],
              tilThresh = args[3],
              cancDir = args[4],
              cancThresh = args[5],
              sampFile = args[6],
              outputFile = args[7],
              outputDir = args[8],
              writePNG = args[9],
              sampInfo = args[10])

## ==== Before anything, flag for missing information ====
if(params$tilDir == "blank"){
   stop("Lymphocyte prediction location not supplied, please supply path to lymph csvs (relative or absolute)")
} else if (!dir.exists(params$tilDir)){
   stop("Lymphocyte prediction location does not exist (or is inaccesible), please check supplied path and permissions")
}

if(params$cancDir == "blank"){
   stop("Cancer prediction location not supplied, please supply path to Cancer csvs (relative or absolute)")
} else if (!dir.exists(params$cancDir)){
   stop("Cancer prediction location does not exist (or is inaccesible), please check supplied path and permissions")
}

## ==== If needed, make output dir ====
if(!dir.exists(params$outputDir)){
   dir.create(params$outputDir)
}

## ==== Check writePNG request ====
if(params$writePNG %in% c("true", "TRUE","t","T",T, TRUE)){
   params$writePNG = TRUE
} else { params$writePNG = FALSE }

## ==== Update and print thresholds ====
params$tilThresh = as.numeric(params$tilThresh)
params$cancThresh = as.numeric(params$cancThresh)
print("TIL Algorithm (Threshold): Frontiers InceptionV4: 0.1")

print("=========== Params after R parsing, if any misalignment please check your flags ===========")
print(params)

### ==== For testing!! comment out ====
# params = list(tilDir = "./scripts_and_docker/data_for_sample_run/tilPreds",
#               tilThresh = 0.1,
#               cancDir = "./scripts_and_docker/data_for_sample_run/cancPreds",
#               cancThresh = 0.5,
#               #sampFile = "/datadrive/shared/image_analysis/ML_output/extdata/flexible_testDir/sampFileShort.csv",
#               outputFile = "Percent_Invasion.csv",
#               outputDir = "outputs",
#               writePNG = T)
# tils = "TCGA-3C-AALI-01Z-00-DX1.F6E9A5DF-D8FB-45CF-B4BD-C6B76294C291.csv"
# canc = tils
### ==== Above is for testing!! comment out ====


### ==== Read in files ====
tils <- sort(list.files(params$tilDir), decreasing = TRUE) # data dirs can be hardcoded or relative
canc <- sort(list.files(params$cancDir),decreasing = TRUE) # data dirs can be hardcoded or relative

## ==== Drop extraneous files ("color-*" and "*.low_res") ====
tils = tils[grep("^prediction", tils)]
writeLines(" . . . Dropping low_res and color- files . . . ")
if(grepl("low_res", tils)){
   tils = tils[-grep("low_res", tils)]
}

canc = canc[grep("^prediction", canc)]
if(grepl("low_res", canc)){
   canc = canc[-grep("low_res", canc)]
}



## ==== Check for missing file pairs (if there is a tumor or lymph prediction but not the other)
sampNames = tils
writeLines(" . . . Checking for tumor/lymph pairs . . . ")
if(length(setdiff(tils,canc)) > 0){
   warning("Some cancer-lymph predictions are missing pairs, skipping non-paired samples")
   writeLines(paste0(length(setdiff(tils,canc)," out of ", length(tils)," are missing pairs. They will be dropped")))
   writeLines(paste0("Samples missing pairs (up to first 6) are: ",
                head(setdiff(tils,canc)))
   )
   tils = intersect(tils,canc)
   canc = intersect(tils,canc)
}

## ==== After removing unequal pairs, make sure that files are a 1 to 1 ordered match. No duplicates or differences whatsoever ====
if(length(tils)==0){
   stop("No predictions had exact pairs. Please ensure lymph and cancer pairs have the exact same name.")
}
if(!all.equal(tils,canc)){
   stop("Supplied files do not match. Must be same names. Please adjust directories and try again.")
}

## ==== This is redundant with the  intersect/setdiff portion above
# if(length(tils) != length(canc)){
#    warning("Unequal number of til prediction files and cancer prediction files, check your inputs")
# }

## ==== Check for read and align format using first entry in both cancer and lymph file lists
lymphFormatCSV = grepl("csv", tils[1])
cancFormatCSV = grepl("csv", canc[1])

## ==== If they are different formats, warn but continue. This is likely not an optimal implementation
if(lymphFormatCSV != cancFormatCSV){
   warning("Prediction formats for cancer and lymphocyte are different, please double check this is intentional")
}

## ==== read in sampInfo csv ahead of time if provided
if(params$sampInfo == 'blank'){
   message("No Sample Information passed, only alignment will be run.")
} else if(!file.exists(params$sampInfo)){
   stop(paste0("Sample information file (",
               params$sampInfo,
               ") not found, please check path and try again"))
} else {
   ## ==== Read in extra sample Info ====
   sampInfo = as.data.frame(
      readr::read_csv(
         params$sampInfo,
         col_names = T,
         col_type = cols()
      )
   )
   
   ## ==== enforce expected formatting on sampInfo
   if(!("slideID" %in% colnames(sampInfo))){ 
      stop("slideID column not found, please check sampInfo.csv formatting")
   }
   if(!("survivalA" %in% colnames(sampInfo))){ # If no survivalA column found
      stop("Time to event (survivalA) not provided, please check sampInfo.csv formatting")
   }
   if(!("censorA.0yes.1no" %in% colnames(sampInfo))){ # If no censor column found
      stop("Event status (censorA.0yes.1no) not provided, please check sampInfo.csv formatting")
   }
}

# ## ==== Align ====
# Initialize empty data frame for output
percent_calls <- data.frame(slideID = "Not Run",
                            n_Canc_patch = numeric(length = length(canc)),
                            n_TIL_patch = numeric(length = length(canc)),
                            n_TIL_patch_overlap = numeric(length = length(canc)),
                            percent_pos = numeric(length = length(canc)),
                            patchRatio = numeric(length = length(canc)),
                            stringsAsFactors = F)
count = 0
for(j in 1:length(canc)){
   count = count + 1
   # =============================================================
   # Load in and reorder Cancer Annotation File
   # =============================================================
   if(cancFormatCSV){
      C1 = as.data.frame(
         readr::read_csv(
            paste(params$cancDir,canc[j], sep = "/"),
            col_names = T,
            col_type = cols()
         )
      )
      
      # Order data as required (y then x)
      C1 = C1[order(C1$miny, C1$minx),]
      C_range = C1$width[1] ## WSIInfer provides patch size, should we include square check? Seems unnecessary
   } else {
      C1 = as.data.frame(
         readr::read_table(
            paste(params$cancDir,canc[j], sep = "/"),
            col_names = F,
            col_type = cols()
         )
      )
      
      # Data is not ordered by position, fourth column is unnecessary
      # Fix: order by Y, then X, remove empty column
      C1 = C1[order(C1$X2, C1$X1),-4]
      names(C1) = c("minx","miny","prob_tumor")
      
      # Identify Cancer patch size
      C_range = (C1$minx[2] - C1$minx[1])
   }
   
   
   # =============================================================
   # Find, load, and reorder corresponding TIL Annotation file
   # =============================================================
   if(lymphFormatCSV){
      T1 = as.data.frame(
         readr::read_csv(
            paste(params$tilDir,tils[j], sep = "/"),
            col_names = T,
            col_type = cols()
         )
      )
      
      # Data is not ordered by position, fourth column is unnecessary
      # Fix: order by Y, then X, remove empty column
      T1 = T1[order(T1$miny,T1$minx),]
      T_range = T1$width[1]
   } else {
      T1 = as.data.frame(
         readr::read_table(
            paste(params$tilDir,tils[j], sep = "/"),
            col_names = F,
            col_type = cols()
         )
      )
      
      # Data is not ordered by position, fourth column is unnecessary
      # Fix: order by Y, then X, remove empty column
      T1 = T1[order(T1$X2,T1$X1),-4]
      names(T1) = c("minx","miny","prob_tils")
      
      # Identify TIL patch size
      T_range = (T1$minx[2] - T1$minx[1])
   }
   
   print(paste0("Sample ", count, ": ", C_range/T_range))
   # ~1.75
   
   ##=== Get LCM for scaling, for now, only get this for first WSI pair ===
   if(count==1){
      asFrac = MASS::fractions(signif(C_range/T_range, digits = 3)) ## for math facilitation, reduce to hundreths place
      numerator = as.integer(strsplit(attr(asFrac,"fracs"),"/")[[1]][1]) ## This can be a messy part if sample doesnt have good scalability
      denom = as.integer(strsplit(attr(asFrac,"fracs"),"/")[[1]][2])
   } else { ## Check for instability
      newSamp = MASS::fractions(signif(C_range/T_range, digits = 3))
      if(asFrac-newSamp > asFrac/2){
         errorName = canc[j]
         warning(paste0("Sample ", errorName, " has patch ratio of significantly different size than previous iterations (varied by > half of initial sample -- check)"))
      }
   }
   
   # =============================================================
   # Condense annotation files to [i,j] format and identify max scaling factor (initial ratio 7/4)
   # =============================================================
   
   ## Divide patch location by size and round up (Note: WSI-Infer does not report background, range may not reach 0 as min)
   T1[,c("minx","miny")] <- ceiling(T1[,c("minx","miny")]/T_range) 
   C1[,c("minx","miny")] <- ceiling(C1[,c("minx","miny")]/C_range)
   
   # Not all images are perfect squares/rectangles, pad empty locations with 0's to prevent shifts
   t1.maxi <-  max(T1$minx)
   t1.maxj <-  max(T1$miny)
   c1.maxi <-  max(C1$minx)
   c1.maxj <-  max(C1$miny)
   
   ## Ensure scaling fits larger samples (put in same feature space)
   big.max.i <- max(t1.maxi*denom,
                    c1.maxi*numerator)
   big.max.i <- (numerator*denom)*ceiling(big.max.i/(numerator*denom))
   
   t1.maxi <- ceiling(big.max.i/denom)
   c1.maxi <- ceiling(big.max.i/numerator)
   
   big.max.j <- max(t1.maxj*denom,c1.maxj*numerator)
   big.max.j <- (numerator*denom)*ceiling(big.max.j/(numerator*denom))
   t1.maxj <- ceiling(big.max.j/denom)
   c1.maxj <- ceiling(big.max.j/numerator)
   
   # =============================================================
   # Print blank padded matrix (account for non-square images)
   # =============================================================
   T2 <- matrix(0, nrow = t1.maxj, ncol = t1.maxi)
   C2 <- matrix(0, nrow = c1.maxj, ncol = c1.maxi)
   
   # =============================================================
   # Fill patches with scaled and padded prediction values - will be different dimensions - only fills in patches with predictions -- then convert to raster for interpolation
   # =============================================================
   for(el in 1:nrow(T1)){
      T2[T1$miny[el],T1$minx[el]] <- T1$prob_tils[el]
   }
   rT <- raster(T2)
   
   for(el in 1:nrow(C1)){
      C2[C1$miny[el],C1$minx[el]] <- C1$prob_tumor[el]
   }
   rC <- raster(C2)
   
   # =============================================================
   # Scale using num/denom factor to full size overlapping images - use raster for nearest-neighbor interpolation
   # =============================================================
   temp_dim <- raster(nrows = nrow(C2)*numerator,
                      ncols = ncol(C2)*numerator)
   
   crs(temp_dim) = NA
   extent(temp_dim) <- extent(c(0, 1, 0, 1))
   C_resized = (raster::resample(x = rC,y = temp_dim, method = 'ngb'))
   
   temp_dim <- raster(nrows = nrow(T2)*denom,
                      ncols = ncol(T2)*denom)
   crs(temp_dim) = NA
   extent(temp_dim) <- extent(c(0, 1, 0, 1))
   T_resized = (raster::resample(x = rT,y = temp_dim, method = 'ngb'))
   
   Cdat <- as.matrix(C_resized)
   #Cdat = round(Cdat)
   Tdat <- as.matrix(T_resized)
   #Tdat = round(Tdat)
   
   # =============================================================
   # Write thresholded images to png if desired
   # =============================================================
   if(params$writePNG == TRUE){
      ## Threshold predictions
      Cdat_thresh = (Cdat >= params$cancThresh)
      Tdat_thresh = (Tdat >= params$tilThresh)
      
      ## Make rgb matrix
      my.rgb <- abind(Tdat_thresh, ## R matrix
                      matrix(0,  ## empty G matrix
                             nrow = nrow(Cdat_thresh),
                             ncol = ncol(Cdat_thresh)),
                      Cdat_thresh, ## B matrix
                      along = 3)
      
      ##write files
      if(!dir.exists(paste0(params$outputDir,"/PNGs"))){
         dir.create(paste0(params$outputDir,"/PNGs"))
      }
      png::writePNG(target = paste0(params$outputDir,"/PNGs/",
                                    sampNames[j],'.thresh.png'),
                    image = my.rgb)
   }
   # =============================================================
   # Extract cancer patches and percent TIL patches
   # =============================================================
   Cancer_patches = sum(Cdat >= params$cancThresh) ## How many predicted canc?
   Til_patches = sum(Tdat >= params$tilThresh) # How many predicted Lymph
   Cancer_patches_with_til = sum(Cdat >= params$cancThresh &
                                    Tdat >= params$tilThresh) # how many predicted both
   
   ## Arrange
   output = data.frame(slideID = sampNames[j],
                       n_Canc_patch = Cancer_patches,
                       n_TIL_patch = Til_patches,
                       n_TIL_patch_overlap = Cancer_patches_with_til,
                       percent_pos = Cancer_patches_with_til / Cancer_patches,
                       patch_ratio = C_range/T_range,
                       stringsAsFactors = F)
   
   ## Append
   percent_calls[j,] = output
   # }
}

## trim .csv from slideID names if it exists
percent_calls$slideID = gsub(pattern = "\\.csv",
                             replacement = "",
                             x = percent_calls$slideID)

## Same for the preceding "prediction-"
percent_calls$slideID = gsub(pattern = "prediction-",
                             replacement = "",
                             x = percent_calls$slideID)

## ==== Save here after loop in case a following line breaks it, don't want to lose progress ====
write.csv(x = percent_calls,
          file = paste(params$outputDir,params$outputFile, sep = "/"),
          row.names = F)

## Calculate additional metrics
percent_calls$scaled_PP = percent_calls$percent_pos/sd(percent_calls$percent_pos, na.rm = T)
percent_calls$TIL_Class = cut(percent_calls$percent_pos,
                              breaks = c(0, mean(percent_calls$percent_pos, na.rm = T), 1),
                              labels = c("Low","High"),
                              include.lowest = T)

# =============================================================
# Save compiled data frame as .csv
write.csv(x = percent_calls,
          file = paste(params$outputDir,params$outputFile, sep = "/"),
          row.names = F)

# =============================================================
# Write file of histogram
tmp = ggplot(data = percent_calls, aes(x = percent_pos, color = TIL_Class, fill = TIL_Class)) + geom_histogram(binwidth = .1)

pdf(file = paste0(params$outputDir,"/", "invasion_histogram.pdf"), width = 8, height = 6)
tmp
suppressMessages(dev.off())

# =============================================================
# If provided, join sample level information to TIL output, replace output file with updated output
if(params$sampInfo != 'blank'){
   percent_calls = dplyr::full_join(percent_calls,sampInfo, by = "slideID")
   write.csv(x = percent_calls,
             file = paste(params$outputDir,params$outputFile, sep = "/"),
             row.names = F)
}
