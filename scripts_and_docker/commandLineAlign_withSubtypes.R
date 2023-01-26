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
              writePNG = args[9])

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
params = list(tilDir = "./example_withSubtypes/tilPreds",
              tilThresh = 0.1,
              cancDir = "./example_withSubtypes/cancPreds",
              cancThresh = 0.5,
              #sampFile = "/datadrive/shared/image_analysis/ML_output/extdata/flexible_testDir/sampFileShort.csv",
              outputFile = "Percent_Invasion.csv",
              outputDir = "outputs",
              writePNG = T)
### ==== Above is for testing!! comment out ====


## ==== Read in files ====
tils <- sort(list.files(params$tilDir), decreasing = TRUE) # data dirs can be hardcoded or relative
canc <- sort(list.files(params$cancDir),decreasing = TRUE) # data dirs can be hardcoded or relative
canc_short = stringr::str_split_fixed(canc, "\\.",2)[,1]
sampNames = tils
if(!all.equal(tils,canc)){
   warning("Supplied files do not match. Must be same names and order")
}
if(length(tils) != length(canc)){
   warning("Unequal number of til prediction files and cancer prediction files, check your inputs")
}

# ## ==== Align ====
count = 0
for(j in 1:length(canc)){
   count = count + 1
   # =============================================================
   # Load in and reorder Cancer Annotation File
   # =============================================================
   ## True file will be csv - comment for now
   C1 = as.data.frame(
      readr::read_csv(
         paste(params$cancDir,canc[j], sep = "/"),
         col_names = T,
         col_type = cols()
      )
   )
   
   ## Old data format
   # C1 = as.data.frame(
   #    readr::read_delim(
   #       paste(params$cancDir,canc[j], sep = "/"),
   #       col_names = T,
   #       col_type = cols()
   #    )
   # )
   
   ## This needs to be flexible for diff tumor types
   predIndices = grep("prob", colnames(C1), ignore.case = T)
   C1$Call = apply(C1[,predIndices], 1,which.max) ## call each patch
   C1$Call[rowSums(C1[,predIndices]) == 0] = 0 ## If all probs are 0, it is background (type 0)
   C1$Call = as.integer(C1$Call)
   
   ## Make a dictionary for subtype to numeric conversion
   callDict = list(background = 0)
   for(tissueType in predIndices){
      trimmed = gsub("prob_","",colnames(C1)[tissueType]) ## extract just tissue type (drop "prob_")
      callDict[[trimmed]] = length(callDict) ## Count number of tissue types and assign pair
   }
   
   # Order data as required (y then x)
   C1 = C1[order(C1$miny, C1$minx),]
   C_range = C1$width[1] ## WSIInfer provides patch size, should we include square check? Seems unnecessary
   
   # =============================================================
   # Find, load, and reorder corresponding TIL Annotation file
   # =============================================================
   # Input will ultiamtely be csv
   T1 = as.data.frame(
      readr::read_csv(
         paste(params$tilDir,tils[j], sep = "/"),
         col_names = T,
         col_type = cols()
      )
   )
   
   ## Old data type, for testing
   # T1 = as.data.frame(
   #    readr::read_delim(
   #       paste(params$tilDir,tils[j], sep = "/"),
   #       col_names = F,
   #       col_type = cols()
   #    )
   # )
   
   # Order data by y then x
   T1 = T1[order(T1$miny,T1$minx),]
   T_range = T1$width[1]
   
   print(paste0("Sample ", count, ": ", C_range/T_range))
   # ~1.75
   
   ## ==== Get LCM for scaling, for now, only get this for first WSI pair ====
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
   T1[,2:3] <- ceiling(T1[,2:3]/T_range) ## col 1 is now slideID, x and y are 2 and 3 respectively
   C1[,2:3] <- ceiling(C1[,2:3]/C_range)
   
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
      T2[T1$miny[el],T1$minx[el]] <- T1$prob_tils[el] ## Same as other align code
   }
   rT <- raster(T2)
   
   for(el in 1:nrow(C1)){
      C2[C1$miny[el],C1$minx[el]] <- C1$Call[el] ## Diff from other code, grabs the maximum subtype prediction
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
   Tdat <- as.matrix(T_resized)
   
   # =============================================================
   # Write overlaid and thresholded (for til) images to png if desired -- NOT READY YET, need to decide on approach
   # =============================================================
   # if(params$writePNG == TRUE){ ## gotta figure out best way to write these
   #    ## Threshold predictions
   #    Tdat_thresh = (Tdat >= params$tilThresh)
   #    
   #    ## Make rgb matrix
   #    my.rgb <- abind(Tdat_thresh, ## R matrix
   #                    matrix(0,  ## empty G matrix
   #                           nrow = nrow(Tdat_thresh),
   #                           ncol = ncol(Tdat_thresh)),
   #                    Cdat, ## B matrix
   #                    along = 3)
   #    
   #    ##write files
   #    if(!dir.exists(paste0(params$outputDir,"/PNGs"))){
   #       dir.create(paste0(params$outputDir,"/PNGs"))
   #    }
   #    png::writePNG(target = paste0(params$outputDir,"/PNGs/",
   #                                  sampNames[j],'.thresh.png'),
   #                  image = my.rgb)
   # }
   # 
   # png::writePNG(target = 'test.png',
   #               image = my.rgb)
   
   
   
   # =============================================================
   # Extract cancer patches and percent TIL patches
   # =============================================================
   ## ==== First, overall calculations ====
   Tissue_patches = sum(Cdat != callDict$background) ## How many predicted to be "NOT BACKGROUND"
   
   Cancer_patches = sum(!(Cdat %in% c(callDict$background, callDict$benign))) ## How many neither background nor benign
   
   Til_patches = sum(Tdat >= params$tilThresh) # How many predicted Lymph
   
   Cancer_patches_with_til = sum(sum(!(Cdat %in% c(callDict$background, callDict$benign))) &
                                    Tdat >= params$tilThresh) # how many predicted both
   
   ## Arrange
   output = data.frame(slideID = sampNames[j],
                       n_Tissue_patch = Tissue_patches,
                       n_Canc_patch = Cancer_patches,
                       n_TIL_patch = Til_patches,
                       n_TIL_patch_overlap = Cancer_patches_with_til,
                       percent_pos = Cancer_patches_with_til / Cancer_patches,
                       patchRatio = C_range/T_range,
                       stringsAsFactors = F)
   
   ## ==== Now, append subtype specific metrics ====
   for(tissueType in names(callDict)){
      ## Calculate metrics
      Cancer_patches = sum(Cdat == callDict[[tissueType]]) ## How many predicted particular cancer type
      Cancer_patches_with_til = sum(Cdat == callDict[[tissueType]] & ## How many of type X have lymph invasion
                                       Tdat >= params$tilThresh)
      percent_pos = Cancer_patches_with_til/Cancer_patches ## Percentage of line above
      
      if(tissueType == "background"){
         percent_patch_subtype = NA
      } else {
         percent_patch_subtype = Cancer_patches / output$n_Canc_patch ## How many cancer patches are type X
      }
      
      ## Append metrics
      output[,paste0(tissueType,"_n_patch")] = Cancer_patches
      output[,paste0(tissueType,"_percent_of_canc_patches")] = percent_patch_subtype
      output[,paste0(tissueType,"_TIL_overlap")] = Cancer_patches_with_til
      output[,paste0(tissueType,"_percent_pos")] = percent_pos
   }
   
   # move outside of loop
   if(count == 1){
      # Initialize empty data frame for output
      percent_calls = output
   } else{
      percent_calls = rbind(percent_calls, output)
   }
}

## ==== Save here after loop in case a following line breaks it, don't want to lose progress ====
write.csv(x = percent_calls,
          file = paste(params$outputDir,params$outputFile, sep = "/"),
          row.names = F)

## Calculate additional metrics for overall classifications
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


