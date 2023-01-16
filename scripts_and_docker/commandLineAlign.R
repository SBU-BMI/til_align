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

## ==== Update Thresholds if algorithm is specified ====
if(params$algorithm %in% c("RESNET","ResNet", "resnet", "rn","r","R")){
   params$tilThresh = 0.5
   print("Lymphocyte ALgorithm listed as ResNet, TIL Thresh 0.5 applied")
} else if(params$algorithm %in% c("VGG","Vgg", "vgg", "V","V")){
   params$tilThresh = 0.5
   print("TIL ALgorithm listed as VGG, TIL Thresh 0.5 applied")
} else if (params$algorithm %in% c("INCEPTION","Inception", "Inception","INCEPT", "Incept","incept","I", "i")){
   params$tilThresh = 0.1
   print("TIL ALgorithm listed as Inception, TIL Thresh 0.1 applied")
} else{
   params$tilThresh = as.numeric(params$tilThresh)
   params$cancThresh = as.numeric(params$cancThresh)
   print("TIL Algorithm unspecified, using supplied thresholds (Default: 0.5)")
}
print("=========== Params after R parsing, if any misalignment please check your flags ===========")
print(params)


## ==== Read in files ====
sampFileLoc = as.character(params$sampFile)
if(params$sampFile == "blank"){ ## Files MUST have exact same names in different directories
   tils <- sort(list.files(params$tilDir), decreasing = TRUE) # data dirs can be hardcoded or relative
   canc <- sort(list.files(params$cancDir),decreasing = TRUE) # data dirs can be hardcoded or relative
   sampNames = tils
   if(!all.equal(tils,canc)){
      warning("Supplied files do not match. Must be same names and order")
   }
} else { ## sampFile must include 3 columns: 1 is sample name, 2 is TIL pred file, 3 is Canc pred file. Files must be in same order
   sampFile = as.matrix(suppressMessages(readr::read_csv(file = sampFileLoc)))
   sampNames = as.character(sampFile[,1])
   tils = as.character(sampFile[,2])
   canc = as.character(sampFile[,3])
}
if(length(tils) != length(canc)){
   warning("Unequal number of til prediction files and cancer prediction files, check your inputs")
}

## For testing
# params = list(tilDir = "/datadrive/shared/image_analysis/ML_output/extdata/flexible_testDir/TIL_preds/",
#               tilThresh = 0.1,
#               cancDir = "/datadrive/shared/image_analysis/ML_output/extdata/flexible_testDir/Canc_preds/",
#               cancThresh = 0.5,
#               sampFile = "/datadrive/shared/image_analysis/ML_output/extdata/flexible_testDir/sampFileShort.csv",
#               outputFile = "Percent_Invasion.csv",
#               outputDir = "outputs",
#               writePNG = T)

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
   names(C1) = c("X_loc","Y_loc","Prediction")

   # Identify patch size (typically 348-352 range)
   C_range = (C1$X_loc[2] - C1$X_loc[1])

   # =============================================================
   # Find, load, and reorder corresponding TIL Annotation file
   # =============================================================
   T1 = as.data.frame(
      readr::read_table(
         paste(params$tilDir,til[j], sep = "/"),
         col_names = F,
         col_type = cols()
      )
   )

   # Data is not ordered by position, fourth column is unnecessary
   # Fix: order by Y, then X, remove empty column
   T1 = T1[order(T1$X2,T1$X1),-4]
   names(T1) = c("X_loc","Y_loc","Prediction")

   # Identify TIL patch size (typically 198-202 range)
   T_range = (T1$X_loc[2] - T1$X_loc[1])

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
   T1[,1:2] <- ceiling(T1[,1:2]/T_range)
   C1[,1:2] <- ceiling(C1[,1:2]/C_range)

   # Not all images are perfect squares/rectangles, pad empty locations with 0's to prevent shifts
   t1.maxi <-  max(T1$X_loc)
   t1.maxj <-  max(T1$Y_loc)
   c1.maxi <-  max(C1$X_loc)
   c1.maxj <-  max(C1$Y_loc)

   ## Ensure scaling fits larger samples
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
   # Fill patches with scaled and padded prediction values - will be different dimensions
   # =============================================================
   for(el in 1:nrow(T1)){
      T2[T1$Y_loc[el],T1$X_loc[el]] <- T1$Prediction[el]
   }
   rT <- raster(T2)

   for(el in 1:nrow(C1)){
      C2[C1$Y_loc[el],C1$X_loc[el]] <- C1$Prediction[el]
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
   # Write thresholded images to png
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
                       patchRatio = C_range/T_range,
                       stringsAsFactors = F)

   ## Append
   percent_calls[j,] = output
   # }
}

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


