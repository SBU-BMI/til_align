#!/usr/bin/env -S Rscript --vanilla

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
              sampInfo = args[6],
              outputFile = args[7],
              outputDir = args[8],
              writePNG = args[9],
              subtypes = args[10])

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

## ==== Check subtypes request ====
if(params$subtypes %in% c("true", "TRUE","t","T",T, TRUE)){
   params$subtypes = TRUE
} else { params$subtypes = FALSE }


## ==== Update and print thresholds ====
params$tilThresh = as.numeric(params$tilThresh)
params$cancThresh = as.numeric(params$cancThresh)
print("TIL Algorithm (Threshold): Frontiers InceptionV4: 0.1")

print("=========== Params after R parsing, if any misalignment please check your flags ===========")
print(params)

### ==== For testing!! comment out ====
# setwd("~/Desktop/tilalign/heathCheck/")
# params = list(tilDir = "./til_inception/",
#               #tilDir = "/datadrive/shared/image_analysis/SEERky_test/til",
#               tilThresh = 0.1,
#               cancDir = "./tumor/heatmap_txt/",
#               #cancDir = "/datadrive/shared/image_analysis/SEERky_test/tumor",
#               cancThresh = 0.5,
#               sampInfo = "data.csv",
#               #sampInfo = "/datadrive/shared/image_analysis/SEERky_test/clinical_data.csv",
#               #sampFile = "/datadrive/shared/image_analysis/ML_output/extdata/flexible_testDir/sampFileShort.csv",
#               outputFile = "Percent_Invasion.csv",
#               outputDir = "results",
#               writePNG = T,
#               subtypes = T)
### ==== Above is for testing!! comment out ====

### ==== Read in files ====
tils <- sort(list.files(params$tilDir), decreasing = TRUE) # data dirs can be hardcoded or relative
canc <- sort(list.files(params$cancDir),decreasing = TRUE) # data dirs can be hardcoded or relative

## ==== Drop extraneous files ("color-*" and "*.low_res") ====
tils = tils[grep("^prediction", tils)]
writeLines(" . . . Dropping low_res and color- files . . . ")
if(any(grepl("low_res", tils))){
   tils = tils[-grep("low_res", tils)]
}

canc = canc[grep("^prediction", canc)]
if(any(grepl("low_res", canc))){
   canc = canc[-grep("low_res", canc)]
}

## ==== Check for missing file pairs (if there is a tumor or lymph prediction but not the other) ====
writeLines(" . . . Checking for tumor/lymph pairs . . . ")

if(length(setdiff(tils,canc)) > 0 | length(setdiff(canc,tils)) > 0){
   warning("Some cancer-lymph predictions are missing pairs, skipping non-paired samples")
   ## Check for Lymph is present but Cancer is missing
   if(length(setdiff(tils,canc)) > 0){
      writeLines(paste0(length(setdiff(tils,canc))," out of ", length(tils)," provided samples have Lymph predictions but no Cancer predictions. They will be dropped"))
      writeLines(paste0("Lymph samples missing Cancer pairs (up to first 6) are: 
                     ",
                        head(setdiff(tils,canc))))
   }
   if(length(setdiff(canc,tils))){
      writeLines(paste0(length(setdiff(canc,tils))," out of ", length(canc)," provided samples have Cancer predictions but no Lymph predictions. They will be dropped"))
      writeLines(paste0("Cancer samples missing Lymph pairs (up to first 6) are: 
                     ",
                        head(setdiff(canc,tils))))
   }
   tils = intersect(tils,canc)
   canc = intersect(canc,tils)
} else {
   writeLines(" . . . All files have pairs . . . ")
}

## ==== After removing unequal pairs, make sure that files are a 1 to 1 ordered match. No duplicates or differences whatsoever ====
if(length(tils)==0){
   stop("No predictions had exact pairs. Please ensure lymph and cancer pairs have the exact same name.")
}
if(!all.equal(tils,canc)){
   stop("Supplied files do not match. Must be same names. Please adjust directories and try again.")
}

## ==== Check for read and align format using first entry in both cancer and lymph file lists ====
lymphFormatCSV = grepl("csv", tils[1])
cancFormatCSV = grepl("csv", canc[1])

## ==== If they are different formats, warn but continue. This is likely not an optimal implementation ====
if(lymphFormatCSV != cancFormatCSV){
   warning("Prediction formats for cancer and lymphocyte are different, please double check this is intentional")
}

## ==== Print alignment trajectory ====
if(lymphFormatCSV & cancFormatCSV){
   if(params$subtypes){
      message(" . . . Inputs detected as WSInfer with Subtypes . . . ")
   } else {
      message(" . . . Inputs detected as WSInfer without Subtypes . . . ")
   }
} else {
   if(params$subtypes){
      message(" . . . Inputs detected as Original with Subtypes . . . ")
   } else {
      message(" . . . Inputs detected as Original without Subtypes . . . ")
   }
}

## ==== read in sampInfo csv ahead of time if provided ====
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

### ==== Define functions to facilitate tryCatch approaches ====
loadAndSort <- function(whichPred){
   ## ============================================
   ## === Run Cancer prediction load and parse ===
   ## ============================================
   if(whichPred == "Canc"){
      out <- tryCatch({
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
            if(params$subtypes){ ## Whitespace sep + w/ subtypes
               C1 = as.data.frame(
                  readr::read_table(
                     paste(params$cancDir,canc[j], sep = "/"),
                     col_names = T,
                     col_type = cols()
                  )
               )
               colnames(C1)[1:2] = c("minx","miny")
               colnames(C1)[3:ncol(C1)] = paste0("prob_",colnames(C1)[3:ncol(C1)])
               C1 = C1[order(C1$miny, C1$minx),]
            } else { ## Whitespace + no subtypes
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
            }
         }
         ## If inputs are flagged to include subtype predictions
         if(params$subtypes){
            ## This needs to be flexible for diff tumor types
            predIndices = grep("prob", colnames(C1), ignore.case = T)
            C1$Call = apply(C1[,predIndices], 1,which.max) ## call each patch
            C1$Call[rowSums(C1[,predIndices]) == 0] = 0 ## If all probs are 0, it is background (type 0)
            C1$Call = as.integer(C1$Call)
         }
         C1
      },
      error=function(cond) {
         message(paste0("Error upon C1 Loading. Original error message:\n",
                        cond,
                        "\nSkipping to next sample")
         )
         # Choose a return value in case of error
         return(NULL)
      })
   } else {
      out <- tryCatch({
         ## ============================================
         ## === Run Lymph prediction load and parse ===
         ## ============================================
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
            T1
         }
      },
      error=function(cond) {
         message(paste0("Error upon T1 Loading. Original error message:\n",
                        cond,
                        "\nSkipping to next sample")
         )
         # Choose a return value in case of error
         return(NULL)
      })
   }
   return(out)
}

makeDict <- function(C1){
   out <- tryCatch({
      ## Make a dictionary for subtype to numeric conversion
      callDict = list(background = 0)
      predIndices = grep("prob", colnames(C1), ignore.case = T)
      for(tissueType in predIndices){
         trimmed = gsub("prob_","",colnames(C1)[tissueType]) ## extract just tissue type (drop "prob_")
         callDict[[trimmed]] = length(callDict) ## Count number of tissue types and assign pair
      }
      callDict
   },
   error=function(cond) {
      return(NULL)
   })
   return(out)
}

rasterAndResize <- function(whichPred){
   if(whichPred == "Canc"){
      out <- tryCatch({
         # =============================================================
         # Print blank padded matrix (account for non-square images)
         # =============================================================
         C2 <- matrix(0, nrow = c1.maxj, ncol = c1.maxi)
         
         # =============================================================
         # Fill patches with scaled and padded prediction values - will be different dimensions - only fills in patches with predictions -- then convert to raster for interpolation
         # =============================================================
         if(params$subtypes){ ## Grab Predicted Call
            for(el in 1:nrow(C1)){
               C2[C1$miny[el],C1$minx[el]] <- C1$Call[el]
            }
         } else {
            for(el in 1:nrow(C1)){ ## Grab tumor probability
               C2[C1$miny[el],C1$minx[el]] <- C1$prob_tumor[el]
            }
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
         as.matrix(C_resized)
      },
      error=function(cond) {
         message(paste0("Error upon Cdat resizing. Original error message:\n",
                        cond,
                        "\n Skipping to next sample")
         )
         
         
         # Choose a return value in case of error
         return(NULL)
      })
   } else{
      out <- tryCatch({
         # =============================================================
         # Print blank padded matrix (account for non-square images)
         # =============================================================
         T2 <- matrix(0, nrow = t1.maxj, ncol = t1.maxi)
         
         # =============================================================
         # Fill patches with scaled and padded prediction values - will be different dimensions - only fills in patches with predictions -- then convert to raster for interpolation
         # =============================================================
         for(el in 1:nrow(T1)){
            T2[T1$miny[el],T1$minx[el]] <- T1$prob_tils[el]
         }
         rT <- raster(T2)
         
         # =============================================================
         # Scale using num/denom factor to full size overlapping images - use raster for nearest-neighbor interpolation
         # =============================================================
         temp_dim <- raster(nrows = nrow(T2)*denom,
                            ncols = ncol(T2)*denom)
         crs(temp_dim) = NA
         extent(temp_dim) <- extent(c(0, 1, 0, 1))
         T_resized = (raster::resample(x = rT,y = temp_dim, method = 'ngb'))
         as.matrix(T_resized)
      },
      error=function(cond) {
         message(paste0("Error upon Tdat resizing. Original error message:\n",
                        cond,
                        "\n Skipping to next sample")
         )
         return(NULL)
      })
   }
   return(out)
}

writePNGs <- function(){
   ## Make dir if needed
   if(!dir.exists(paste0(params$outputDir,"/PNGs"))){
      dir.create(paste0(params$outputDir,"/PNGs"))
   }
   if(params$subtypes){
      plots2 = list()
      ## Print a side-by-side PNG of thresholded overlay and subtyped patches
      ## Threshold predictions
      Cdat_thresh = Cdat != callDict$background & Cdat != callDict$benign
      Tdat_thresh = (Tdat >= params$tilThresh)
      
      Cdat_thresh_melted = reshape2::melt(Cdat)
      colnames(Cdat_thresh_melted) = c("x","y","Canc_type")
      
      ## Parse Canc Type 
      Cdat_thresh_melted$Canc_type = factor(Cdat_thresh_melted$Canc_type,
                                            levels = 0:6,
                                            labels = names(callDict))
      Cdat_thresh_melted$Canc_presence = !(Cdat_thresh_melted$Canc_type %in% c("benign","background"))
      
      Tdat_thresh_melted = reshape2::melt(Tdat_thresh)
      
      to_plot = data.frame(x = Cdat_thresh_melted$x, 
                           y = Cdat_thresh_melted$y, 
                           presence = factor("No Tissue", levels = c("No Tissue", "Benign","Cancer","TIL","Both")))
      
      to_plot$presence[Cdat_thresh_melted$Canc_type == "benign"] = "Benign"
      to_plot$presence[Cdat_thresh_melted$Canc_presence & !Tdat_thresh_melted$value] = "Cancer"
      to_plot$presence[!Cdat_thresh_melted$Canc_presence & Tdat_thresh_melted$value] = "TIL"
      to_plot$presence[Cdat_thresh_melted$Canc_presence & Tdat_thresh_melted$value] = "Both"
      
      plots2[["invasion"]] = ggplot(to_plot, aes(x=y,y=x,color=presence))+geom_point() +
         scale_color_manual(values = c("black", "grey","yellow","red","red")) +
         theme(axis.title = element_blank(),
               axis.text = element_blank(),
               legend.title = element_text(size = 14)) +
         guides(color=guide_legend(title="Patch Type"))
      
      ## Also print just a color coded map of subtypes
      plots2[["types"]] = ggplot(C1, aes(x = minx, y = miny, 
                                         fill = as.character(Call)), color = as.character(Call)) + 
         geom_point(shape = 22, size = 3) + #theme_pubr() +
         scale_fill_manual(values=palette(),
                           labels = names(callDict)) +
         theme(axis.title = element_blank(),
               axis.text = element_blank(),
               legend.title = element_text(size = 10)) +
         guides(fill=guide_legend(title="Tissue Type"))
      
      tmp = ggpubr::ggarrange(plotlist = plots2, ncol = 1)
      png( paste0(params$outputDir,"/PNGs/", tils[j],'subtypes.png') )
      print(tmp)
      dev.off()
   } else {
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
      png::writePNG(target = paste0(params$outputDir,"/PNGs/",
                                    tils[j],'.thresh.png'),
                    image = my.rgb)
   }
}

calculateAlignment <- function(){
   out <- tryCatch({
      if(params$subtypes){
         ## ==== First, overall calculations ====
         Tissue_patches = Cdat != callDict$background ## Boolean: patch locations predicted to be "NOT BACKGROUND"
         
         nonCancVals = c(grep("background",names(callDict)),
                         grep("benign",names(callDict)),
                         grep("stroma",names(callDict)),
                         grep("necrosis",names(callDict)))-1
         
         Cancer_patches = (!(Cdat %in% nonCancVals)) ## Boolean: patch locations NEITHER background nor any non-cancer (varies by organ)
         
         Til_patches = Tdat >= params$tilThresh # Boolean: Patch locations predicted to contain Lymphocytes
         
         Cancer_patches_with_til = Cancer_patches & Til_patches # Boolean: Patch locations predicted Cancer & Lymphocyte
         
         ## Arrange
         output = data.frame(slideID = tils[j],
                             n_Tissue_patch = sum(Tissue_patches),
                             n_Canc_patch = sum(Cancer_patches),
                             n_TIL_patch = sum(Til_patches),
                             n_TIL_patch_overlap = sum(Cancer_patches_with_til),
                             percent_pos = sum(Cancer_patches_with_til) / sum(Cancer_patches),
                             benign_percent_of_tissue_patches = (sum(Tissue_patches)-sum(Cancer_patches))/sum(Tissue_patches),
                             patchRatio = C_range/T_range,
                             stringsAsFactors = F)
         
         ## ==== Now, append subtype specific metrics ====
         for(tissueType in names(callDict)){
            if(tissueType != "background"){
               ## Calculate metrics
               subtype_patches = sum(Cdat == callDict[[tissueType]]) ## How many predicted particular cancer type
               subtype_patches_with_til = sum(Cdat == callDict[[tissueType]] & ## How many of type X have lymph invasion
                                                 Tdat >= params$tilThresh)
               percent_pos = subtype_patches_with_til/subtype_patches ## Percentage of line above
               
               if(tissueType %in% names(callDict)[callDict %in% nonCancVals]){
                  percent_patch_subtype = NA
               } else {
                  percent_patch_subtype = subtype_patches / output$n_Canc_patch ## How many cancer patches are type X
               }
               
               ## Append metrics
               output[,paste0(tissueType,"_n_patch")] = subtype_patches
               output[,paste0(tissueType,"_percent_of_canc_patches")] = percent_patch_subtype
               output[,paste0(tissueType,"_TIL_overlap")] = subtype_patches_with_til
               output[,paste0(tissueType,"_percent_pos")] = percent_pos
            }
         }
         output$slide_level_call = names(which.max(output[,grep("percent_of_canc",colnames(output))]))
         output$slide_level_call = stringr::str_split_fixed(output$slide_level_call,"_",2)[[1]]
         output
      } else {
         Cancer_patches = sum(Cdat >= params$cancThresh) ## How many predicted canc?
         Til_patches = sum(Tdat >= params$tilThresh) # How many predicted Lymph
         Cancer_patches_with_til = sum(Cdat >= params$cancThresh &
                                          Tdat >= params$tilThresh) # how many predicted both
         
         ## Arrange
         data.frame(slideID = tils[j],
                    n_Canc_patch = Cancer_patches,
                    n_TIL_patch = Til_patches,
                    n_TIL_patch_overlap = Cancer_patches_with_til,
                    percent_pos = Cancer_patches_with_til / Cancer_patches,
                    patch_ratio = C_range/T_range,
                    stringsAsFactors = F)
      }
   },
   error=function(cond) {
      message(paste0("Error upon alignment calculation. Original error message:\n",
                     cond,
                     "\n Returning NA and starting next sample")
      )
      # Choose a return value in case of error
      return(data.frame(slideID = tils[j],
                        n_Canc_patch = NA,
                        n_TIL_patch = NA,
                        n_TIL_patch_overlap = NA,
                        percent_pos = NA,
                        patch_ratio = NA,
                        stringsAsFactors = F)
      )
   })
   return(out)
}

## ==== Begin Alignment Loop ====
writeLines(paste0(" . . . Running alignment and analyses for ", length(canc), " samples . . . "))

count = 0
for(j in 1:length(canc)){
   writeLines(paste0("============================================= \n",
                     "# Processing ", tils[j], " # \n", 
                     "============================================="))
   count = count + 1
   # =============================================================
   # Load in and reorder Cancer Annotation File
   # =============================================================
   C1 <- loadAndSort("Canc")
   if(is.null(C1)){
      percent_calls[j,2:ncol(percent_calls)] = NA
      next
   }
   
   # Identify and log Cancer patch size
   C_range = (C1$minx[2] - C1$minx[1])
   
   ## =============================================================
   # Make a dictionary
   # =============================================================
   if(params$subtypes & count == 1){ ## Only need to make this once
      callDict = makeDict(C1)   
   }
   
   # =============================================================
   # Find, load, and reorder corresponding TIL Annotation file
   # =============================================================
   T1 <- loadAndSort("TIL")
   if(is.null(T1)){
      percent_calls[j,2:ncol(percent_calls)] = NA
      next
   }
   
   # Identify and log TIL patch size
   T_range = (T1$minx[2] - T1$minx[1])
   
   writeLines(paste0("Patch Ratio (Canc/Til): ", C_range/T_range, "\n"))
   # ===============================================================================
   ## == Get LCM for scaling (Assign with sample 1, check for stability with rest ==
   # ===============================================================================
   if(count==1){ ## Assign
      firstFrac = MASS::fractions(signif(C_range/T_range, digits = 3))
      # =================================================================================
      ## == Explicitly define based on proximity to best LCM for memory considerations ==
      # =================================================================================
      if(abs(C_range/T_range - 1.75) < .02){ ## BRCA, PRAD, READ, COAD are roughly 1.75 (7/4)
         numerator = as.integer(7)
         denom = as.integer(4)
      }  else if(abs(C_range/T_range - 3.5) < .02){ ## LUAD is roughly 3.5 (7/2)
         numerator = as.integer(7)
         denom = as.integer(2)
      } else if(abs(C_range/T_range - 10.5) < .02){ ## PAAD is roughly 10.5 (21/2)
         numerator = as.integer(21)
         denom = as.integer(2)
      } else { ## If ratio is much different from the above, do it manually.
         firstFrac = MASS::fractions(signif(C_range/T_range, digits = 3)) ## Calc fraction (for math facilitation, reduce to hundreths place)
         numerator = as.integer(strsplit(attr(firstFrac,"fracs"),"/")[[1]][1]) ## This can be a messy part if sample doesnt have good scalability
         denom = as.integer(strsplit(attr(firstFrac,"fracs"),"/")[[1]][2])
      }
      # ======================
      ## == Stability check ==
      # ======================
   } else { 
      ## Possible future change, assign ratio each time, allows for single run of multiple tumor types.
      newSamp = MASS::fractions(signif(C_range/T_range, digits = 3))
      if(firstFrac-newSamp > firstFrac/2){
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
   # Resize and rescale predictions. If any errors, skip to next samp
   # =============================================================
   Tdat <- rasterAndResize("TIL")
   if(is.null(Tdat)){
      percent_calls[j,2:ncol(percent_calls)] = NA
      next
   }
   
   Cdat <- rasterAndResize("Canc")
   if(is.null(Cdat)){
      percent_calls[j,2:ncol(percent_calls)] = NA
      next
   }
   
   # =============================================================
   # Write thresholded images to png if desired
   # =============================================================
   if(params$writePNG){
      try(writePNGs())
   }
   
   # =============================================================
   # Extract cancer patches and percent TIL patches
   # =============================================================
   if(count==1){
      percent_calls <- calculateAlignment()
   } else {
      percent_calls[j,] <- calculateAlignment()
   }
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
# If provided, join sample level information to TIL output, replace output file with updated output
if(params$sampInfo != 'blank'){
   percent_calls = dplyr::full_join(percent_calls,sampInfo, by = "slideID")
   write.csv(x = percent_calls,
             file = paste(params$outputDir,params$outputFile, sep = "/"),
             row.names = F)
}

# =============================================================
# Write file of histogram
tmp = ggplot(data = percent_calls, aes(x = percent_pos, color = TIL_Class, fill = TIL_Class)) + geom_histogram(binwidth = .1)

pdf(file = paste0(params$outputDir,"/", "invasion_histogram.pdf"), width = 8, height = 6)
tmp
suppressMessages(dev.off())
