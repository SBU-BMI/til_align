#! /bin/bash

################################################################################
Help()
{
   # Display Help
   echo "This function takes as input two directories: one of cancer predications and one of lymphocyte predictions. It will align the predictions and generate percent invasion for all included samples, as well as call class as high or low."
   echo
   echo "Syntax: callAlign [-a|t|T|c|C|s|o|O|w|h]"
   echo "options:"
   echo "a     If you know the lymph prediction version you used, listing here will automatically define the Lymph threshold. Acceptable options are r (resnet), v (VGG16), or i (Inception)"
   echo "t     (Default: tilPreds) path to lymph predictions"
   echo "T     (Default: 0.5; Overwritten if -a is provided) Threshold for lymph calling at a patch level"
   echo "c     (Default: cancPreds) path to cancer predictions"
   echo "C     (Default: 0.5) Threshold for cancer calling at a patch level "
   echo "s     (OPTIONAL) Only necessary if file names for Lymph and Cancer preds do not exactly match. Path to dataset csv with 3 columns: 1) Sample name, 2) TIL file name, 3) Canc file name"
   echo "o     (Default: Percent_Invasion.csv) Name for csv output"
   echo "O     (Default: outputs) Path to output directory."
   echo "w     Flag to write PNGs of Tumor-TIL maps. Makes dir PNGs/ within output directory"
   echo "h     Prints this help"
}
################################################################################
## Set Defaults
algorithm='blank'
tilDir='/data/tilPreds'
tilThresh=0.5 
cancDir='/data/cancPreds'
cancThresh=0.5
sampFile='blank'
outputFile='Percent_Invasion.csv'
outputDir='/data/outputs'
writePNG=false

while getopts :ha:t:T:c:C:s:o:O:w flag
do
    case "${flag}" in
        h) Help
        exit;;
        a) algorithm=${OPTARG};;
        t) tilDir=${OPTARG};;
        T) tilThresh=${OPTARG};;
        c) cancDir=${OPTARG};;
        C) cancThresh=${OPTARG};;
        s) sampFile=${OPTARG};;
        o) outputFile=${OPTARG};;
        O) outputDir=${OPTARG};;
        w) writePNG=true;;
    esac
done
echo " ======================================================================================== "
echo "Flags specified as the following (BASH parsed): "
echo "algorithm: $algorithm";
echo "tilDir: $tilDir";
echo "tilThresh: $tilThresh";
echo "cancDir: $cancDir";
echo "cancThresh: $cancThresh";
echo "sampFile: $sampFile";
echo "outputFile: $outputFile";
echo "outputDir: $outputDir";
echo "writePNG: $writePNG";
echo " ======================================================================================== "
## Pass these flags to R script for alignment
Rscript commandLineAlign.R $algorithm $tilDir $tilThresh $cancDir $cancThresh $sampFile $outputFile $outputDir $writePNG
