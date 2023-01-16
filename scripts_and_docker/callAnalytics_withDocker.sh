#! /bin/bash

################################################################################
Help()
{
   # Display Help
   echo "This function takes as input a csv file (output from alignment pipeline merged with variables of interest) and the column names of those columns of interest"
   echo
   echo "Syntax: callAlign [-p|s|t|c|v|h]"
   echo "options:"
   echo "p     (Default: Percent_Invasion) Path to csv file to analyze"
   echo "s     Flag to exclude survival analysis"
   echo "t     (Default: survTime) Column name in csv containing time to event"
   echo "c     (Default: survCensor) Column name in csv containing if patient is censored or not (0 = censored, 1 = not censored)"
   echo "h     Prints this help"
}
################################################################################
## Set Defaults
csvPath='Percent_Invasion.csv'
excludeSurv=false
survTime='survivalA'
survCensor='censorA.0yes.1no'
outputFormat='pdf_document'

while getopts :hp:st:c:f: flag
do
    case "${flag}" in
        h) Help
        exit;;
        p) csvPath=${OPTARG};;
        s) excludeSurv = true;;
        t) survTime=${OPTARG};;
        c) survCensor=${OPTARG};;
        f) outputFormat=${OPTARG};;
    esac
done
echo " ======================================================================================== "
echo "Flags specified as the following (BASH parsed): "
echo "csvPath: $csvPath";
echo "excludeSurv: $excludeSurv";
echo "survTime: $survTime";
echo "survCensor: $survCensor";
echo "outputFormat: $outputFormat";
echo " ======================================================================================== "

## Pass these flags to R script for alignment
Rscript renderWrapper.R $csvPath $excludeSurv $survTime $survCensor $outputFormat
