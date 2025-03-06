
# Eva Caamano Gutierrez and Arturas Grauslys, 2017. 
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
# https://creativecommons.org/licenses/by-nc-sa/4.0/

# R script for data normalisation and scaling
# Procedures are performed by calling the do_norm_scale() function with parameters:
#
# - data: a data frame with samples in the rows and variables (bins) in the columns.
#   The first column in the data (not counting sample names) should be a grouping variable
#
# - normalisation: a normalisation method to be used. 
#     Available methods: "PQN" - Probabilistic quotient normalisation
#                       "TotArea" - normalisation by the total area under the curve
#                       "Bin" - normalisation by 1 bin in the dataset
#
# - bin: the number of column in the dataset by which normalisation has to be performed (only used for "Bin" normalisation)
#
# - scaling: a scaling method to be used.
#     Available methods: "Auto" - mean centering and scaling by the standard deviation
#                       "Pareto" - mean centering and scaling by the square root of the standard deviation
#                       "Range" - mean centering and scaling such that the data ranges from 0 to 1.
#                       "Mean" - mean centering

NMRMetab_norm_scale = function(data, normalisation = 'None', bin = NA, scaling = 'None', writeToFile = F){
  
  #separate data from groups
  data_ = as.matrix(data[,3:ncol(data)])
  grp = as.factor(data[,2])
  
  # apply normalisation
  if (normalisation == 'None'){
    cat('Normalisation: None\n')
  } else if (normalisation == 'PQN'){
    data_ = PQN(data_)
    cat('Normalisation: PQN\n')
  } else if (normalisation == 'TotArea'){
    data_ = TotArea(data_)
    cat('Normalisation: Total area\n')
  } else if (normalisation == 'Bin'){
    if (!is.na(bin)){
      data_ = NormByBin(data_, bin)
    cat('Normalisation: Bin\n')
    cat(sprintf('Selected bin: %s \n', as.character(bin)))
    } else {
      print(paste('Method does not exist: ', normalisation, sep=''))
    }
  }
  
  #apply scaling
  if (scaling == 'None'){
    cat('Scaling: None\n')
  } else if (scaling == 'Auto'){
    data_ = apply(data_, 2, AutoScale)
    cat('Scaling: Auto\n')
  } else if (scaling == 'Pareto'){
    data_ = apply(data_, 2, ParetoScale)
    cat('Scaling: Pareto\n')
  } else if (scaling == 'Range') {
    data_ = apply(data_, 2, RangeScale)
    cat('Scaling: Range\n')
  } else if (scaling == 'Mean') {
    data_ = scale(data_, center=T, scale=F)
    cat('Scaling: Mean\n')
  } else {
    cat(sprint('Method does not exist: %s \n', scaling))
  }
 
  dataLabs = data[,1:2]
  out_data = cbind(dataLabs, as.data.frame(data_))
  
  if(writeToFile){
    # make output folder and write the data to file
    outDir = makeTSFolder('PLSDA')
    outPath = file.path(outDir, "NormScale_data.csv")
    write.csv(out_data, outPath,row.names = F)
    print(paste('Data written to ', outPath, sep=''))
  }
  
  return(out_data)
}

# scaling funcions
AutoScale<-function(x){
  (x - mean(x))/sd(x, na.rm=T)
}

ParetoScale<-function(x){
  (x - mean(x))/sqrt(sd(x, na.rm=T))
}

RangeScale<-function(x){
  if(max(x) == min(x)){
    x
  }else{
    (x - mean(x))/(max(x)-min(x))
  }
}

# normalisation

PQN <- function(data, loc = "median"){
  
  if (loc == "mean") {
    locFunc <- mean
  } else if (loc == 'median'){
    locFunc <- median
  } else {
    cat(sprintf("non such location metric %d", loc))
  }

  #if(ncol(data)>nrow(data)) data <- t(data)
  #data = abs(data)
  data_ = t(data)
  reference <- apply(data_,1,locFunc)
  # sometimes reference produces 0s so we turn them into 1s before division
  # so spectrum stays unchanged
  reference[reference==0] <- 1

  quotient <- data_/reference
  quotient.withLocFunc <- apply(quotient,2,locFunc)

  pqn.data <- t(data_)/quotient.withLocFunc
  pqn.data
}

# normalisation to total integral
TotArea <- function(data) {
  data_ = t(data)
  meanInt = sum(apply(data_,1,mean))
  scalingFactor = apply(data_,2,sum) / meanInt
  data_ = t(t(data_) / scalingFactor)
  t(data_)
}

# normalisation to reference peak
NormByBin <- function(data, bin) {
  data_ = t(data)
  refPeakInt = data[,bin]
  # adjust the integral for mean peak to preserve scale of spectra in the dataset
  refPeaksAdj = refPeakInt/mean(refPeakInt)
  data_ = t(t(data_)/refPeaksAdj)
  t(data_)
}

makeTSFolder = function(prefix){
  ts = format(Sys.Date(), "%b_%d_%Y")
  ts<-gsub(":","-",ts)
  tsDir = paste(prefix, ts, sep='_')
  if(!file.exists(file.path(getwd(),tsDir))) dir.create(file.path(getwd(),tsDir))
  return(file.path(getwd(),tsDir))
}