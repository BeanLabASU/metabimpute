#' iccEval
#'
#' a function to compare ICC of original data to imputed data. MAKE SURE THAT THE LAST IN METHODS VECTOR is the comparison imputation method (zero usually)
#'
#' @param origData
#' @param reps number of replicates
#' @param imputed this is a list of imputed dataframes
#' @param methods this is a vector of all the imputation methods used

#'
#' @return a list of various measures
#' @export

iccEval<-function(origData, reps, imputed, methods){

  require(Rmisc)
  require(irr)
  require(ggplot2)
  require(ICC)

  results<-list()

  groups<-as.factor(c(rep(1:(nrow(data)/reps), times=1, each=reps)))

  vector<-vector(mode="numeric")
  iccDF<-data.frame(vector)


  #calculating the ICC of the imputed matrices
  for(j in 1:length(methods)){

    for(i in 1:ncol(imputed[[j]])){
      #temp<-as.matrix(imputed[[j]][groups==i,])
      #tempT<-t(temp)

      iccDF[i,j]<-ICCest(x=groups,y=imputed[[j]][,i])[[1]]


    }

    colnames(iccDF)[j]<-methods[j]
  }

  #calculating the ICC of orig data
  for (i in 1:ncol(origData)){
    iccDF[i,(length(methods)+1)]<-(ICCest(x=groups,y=origData[,i]))
    print(i)
  }

  colnames(iccDF)[(length(methods)+1)]<-"original data"

  #replacing NaNs with 0, since the reason for NaNs will be variances of 0, or imputation of the whole vector to 0
  iccDF[is.na(iccDF)]<-0


  results[[1]]<-iccDF
  names(results)[1]<-"ICC dataframe"

  #creating dataframe of absolute value measures
  df1<-iccDF[iccDF[,length(methods)]>=0,1:(length(methods)-1)]
  Abs_Measure<-iccDF[iccDF[,length(methods)]>=0,1:(length(methods))]
  for(i in 1:(length(methods)-1)){
    df1[,i]<-abs(Abs_Measure[,length(methods)]-Abs_Measure[,i])
  }
  iccMeasure1<-apply(df1, 2, function(x) CI(x))

  results[[2]]<-iccMeasure1
  names(results)[2]<-"Absolute Value measure"


  #creating a data frame of difference measure
  df2<-iccDF[iccDF[,length(methods)]>=0,1:(length(methods)-1)]
  diff_Measure<-iccDF[iccDF[,length(methods)]>=0,1:(length(methods))]
  for(i in 1:(length(methods)-1)){
    df2[,i]<-diff_Measure[,i]-diff_Measure[,length(methods)]
  }

  iccMeasure2<-apply(df2,2, function(x) CI(x))

  results[[3]]<-iccMeasure2
  names(results)[3]<- "Difference measure"



  df3<-iccDF[iccDF[,length(methods)]>=0,1:(length(methods)-1)]
  SOS_Measure<-iccDF[iccDF[,length(methods)]>=0,1:(length(methods))]
  for(i in 1:(length(methods)-1)){
    df3[,i]<-(SOS_Measure[,length(methods)]-SOS_Measure[,i])^2
  }
  iccMeasure3<-apply(df3, 2, function(x) CI(x))

  results[[4]]<-iccMeasure3
  names(results)[4]<-"Sum of Squares measure"



  return(results)

}

#' ICC Change Plot
#'
#'@param iccMeasure matrix of a measure of ICC change from ICC Eval
#'@param methods list of methods of imputation with baseline comparison matrix last
#'
#'@return a ggplot of measure of icc change with CIs
#'@export

ICC_Change_Plot<-function(iccMeasure, methods, title){
  dfPlot<-data.frame(x =methods[-c(length(methods))],
                      F =iccMeasure[2,],
                      L =iccMeasure[3,],
                      U =iccMeasure[1,])

  plot<- ggplot(dfPlot, aes(x = x, y = F)) + geom_point(size = 4) + geom_errorbar(aes(ymax = U, ymin = L))+ ggtitle(title)
  return(plot)

}

#'ICC Scatter Plot
#'
#'@param rawData original dataset
#'@param reps number of replicate groups
#'@param iccImputed ICCs of the data after an imputation method
#'@param iccComparison ICCs of the comparison data (eg zero imputed matrix)
#'@param plotTitle string
#'
#'@return result list of the scatter plot data and the ggplot
#'@export

ICC_Scatter_Plot<-function(data, reps, iccImputed, iccComparison, plotTitle){
  result<-list()
  rep_groups <- c(rep(1:nrow(data)/reps, times=1, each=reps))

  #filtering data where ICC cannot be calculated due to missingness
  rawData_Filtered<-ICCformatting(data, reps=reps)

  #proportion of replicate groups permuted to zero per feature
  Proportion<-rawData_Filtered[1,]
  Proportion[1,]<-NA
  Proportion<-t(Proportion)

  for(i in 1:ncol(rawData_Filtered)){
    count<-0
    for(j in 1:length(unique(rep_groups))){
      if(sum(is.na(rawData_Filtered[rep_groups==j,i]))>=0.5*reps){
        count<-count+1
      }
    }
    Proportion[i,]<-count/(length(unique(rep_groups)))

  }

  ScatterplotData<-cbind(Proportion,iccImputed-iccComparison, iccComparison)
  #removing ICC estimates less than 0 from the comparison matrix
  ScatterplotData<-ScatterplotData[ScatterplotData[,3]>=0,]

  result[[1]]<-ScatterplotData

  plot <- ggplot(as.data.frame(ScatterplotData), aes(x=ScatterplotData[,1], y=ScatterplotData[,2])) +
    geom_point() + xlab("Proportion of Groups Permuted to Zero")+ylab("ICC Change")+ggtitle(plotTitle)+
    geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE)+theme_classic()

  result[[2]]<-plot

  return(result)

}


#' ICC Counts
#'
#' Function to count ICCs in a certain classification by Koo et al. Excellent is >.9, good is from .75 to .9, moderate from .5 to .75, poor from 0 to 0.5 and inconclusive
#' below 0
#'
#' @param iccMatrix matrix of ICCs for all imputation methods
#'
#' @return
#' @export

ICC_Counts<-function(iccMatrix){

  icc_Counts<-iccMatrix[1,]
  icc_Counts[1,]<-0
  icc_Counts[2,]<-0
  icc_Counts[3,]<-0
  icc_Counts[4,]<-0
  icc_Counts[5,]<-0
  rownames(icc_Counts)<-c("Excellent", "Good", "Moderate", "Poor", "Inconclusive")

  for(j in 1:ncol(icc_Counts)){
    for(i in 1:nrow(iccMatrix)){
      if(iccMatrix[i,j]>.9){
        icc_Counts[1,j]<-icc_Counts[1,j]+1
      }else if (iccMatrix[i,j]>.75){
        icc_Counts[2,j]<-icc_Counts[2,j]+1
      }else if (iccMatrix[i,j]>.5){
        icc_Counts[3,j]<-icc_Counts[3,j]+1
      }else if(iccMatrix[i,j]>=0){
        icc_Counts[4,j]<-icc_Counts[4,j]+1
      }else{
        icc_Counts[5,j]<-icc_Counts[5,j]+1
      }
    }
  }

  return(icc_Counts)

}

#'ICC_Change_Counts
#'
#'Counts of ICC changes from one class of Koo et al. to another after imputation
#'
#'@param iccImputed imputed ICCs vector
#'@param iccComparison vector of ICCs of the comparison method (eg zero imputation)
#'
#'@return iccChange a matrix displaying the icc changes from iccImputed to iccComparison
#'@export

ICC_Change_Counts<-function (iccImputed, iccComparison){

  names<-  c("Excellent to Good", "Excellent to Moderate", "Excellent to Poor", "Excellent to Inconclusive",
                              "Good to Excellent", "Good to Moderate", "Good to Poor", "Good to Inconclusive",
                              "Moderate to Excellent", "Moderate to Good", "Moderate to Poor", "Moderate to Inconclusive",
                              "Poor to Excellent", "Poor to Good", "Poor to Moderate", "Poor to Inconclusive",
                              "Inconclusive to Excellent", "Inconclusive to Good", "Inconclusive to Moderate", "Inconclusive to Poor")

  iccChange<-matrix(data=0.0,nrow=20,ncol=1)
  rownames(iccChange)<-names


    for (i in 1:length(iccImputed)){
      #Ex
      if (iccComparison[i]>.9){
        if (iccImputed[i]>.9){

        }else if(iccImputed[i]>.75){
          iccChange[1,1]<-iccChange[1,1]+1
        }else if(iccImputed[i]>.5){
          iccChange[2,1]<-iccChange[2,1]+1
        }else if (iccImputed[i]>=0){
          iccChange[3,1]<-iccChange[3,1]+1
        }else{
          iccChange[4,1]<-iccChange[4,1]+1
        }
      }
      #Good
      else if(iccComparison[i]>.75){
        if(iccImputed[i]>0.9){
          iccChange[5,1]<-iccChange[5,1]+1
        }else if(iccImputed[i]>.75){

        }else if(iccImputed[i]>.5){
          iccChange[6,1]<-iccChange[6,1]+1
        }else if(iccImputed[i]>=0){
          iccChange[7,1]<-iccChange[7,1]+1
        }else{
          iccChange[8,1]<-iccChange[8,1]+1
        }
      }
      #Mod
      else if (iccComparison[i]>.5){
        if(iccImputed[i]>.9){
          iccChange[9,1]<-iccChange[9,1]+1
        }else if (iccImputed[i]>.75){
          iccChange[10,1]<-iccChange[10,1]+1
        }else if (iccImputed[i]>.5){

        }else if (iccImputed[i]>=0){
          iccChange[11,1]<-iccChange[11,1]+1
        }else{
          iccChange[12,1]<-iccChange[12,1]+1
        }
      }
      #Poor
      else if (iccComparison[i]>=0){
        if(iccImputed[i]>.9){
          iccChange[13,1]<-iccChange[13,1]+1
        }else if(iccImputed[i]>.75){
          iccChange[14,1]<-iccChange[14,1]+1
        }else if (iccImputed[i]>.5){
          iccChange[15,1]<-iccChange[15,1]+1
        }else if(iccImputed[i]>=0){

        }else{
          iccChange[16,1]<-iccChange[16,1]+1
        }
      }
      # Inconclusive
      else if(iccComparison[i]<0){
        if(iccImputed[i]>.9){
          iccChange[17,1]<-iccChange[17,1]+1
        }else if(iccImputed[i]>.75){
          iccChange[18,1]<-iccChange[18,1]+1
        }else if(iccImputed[i]>.5){
          iccChange[19,1]<-iccChange[19,1]+1
        }else if(iccImputed[i]>=0){
          iccChange[20,1]<-iccChange[20,1]+1
        }
      }


  }

  return(iccChange)

}

#' ICCformatting
#'
#' Function that formats a dataframe or matrix into a dataframe where there are at least (replicate # - 1) factor levels that have non
#' missing data.
#'
#' @param data
#' @param reps replicate number
#' @return matrix
#' @export

ICCformatting<-function(data, reps){
  keep<-vector()
  rep_groups <- c(rep(1:nrow(data)/reps, times=1, each=reps))

  for (j in 1:ncol(data)){
    count<-0
    for(i in 1:length(unique(rep_groups))){

      if (sum(is.na(data[rep_groups==i,j]))<reps){
        count<-count+1
      }

    }
    if(count>1){
      keep[j]<-TRUE
    }else{
      keep[j]<-FALSE
    }

  }

  return(data[,keep==TRUE])

}
