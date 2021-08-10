#' errorEvals
#'
#' this will be a function that utilizes all of the different error analysis tools into one nice function
#' @param origData
#' @param missData
#' @param method string "NRMSE" "NRMSE-SOR", "PCA-P"
#' @param imputationResults is specifically a list of imputed dataframes from all the chosen methods, used in SOR
#' @results a list of error measurements across imputations methods
#' @export

errorEvals<- function(origData, missData, method, imputationResults, simulate_Data=T){

  require(vegan)
  require(foreach)
  require(doParallel)
  require(reshape2)
  require(ggplot2)
  require(pls)

  score<-vector(mode="numeric", length=length(imputationResults))
  for (i in seq_along(imputationResults)){
    names(score)[i]<-names(imputationResults)[i]
  }

  for (i in seq_along(imputationResults)){
    imputationResults[[i]][is.na(imputationResults[[i]])]<-0.0
  }

  if (method=="NRMSE"){

    for (i in seq_along(imputationResults)){
      score[i]<-nrmse(ximp=imputationResults[[i]], xmis=missData, xtrue=origData)
    }


  }


  if (method=="NRMSE-SOR"){


    for (i in 1:ncol(origData)){
      temp_NRMSE<- vector(mode="numeric", length=length(imputationResults))
      rank_NRMSE<-temp_NRMSE
      for (k in seq_along(imputationResults)){
        temp_NRMSE[k]<-nrmse(ximp=imputationResults[[k]][,i], xmis=missData[,i],xtrue=origData[,i])

      }
      rank_NRMSE<-rank(temp_NRMSE, ties.method = "min")
      score<-score+rank_NRMSE

    }


  }

  if (method=="PCA-P"){

    if(simulate_Data==T){

    for (i in seq_along(imputationResults)){


      pca_Orig<-prcomp(origData, scale. = F, center = F)$x[,1:2]
      score[i]<-procrustes(pca_Orig, prcomp(imputationResults[[i]], scale. =F, center=F)$x[,1:2], symmetric = T)$ss
    }
    }else{
      for (i in seq_along(imputationResults)){


      pca_Orig<-prcomp(origData, scale. = T, center = T)$x[,1:2]
      score[i]<-procrustes(pca_Orig, prcomp(imputationResults[[i]], scale. =T, center=T)$x[,1:2], symmetric = T)$ss
    }
    }

  }


  return(score)

}


#' simEval
#'
#' @param origData the dataset prior to simulating missingness
#' @param missData the dataset with missingness applied
#' @param impData the imputed dataset
#' @param methods a list of evaluation methods to be utilized
#' @param outcome is a grouping list for supervized learning to use PLS or pearson
#' @param imputationResults
#'
#' @return results a dataframe displaying NRMSE, PCA-Procustes, PLS Procrustes and Student's T test and Pearson Correl
#' @export

simEval<- function(origData, missData, impData, methods, imputationResults, simulate_Data){

  results<-data.frame(row.names=names(imputationResults))

  for (i in seq_along(methods)){
    results[,i]<-errorEvals(origData = origData, missData = missData, method = methods[i], imputationResults = imputationResults, simulate_Data=simulate_Data)
    colnames(results)[i]<-methods[i]

}

return(results)

}


#' graphEval
#'
#' Function that takes in the error result list from rearrangeList which is from simulateEngine followed by aggregateDF and outputs multiplots for each missingness
#' ratio and makes graphs for each error type in the multiplot
#'
#' @param ggplotlist this is a list by missingness ratio and then percents
#' @return a list of plot collections by error evaluation type.
#' @export

graphEval<- function( ggplotlist){

require(ggplot2)
require(tidyr)
require(gridExtra)

  plotList<-list()
  for (i in 1:length(ggplotlist)){
    myPlots<-list()
    for (j in 1:length(ggplotlist[[i]])){
      df<-ggplotlist[[i]][[j]]
      s<- names(ggplotlist[[i]])[j]
      data<- df %>% gather(key="variable", value= "value", -Percent)
      plot<- ggplot(data, aes(x = Percent, y = value)) +
        geom_line(aes(color = variable))+ labs(y= s, x = "Missing Proportion")+
        geom_point(aes(color = variable))
      myPlots[[j]]<-plot
    }

    plotList[[i]]<-myPlots
    names(plotList)[i]<-names(ggplotlist)[i]


  }

  #grid.arrange(grobs=plotList[[1]], top=names(plotList)[1], ncol=2) Use this to display the plots!!!!

return(plotList)

}


#' iccEval
#'
#' a function to compare ICC of original data to imputed data. MAKE SURE THAT THE LAST IN METHODS VECTOR is zero imputation
#'
#' @param origData
#' @param groups this is a factor containing group levels
#' @param imputed this is a list of imputed dataframes
#' @param methods this is a vector of all the imputation methods used

#'
#' @return a list of various measures
#' @export

iccEval<-function(origData, groups, imputed, methods){

  require(Rmisc)
  require(irr)
  require(ggplot2)
  require(ICC)

  results<-list()


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
  names(results)[2]<-"Absolute Value measures"

  dfPlot1 <- data.frame(x =methods[-c(length(methods))],
                        F =iccMeasure1[2,],
                        L =iccMeasure1[3,],
                        U =iccMeasure1[1,])



  plot1<- ggplot(dfPlot1, aes(x = x, y = F)) +
    geom_point(size = 4) +
    geom_errorbar(aes(ymax = U, ymin = L))+ ggtitle("Absolute Value Measure")

  results[[3]]<- plot1
  names(results)[3]<- "Abs value measure chart"


  #creating a data frame of difference measure
  df2<-iccDF[iccDF[,length(methods)]>=0,1:(length(methods)-1)]
  diff_Measure<-iccDF[iccDF[,length(methods)]>=0,1:(length(methods))]
  for(i in 1:(length(methods)-1)){
    df2[,i]<-diff_Measure[,i]-diff_Measure[,length(methods)]
  }

  iccMeasure2<-apply(df2,2, function(x) CI(x))

  results[[4]]<-iccMeasure2
  names(results)[4]<- "Diff measure matrix"

  dfPlot2 <- data.frame(x =methods[-c(length(methods))],
                        F =iccMeasure2[2,],
                        L =iccMeasure2[3,],
                        U =iccMeasure2[1,])

  plot2<- ggplot(dfPlot2, aes(x = x, y = F)) +
    geom_point(size = 4) +
    geom_errorbar(aes(ymax = U, ymin = L))+ ggtitle("Difference Measure")

  results[[5]]<-plot2
  names(results)[5]<-"Diff measure plot"

  df3<-iccDF[iccDF[,length(methods)]>=0,1:(length(methods)-1)]
  SOS_Measure<-iccDF[iccDF[,length(methods)]>=0,1:(length(methods))]
  for(i in 1:(length(methods)-1)){
    df3[,i]<-(SOS_Measure[,length(methods)]-SOS_Measure[,i])^2
  }
  iccMeasure3<-apply(df3, 2, function(x) CI(x))

  results[[6]]<-iccMeasure3
  names(results)[6]<-"Absolute Value measures"

  dfPlot3 <- data.frame(x =methods[-c(length(methods))],
                        F =iccMeasure3[2,],
                        L =iccMeasure3[3,],
                        U =iccMeasure3[1,])



  plot3<- ggplot(dfPlot3, aes(x = x, y = F)) +
    geom_point(size = 4) +
    geom_errorbar(aes(ymax = U, ymin = L))+ ggtitle("SOS Value Measure")

  results[[7]]<- plot3
  names(results)[7]<- "SOS value measure chart"

  return(results)

}


#' ICCformatting
#'
#' Function that formats a dataframe or matrix into a dataframe where there are at least 2 factor levels that have non
#' missing data.
#'
#' @param data
#' @param groups
#' @return matrix
#' @export

ICCformatting<-function(data, groups){
  keep<-vector()
  for (j in 1:ncol(data)){
    count<-0
    for(i in 1:length(unique(groups))){

      if (sum(is.na(data[groups==i,j]))<3){
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
