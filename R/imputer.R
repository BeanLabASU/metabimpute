#' Title impute
#' This functions contains different imputation methods and imputes the data with all
#' the different imputation methods
#' @param data data matrix with simulated data
#' @param methods vector containing the method names. Note GSimp_Real is for real non-neg data. GSimp_Sim is for already
#' @param local a boolean to determine if local rep_impute method is to be used, default to true.
#' @param reps the number of replicate groups
#' simulated data that doesn't require log preprocessing.
#'
#' @return results_data the matrix containing the imputed data
#' @export
#'
#' @examples
#' imputed_data <- impute(data=miss_data, methods=imputation_methods)
#' #'####################################################

impute <- function(data, methods, local=TRUE, reps) {
  require(doParallel)
  require(missForest)
  require(missRanger)
  require(pcaMethods)
  require(impute)
  require(PEMM)
  require(imputeLCMD)
  require(magrittr)
  require(matrixStats)

  #source('MetabImpute/R/MVI_global.R')
  #source('MetabImpute/R/GSimp.R')

  if (length(methods) != 1 & length(methods) != ncol(data)) {
    stop("Methods needs to be either one value or of the same length as number of columns in data.")
  }
  imputed_data <- matrix(NA,nrow = nrow(data),ncol = ncol(data))
  results_data <- data

  if (length(methods) == 1) {
    methods <- rep(methods, times=ncol(data))
  }

  if ("RF_P" %in% methods) {
    cl<-makeCluster(4)
    registerDoParallel(cl)
    imputed_data <- missForest::missForest(xmis = data,maxiter = 10,verbose = FALSE, parallelize = 'variables')$ximp
    index <- which(methods == "RF")
    results_data <- imputed_data
    stopCluster(cl=cl)

  }

  if ("RF" %in% methods) {

    #imputed_data <- missForest::missForest(xmis = data,maxiter = 10,verbose = FALSE, parallelize = 'no')$ximp
    imputed_data<-missRanger(data=as.data.frame(data), num.trees=100)
    index <- which(methods == "RF")
    results_data <- imputed_data


  }

  if ("PPCA" %in% methods) {
    # Do cross validation with ppca for component 2:10
    #esti <- kEstimate(Matrix= data, method = "ppca", evalPcs = 2:10, nruncv = 1 , em="nrmsep")
    esti<-kEstimateFast(Matrix=data, method="ppca", evalPcs=1:10, em='q2')
    # The best result was obtained for this number of PCs:esti$bestNPcs
    pc <- pcaMethods::pca(object = data,nPcs=esti$bestNPcs, method="ppca")
    #index <- which(methods == "PPCA")
    imputed_data <- completeObs(pc)

    index <- which(methods == "PPCA")
    results_data[,index] <- imputed_data[,index]

  }


  if ("BPCA" %in% methods){
    # bayesian principal component analysis
    pc <-pcaMethods:: pca(object = data, method="bpca", nPcs=10)

    ## Get the estimated complete observations
    imputed_data <- completeObs(pc)

    index <- which(methods == "BPCA")
    results_data[,index] <- imputed_data[,index]

  }

  if ("QRILC" %in% methods){


    qrilc<- QRILC_Prime(dataSet.mvs = data)
    imputed_data<-qrilc

    index <- which(methods=="QRILC")
    results_data[,index]<- imputed_data[,index]
  }

  if ("GSimp_Log" %in% methods){

    data_raw <- as.matrix(data)
    ## log transformation ##
    data_raw_log <- data_raw %>% log()
    data_raw_log[data_raw_log==-Inf]<-0
    ## Initialization ##
    data_raw_log_qrilc <- QRILC_Prime(data_raw_log)
    ## Centralization and scaling ##
    data_raw_log_qrilc_sc <- scale_recover(data_raw_log_qrilc, method = 'scale')
    ## Data after centralization and scaling ##
    data_raw_log_qrilc_sc_df <- data_raw_log_qrilc_sc[[1]]
    ## Parameters for centralization and scaling ##
    ## For scaling recovery ##
    data_raw_log_qrilc_sc_df_param <- data_raw_log_qrilc_sc[[2]]
    ## NA position ##
    NA_pos <- which(is.na(data_raw), arr.ind = T)
    ## bala bala bala ##
    data_raw_log_sc <- data_raw_log_qrilc_sc_df
    data_raw_log_sc[NA_pos] <- NA
    ## GSimp imputation with initialized data and missing data ##
    result <- data_raw_log_sc %>% GS_impute(., iters_each=30, iters_all=5,
                                            initial = as.data.frame(data_raw_log_qrilc_sc_df),
                                            lo=-Inf, hi= 'min', n_cores=4,
                                            imp_model='glmnet_pred')
    data_imp_log_sc <- result$data_imp
    ## Data recovery ##
    data_imp <- data_imp_log_sc %>%
      scale_recover(., method = 'recover',
                    param_df = data_raw_log_qrilc_sc_df_param) %>%
      extract2(1) %>% exp()

    #imputed_data<-pre_processing_GS_wrapper(data)
    #index<- which(methods=="GSimp_Real")
    results_data<- data_imp



  }

  if ("GSimp_NoLog" %in% methods){
    data_raw <- data

    data_raw_qrilc <- as.data.frame(QRILC_Prime(data_raw))

    result <- data_raw%>% GS_impute(., iters_each=50, iters_all=10,
                                    initial = data_raw_qrilc,
                                    lo=-Inf, hi= 'min', n_cores=4,
                                    imp_model='glmnet_pred')
    imputed_data <- result$data_imp
    index<-which(methods=="GSimp_Sim")
    results_data<-imputed_data
  }

  if("Rep_Imp_HM" %in% methods){
    rep_groups <- c(rep(1:(nrow(data)/reps), times=1, each=reps))
    data[data==0]<-NA
    data<- cbind(data,rep_groups)

      newData <-data.frame(matrix(nrow=nrow(data),ncol=ncol(data)-1))
      for (i in 1:length(unique(data$rep_groups))){
        tempData<-data[data[,ncol(data)]==i,]
        for (j in 1:(ncol(data)-1)){
          halfmin<-NA
          if(local==T){
            halfmin<- min(tempData[,j],na.rm = TRUE)/2
          } else{
            halfmin<-min(data[,j],na.rm=TRUE)/2
          }

          if (sum(is.na(tempData[,j]))>=(0.5*reps) ){
            tempData[,j]<-0
          }else {
            tempData[is.na(tempData[,j]),j]<- halfmin
          }
          newData[((reps*i)-(reps-1)):(reps*i),]<-tempData[,1:(ncol(tempData)-1)]
        }
      }
      rownames(newData)<-rownames(data)
      colnames(newData)<-colnames(data)[1:(ncol(data)-1)]

      results_data<newData



  }

  if("Rep_Imp_mean" %in% methods){
    rep_groups <- c(rep(1:(nrow(data)/reps), times=1, each=reps))
    data[data==0]<-NA
    data<- cbind(data,rep_groups)

    newData <- data.frame(matrix(nrow=nrow(data),ncol=ncol(data)-1))
    for (i in 1:length(unique(data$rep_groups))){
      tempData<-data[data[,ncol(data)]==i,]
      for (j in 1:(ncol(data)-1)){
        mean<-NA
        if(local==T){
          mean<- mean(tempData[,j],na.rm = TRUE)
        } else{
          mean<-mean(data[,j],na.rm=TRUE)
        }

        if (sum(is.na(tempData[,j]))>=(0.5*reps)){
          tempData[,j]<-0
        }else {
          tempData[is.na(tempData[,j]),j]<- mean
        }
        newData[((reps*i)-(reps-1)):(reps*i),]<-tempData[,1:(ncol(tempData)-1)]
      }
    }
    rownames(newData)<-rownames(data)
    colnames(newData)<-colnames(data)[1:(ncol(data)-1)]

    results_data<newData


  }

  if("Rep_Imp_median" %in% methods){
    rep_groups <- c(rep(1:(nrow(data)/reps), times=1, each=reps))
    data[data==0]<-NA
    data<- cbind(data,rep_groups)

    newData <- data.frame(matrix(nrow=nrow(data),ncol=ncol(data)-1))
    for (i in 1:length(unique(data$rep_groups))){
      tempData<-data[data[,ncol(data)]==i,]
      for (j in 1:(ncol(data)-1)){
        median<-NA
        if(local==T){
          median<- median(tempData[,j],na.rm = TRUE)
        } else{
          median<-median(data[,j],na.rm=TRUE)
        }

        if (sum(is.na(tempData[,j]))>=(0.5*reps)){
          tempData[,j]<-0
        }else {
          tempData[is.na(tempData[,j]),j]<- median
        }
        newData[((reps*i)-(reps-1)):(reps*i),]<-tempData[,1:(ncol(tempData)-1)]
      }
    }
    rownames(newData)<-rownames(data)
    colnames(newData)<-colnames(data)[1:(ncol(data)-1)]

    results_data<newData



  }

  if("Rep_Imp_min" %in% methods){
    rep_groups <- c(rep(1:(nrow(data)/reps), times=1, each=reps))
    data[data==0]<-NA
    data<- cbind(data,rep_groups)

    newData <- data.frame(matrix(nrow=nrow(data),ncol=ncol(data)-1))
    for (i in 1:length(unique(data$rep_groups))){
      tempData<-data[data[,ncol(data)]==i,]
      for (j in 1:(ncol(data)-1)){
        min<-NA
        if(local==T){
          min<- min(tempData[,j],na.rm = TRUE)
        } else{
          min<-min(data[,j],na.rm=TRUE)
        }

        if (sum(is.na(tempData[,j]))>=(0.5*reps)){
          tempData[,j]<-0
        }else {
          tempData[is.na(tempData[,j]),j]<- min
        }
        newData[((reps*i)-(reps-1)):(reps*i),]<-tempData[,1:(ncol(tempData)-1)]
      }
    }
    rownames(newData)<-rownames(data)
    colnames(newData)<-colnames(data)[1:(ncol(data)-1)]

    results_data<newData





  }

  if("Rep_Zero" %in% methods){
    rep_groups <- c(rep(1:(nrow(data)/reps), times=1, each=reps))
    data[data==0]<-NA
    data<- cbind(data,rep_groups)

    newData <- data.frame(matrix(nrow=nrow(data),ncol=ncol(data)-1))
    for (i in 1:length(unique(data$rep_groups))){
      tempData<-data[data[,ncol(data)]==i,]
      for (j in 1:(ncol(data)-1)){

        if (sum(is.na(tempData[,j]))>=(0.5*reps)){
          tempData[,j]<-0
        }else {
          tempData[is.na(tempData[,j]),j]<-0
        }
        newData[((reps*i)-(reps-1)):(reps*i),]<-tempData[,1:(ncol(tempData)-1)]
      }
    }
    rownames(newData)<-rownames(data)
    colnames(newData)<-colnames(data)[1:(ncol(data)-1)]

    results_data<newData
  }

  if("Rep_Imp_RF" %in% methods){
    rep_groups <- c(rep(1:(nrow(data)/reps), times=1, each=reps))
    data[data==0]<-NA
    data<- cbind(data,rep_groups)

    newData <-data.frame(matrix(nrow=nrow(data),ncol=ncol(data)-1))

    for (i in 1:length(unique(rep_groups))){
      tempData<-data[data[,ncol(data)]==i,]
      for (j in 1:(ncol(data)-1)){
        if (sum(is.na(tempData[,j]))>=(0.5*reps)){
          tempData[,j]<-0
        }
        newData[((reps*i)-(reps-1)):(reps*i),]<-tempData[,1:(ncol(tempData)-1)]
      }
    }

    results_data<-missRanger(data=as.data.frame(newData), num.trees=100)



  }

  if("Rep_Imp_GSimp" %in% methods){
    rep_groups <- c(rep(1:(nrow(data)/reps), times=1, each=reps))
    data[data==0]<-NA
    data<- cbind(data,rep_groups)

    newData <-data.frame(matrix(nrow=nrow(data),ncol=ncol(data)-1))

    for (i in 1:length(unique(rep_groups))){
      tempData<-data[data[,ncol(data)]==i,]
      for (j in 1:(ncol(data)-1)){
        if (sum(is.na(tempData[,j]))>=(0.5*reps)){
          tempData[,j]<-0
        }
        newData[((reps*i)-(reps-1)):(reps*i),]<-tempData[,1:(ncol(tempData)-1)]
      }
    }

    data_raw <- as.matrix(newData)
    ## log transformation ##
    data_raw_log <- data_raw %>% log()
    data_raw_log[data_raw_log==-Inf]<-0
    ## Initialization ##
    data_raw_log_qrilc <- QRILC_Prime(data_raw_log)
    ## Centralization and scaling ##
    data_raw_log_qrilc_sc <- scale_recover(data_raw_log_qrilc, method = 'scale')
    ## Data after centralization and scaling ##
    data_raw_log_qrilc_sc_df <- data_raw_log_qrilc_sc[[1]]
    ## Parameters for centralization and scaling ##
    ## For scaling recovery ##
    data_raw_log_qrilc_sc_df_param <- data_raw_log_qrilc_sc[[2]]
    ## NA position ##
    NA_pos <- which(is.na(data_raw), arr.ind = T)
    ## bala bala bala ##
    data_raw_log_sc <- data_raw_log_qrilc_sc_df
    data_raw_log_sc[NA_pos] <- NA
    ## GSimp imputation with initialized data and missing data ##
    result <- data_raw_log_sc %>% GS_impute(., iters_each=30, iters_all=5,
                                            initial = as.data.frame(data_raw_log_qrilc_sc_df),
                                            lo=-Inf, hi= 'min', n_cores=4,
                                            imp_model='glmnet_pred')
    data_imp_log_sc <- result$data_imp
    ## Data recovery ##
    data_imp <- data_imp_log_sc %>%
      scale_recover(., method = 'recover',
                    param_df = data_raw_log_qrilc_sc_df_param) %>%
      extract2(1) %>% exp()

    #imputed_data<-pre_processing_GS_wrapper(data)
    #index<- which(methods=="GSimp_Real")
    data_imp[data_imp==0.001]<-0
    results_data<- data_imp



  }


  if("Rep_Imp_QRILC" %in% methods){
    rep_groups <- c(rep(1:(nrow(data)/reps), times=1, each=reps))
    data[data==0]<-NA
    data<- cbind(data,rep_groups)

    newData <- data.frame(matrix(nrow=nrow(data),ncol=ncol(data)-1))

    for (i in 1:length(unique(rep_groups))){
      tempData<-data[data[,ncol(data)]==i,]
      for (j in 1:(ncol(data)-1)){
        if (sum(is.na(tempData[,j]))>=(0.5*reps)){
          tempData[,j]<-0.001
        }
        newData[((reps*i)-(reps-1)):(reps*i),]<-tempData[,1:(ncol(tempData)-1)]
      }
    }

    results_data<-QRILC_Prime(dataSet.mvs = data)
    results_data[results_data<=0.001]<-0



  }

  if("Rep_Imp_BPCA" %in% methods){
    rep_groups <- c(rep(1:(nrow(data)/reps), times=1, each=reps))
    data[data==0]<-NA
    data<- cbind(data,rep_groups)

    newData <- data.frame(matrix(nrow=nrow(data),ncol=ncol(data)-1))

    for (i in 1:length(unique(rep_groups))){
      tempData<-data[data[,ncol(data)]==i,]
      for (j in 1:(ncol(data)-1)){
        if (sum(is.na(tempData[,j]))>=(0.5*reps)){
          tempData[,j]<-0
        }
        newData[((reps*i)-(reps-1)):(reps*i),]<-tempData[,1:(ncol(tempData)-1)]
      }
    }

    pc <-pcaMethods:: pca(object = newData, method="bpca", nPcs=10)


    results_data<-completeObs(pc)



  }


  if ("min" %in% methods || "MIN" %in% methods) {
    imputed_data <- data
    foreach (data_column=which(methods == "min")) %do% {
      method <- methods[data_column]
      impu_value <- min(data[,data_column], na.rm=TRUE)
      results_data[is.na(imputed_data[,data_column]), data_column] <- impu_value
    }

  }





  if (methods=="halfmin") {
    imputed_data <- data

    if(local==T){
      rep_groups <- c(rep(1:(nrow(data)/reps), times=1, each=reps))
      data<- cbind(data,rep_groups)

      newData <-data.frame(matrix(nrow=nrow(data),ncol=ncol(data)-1))
      for (i in 1:length(unique(data$rep_groups))){
        tempData<-data[data[,ncol(data)]==i,]
        for (j in 1:(ncol(data)-1)){
          halfmin<- min(tempData[,j],na.rm = TRUE)/2
          tempData[is.na(tempData[,j]),j]<- halfmin

          newData[((reps*i)-(reps-1)):(reps*i),]<-tempData[,1:(ncol(tempData)-1)]
        }
      }
      rownames(newData)<-rownames(data)
      colnames(newData)<-colnames(data)[1:(ncol(data)-1)]

      results_data<newData


    }else{
    foreach (data_column=which(methods == "halfmin")) %do% {
      method <- methods[data_column]
      impu_value <- (min(data[,data_column], na.rm=TRUE)/2)
      results_data[is.na(imputed_data[,data_column]), data_column] <- impu_value
    }

    }
    return(results_data)
  }
  if ("mean" %in% methods || "MEAN" %in% methods){
    imputed_data <- data
    foreach (data_column=which(methods == "mean")) %do% {
      method <- methods[data_column]
      impu_value <- round(mean(data[,data_column], na.rm=TRUE),digits = 2)
      results_data[is.na(imputed_data[,data_column]), data_column] <- impu_value
    }

  }

  if ("median" %in% methods || "MEDIAN" %in% methods){
    imputed_data <- data
    foreach (data_column=which(methods == "median")) %do% {
      method <- methods[data_column]
      impu_value <- round(median(data[,data_column], na.rm=TRUE),digits = 2)
      results_data[is.na(imputed_data[,data_column]), data_column] <- impu_value
    }

  }

  if ("zero" %in% methods || "ZERO" %in% methods){
    imputed_data <- data
    foreach (data_column=which(methods == "zero")) %do% {
      method <- methods[data_column]
      impu_value <- 0
      results_data[is.na(imputed_data[,data_column]), data_column] <- impu_value
    }

  }
  if ("EX" %in% methods){


    index <- which(methods == "EX")
    results_data[,index] <-  data[,index]
  }

  if ("NONE" %in% methods){


    index <- which(methods == "NONE")
    results_data[,index] <-  data[,index]

  }


  return(results_data)
}





#' GSimp function
#'
#' This does the imputation from WandeRum/GSimp. It does rely on MVI_Global.R and Prediction_funcs.R
#' @param data
#' @return imputed dataframe.
pre_processing_GS_wrapper <- function(data) {
  data_raw <- data
  ## log transformation ##
  data_raw_log <- data_raw %>% log()

  data_raw_log=do.call(data.frame,lapply(data_raw_log, function(x) replace(x, is.infinite(x),NA)))

  ## Initialization ##
  data_raw_log_qrilc <- QRILC_Prime(data_raw_log)
  ## Centralization and scaling ##
  data_raw_log_qrilc_sc <- scale_recover(as.data.frame(data_raw_log_qrilc), method = 'scale')
  ## Data after centralization and scaling ##
  data_raw_log_qrilc_sc_df <- data_raw_log_qrilc_sc[[1]]
  ## Parameters for centralization and scaling ##
  ## For scaling recovery ##
  data_raw_log_qrilc_sc_df_param <- data_raw_log_qrilc_sc[[2]]
  ## NA position ##
  NA_pos <- which(is.na(data_raw), arr.ind = T)
  ## bala bala bala ##
  data_raw_log_sc <- data_raw_log_qrilc_sc_df
  data_raw_log_sc[NA_pos] <- NA
  ## GSimp imputation with initialized data and missing data ##
  result <- data_raw_log_sc %>% GS_impute(., iters_each=50, iters_all=10,
                                          initial = data_raw_log_qrilc_sc_df,
                                          lo=-Inf, hi= 'min', n_cores=2,
                                          imp_model='glmnet_pred')
  data_imp_log_sc <- result$data_imp
  ## Data recovery ##
  data_imp <- data_imp_log_sc %>%
    scale_recover(., method = 'recover',
                  param_df = data_raw_log_qrilc_sc_df_param) %>%
    extract2(1) %>% exp()
  return(data_imp)
}


#' imputeMulti
#'
#' Uses the impute function to do multiple specified imputation methods and contains them within a list
#'
#' @param methods a list of imputation methods
#' @param data the dataset with missing values contained
#' @param reps number of replicates, if applicable
#' @return a list containing the dataframes of the imputed data
#' @export

imputeMulti<- function(methods, data, reps=NULL){
  results<-list()
  for (i in 1:length(methods)){
    results[[i]]<-impute(data, methods = methods[i], reps = reps)
    print(paste("imputed using", methods[i], sep=" "))
    names(results)[i]<-methods[i]

  }

  return(results)

}

#'QRILC_Prime
#'
#'Modified from source code written by Cosmin Lazar at https://www.rdocumentation.org/packages/imputeLCMD/versions/2.0/topics/impute.QRILC
#' ammended to handle cases in which there are columns without any missing values. Accessed from git
#'
#'@param dataSet.mvs the data to be imputed which can include non-missing rows which will be left untouched
#'@param tune.sigma
#'@return imputed data frame
#'@export

QRILC_Prime<-function (dataSet.mvs, tune.sigma = 1)
{
  nFeatures = dim(dataSet.mvs)[1]
  nSamples = dim(dataSet.mvs)[2]
  dataSet.imputed = dataSet.mvs
  QR.obj = list()
  for (i in 1:nSamples) {
    curr.sample = dataSet.mvs[, i]
    pNAs = length(which(is.na(curr.sample)))/length(curr.sample)

    if(pNAs!=0){
      upper.q = 0.95
      q.normal = qnorm(seq(pNAs, upper.q, (upper.q - pNAs)/(upper.q *
                                                            10000)), mean = 0, sd = 1)
      q.curr.sample = quantile(curr.sample, probs = seq(0,
                                                      upper.q, 1e-04), na.rm = T)
      q.curr.sample[q.curr.sample==-Inf]<-0
      temp.QR = lm(q.curr.sample ~ q.normal)
      QR.obj[[i]] = temp.QR
      mean.CDD = temp.QR$coefficients[1]
      sd.CDD = abs(as.numeric(temp.QR$coefficients[2]))
      data.to.imp = rtmvnorm(n = nFeatures, mean = mean.CDD,
                           sigma = sd.CDD * tune.sigma, upper = qnorm(pNAs,
                                                                      mean = mean.CDD, sd = sd.CDD), algorithm = c("gibbs"))
      curr.sample.imputed = curr.sample
      curr.sample.imputed[which(is.na(curr.sample))] = data.to.imp[which(is.na(curr.sample))]
      dataSet.imputed[, i] = curr.sample.imputed
    }

  }
  results = dataSet.imputed
  return(results)
}
