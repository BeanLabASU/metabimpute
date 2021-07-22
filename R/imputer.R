#' Title impute
#' This functions contains different imputation methods and imputes the data with all
#' the different imputation methods
#' @param data data matrix with simulated data
#' @param method the imputation method Note within replicate methods are preceded by 'R', eg 'RMIN'
#' method=c('RF', 'BPCA', 'QRILC', 'GSIMP', 'RHM','RMEAN', 'RMEDIAN', 'RMIN','RZERO', 'RRF',
#' 'RGSIMP', 'RQRILC','RBPCA','min','halfmin', 'mean', 'median', 'zero')
#' @param local a boolean to determine if local rep_impute method is to be used, default to true.
#' @param reps the number of replicate groups
#' simulated data that doesn't require log preprocessing.
#'
#' @return results_data the matrix containing the imputed data
#' @export
#'
#' @examples
#' imputed_data <- impute(data=miss_data, methods=imputation_methods, local=T, reps=3)
#' #'####################################################
impute <- function(data, method, local=TRUE, reps) {
  require(doParallel)
  require(missForest)
  require(missRanger)
  require(pcaMethods)
  require(impute)
  require(PEMM)
  require(imputeLCMD)
  require(magrittr)
  require(matrixStats)
  require(foreach)
  require(MASS)
  require(abind)
  require(randomForest)
  require(glmnet)
  require(rpart)
  require(FNN)


  imputed_data <- matrix(NA,nrow = nrow(data),ncol = ncol(data))
  results_data <- data



  if (method == 'RF') {

    #imputed_data <- missForest::missForest(xmis = data,maxiter = 10,verbose = FALSE, parallelize = 'no')$ximp
    imputed_data<-missRanger(data=as.data.frame(data), num.trees=100)
    rownames(imputed_data)<-rownames(data)
    colnames(imputed_data)<-colnames(data)
    results_data <- imputed_data



  }

  if (method == "BPCA"){
    # bayesian principal component analysis
    pc <-pcaMethods:: pca(object = data, method="bpca", nPcs=min(c(10,ncol(data)/2)))

    ## Get the estimated complete observations
    imputed_data <- completeObs(pc)
    rownames(imputed_data)<-rownames(data)
    colnames(imputed_data)<-colnames(data)
    results_data <- imputed_data

  }

  if (method=='QRILC'){


    qrilc<- impute.QRILC(dataSet.mvs = data)
    imputed_data<-qrilc[[1]]
    rownames(imputed_data)<-rownames(data)
    colnames(imputed_data)<-colnames(data)
    results_data<- imputed_data
  }

  if (method=='GSIMP'){

    data_raw <- as.matrix(data)
    ## log transformation ##
    data_raw_log <- data_raw %>% log()
    data_raw_log[data_raw_log==-Inf]<-0
    ## Initialization ##
    data_raw_log_qrilc <- impute.QRILC(data_raw_log)[[1]]
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
    rownames(results_data)<-rownames(data)
    colnames(results_data)<-colnames(data)



  }

  if(method=='RHM'){
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

      results_data<-newData




  }

  if(method=="RMEAN"){
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

    results_data<-newData



  }

  if(method=="RMEDIAN"){
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

    results_data<-newData




  }

  if(method=="RMIN"){
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

    results_data<-newData






  }

  if(method=="RZERO"){
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

    results_data<-newData

  }

  if(method=="RRF"){
    rep_groups <- c(rep(1:(nrow(data)/reps), times=1, each=reps))
    rownames<-rownames(data)
    colnames<-colnames(data)
    newData <-data.frame(matrix(nrow=nrow(data),ncol=ncol(data)))
    data<- cbind(data,rep_groups)



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
    rownames(results_data)<-rownames
    colnames(results_data)<-colnames


  }

  if(method=="RGSIMP"){
    rep_groups <- c(rep(1:(nrow(data)/reps), times=1, each=reps))
    rownames<-rownames(data)
    colnames<-colnames(data)
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
    miss<-as.vector(lapply(newData, function(x) sd(x,na.rm=T)==0)==T)

    data_raw <- as.data.frame(newData)
    ## log transformation ##
    data_raw_log <- data_raw %>% log()

    #added to deal with fully present columns
    data_raw_log[data_raw_log==-Inf]<- -1
    data_keep<-data_raw_log[,miss==F]
    data_rm<-newData[,miss==T]
    ## Initialization ##

    data_raw_log_qrilc<- impute.QRILC(as.matrix(data_keep))[[1]]

    ## Centralization and scaling ##
    data_raw_log_qrilc_sc <- scale_recover(as.matrix(data_raw_log_qrilc), method = 'scale')
    ## Data after centralization and scaling ##
    data_raw_log_qrilc_sc_df <- data_raw_log_qrilc_sc[[1]]


    ## Parameters for centralization and scaling ##
    ## For scaling recovery ##
    data_raw_log_qrilc_sc_df_param <- data_raw_log_qrilc_sc[[2]]
    ## NA position ##
    NA_pos <- which(is.na(data_keep), arr.ind = T)
    ## bala bala bala ##
    data_raw_log_sc <- as.matrix(data_raw_log_qrilc_sc_df)
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

    data_imp[newData[,miss==F]==0]<-0

    newData[,miss==F]<-data_imp
    newData[,miss==T]<-data_rm


    #imputed_data<-pre_processing_GS_wrapper(data)
    #index<- which(methods=="GSimp_Real")

    results_data<- newData
    rownames(results_data)<-rownames
    colnames(results_data)<-colnames



  }


  if(method=="RQRILC"){
    rep_groups <- c(rep(1:(nrow(data)/reps), times=1, each=reps))
    rownames<-rownames(data)
    colnames<-colnames(data)
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

    results_data<-impute.QRILC(dataSet.mvs = newData)[[1]]
    results_data[results_data<=0.001]<-0
    rownames(results_data)<-rownames
    colnames(results_data)<-colnames


  }

  if(method=='RBPCA'){
    rep_groups <- c(rep(1:(nrow(data)/reps), times=1, each=reps))
    rownames<-rownames(data)
    colnames<-colnames(data)
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

    pc <-pcaMethods:: pca(object = newData, method="bpca", nPcs=min(c(10,ncol(data)/2)))


    results_data<-completeObs(pc)
    rownames(results_data)<-rownames
    colnames(results_data)<-colnames


  }


  if (method=='min' | method=='MIN') {
    imputed_data <- data
    if(local==T){
      rep_groups <- c(rep(1:(nrow(data)/reps), times=1, each=reps))
      data<- cbind(data,rep_groups)

      newData <-data.frame(matrix(nrow=nrow(data),ncol=ncol(data)-1))
      for (i in 1:length(unique(data$rep_groups))){
        tempData<-data[data[,ncol(data)]==i,]
        for (j in 1:(ncol(data)-1)){
          if(sum(is.na(tempData[,j]))==reps){
            tempData[,j]<-0
          }else{
           min<- min(tempData[,j],na.rm = TRUE)
            tempData[is.na(tempData[,j]),j]<- min

            newData[((reps*i)-(reps-1)):(reps*i),]<-tempData[,1:(ncol(tempData)-1)]
          }

        }
      }
      rownames(newData)<-rownames(data)
      colnames(newData)<-colnames(data)[1:(ncol(data)-1)]

      results_data<-newData


    }else{
      for(i in 1:ncol(data)){
        min<-min(data[,i],na.rm=TRUE)
        imputed_data[is.na(data[,i]),i]<-min

      }
      results_data<-imputed_data
    }

  }

  if (method=="halfmin") {
    imputed_data <- data

    if(local==T){
      rep_groups <- c(rep(1:(nrow(data)/reps), times=1, each=reps))
      data<- cbind(data,rep_groups)

      newData <-data.frame(matrix(nrow=nrow(data),ncol=ncol(data)-1))
      for (i in 1:length(unique(data$rep_groups))){
        tempData<-data[data[,ncol(data)]==i,]
        for (j in 1:(ncol(data)-1)){
          if(sum(is.na(tempData[,j]))==reps){
            tempData[,j]<-0
          }else{
          halfmin<- min(tempData[,j],na.rm = TRUE)/2
          tempData[is.na(tempData[,j]),j]<- halfmin

          newData[((reps*i)-(reps-1)):(reps*i),]<-tempData[,1:(ncol(tempData)-1)]
        }}
      }
      rownames(newData)<-rownames(data)
      colnames(newData)<-colnames(data)[1:(ncol(data)-1)]

      results_data<-newData


    }else{
      for(i in 1:ncol(data)){
        halfmin<-min(data[,i],na.rm=TRUE)/2
        imputed_data[is.na(data[,i]),i]<-halfmin

      }
      results_data<-imputed_data

    }

  }
  if (method=='mean' | method=='MEAN'){
    imputed_data <- data

    if(local==T){
      rep_groups <- c(rep(1:(nrow(data)/reps), times=1, each=reps))
      data<- cbind(data,rep_groups)

      newData <-data.frame(matrix(nrow=nrow(data),ncol=ncol(data)-1))
      for (i in 1:length(unique(data$rep_groups))){
        tempData<-data[data[,ncol(data)]==i,]
        for (j in 1:(ncol(data)-1)){
          if(sum(is.na(tempData[,j]))==reps){
            tempData[,j]<-0
          }else{
          mean<- mean(tempData[,j],na.rm = TRUE)
          tempData[is.na(tempData[,j]),j]<- mean

          newData[((reps*i)-(reps-1)):(reps*i),]<-tempData[,1:(ncol(tempData)-1)]
        }}
      }
      rownames(newData)<-rownames(data)
      colnames(newData)<-colnames(data)[1:(ncol(data)-1)]

      results_data<-newData


    }else{
      for(i in 1:ncol(data)){
        mean<-mean(data[,i],na.rm=TRUE)
        imputed_data[is.na(data[,i]),i]<-mean

      }
      results_data<-imputed_data

    }
  }

  if (method=='median' | method=='MEDIAN'){
    imputed_data <- data

    if(local==T){
      rep_groups <- c(rep(1:(nrow(data)/reps), times=1, each=reps))
      data<- cbind(data,rep_groups)

      newData <-data.frame(matrix(nrow=nrow(data),ncol=ncol(data)-1))
      for (i in 1:length(unique(data$rep_groups))){
        tempData<-data[data[,ncol(data)]==i,]
        for (j in 1:(ncol(data)-1)){
          if(sum(is.na(tempData[,j]))==reps){
            tempData[,j]<-0
          }else{
          median<- median(tempData[,j],na.rm = TRUE)
          tempData[is.na(tempData[,j]),j]<- median

          newData[((reps*i)-(reps-1)):(reps*i),]<-tempData[,1:(ncol(tempData)-1)]
        }
      }}
      rownames(newData)<-rownames(data)
      colnames(newData)<-colnames(data)[1:(ncol(data)-1)]

      results_data<-newData


    }else{
      for(i in 1:ncol(data)){
        median<-median(data[,i],na.rm=TRUE)
        imputed_data[is.na(data[,i]),i]<-median

      }
      results_data<-imputed_data

  }}

  if (method=='zero' |method=='ZERO'){
    imputed_data <- data

    for(i in 1:ncol(data)){
      imputed_data[is.na(data[,i]),i]<-0

    }

    results_data<-imputed_data

  }

  results_data[results_data<0]<-0.0

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
  data_raw_log_qrilc <-imputeLCMD::impute.QRILC(data_raw_log)[[1]]
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
#' @param local boolean whether to use a local imputation method or global, only applicable to non-replicate methods that do single value imputation
#' @return a list containing the dataframes of the imputed data
#' @export

imputeMulti<- function(methods, data, local=T, reps=NULL){
  results<-list()
  for (i in 1:length(methods)){
    results[[i]]<-impute(data, method = methods[i], local=local,reps = reps)
    print(paste("imputed using", methods[i], sep=" "))
    names(results)[i]<-methods[i]

  }

  return(results)

}


