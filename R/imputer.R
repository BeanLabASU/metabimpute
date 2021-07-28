


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
          tempData[,j]<-0.0
        }
        newData[((reps*i)-(reps-1)):(reps*i),]<-tempData[,1:(ncol(tempData)-1)]
      }
    }

    results_data<-impute.QRILC(dataSet.mvs = newData)[[1]]
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


#' Scale and recover -------------------------------------------------------
#' function copied from GSimp by WandeRum/GSimp
#'
scale_recover <- function(data, method='scale', param_df = NULL) {

  results <- list()
  data_res <- data
  if (!is.null(param_df)) {
    if (method=='scale') {
      data_res[] <- scale(data, center=param_df$mean, scale=param_df$std)
    } else if (method=='recover') {
      data_res[] <- t(t(data)*param_df$std+param_df$mean)
    }
  } else {
    if (method=='scale') {
      param_df<- data.frame(mean=colMeans(data), std=colSds(data))
      data_res[] <- scale(data, center=param_df$mean, scale=param_df$std)
    } else {stop('no param_df found for recover...')}
  }
  results[[1]] <- data_res
  results[[2]] <- param_df
  return(results)
}


#' Parallel combination ----------------------------------------------------
#' function copied from GSimp by WandeRum/GSimp
cbind_abind <- function(a, b) {
  res <- list()
  res$y_imp <- cbind(a$y_imp, b$y_imp)
  res$gibbs_res <- abind(a$gibbs_res, b$gibbs_res, along=2)
  return(res)
}


#'copied from WandeRum/GSimp
glmnet_pred <- function(x, y, alpha=.5, lambda=.01) {
  require(glmnet)
  x_mat <- as.matrix(x)
  model <- glmnet(x=x_mat, y=y, alpha=alpha, lambda=lambda)
  y_hat <- predict(model, newx=x_mat)[, 1]
  return(y_hat)
}

## Draw n samples from a truncated normal distribution N(mu, std^2|[lo, hi]) ##
#' copied from WandeRum/GSimp
rnorm_trunc <- function (n, mu, std, lo=-Inf, hi=Inf) {

  p_lo <- pnorm(lo, mu, std)
  p_hi <- pnorm(hi, mu, std)
  p_hi[p_hi < .01] <- .01
  u <- runif(n, p_lo, p_hi)
  return(qnorm(u, mu, std))
}



## Single missing variable imputation based on Gibbs sampler ##
#' copied from WandeRum/GSimp
single_impute_iters <- function(x, y, y_miss, y_real=NULL, imp_model='glmnet_pred', lo=-Inf, hi=Inf, iters_each=100, gibbs=c()) {

  y_res <- y
  x <- as.matrix(x)
  na_idx <- which(is.na(y_miss))
  imp_model_func <- getFunction(imp_model)
  nrmse_vec <- c()
  gibbs_res <- array(NA, dim=c(3, length(gibbs), iters_each))
  dimnames(gibbs_res) <- list(c('std', 'yhat', 'yres'), NULL, NULL)

  for (i in 1:iters_each) {
    y_hat <- imp_model_func(x, y_res)
    std <- sqrt(sum((y_hat[na_idx]-y_res[na_idx])^2)/length(na_idx))
    y_res[na_idx] <- rnorm_trunc(length(na_idx), y_hat[na_idx], std, lo, hi)
    if (length(gibbs)>0) {
      gibbs_res[1, , i] <- std
      gibbs_res[2, , i] <- y_hat[gibbs]
      gibbs_res[3, , i] <- y_res[gibbs]
    }
    ## The following code is for prediction function testing when y_real availabe ##
    if (!is.null(y_real)) {
      Sys.sleep(.5)
      par(mfrow=c(2, 2))
      nrmse_vec <- c(nrmse_vec, nrmse(y_res, y_miss, y_real))
      plot(y_real~y_res)
      plot(y_real~y_hat)
      plot(y_hat~y_res)
      plot(nrmse_vec)
    }
  }
  return(list(y_imp=y_res, gibbs_res=gibbs_res))
}

## Multiple missing variables imputation ##
## iters_each=number (100); vector of numbers, e.g. rep(100, 20) while iters_all=20
## lo/hi=numer; vector; functions like min/max/median/mean...
## initial=character ('qrilc'/'lysm'); initialized data maatrix
## n_cores=1 is sequentially (non-parallel) computing
#' copied from WandeRum/GSimp
multi_impute <- function(data_miss, iters_each=100, iters_all=20, initial='qrilc', lo=-Inf, hi='min',
                         n_cores=1, imp_model='glmnet_pred', gibbs=data.frame(row=integer(), col=integer())) {
  ## Convert to data.frame ##

  data_miss %<>% data.frame()

  ## Make vector for iters_each ##
  if (length(iters_each)==1) {
    iters_each <- rep(iters_each, iters_all)
  } else if (length(iters_each)==iters_all) {
    iters_each <- iters_each
  } else {stop('improper argument: iters_each')}


  ## Missing count in each column ##
  miss_count <- data_miss %>% apply(., 2, function(x) sum(is.na(x)))
  ## Index of missing variables, sorted (increasing) by the number of missings


  #miss_col_idx <- order(miss_count, decreasing = T) %>% extract(1:sum(miss_count!=0)) %>% rev()
  miss_col_idx <- order(miss_count, decreasing = T)[1:sum(miss_count!=0)] %>% rev()


  if (!all(gibbs$col %in% miss_col_idx)) {stop('improper argument: gibbs')}
  gibbs_sort <- gibbs
  if (nrow(gibbs_sort)>0) {
    gibbs_sort$order <- c(1:nrow(gibbs_sort))
    gibbs_sort <- gibbs_sort[order(gibbs_sort$row), ]
    gibbs_sort <- gibbs_sort[order(match(gibbs_sort$col, miss_col_idx)), ]
  } else {gibbs_sort$order <- integer()}

  ## Make vectors for lo and hi ##
  if (length(lo)>1) {
    if (length(lo)!=ncol(data_miss)) {stop('Length of lo should equal to one or the number of variables')}
    else {lo_vec <- lo}
  } else if (is.numeric(lo)) {
    lo_vec <- rep(lo, ncol(data_miss))
  } else if (is.character(lo)) {
    lo_fun <- getFunction(lo)
    lo_vec <- apply(data_miss, 2, function(x) x %>% na.omit %>% lo_fun)
  }

  if (length(hi)>1) {
    if (length(hi)!=ncol(data_miss)) {stop('Length of hi should equal to one or the number of variables')}
    else {hi_vec <- hi}
  } else if (is.numeric(hi)) {
    hi_vec <- rep(hi, ncol(data_miss))
  } else if (is.character(hi)) {
    hi_fun <- getFunction(hi)
    hi_vec <- apply(data_miss, 2, function(x) x %>% na.omit %>% hi_fun)
  }

  # Check whether lo is lower than hi
  if(!all(lo_vec < hi_vec)) {stop('lo should be lower than hi')}

  ## Initialization using build-in method or input initial matrix ##
  if(is.character(initial)) {
    data_init <- miss_init(data_miss, method=initial)
  } else if(is.data.frame(initial) & identical(data_miss[!is.na(data_miss)], initial[!is.na(data_miss)])) {
    data_init <- initial
  } else {stop('improper argument: initial')}

  data_imp <- data_init
  gibbs_res_final <- array(NA, dim=c(3, nrow(gibbs), 0))

  ## Iterations for the whole data matrix ##
  for (i in 1:iters_all) {
    cat('Iteration', i, 'start...')

    ## Parallel computing ##
    if (n_cores>1) {
      cat(paste0('Parallel computing (n_cores=', n_cores, ')...'))
      ## Parallel on missing variables
      cl <- makeCluster(n_cores)
      registerDoParallel(cl)
      core_res <- foreach (k=miss_col_idx, .combine='cbind_abind', .export=c('single_impute_iters', 'rnorm_trunc'), .packages=c('magrittr')) %dopar% {
        #source('Prediction_funcs.R')
        gibbs_sort_temp <- gibbs_sort[gibbs_sort$col==k, ]
        y_imp_res <- single_impute_iters(data_imp[, -k], data_imp[, k], data_miss[, k], imp_model=imp_model,
                                         lo=lo_vec[k], hi=hi_vec[k], iters_each=iters_each[i], gibbs=gibbs_sort_temp$row)
        y_imp_df <- y_imp_res$y_imp %>% data.frame
        colnames(y_imp_df) <- colnames(data_miss)[k]
        gibbs_res <- y_imp_res$gibbs_res
        list(y_imp=y_imp_df, gibbs_res=gibbs_res)
      }
      stopCluster(cl)
      y_imp_df <- core_res$y_imp
      gibbs_res_final <- abind(gibbs_res_final, core_res$gibbs_res, along=3)
      miss_col_idx_match <- match(colnames(y_imp_df), colnames(data_miss))
      data_imp[, miss_col_idx_match] <- y_imp_df
    } else {
      ## Sequential computing ##
      gibbs_res_j <- array(NA, dim=c(3, 0, iters_each[i]))
      for (j in miss_col_idx) {
        gibbs_sort_temp <- gibbs_sort[gibbs_sort$col==j, ]
        y_miss <- data_miss[, j]
        y_imp_res <- single_impute_iters(data_imp[, -j], data_imp[, j], y_miss, imp_model=imp_model, lo=lo_vec[j], hi=hi_vec[j],
                                         iters_each=iters_each[i], gibbs=gibbs_sort_temp$row)
        y_imp <- y_imp_res$y_imp
        gibbs_res_j <- abind(gibbs_res_j, y_imp_res$gibbs_res, along=2)
        data_imp[is.na(y_miss), j] <- y_imp[is.na(y_miss)]
      }
      gibbs_res_final <- abind(gibbs_res_final, gibbs_res_j, along=3)
    }
    cat('end!\n')
  }
  gibbs_res_final_reorder <- gibbs_res_final[, gibbs_sort$order, ]
  return(list(data_imp=data_imp, gibbs_res=gibbs_res_final_reorder))
}

# GS_impute ---------------------------------------------------------------
#' copied from WandeRum/GSimp
GS_impute <- multi_impute
