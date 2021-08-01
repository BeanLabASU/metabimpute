#' correlationMatrix
#' function returning a correlation matrix and NA_correlation matrix
#'
#' @param data a dataframe
#' @return a list containing a correlation matrix, NA_correlation matrix
#' @export

correlationMatrix<-function(data){
  require(stats)
  require(ltm)
  out<-list()
  mat <- cor(data, use = "pairwise.complete.obs", method = "pearson")

  na_cor <- matrix(nrow = ncol(data), ncol = ncol(data))
  colnames(na_cor) <- colnames(data)
  row.names(na_cor) <- paste0(colnames(data), "_is.na")
  for (i in 1:ncol(data)) {
    for (j in 1:ncol(data)) {
      na_cor[i,j] <- biserial.cor(data[,j], is.na(data)[,i], use="complete.obs")
    }
  }
  mat[is.na(mat)]<-0
  out$corMat<-mat
  out$NA_corMat<-na_cor
  return(out)
}


#' simulate
#' function to simulate a similar data set with similar correlation matrix.
#' modified from source code written by Tibor V. Varga from https://github.com/Tirgit/missCompare
#'
#' @param corMat correlation matrix.
#' @param rownum number of rows of the original data matrix
#' @param colnum number of columns of the original data matrix
#' @return a matrix of the original size of your data matrix with simulated values (normalized)
#' @export

simulate <- function(rownum, colnum, corMat) {

  require(Matrix)
  require(MASS)
  require(stats)

  pd_corr_matrix <- nearPD(corMat, keepDiag = T, conv.tol = 1e-07, corr = T, maxit=1000, conv.norm.type = "F")
  mu <- rep(2, colnum)
  stddev <- rep(1, colnum)
  covMat <- stddev %*% t(stddev) * pd_corr_matrix$mat
  simMat <- mvrnorm(n = rownum, mu = mu, Sigma = covMat, empirical = FALSE)  # Simulated values

  rownames(simMat) <- 1:nrow(simMat)

  return(simMat)
}




#' Title simulate_missingness
#' Uses the simulated data to create the different types of missingness (MCAR,MNAR,MAR)
#' using different percentages of missigness
#' @param data matrix with simulated data
#' @param mcar percenatge of missigness in Missing Completly at Random
#' @param mar percenatge of missigness in Missing at Random
#' @param mnar percenatge of missignessMisiing not at Random
#' @param mnar.type type of truncation 'left' or right
#' @param mar.type type of truncation
#'
#' @return simulated_data
#' @export

simulate_missingness <- function(data, mcar=0, mar=0, mnar=0, mnar.type="left", mar.type="left") {

  if(class(data) != "matrix") {
    stop("Variable data should be a matrix.")
  }
  if ( ((mcar + mar + mnar) > 1) || ((mcar + mar + mnar) < 0)) {
    stop("Sum of mcar, mar and mnar should be between 0 and 1.")
  }
  print(1)
  simulated_data <- data
  if (mcar > 0){
    mcar_distribution   = runif(nrow(data)*ncol(data), min=0, max=1)
    simulated_data = matrix(ifelse(mcar_distribution<mcar, NA, data), nrow=nrow(data), ncol=ncol(data))
  }
  print(2)
  if (mnar > 0) {
    added_mnar <- 0
    initial_nas <- sum(colSums(is.na(simulated_data)))
    current_nas <- initial_nas
    simulated_data_n <- ncol(simulated_data) * nrow(simulated_data)
    # condition of porpotion of missigness
    while (((current_nas - initial_nas) / simulated_data_n) < mnar) {

      # Select random variable
      variable_index <- sample(1:ncol(simulated_data), 1)

      # What percentage of variable to set missing
      cut_percentage <- 0
      while (cut_percentage <= 0) {
        cut_percentage <- rchisq(1, df=1) / 30
      }
      cut_percentage <- min(cut_percentage, 1)

      # How many values to set missing
      cut_index <- min(floor(cut_percentage * nrow(simulated_data)), nrow(data)-3)

      sorted_variable <- sort(simulated_data[,variable_index])

      # Set values to missing
      if (mnar.type == "right") {
        # Corresponding cut-off point for values
        cut_point <- sorted_variable[length(sorted_variable) - cut_index]

        #making sure that no columns end up with one or no data points
        if ((sum(!is.na(simulated_data[,variable_index]))- cut_index)>2){
          simulated_data[simulated_data[,variable_index] < cut_point, variable_index] <- NA
        }

      } else {
        # Corresponding cut-off point for values
        cut_point <- sorted_variable[cut_index]

        #making sure that no columns end up with one or no data points
        if ((sum(!is.na(simulated_data[,variable_index]))- cut_index)>2){
          simulated_data[simulated_data[,variable_index] < cut_point, variable_index] <- NA
        }

      }

      # Counter to check how much MNAR missingness has been added to data
      current_nas <- sum(colSums(is.na(simulated_data)))
    }
  }
  print(3)
  if (mar > 0) {
    added_mar <- 0
    initial_nas <- sum(colSums(is.na(simulated_data)))
    current_nas <- initial_nas
    simulated_data_n <- ncol(simulated_data) * nrow(simulated_data)
    # condition of porpotion of missigness
    while (((current_nas - initial_nas) / simulated_data_n) < mar) {

      # Select random variable
      variable_index <- sample(1:ncol(simulated_data), 1)
      variable_index2 <- sample(setdiff(1:ncol(simulated_data), variable_index), 1)


      # What percentage of variable to set missing
      cut_percentage <- 0
      while (cut_percentage <= 0) {
        cut_percentage <- rchisq(1, df=1) / 30
      }
      cut_percentage <- min(cut_percentage, 1)

      # How many values to set missing
      cut_index <- min(floor(cut_percentage * nrow(simulated_data)), nrow(data)-3)

      sorted_variable <- sort(simulated_data[,variable_index])

      # Set values to missing
      if (mar.type == "right") {
        # Corresponding cut-off point for values
        cut_point <- sorted_variable[length(sorted_variable) - cut_index]
        #only allows the code to remove data if there are more than 3 values present
        if ((sum(!is.na(simulated_data[,variable_index2]))- cut_index)>2){
          simulated_data[simulated_data[,variable_index] < cut_point, variable_index2] <- NA
        }

      } else {
        # Corresponding cut-off point for values
        cut_point <- sorted_variable[cut_index]
        #only allows the code to remove data if there are more than 3 values present
        if ((sum(!is.na(simulated_data[,variable_index2]))- cut_index)>2){
          simulated_data[simulated_data[,variable_index] < cut_point, variable_index2] <- NA
        }

      }

      # Counter to check how much MNAR missingness has been added to data
      current_nas <- sum(colSums(is.na(simulated_data)))
    }
  }

  simulated_data

}


#' simulateEngine
#'
#' This function is going to be the main engine that actually runs the simulations and tracks the data and the
#' error calculations at differing levels of missingness, ratios of missingness and iterations of the sim.
#' @param data the original data set
#' @param simIter this is the number of unique simulation matrices to be made from the covMat of the origData
#' @param simMissIter this is the number of simulated missing data matrices to be created for a given missingness
#' proportion
#' @param missMax the highest amount of missingness wanted in the sim
#' @param missMin the min amount of missingness wanted in the sim
#' @param missInc the amount of missingness to increment by
#' @param missRatios this is a list of triplets defining the ratio of MCAR, MAR, MNAR in form c(0,0.5, 0.5) would
#' specify that we want a split of 50-50 MAR MNAR for each missingness level
#' @param methodsImp the imputation methods to be utilized.
#' @param medthodsEval the methods of error evaluation to be utilized.
#'
#' @returns a collection of dataframes containing the errors of imp methods for each eval method listed by
#' proportion of missingness and total percent of missingness.
#' @export

simulateEngine<- function (data, simIter, simMissIter, missMax, missMin, missInc, missRatios, methodsImp,
                           methodsEval, simulate_Data=T, reps){

  rep_groups <- c(rep(1:(nrow(data)/reps), times=1, each=reps))
  startTime<-Sys.time()
  resultList<- list()
  simList<-list()
  missingMatrixList<-list()
  counter<-1
  for (i in 1:simIter){
    if(simulate_Data==T){
      simData<-simulate(rownum = nrow(data), colnum = ncol(data), corMat = correlationMatrix(data)$corMat)
      simList[[i]]<-simData

    }else{
      simData<-as.matrix(data)
      simList[[i]]<-simData
    }

    print(paste("Simulated Data Matrix: ", i, "created", sep=" "))

    missLevel<-missMin
    percentList<-list()
    for(z in 1:(((missMax-missMin)/missInc)+1)){
      ratioList<-list()
      for (k in 1:length(missRatios)){

        mcar<- missRatios[[k]][1]*missLevel
        mar<- missRatios[[k]][2]*missLevel
        mnar<-missRatios[[k]][3]*missLevel

        missingDataList<-list()

        for (m in 1:simMissIter){

          if(simulate_Data==T){
            missData<- simulate_missingness(data=simData, mcar=mcar, mar=mar, mnar=mnar)
          }else{
            missData<- simulate_missingness(data=simData, mcar=mcar, mar=mar, mnar=mnar)
          }


          missingMatrixList[[counter]]<-missData
          names(missingMatrixList)[counter]<-paste("Sim matrix: " , i,"Missing matrix", m, "of",
                                                   missLevel, "percent with ratios", mcar, mar, mnar, sep=" " )
          counter<-counter+1
          print(paste("Missing matrix", m, "of", missLevel, "percent with ratios", mcar, mar, mnar, sep=" "))

          imputeResults<-imputeMulti(methods=methodsImp, data=missData, reps=reps)
          errors<-simEval(origData=simData, missData=missData, methods=methodsEval, imputationResults = imputeResults, simulate_Data=simulate_Data)
          print("evaluted errors, added to list")
          missingDataList[[m]]<-errors

        }

        #initialize and use helper function
        ratioList[[k]]<-aggregateDFList(missingDataList)
        print("aggregated errors for all copies of missing matrix")
        names(ratioList)[k]<- paste('MCAR:', mcar, 'MAR:',
                                    mar, 'MNAR:', mnar, sep=' ' )

      }

      #add ratio list (which were aggregated) to percent list
      #z<-((missLevel-missMin)/missInc)+1
      print(z)
      percentList[[z]]<-ratioList
      names(percentList)[z]<- paste('percent:', missLevel, sep=' ')

      missLevel<-missLevel+missInc

    }

    # add percent list to result list (remember this will have multiple levels, all corresponding ones will need
    # aggregation)
    resultList[[i]]<-percentList
    names(resultList)[i]<- paste('sim: ', i, sep=' ')


  }


  endTime<-Sys.time()

  return(list(resultList, (endTime-startTime), simList, missingMatrixList))

}


#' aggregateDFList
#'
#' Helper function to aggregate and give means of a list of same sized data frames
#'
#' @param list a list of same sized dataframes
#' @result a data frame

aggregateDFList<- function (list ){
  result<-list[[1]]
  if (length(list)>1){
    for (i in 1:nrow(list[[1]])){
    for (j in 1:ncol(list[[1]])){
      sum<-0
      for (k in seq_along(list)){
        sum<- sum+list[[k]][i,j]
      }
      result[i,j]<-sum/length(list)

    }


  }
  }


  return(result)
}


#' rearrangeList
#'
#' The purpose of this functino is to take the results from simulateEngine and output a list of data frames
#' that are useable by ggPlots
#' @param resultList
#' @param missMax from simulate engine
#' @param missMin from simulate engine
#' @param missInc from simulae engine
#' @param simIter from simulate engine
#' @param missRatios from sim engine
#' @return a ggplot formatted list
#' @export

rearrangeList<- function (resultList, missRatios, missMax, missMin, missInc, simIter){
  #reorganizing the results of resultsList so that we list by ratio first, percent second and aggregate all
  #together
  results<-list()



  for (i in 1:length(missRatios)){
    perPercentList<-list()
    for (j in 1:(((missMax-missMin)/missInc)+1)){
      perSimList<-list()
      for (k in 1:simIter){
        temp_simList<-resultList[[k]]
        temp_percentList<-temp_simList[[j]]
        perSimList[[k]]<-temp_percentList[[i]]



      }

      #aggregate the per simulation per missing% and ratio
      perPercentList[[j]]<-aggregateDFList(perSimList)
      names(perPercentList)[j]<- paste('Percent: ', (missMin+((j-1)* missInc)), sep=' ')
    }

    results[[i]]<- perPercentList
    names(results)[i]<- paste('MCAR:',  missRatios[[i]][1], 'MAR: ', missRatios[[i]][2],
                              'MNAR: ', missRatios[[i]][3], sep=' ')

  }


  #once more refactoring into a graphable usefulness
  multiplotList<-list()
  for (i in 1:length(results)){ # each missing ratio
    errorEvalList<-list()
    for (j in 1:ncol(results[[i]][[1]])){ #each error eval type
      df<-data.frame(row.names =rownames(results[[i]][[1]]))
      for (k in 1:length(results[[i]])){ #each percent
        df[,k]<- results[[i]][[k]][,j]
        colnames(df)[k]<-substring(names(results[[i]][k]),11)
      }
      errorEvalList[[j]]<-df
      names(errorEvalList)[j]<-colnames(results[[i]][[1]])[j]
    }
    multiplotList[[i]]<-errorEvalList
    names(multiplotList)[i]<-names(results[i])
  }


  ggplotlist<-list()
  for (i in 1:length(multiplotList)){
    tempList<-list()
    for (j in 1:length(multiplotList[[i]])){
      df<-data.frame()

      for (k in 1:nrow(multiplotList[[i]][[j]])){
        for (m in 1:ncol(multiplotList[[i]][[j]])){
          df[m, 1]<- as.numeric(colnames(multiplotList[[i]][[j]])[m])
          df[m, 1+k]<- multiplotList[[i]][[j]][k,m]
          colnames(df)[1+k]<- rownames(multiplotList[[i]][[j]])[k]

        }
      }
      colnames(df)[1]<-"Percent"

      tempList[[j]]<-df
      names(tempList)[j]<-names(multiplotList[[i]])[j]
    }

    ggplotlist[[i]]<-tempList
    names(ggplotlist)[i]<-names(multiplotList)[i]

  }
  return (list(ggplotlist, multiplotList, results))
}
