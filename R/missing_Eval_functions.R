#' binaryMatrix
#'
#' A function returning a binary matrix with 0s for missing and 1s for present data
#' @param data A data frame
#' @return result A data frame with 0s for missing data and 1s for present data
#' @export

binaryMatrix<- function(data){
  result<-data
  for (i in 1:ncol(data)) {

    for (j in 1:nrow(data)){
      x<-d[c(j),c(i)]
      if (is.na(x)==TRUE){
        result[c(j),c(i)]<-0
      }
      else {
        result[c(j),c(i)]<-1
      }
    }
  }
  return(result)
}

#' variableStatistics
#'
#' A function that does various statistics on missingness per variables
#' @param data A data frame (rows as samples, columns as variables)
#' @param gof_method the method to compute goodness of fit with left-censored data 'kolmogorov' (default), 'cucconi'
#' @param correlation_method the method to compute correlations 'pearson' (default), 'spearman', 'kendall'
#' @return result A matrix containing: percent missingness, quartile cutoffs, mean, missingness type, variance, distribution
#' @export

variableStatistics<- function (data, correlation_method='pearson', gof_method='kolmogorov'){

  require(matrixStats)
  require(stats)
  require(ADGofTest)
  require(gamlss)
  options(scipen=999)

  out<-list()

  result<-data.frame(row.names=c("Missing Number", "Missing %:", "Min:", "Mean:", "Max:", "Variance:", "25th Percentile:",
                                 "50th Percentile:", "75th Percentile:", "Missingness Type:", "Distribution (by GAMLSS)"))


  missingNum<-0
  #Fills in Missing number and Missing %
  for (i in 1:ncol(data)){
    missList<-missingCounter(data[,i])
    result[1,i]<-missList[[1]]
    result[2,i]<-missList[[3]]
    result[3,i]<- min(data[,i], na.rm=TRUE)
    result[4,i]<- mean(data[,i], na.rm=TRUE)
    result[5,i]<- max(data[,i], na.rm=TRUE)
    result[6,i]<- var(data[,i], na.rm=TRUE)
    quant<-quantile(data[,i],na.rm=TRUE)

    result[7,i]<-quant[2]
    result[8,i]<-quant[3]
    result[9,i]<-quant[4]

    missingNum<-missingNum+ missList[[1]]

  }

  detect_MNAR<- detect.miss.MNAR.MAR(data, alpha=0.05, correlation_method = correlation_method)
  print('check1')
  MissingVar<-detect_MNAR[[1]]

  MAR_MNAR<- detect_MNAR[[3]]
  missigness <- detect.MCAR.MNAR.MAR (data,MissingVar,MAR_MNAR,alpha=0.05, percentage=.6, gof_Method=gof_method)
  print('check2')
  print(missigness)
  listMiss <- detect_missingness_type(missigness)
  print('check3')





  #fit_matrix<- lapply(data, function(x) fitDist(x, type="realplus")[[1]])

  #result[11:12,]<-data.frame(fit_matrix, stringsAsFactors = FALSE)
  result[10,]<-listMiss


  colnames(result)<-colnames(data)


  result2<-missingProportions(result)



  out$varStats<-result
  out$proportions<-result2
  out$percentMissing<-100*missingNum/(nrow(data)*ncol(data))

  return(out)
}


#' missingProportions
#'
#' A function that estimates missingness proportions across a data set
#' @param varStats a datafram of variable statistics as given by variableStatistics function
#' @return result a data frome of the results
#'

missingProportions<- function(varStats){

  result<-data.frame(row.names=c("MCAR","MAR","MNAR","EX","NONE"))
  result[,1]<-0
  result[,2]<-0
  colnames(result)<- c("Missingness Type By Variable", "Proportion of Overall Missingness")
  missingNum<-0
  tot_MNAR_Vars<-0
  tot_MAR_Vars<-0
  tot_MCAR_Vars<-0
  tot_EX_Vars<-0
  tot_NONE_Vars<-0
  tot_MNAR<-0
  tot_MAR<-0
  tot_MCAR<-0
  tot_EX<-0
  tot_NONE<-0

  for (i in 1:ncol(varStats)){
    print(i)
    if (varStats[10,i]=="MCAR"){
      tot_MCAR_Vars<-tot_MCAR_Vars+1
      tot_MCAR<-tot_MCAR+as.numeric(varStats[1,i])
    } else if (varStats[10,i]=="MAR"){
      tot_MAR_Vars<-tot_MAR_Vars+1
      tot_MAR<-tot_MAR+as.numeric(varStats[1,i])
    }else if (varStats[10,i]=="MNAR"){
      tot_MNAR_Vars<-tot_MNAR_Vars+1
      tot_MNAR<-tot_MNAR+as.numeric(varStats[1,i])
    }else if (varStats[10,i]=="EX"){
      tot_EX_Vars<-tot_EX_Vars+1
      tot_EX<-tot_EX+as.numeric(varStats[1,i])
    }else if (varStats[10,i]=="NONE"){
      tot_NONE_Vars<-tot_NONE_Vars+1
      tot_NONE<-tot_NONE+as.numeric(varStats[1,i])
    }

    missingNum<-missingNum+as.numeric(varStats[1,i])
  }
  cols<-ncol(varStats)

  result[1,1]<-100*tot_MCAR_Vars/cols
  result[2,1]<-100*tot_MAR_Vars/cols
  result[3,1]<-100*tot_MNAR_Vars/cols
  result[4,1]<-100*tot_EX_Vars/cols
  result[5,1]<-100*tot_NONE_Vars/cols

  result[1,2]<-100*tot_MCAR/missingNum
  result[2,2]<-100*tot_MAR/missingNum
  result[3,2]<-100*tot_MNAR/missingNum
  result[4,2]<-100*tot_EX/missingNum
  result[5,2]<-100*tot_NONE/missingNum

  return(result)

}


#'CullenFrey
#'
#'Function to produce a cullen and frey plot with every feature plotted
#'
#'@param data dataframe
#'
#'@export

CullenFrey<-function(data){

  a<-lapply(data, function(x) sum(is.na(x)))
  skew_list<-list()
  kurt_list<-list()
  data<-data[,(a<=(nrow(data)-3))]
  count<-1

  for (i in 1:ncol(data)){
    if (sum(!is.na(data[,i]))>4){


      skew_list[[count]]<-skewness(data[,i], na.rm = T)
      kurt_list[[count]]<-kurtosis(data[,i], na.rm = T)
      count<-count+1
    }


  }

  skew_DF<-data.frame(skew_list)
  kurt_DF<-data.frame(kurt_list)+3

  kurtmax<-max(10,ceiling(max(kurt_DF)))
  ymax<-kurtmax-1
  xmax<-max(4,ceiling(max(skew_DF)^2))

  skew_DF<-skew_DF^2
  kurt_DF<-kurtmax-kurt_DF

  CF_DF<-data.frame(skew_DF[1,])
  CF_DF[2,]<-kurt_DF[1,]
  CF_DF<-t(CF_DF)
  plot(x=CF_DF[,1],y=CF_DF[,2], yaxt="n",xlim=c(0,xmax),ylim=c(0,ymax), main="Cullen and Frey Plot", xlab="skewness squared", ylab="kurtosis")
  yax<-as.character(kurtmax-0:ymax)
  axis(side=2,at=0:ymax,labels=yax)

  # beta dist
  p<-exp(-100)
  lq<-seq(-100,100,0.1)
  q<-exp(lq)
  s2a<-(4*(q-p)^2*(p+q+1))/((p+q+2)^2*p*q)
  ya<-kurtmax-(3*(p+q+1)*(p*q*(p+q-6)+2*(p+q)^2)/(p*q*(p+q+2)*(p+q+3)))
  p<-exp(100)
  lq<-seq(-100,100,0.1)
  q<-exp(lq)
  s2b<-(4*(q-p)^2*(p+q+1))/((p+q+2)^2*p*q)
  yb<-kurtmax-(3*(p+q+1)*(p*q*(p+q-6)+2*(p+q)^2)/(p*q*(p+q+2)*(p+q+3)))
  s2<-c(s2a,s2b)
  y<-c(ya,yb)
  polygon(s2,y,col="lightgrey",border="lightgrey")
  # gamma dist
  lshape<-seq(-100,100,0.1)
  shape<-exp(lshape)
  s2<-4/shape
  y<-kurtmax-(3+6/shape)
  lines(s2,y,lty=2)
  # lnorm dist
  lshape<-seq(-100,100,0.1)
  shape<-exp(lshape)
  es2<-exp(shape^2)
  s2<-(es2+2)^2*(es2-1)
  y<-kurtmax-(es2^4+2*es2^3+3*es2^2-3)
  lines(s2,y,lty=3)

  points(x=CF_DF_T)


  obs.pch<-16
  obs.col<-"darkblue"


  legend(xmax*0.55,ymax*1.03,legend="Theoretical distributions",bty="n",cex=0.8)
  legend(xmax*0.6,0.98*ymax,pch=8,legend="normal",bty="n",cex=0.8,col="red")
  legend(xmax*0.6,0.94*ymax,pch=2,legend="uniform",bty="n",cex=0.8, col="red")
  legend(xmax*0.6,0.90*ymax,pch=7,legend="exponential",bty="n",cex=0.8, col="red")
  legend(xmax*0.6,0.86*ymax,pch=3,legend="logistic",bty="n",cex=0.8, col="red")
  legend(xmax*0.6,0.82*ymax,fill="grey80",legend="beta",bty="n",cex=0.8)
  legend(xmax*0.6,0.78*ymax,lty=3,legend="lognormal",bty="n",cex=0.8)
  legend(xmax*0.6,0.74*ymax,lty=2,legend="gamma",bty="n",cex=0.8)
  legend(xmax*0.58,0.69*ymax,legend=c("(Weibull is close to gamma and lognormal)"),
         bty="n",cex=0.6)

  points(0,kurtmax-3,pch=8,cex=1.5,lwd=2, col="red")
  points(0,kurtmax-9/5,pch=2,cex=1.5,lwd=2, col="red")
  # exp dist
  points(2^2,kurtmax-9,pch=7,cex=1.5,lwd=2, col="red")
  # logistic dist
  points(0,kurtmax-4.2,pch=3,cex=1.5,lwd=2, col="red")


}


#' missingCounter
#'
#' A function that counts missingness in a vector
#' @param vector A vector of values
#' @result list A list containing number of missing elements, percent missingness

missingCounter<- function(vector){
  missing<-0
  total<-0
  for (i in 1:length(vector)){
    if (is.na(vector[i])==TRUE){
      missing<-missing+1
    }
    total<-total+1
  }
  return(list(missing, total, 100*(missing/total)))
}



#' Title check.miss
#' check the percent missingness. From juuussi/impute-metabo
#' @param data
#' data matrix
#'
#' @return listVar
#'  list of 3 vectors
#'  MissingVar     = vector containing the columns numbers of the data matrix with missing values that are less than 80 % NAs
#'  ExcludedVar    = vector containing the columns numbers of the data matrix with missing values that are more than 80% NAs
#'  CompleteVar    = vector containing the columns numbers of the data matrix without NAs
#' @export

check.miss <- function(data,percentage = 0.80){
  if (is.data.frame(data)) {
    data <- as.matrix(data)

  }

  AllVar      <- as.numeric(1:ncol(data))
  MissingVar  <- numeric(0)
  ExcludedVar <- numeric(0)
  CompleteVar <- numeric(0)
  listVar <- list()



  perc.col <- round(colMeans(is.na(data)),digits = 2)


  for (j in 1:length(perc.col)){

    if(perc.col [j] == 0){
      CompleteVar <- c(CompleteVar,j)
    }
    if(perc.col [j] >= percentage){
      ExcludedVar <- c(ExcludedVar,j)
      cat('The variable: ',j ,' has a % miss  value:' , perc.col [j] ,"\n" )

    }else {
      MissingVar <- c(MissingVar,j)
    }

  }

  listVar <- list( MissingVar = MissingVar , ExcludedVar = ExcludedVar, CompleteVar = CompleteVar )

  return(listVar)
}




#' detect.miss.MNAR.MAR
#' It detects if the missigness depends on other varibales
#' Correlates an indicator matrix (0 not missing , 1 missing) with  the original data that can help determine
#' if variables tend to be missing together (MAR) or not (MCAR).
#' Kabacoff, Robert I. R in Action. manning, 2010.
#' From juuussi/impute-metabo
#'
#'
#' @param data ,
#' matrix with missing values
#' @param alpha,
#' significance level ,  default is 0.05
#'
#' @return results
#' list of vectors :
#' MissingVar  = vector containing the columns numbers of the data matrix with missing values that are less than 80 % NAs,
#' PairsCorVar = data frame containing the pairs of correlated variables in the data matrix,
#' MAR_MNAR    = vector containing the columns numbers of the data matrix with MAR or MNAR missingness
#' ExcludedVar = vector containing the columns numbers of the data matrix with missing values that are more than 80% NAs
#' @export
#'
#' @examples marietta <-detect.miss.MNAR.MAR (miss_data)


detect.miss.MNAR.MAR <- function(data,alpha=0.05, correlation_method) {
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }

  # check 80% missigness, exclude variables
  Miss_Var <- check.miss(data)

  data     <- data[,Miss_Var[[1]]]

  colnames(data) <- as.character(Miss_Var[[1]])

  #  Elements of x are 1 if a value in the data is missing and 0 if non-missing.

  x <- data.frame(abs(is.na(data)))
  colnames(x) <- colnames(data)

  if(any(is.na(data))){


    #  Extracting variables that have some missing values.
    #  if the sd is zero all elements in column are missing
    cols <- which(sapply(x, sd) > 0)
    y    <- x[,cols]
    y    <- as.matrix(y)

    #  the list of the all missing variables

    MissingVar <- as.numeric(colnames(y))

    #  Now, looking at the relationship between the
    #  presence of missing values in each variable and the observed values
    #  in other variables:


    #  correlation matrix
    corr_matrix <- stats:: cor(data, y, use="pairwise.complete.obs" ,method = correlation_method)

    #   calculate the probability values from the correlations
    CorTest     <- psych::corr.p(r=corr_matrix,n=nrow(data), adjust="fdr",alpha=.05)

    #THIS IS WHERE IT STOPS WORKING

    #   Correlations between variables in data and y together with confidence interval and pvalues
    table_P_CI  <- CorTest$ci

    #   find which variables are significally correlated

    SigPvalues  <- which(abs(table_P_CI$p) <=alpha)



    #   correlated Varibales pairs
    PairsCorVar <- data.frame(PairVar = rownames(table_P_CI)[SigPvalues])

    #   checking which of the columns have missing values from the pairs of correlated variables
    tmp         <- data.frame(do.call('rbind', strsplit(as.character(PairsCorVar$PairVar),'-',fixed=TRUE)))
    CorVar <- sort(union(as.numeric(as.character(tmp$X1)),as.numeric(as.character(tmp$X2))),decreasing = F)
    #   check which variables are trully missing
    CorVarTF    <-  is.element(CorVar,MissingVar)

    # MAR and MNAR variables
    MAR.MNAR    <- CorVar[CorVarTF]

    #  final results : all missing variables, the pairs of correlated variables and the MNAR_MAR Variables
    results     <- list(MissingVar = MissingVar,PairsCorVar = PairsCorVar,MAR_MNAR = MAR.MNAR ,ExcludedVar = Miss_Var[[2]])

    return(results)

  }else{
    print("matrix does not contain any missing values")
  }
}


#' detect.MCAR.MNAR.MAR
#' from juuussi/impute-metabo
#' @param data ,
#'  data matrix
#' @param MissingVar
#' vector containing the columns numbers of the data matrix with missing values that are less than 80 % (use the detect.miss.MNAR.MAR function)
#' @param MAR_MNAR
#'  vector containing the columns with missing values that have MAR or MNAR missigness (use the detect.miss.MNAR.MAR function)
#' @param gof_Method
#' method of goodness of fit testing, 'kolmogorov' or 'cucconi'
#' @param alpha ,
#' Significance level ,  default is 0.05
#' @param percentage
#  the percenatge of NAs in the MAR_MNAR variables, default 0.6 (above that threshold the ks.test doesnt work)
#'
#' @return results
#' list of  vectors :
#' MCAR = vector conatining the column numbers of MCAR variables,
#' MNAR = vector conatining the column numbers of MNAR variables ,
#' MAR  = vector conatining the column numbers of MAR variables ,
#' Excluded_marmnar  = vector conatining the column numbers of excluded MAR or MNAR variables that have more than 60% NAs  ,
#' ExcludedVar       = vector conatining the column numbers of excluded variables that have more than 80% NAs  ,
#' CompleteVar       = conatining the column numbers of variables without missing values.
#' @export
#'
#' @examples
detect.MCAR.MNAR.MAR <- function(data ,MissingVar, MAR_MNAR ,alpha = 0.05, percentage = 0.6, gof_Method='kolmogorov'){

  MAR  <- numeric(0)
  MNAR <- numeric(0)
  MAR  <- numeric(0)
  NewMAR_MNAR <- numeric(0)
  rm_MAR_MNAR <- numeric(0)

  if (is.data.frame(data)) {
    data <- as.matrix(data)

  }
  if(length(MissingVar) > 0){

    if (length(MAR_MNAR) == 0) {

      MCAR <- MissingVar

    }else{

      MCAR <- setdiff(MissingVar,MAR_MNAR)
      ## left truncation MNAR
      #  Kolmogorov-Smirnov test providing a comparison of a fitted distribution with the empirical distribution
      #  if the distributions are the same then p-values are high and that means they are left trancated if
      #  they are different then they are MAR

      #  Goodness of fit for left truncated data
      models <- list()
      Pval   <- list()
      Padj   <-list()
      perc.col    <- round(colMeans(is.na(data[,MAR_MNAR])),digits = 2)
      marmnar_mat  <- data[,MAR_MNAR]
      colnames(marmnar_mat) <- as.character(MAR_MNAR)
    if(gof_Method=='kolmogorov'){
      for (i in 1:length(MAR_MNAR)){
        # the KS test doesnt work for the 60 % NAs per column
        if (perc.col[i] < percentage){

          xt     <- na.omit( marmnar_mat[,i])
          threshold <- min(na.omit(marmnar_mat[,i]))
          #  truncgof::dplot(xt, "pnorm", list(mean(simulated_data),  sd(simulated_data)), H = threshold, vertical = TRUE)
          models <- c(models, list(truncgof::ks.test(xt, "pnorm",list(mean(data, na.rm = T),  sd(data,na.rm = T)), H = threshold,  alternative ="two.sided")))
          # pvalues from the KS test
          Pval <- c(Pval, tail(models,1)[[1]]$p.value)

        }else{
          cat("Variable the MAR_MNAR is MNAR: ",MAR_MNAR [i],"percentage % of NAs:",perc.col[i], "\n")
          marmnar_mat[is.na(marmnar_mat[,i]),i] <- 0
          # list of exluded variables above 60% missigness
          rm_MAR_MNAR <- c(rm_MAR_MNAR, MAR_MNAR [i])
        }
      }
      if  (length(rm_MAR_MNAR )== 0){
        # adjust p values obtained from KS test using fdr correction
        Padj <- p.adjust(Pval, method = "fdr")
        MAR  <- MAR_MNAR[which(Padj <= alpha)]
        MNAR <- setdiff(MAR_MNAR,MAR)
      }else {
        #  p values adjusted excluding the variables have MNAR or MAR that have been excluded because more than 60% NAs
        NewMAR_MNAR <- setdiff(MAR_MNAR,rm_MAR_MNAR)
        Padj <- p.adjust(Pval, method = "fdr")
        MAR  <- NewMAR_MNAR[which(Padj <= alpha)]
        MNAR <- setdiff(NewMAR_MNAR,MAR)
      }
    }else if(gof_Method=='cucconi'){
    print('inELSE')
    for (i in 1:length(MAR_MNAR)){
      # the KS test doesnt work for the 60 % NAs per column
      if (perc.col[i] < 0.6){

        sd<-sd(marmnar_mat[,i], na.rm=T)
        mean<-mean(marmnar_mat[,i], na.rm=T)
        min<-min(marmnar_mat[,i],na.rm = T)

        x <- seq(-(mean/sd), (mean/sd), length =nrow(marmnar_mat)) * sd + mean
        y<-x[x>=min]
        z<-cucconi.test(marmnar_mat[,i],y)
        #  truncgof::dplot(xt, "pnorm", list(mean(simulated_data),  sd(simulated_data)), H = threshold, vertical = TRUE)
        ##models <- c(models, list(truncgof::ks.test(xt, "pnorm",list(mean(data, na.rm = T),  sd(data,na.rm = T)), H = threshold,  alternative ="two.sided")))
        # pvalues from the KS test
        Pval <- c(Pval, z[[3]])

      }else{
        cat("Variable the MAR_MNAR is MNAR: ",MAR_MNAR [i],"percentage % of NAs:",perc.col[i], "\n")
        marmnar_mat[is.na(marmnar_mat[,i]),i] <- 0
        # list of exluded variables above 60% missigness
        rm_MAR_MNAR <- c(rm_MAR_MNAR, MAR_MNAR [i])
      }
    }
    if  (length(rm_MAR_MNAR )== 0){
      # adjust p values obtained from KS test using fdr correction
      Padj <- p.adjust(Pval, method = "fdr")
      MAR  <- MAR_MNAR[which(Padj <= alpha)]
      MNAR <- setdiff(MAR_MNAR,MAR)
    }else {
      #  p values adjusted excluding the variables have MNAR or MAR that have been excluded because more than 60% NAs
      NewMAR_MNAR <- setdiff(MAR_MNAR,rm_MAR_MNAR)
      Padj <- p.adjust(Pval, method = "fdr")
      MAR  <- NewMAR_MNAR[which(Padj <= alpha)]
      MNAR <- setdiff(NewMAR_MNAR,MAR)
    }
  }






  Miss_Var <-  check.miss(data)
  if (percentage >= 0.6){
    ExcludedVar <- Miss_Var[[2]]
    MNARnew     <- sort(c(MNAR,rm_MAR_MNAR),decreasing = F)
    results     <- list(MCAR = MCAR,MNAR = MNARnew ,MAR = MAR, Excluded_marmnar = numeric(0),ExcludedVar = Miss_Var[[2]] ,CompleteVar = Miss_Var[[3]])
    return(results)


  }else {

    results <- list(MCAR = MCAR,MNAR = MNAR ,MAR = MAR, Excluded_marmnar = rm_MAR_MNAR , ExcludedVar = Miss_Var[[2]], CompleteVar = Miss_Var[[3]])
    return(results)

  }

}
  }
}


#' detect_missingness_type
#' Vector conatining the the columns type of missingness. From juuussi/impute-metabo
#' @param missigness list of vectors obtained by the detect.MCAR.MNAR.MAR function
#'
#' @return miss_type =  vector conatining the the columns type of missingness
#' @export
#'
#' @examples
detect_missingness_type <- function(missigness){

  miss_type <- numeric(0)
  AllVar    <- numeric(0)
  MCAR <- missigness[[1]]
  MNAR <- missigness[[2]]
  MAR  <- missigness[[3]]
  ExcludedVar <- union(missigness[[4]],missigness[[5]])
  CompleteVar <- missigness [[6]]





  miss_type <- character(length(MCAR) + length(MNAR) + length(MAR) + length(ExcludedVar) + length(CompleteVar))
  miss_type[MCAR]        <- "MCAR"
  miss_type[MAR]         <- "MAR"
  miss_type[MNAR]        <- "MNAR"
  miss_type[ExcludedVar] <- "EX"
  miss_type[CompleteVar] <- "NONE"


  return(miss_type)
}

