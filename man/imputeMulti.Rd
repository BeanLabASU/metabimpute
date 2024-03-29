% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imputer.R
\name{imputeMulti}
\alias{imputeMulti}
\title{imputeMulti}
\usage{
imputeMulti(methods, data, local = T, reps = NULL)
}
\arguments{
\item{methods}{a list of imputation methods}

\item{data}{the dataset with missing values contained}

\item{local}{boolean whether to use a local imputation method or global, only applicable to non-replicate methods that do single value imputation}

\item{reps}{number of replicates, if applicable}
}
\value{
a list containing the dataframes of the imputed data
}
\description{
Uses the impute function to do multiple specified imputation methods and contains them within a list
}
