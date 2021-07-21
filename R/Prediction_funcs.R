# Packages ----------------------------------------------------------------
#This is from GSimp
require(randomForest)
require(glmnet)
require(rpart)
require(FNN)

#'@export
lm_pred <- function(x, y) {
  require(randomForest)
  require(glmnet)
  require(rpart)
  require(FNN)
  data <- data.frame(y=y, x)
  model <- lm(y ~ ., data=data)
  y_hat <- predict(model, newdata=data)
  return(y_hat)
}

#'@export
rlm_pred <- function(x, y) {
  require(randomForest)
  require(glmnet)
  require(rpart)
  require(FNN)
  data <- data.frame(y=y, x)
  model <- rlm(y ~ ., data=data)
  y_hat <- predict(model, newdata=data)
  return(y_hat)
}

#'@export
rf_pred <- function(x, y, ntree=200, ...) {
  require(randomForest)
  require(glmnet)
  require(rpart)
  require(FNN)
  model <- randomForest(x=x, y=y, ntree=ntree, ...)
  y_hat <- predict(model, newdata=x)
  return(y_hat)
}

#'@export
glmnet_pred <- function(x, y, alpha=.5, lambda=.01) {
  require(randomForest)
  require(glmnet)
  require(rpart)
  require(FNN)
  x_mat <- as.matrix(x)
  model <- glmnet(x=x_mat, y=y, alpha=alpha, lambda=lambda)
  y_hat <- predict(model, newx=x_mat)[, 1]
  return(y_hat)
}

#'@export
rpart_pred <- function(x, y) {
  require(randomForest)
  require(glmnet)
  require(rpart)
  require(FNN)
  data <- data.frame(y=y, x)
  model <- rpart(y ~ ., data=data)
  y_hat <- predict(model, newdata=data)
  return(y_hat)
}

#'@export
knn_pred <- function(x, y) {
  require(randomForest)
  require(glmnet)
  require(rpart)
  require(FNN)
  model <- knn.reg(train=x, y=y, k=5)
  y_hat <- model$pred
  return(y_hat)
}
