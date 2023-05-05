library(nlme)
library(MASS)

make_Z_noPcov <- function(all_data,data_frame){
  ## This function constructs the Z design matrix for a new dataset.
  
  cov_data_used_only <- all_data
  new_country <- cov_data_used_only$country
  unq_cnt <- unique(new_country)
  ni <- length(new_country[new_country==unq_cnt[1]])
  Z_b <- matrix(1,ni,1)
  unq_cnt <- unique(new_country)
  for(k in 2:length(unq_cnt)){
    ni <- length(new_country[new_country==unq_cnt[k]])
    Z_b2 <- rbind(Z_b,matrix(0,ni,dim(Z_b)[2]))
    Z_b <- cbind(Z_b2,c(rep(0,dim(Z_b)[1]),rep(1,ni)))
  }
  return(Z_b)
}


not_in <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))



pred_fun_noPcov <- function(X,Z,Y,fitted_model,N,SE_var=1,full_b=0,S=0){  
  ## This function predicts the EBLUPs, predicted 
  # values and their variances for a particular X 
  # and Z design matrices.
  ##ARGUMENTS:
  # X: fixed effect design matrix.
  # Z: random effect design matrix.
  # fitted_model: fitted lme object.
  # N: number of levels used in the orginal analysis
  # SE_var: the residual variance covariate (constant 
  #         plus SE_var only)
  # full_b: if X and Z are from new data the original 
  #         b must be given
  
  # Export variance components from the fitted model.
  
  # Calculate the marginal covariance, fixed effect mean, and random effecs.
  gamma <- exp(as.vector(fitted_model$modelStruct$varStruct))-1
  sigma2 <- fitted_model$sigma^2
  
  C <- as.matrix(cbind(X,Z))
  R <- diag(length(sigma2*(1+SE_var%*%gamma)))*sigma2*(c(1+SE_var%*%gamma)^2)
  if(length(full_b)==1){
    b_var <- as.matrix(fitted_model$modelStruct$reStruct$country)*sigma2
    D <- diag(N)*b_var[1]
    R_inv <- diag(length(sigma2*(1+SE_var%*%gamma)))*1/(sigma2*(c(1+SE_var%*%gamma)^2))
    D_inv <- diag(N)*1/b_var[1]
    G <- t(C)%*%R_inv%*%C
    B <- blockMatrixDiagonal(0*diag(dim(X)[2]),D_inv)
    S <- ginv(G + B)
    Sigma <- Z%*%D%*%t(Z) + R
    inv_test <- try(Sigma_inv <- ginv(Sigma),silent = TRUE)
    if(!is.null(attr(inv_test,"class"))){Sigma_inv <- diag(1/diag(Sigma))}
    beta_hat <- ginv(t(X)%*%Sigma_inv%*%X)%*%t(X)%*%Sigma_inv%*%Y
    full_b <- c(beta_hat,D%*%t(Z)%*%Sigma_inv%*%(Y - X%*%beta_hat))
  }
  
  mu_T <- C%*%full_b
  sigma_T <- C%*%S%*%t(C)
  sigma_Y <- R+sigma_T 
  
  # Export the mean, covariance matrix, and random effects 
  preds <- list(mu_T=mu_T,sigma_T=sigma_T,sigma_Y=sigma_Y,b=full_b,D=D,S=S)
  return(preds)
}


pred_fun_noPcov_noC <- function(X,Y,fitted_model,N,SE_var=1,full_b=0,S=0){  
  ## This function predicts the EBLUPs, predicted 
  # values and their variances for a particular X 
  # maxtrix for observations with no Z matrix.
  ##ARGUMENTS:
  # X: fixed effect design matrix.
  # fitted_model: fitted lme object.
  # N: number of levels used in the orginal analysis
  # SE_var: the residual variance covariate (constant 
  #         plus SE_var only)
  # full_b: if X and Z are from new data the original 
  #         b must be given
  
  # Export variance components from the fitted model.
  
  # Calculate the marginal covariance, fixed effect mean, and random effecs.
  gamma <- exp(as.vector(fitted_model$modelStruct$varStruct))-1
  sigma2 <- fitted_model$sigma^2
  
  R <- diag(length(sigma2*(1+SE_var%*%gamma)))*sigma2*(c(1+SE_var%*%gamma)^2)
  mu_T <- X%*%full_b[1:(dim(X)[2])]
  
  b_var <- as.matrix(fitted_model$modelStruct$reStruct$country)*sigma2
  D <- diag(N)*b_var[1]
  full_b <- c(full_b[1:(dim(X)[2])],0*full_b[(dim(X)[2]+1):length(full_b)])
  
  sigma_T <- X%*%S[1:(dim(X)[2]),1:(dim(X)[2])]%*%t(X) + R + b_var[1]
  sigma_Y <- R+sigma_T 
  
  # Export the mean, covariance matrix, and random effects 
  preds <- list(mu_T=mu_T,sigma_T=sigma_T,sigma_Y=sigma_Y)
  return(preds)
}


pred_func_NoN <- function(all_data,data_frame,fitted_model){
  # Make the X and Z matrices for the original fitted model
  Xi.matrix <- cbind(1,log(data_frame$r1mr),data_frame$source_cat1,data_frame$source_cat2)
  SE.matrix <- cbind(data_frame$source_cat1,data_frame$source_cat2)
  N <- length(unique(data_frame$country))
  X_orig <- Xi.matrix
  Z_orig <- make_Z_noPcov(data_frame,data_frame)
  
  # Estimate EBLUPs, predicted means, and variances of the predicted means.
  pred_mat <- pred_fun_noPcov(X=X_orig,Z=Z_orig,Y=data_frame$Y,fitted_model,N,SE_var=SE.matrix,full_b = 1)
  S <- pred_mat$S
  b <- pred_mat$b
  
  
  # Extracting only those countries that were used in the original analysis
  cov_data_used_only <- all_data[all_data$country %in% data_frame$country,]
  # Ording the data so the random effects align.
  div=c(cov_data_used_only$country)
  cov_data_used_only <- cov_data_used_only[order(div,cov_data_used_only$year),]
  # Making the X and Z matrices for the new data.
  
  X <- cbind(1,log(cov_data_used_only$r1mr),cov_data_used_only$source_cat1,cov_data_used_only$source_cat2)
  SE<- cbind(cov_data_used_only$source_cat1,cov_data_used_only$source_cat2)
  Z <- make_Z_noPcov(cov_data_used_only,data_frame)
  
  # Getting rid of NA's then making X and Z again
  cov_data_used_only <- cov_data_used_only[!c(apply((apply(cbind(X,Z,SE),1,is.na)),2,any)),]
  # Making the X and Z matrices for the new data.
  X <- cbind(1,log(cov_data_used_only$r1mr),cov_data_used_only$source_cat1,cov_data_used_only$source_cat2)
  SE<- cbind(cov_data_used_only$source_cat1,cov_data_used_only$source_cat2)
  Z <- make_Z_noPcov(cov_data_used_only,data_frame)
  
  # Predicting the mean and variance for each new observation.
  pred_mat_new <- pred_fun_noPcov(X,Z,Y=0,fitted_model,N,SE_var=SE,full_b = b,S=S)
  pred_all  <- pred_mat_new$mu_T
  sigma_Y_all <- diag(pred_mat_new$sigma_Y)
  pred_new_data <- data.frame(cov_data_used_only,pred_all,sigma_Y_all)
  
  all_data <- merge(all_data,pred_new_data,all.x=TRUE)
  
  # Extracting countries that were not used in the original analysis
  cov_data_not_used <- all_data[not_in(all_data$country,data_frame$country),]
  if(dim(cov_data_not_used)[1]>0){
  # Ording the data so the random effects align.
  div=c(cov_data_not_used$country)
  cov_data_not_used <- cov_data_not_used[order(div,cov_data_not_used$year),]
  # Making the X and SE matrices for the new data.
  X <- cbind(1,log(cov_data_not_used$r1mr),cov_data_not_used$source_cat1,cov_data_not_used$source_cat2)
  SE<- cbind(cov_data_not_used$source_cat1,cov_data_not_used$source_cat2)
  
  # Getting rid of NA's then making X and SE again
  cov_data_not_used <- cov_data_not_used[!c(apply((apply(cbind(X,SE),1,is.na)),2,any)),]
  # Making the X and SE matrices for the new data.
  X <- cbind(1,log(cov_data_not_used$r1mr),cov_data_not_used$source_cat1,cov_data_not_used$source_cat2)
  SE<- cbind(cov_data_not_used$source_cat1,cov_data_not_used$source_cat2)
  N <- length(unique(cov_data_not_used$country))
  pred_mat_new_noC <- pred_fun_noPcov_noC(X,Y,fitted_model,N,SE_var=SE,full_b=b,S=S)
  pred_all_noC  <- pred_mat_new_noC$mu_T
  sigma_Y_all_noC <- diag(pred_mat_new_noC$sigma_Y)
  pred_new_data_noC <- data.frame(cov_data_not_used,pred_all_noC,sigma_Y_all_noC)
  
  all_data <- merge(all_data,pred_new_data_noC,all.x=TRUE)
  all_data$pred_all[is.na(all_data$pred_all)] <- all_data$pred_all_noC[is.na(all_data$pred_all)]
  all_data$sigma_Y_all[is.na(all_data$sigma_Y_all)] <- all_data$sigma_Y_all_noC[is.na(all_data$sigma_Y_all)]
  }
  div=c(all_data$country)
  all_data <- all_data[order(div,all_data$year),]
  if(dim(cov_data_not_used)[1]>0){all_data <- subset(all_data, select = -c(pred_all_noC,sigma_Y_all_noC) )}
  
  cov_w_preds <- data.frame(all_data)
  return(cov_w_preds)
}


pred_func <- function(all_data,data_frame,fitted_model){
  # Make the X and Z matrices for the original fitted model
  Xi.matrix <- cbind(1,log(data_frame$r1mr),log(data_frame$hazB_N),data_frame$source_cat1,data_frame$source_cat2)
  SE.matrix <- cbind(data_frame$source_cat1,data_frame$source_cat2)
  N <- length(unique(data_frame$country))
  X_orig <- Xi.matrix
  Z_orig <- make_Z_noPcov(data_frame,data_frame)
  
  # Estimate EBLUPs, predicted means, and variances of the predicted means.
  pred_mat <- pred_fun_noPcov(X=X_orig,Z=Z_orig,Y=data_frame$Y,fitted_model,N,SE_var=SE.matrix,full_b = 1)
  S <- pred_mat$S
  b <- pred_mat$b
  
  
  # Extracting only those countries that were used in the original analysis
  cov_data_used_only <- all_data[all_data$country %in% data_frame$country,]
  cov_data_not_used <- all_data[not_in(all_data$country,data_frame$country),]
  # Ording the data so the random effects align.
  div=c(cov_data_used_only$country)
  cov_data_used_only <- cov_data_used_only[order(div,cov_data_used_only$year),]
  # Making the X and Z matrices for the new data.
  
  X <- cbind(1,log(cov_data_used_only$r1mr),log(cov_data_used_only$hazB_N),cov_data_used_only$source_cat1,cov_data_used_only$source_cat2)
  SE<- cbind(cov_data_used_only$source_cat1,cov_data_used_only$source_cat2)
  Z <- make_Z_noPcov(cov_data_used_only,data_frame)
  
  # Getting rid of NA's then making X and Z again
  cov_data_used_only <- cov_data_used_only[!c(apply((apply(cbind(X,Z,SE),1,is.na)),2,any)),]
  # Making the X and Z matrices for the new data.
  X <- cbind(1,log(cov_data_used_only$r1mr),log(cov_data_used_only$hazB_N),cov_data_used_only$source_cat1,cov_data_used_only$source_cat2)
  SE<- cbind(cov_data_used_only$source_cat1,cov_data_used_only$source_cat2)
  Z <- make_Z_noPcov(cov_data_used_only,data_frame)
  
  # Predicting the mean and variance for each new observation.
  pred_mat_new <- pred_fun_noPcov(X,Z,Y=0,fitted_model,N,SE_var=SE,full_b = b,S=S)
  pred_all  <- pred_mat_new$mu_T
  sigma_Y_all <- diag(pred_mat_new$sigma_Y)
  pred_new_data <- data.frame(cov_data_used_only,pred_all,sigma_Y_all)
  
  all_data <- merge(all_data,pred_new_data,all.x=TRUE)
  
  cov_data_not_used <- all_data[not_in(all_data$country,data_frame$country),]
  if(dim(cov_data_not_used)[1]>0){
    # Ording the data so the random effects align.
    div=c(cov_data_not_used$country)
    cov_data_not_used <- cov_data_not_used[order(div,cov_data_not_used$year),]
    # Making the X and SE matrices for the new data.
    X <- cbind(1,log(cov_data_not_used$r1mr),log(cov_data_not_used$hazB_N),cov_data_not_used$source_cat1,cov_data_not_used$source_cat2)
    SE<- cbind(cov_data_not_used$source_cat1,cov_data_not_used$source_cat2)
    
    # Getting rid of NA's then making X and SE again
    cov_data_not_used <- cov_data_not_used[!c(apply((apply(cbind(X,SE),1,is.na)),2,any)),]
    # Making the X and SE matrices for the new data.
    X <- cbind(1,log(cov_data_not_used$r1mr),log(cov_data_not_used$hazB_N),cov_data_not_used$source_cat1,cov_data_not_used$source_cat2)
    SE<- cbind(cov_data_not_used$source_cat1,cov_data_not_used$source_cat2)
    N <- length(unique(cov_data_not_used$country))
    pred_mat_new_noC <- pred_fun_noPcov_noC(X,Y,fitted_model,N,SE_var=SE,full_b=b,S=S)
    pred_all_noC  <- pred_mat_new_noC$mu_T
    sigma_Y_all_noC <- diag(pred_mat_new_noC$sigma_Y)
    pred_new_data_noC <- data.frame(cov_data_not_used,pred_all_noC,sigma_Y_all_noC)
    
    all_data <- merge(all_data,pred_new_data_noC,all.x=TRUE)
    all_data$pred_all[is.na(all_data$pred_all)] <- all_data$pred_all_noC[is.na(all_data$pred_all)]
    all_data$sigma_Y_all[is.na(all_data$sigma_Y_all)] <- all_data$sigma_Y_all_noC[is.na(all_data$sigma_Y_all)]
  }
  div=c(all_data$country)
  all_data <- all_data[order(div,all_data$year),]
  if(dim(cov_data_not_used)[1]>0){all_data <- subset(all_data, select = -c(pred_all_noC,sigma_Y_all_noC) )}
  
  cov_w_preds <- data.frame(all_data)
  return(cov_w_preds)
}




blockMatrixDiagonal<-function(...){  
  ## This is a function for making block diagonal matrices.
  matrixList<-list(...)
  if(is.list(matrixList[[1]])) matrixList<-matrixList[[1]]
  
  dimensions<-sapply(matrixList,FUN=function(x) dim(x)[1])
  finalDimension<-sum(dimensions)
  finalMatrix<-matrix(0,nrow=finalDimension,ncol=finalDimension)
  index<-1
  for(k in 1:length(dimensions)){
    finalMatrix[index:(index+dimensions[k]-1),index:(index+dimensions[k]-1)]<-matrixList[[k]]
    index<-index+dimensions[k]
  }
  finalMatrix
}


