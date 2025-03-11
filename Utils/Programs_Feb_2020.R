require(ggplot2)
require(nlme)
require(splines)
require(Matrix)






cmnpe <- function(data_w_out, DF_P, DF_R, B.knots, 
                  q.order = 2, cov_data=matrix(0,1,1), 
                  Pcov_data = NULL, cov_mat="CS", 
                  plots=TRUE, boot=FALSE, B=1000, seed.val = 234, TRANS=FALSE,full_sig=FALSE,zero_covs=NULL,slope=FALSE){
  
  if(!is.matrix(cov_data)){stop("Covariate data must be a matrix even if it has only 1 dimension.")}
  #### Fix all the SE variables if missing ####
  if(is.null(data_w_out$SE_pred)){
    data_w_out$SE_pred <- data_w_out$SE_var
  }
  options(warn = -1)
  if(is.null(data_w_out$SE_pred_noN)){
    data_w_out$SE_pred_noN <- rep(median(data_w_out$SE_var,na.rm = TRUE),length(data_w_out$SE_var))
  }
  if(is.null(data_w_out$SE_pred_SE)){
    data_w_out$SE_pred_SE <- rep(0,length(data_w_out$SE_var))
  }
  options(warn = 0)
  data_w_out$noN_ind <- 1*I(is.na(data_w_out$SE_pred) & is.na(data_w_out$SE_var))
  data_w_out$no_SE_w_N <- 1*I(!is.na(data_w_out$SE_pred) & is.na(data_w_out$SE_var))
  data_w_out$SE_var[is.na(data_w_out$SE_var)] <- data_w_out$SE_pred[is.na(data_w_out$SE_var)]
  data_w_out$SE_var[is.na(data_w_out$SE_var)] <- data_w_out$SE_pred_noN[is.na(data_w_out$SE_var)]
  
  #### Transform the data ####
  if(TRANS){
    data_w_out$SE_var <- data_w_out$SE_var*(data_w_out$Y)^(-1)
    data_w_out$Y <- log(data_w_out$Y)
  }
  
  all_data <- cbind(data_w_out)
  ####### Prepping data for analyses #########
  if(length(cov_data)>1){
    all_data <- cbind(all_data,cov_data)
    cnames = colnames(cov_data)
    cov_data <- as.matrix(cov_data[!is.na(data_w_out$Y),])
    colnames(cov_data) = cnames
  }
  if(!is.null(Pcov_data)){
    all_data <- cbind(all_data,Pcov_data)
    Pcov_data<- Pcov_data[!is.na(data_w_out$Y),]
    all_data <- all_data[,!duplicated(colnames(all_data))]
  }
  data_w_out <- data_w_out[!is.na(data_w_out$Y),]
  
  
  ###### Order the data by country and year (required) #######
  if(length(cov_data)>1){
    if(dim(cov_data)[2]==1){
      cnames = colnames(cov_data)
      cov_data <- as.matrix(cov_data[order(data_w_out$country,data_w_out$year)])
      colnames(cov_data) = cnames
    }
    if(dim(cov_data)[2]>1){cov_data <- cov_data[order(data_w_out$country,data_w_out$year),]}
  }
  if(!is.null(Pcov_data)){
    Pcov_data <- Pcov_data[order(data_w_out$country,data_w_out$year),]
  }
  data_w_out <- data_w_out[order(data_w_out$country,data_w_out$year),]
  all_data <- all_data[order(all_data$country,all_data$year),]
  
  ##### Specifying the covariance matrix ######
  symm = FALSE 
  diag = FALSE
  if(cov_mat=="UN"){symm = TRUE}
  if(cov_mat=="VC"){diag = TRUE}
  
  #Set some covariates equal to zero for plotting.
  if(!is.null(zero_covs)){
    all_data[,colnames(all_data) %in% zero_covs] = 0
  }
  
  ##### Run estimation. #####  
  if(!slope){
    if(is.null(Pcov_data)){  Estimation <- model_fit_noPcov(data_w_out,all_data,DF_P,DF_R,B.knots,q.order,cov_data=cov_data,symm=symm,diag=diag,plots=plots,boot=boot,B=B,seed.val=seed.val,TRANS=TRANS,full_sig=full_sig)}
    if(!is.null(Pcov_data)){ Estimation <- model_fit(data_w_out,all_data,DF_P,DF_R,B.knots,q.order,cov_data=cov_data,Pcov_data=Pcov_data,symm=symm,diag=diag,plots=plots,boot=boot,B=B,seed.val=seed.val,TRANS=TRANS,full_sig=full_sig)}
  }
  
  if(slope){
    if(is.null(Pcov_data)){  Estimation <- model_fit_noPcov_slope(data_w_out,all_data,DF_P,DF_R,B.knots,q.order,cov_data=cov_data,symm=symm,diag=diag,plots=plots,boot=boot,B=B,seed.val=seed.val,TRANS=TRANS,full_sig=full_sig)}
    if(!is.null(Pcov_data)){ Estimation <- model_fit_slope(data_w_out,all_data,DF_P,DF_R,B.knots,q.order,cov_data=cov_data,Pcov_data=Pcov_data,symm=symm,diag=diag,plots=plots,TRANS=TRANS,full_sig=full_sig)}
  }
  
  plot_data <- Estimation$plot_data
  
  if(plots & boot){
    p <- ggplot(data=plot_data, aes(x=year,y=(Y))) + geom_point() + facet_wrap(~country,scales = "free")+ ylab("Malnutrition Prevalence")
    q <- p + geom_line(data=plot_data,aes(x=year,y=(pred)),color="blue")  + geom_line(data=plot_data,aes(x=year,y=lower_CI),color="blue",linetype = 2)+ geom_line(data=plot_data,aes(x=year,y=upper_CI),color="blue",linetype = 2) #+ geom_line(data=plot_data,aes(x=year,y=lower_CI2),color="green",linetype = 2)+ geom_line(data=plot_data,aes(x=year,y=upper_CI2),color="green",linetype = 2) 
    options(warn=-1)
    print(q)
    options(warn=0)
  }
  
  return(Estimation)
}

mis_data_func <- function(data_F,data_R){
  mis_data <- NULL
  for(j in 1:length(data_F[,1])){
    if(all(data_F$Y[j] != data_R$Y)){
      mis_data <- rbind(mis_data,data_F[j,])
    }
  }
  mis_data
}


model_fit <- function(data_w_out,all_data,DF_P,DF_R,B.knots,q.order,cov_data=0,Pcov_data,symm=FALSE,diag=FALSE,plots=TRUE,boot=FALSE,B=2000,seed.val=234,TRANS=FALSE,full_sig=FALSE){
  
  ##ARGUEMENTS:
  # data_w_covs: Data set that contains:
  #     -Y: Vector of the outcome of interest {length=sum(ni)}
  #     -COV: Matrix of covariate data not including year {dimension=sum(ni) x p}.
  #     -year: vector of year {length=sum(ni)}
  #     -SE_var: variable under consideration to effect the residual variability {length=sum(ni)}
  #     -country: country level ID variable {length=sum(ni)}
  # knots: the knots for the PT-splines.
  # knots.r: additional random knots.
  # threshold: cutoff to leave a variable in the model.
  ##This function will do backwards elimination of all the elements in "cov_data" with the proposed LMM
  COV = 0
  data_w_rm_NA <- data_w_out
  P_COV <- as.matrix(Pcov_data)
  if(length(cov_data)>1){if(dim(data_w_out)[1] != dim(cov_data)[1]){stop("Outcome and Covariate data have different n.")}
    
    if(dim(cov_data)[2]==1){
      cnames = colnames(cov_data)
      COV <- as.matrix(cov_data[!is.na(cov_data)])
      P_COV <- as.matrix(Pcov_data[!is.na(cov_data),])
      colnames(COV) = cnames
    }
    if(dim(cov_data)[2]>1){
      data_w_rm_NA <- data_w_out[!c(apply((apply(cov_data,1,is.na)),2,any)),]
      P_COV <- as.matrix(Pcov_data[!c(apply((apply(cov_data,1,is.na)),2,any)),])
      COV <- as.matrix(cov_data[!c(apply((apply(cov_data,1,is.na)),2,any)),])
    }
  }
  
  Y <- data_w_rm_NA$Y
  SE_var <- (data_w_rm_NA$SE_var)
  year <- data_w_rm_NA$year
  country <- data_w_rm_NA$country
  country2 <- country
  
  bsplines <- bs(year,knots=DF_P,Boundary.knots=B.knots,intercept= F)
  B.matrix <- as.matrix(bsplines)
  b.col <- ncol(B.matrix)
  b.row <- nrow(B.matrix)
  D.matrix <- diff(diag(b.col),differences = q.order)
  P.matrix <- crossprod(D.matrix)
  Pi.matrix.svd <- svd(P.matrix)
  
  Ui.matrix <-(Pi.matrix.svd$u)[,1:(b.col-q.order)]  		# matrix of eigenvectors
  eigen.vec <-(Pi.matrix.svd$d)[1:(b.col-q.order)]   		# vector of eigenvalues
  Sigmai.inv <- diag(1/sqrt(eigen.vec))  	  				    # diagonal matrix of eigenvalues
  Zi.matrix <- B.matrix%*%Ui.matrix%*%Sigmai.inv        # instead of Z=BU
  Xi.matrix <- cbind(1,year)
  rand_spline <- as.matrix(bs(year,knots=DF_R,Boundary.knots = B.knots,intercept = F) )
  
  # The "id" for the PT-splines
  ind <- rep(1,length(country))
  ind2<- ind
  
  ctrl <- lmeControl(opt = c("nlminb"),maxIter = 1000, msMaxIter = 1000, niterEM = 500,  msMaxEval = 2000,tolerance = 1e-6)
  
  ## Fit the model
  if(length(cov_data)==1){
    if(!symm & diag){fitted_model <- lme(Y ~ 0+Xi.matrix,random=(list(ind=pdIdent(~0+Zi.matrix+P_COV[,-1]:Zi.matrix),ind2=pdIdent(~0+P_COV),country=pdSymm(~1),country2=pdDiag(~0+rand_spline))) ,weights = varSum(form=~(SE_var^2)),control = ctrl)}
    if(symm){         fitted_model <- lme(Y ~ 0+Xi.matrix,random=(list(ind=pdIdent(~0+Zi.matrix+P_COV[,-1]:Zi.matrix),ind2=pdIdent(~0+P_COV),country=pdSymm(~1),country2=pdSymm(~0+rand_spline))) ,weights = varSum(form=~(SE_var^2)),control = ctrl)}
    if(!symm & !diag){ fitted_model <- lme(Y ~ 0+Xi.matrix,random=(list(ind=pdIdent(~0+Zi.matrix+P_COV[,-1]:Zi.matrix),ind2=pdIdent(~0+P_COV),country=pdSymm(~1),country2=pdCompSymm(~0+rand_spline))) ,weights = varSum(form=~(SE_var^2)),control = ctrl)}
  }
  if(length(cov_data)>1){
    if(!symm & diag){fitted_model <- lme(Y ~ 0+Xi.matrix+COV,random=(list(ind=pdIdent(~0+Zi.matrix+P_COV[,-1]:Zi.matrix),ind2=pdIdent(~0+P_COV),country=pdSymm(~1),country2=pdDiag(~0+rand_spline))) ,weights = varSum(form=~(SE_var^2)),control = ctrl)}
    if(symm){         fitted_model <- lme(Y ~ 0+Xi.matrix+COV,random=(list(ind=pdIdent(~0+Zi.matrix+P_COV[,-1]:Zi.matrix),ind2=pdIdent(~0+P_COV),country=pdSymm(~1),country2=pdSymm(~0+rand_spline))) ,weights = varSum(form=~(SE_var^2)),control = ctrl)}
    if(!symm & !diag){ fitted_model <- lme(Y ~ 0+Xi.matrix+COV,random=(list(ind=pdIdent(~0+Zi.matrix+P_COV[,-1]:Zi.matrix),ind2=pdIdent(~0+P_COV),country=pdSymm(~1),country2=pdCompSymm(~0+rand_spline))) ,weights = varSum(form=~(SE_var^2)),control = ctrl)}
  }
  
  
  ## Exporting the results.
  data_frame <- data.frame(Y=Y,year=year,ind=ind,country=country,SE_var=SE_var,country2=country)
  result <- list(model=fitted_model,data_frame=data_frame,Xi.matrix=Xi.matrix,DF_P=DF_P,DF_R=DF_R,B.knots=B.knots,q.order=q.order,COV = COV,P_COV = P_COV)
  pred_results <- prediction(result,all_data,plots,boot)
  
  orig_data <- pred_results$newdata[!is.na(pred_results$newdata$Y),]
  sub_2_orig_ind <- (orig_data$SE_var %in% data_w_out$SE_var & orig_data$Y %in% data_w_out$Y)
  pred_newdata <- pred_results$pred_new
  df <- pred_results$df
  plot_data <- NULL
  Boot_data <- NULL
  if(!boot){
    mis_ordered <- cbind(pred_results$newdata)
    SIGMA_T_mis_ordered <- pred_results$SIGMA_T
    SIGMA_Y_mis_ordered <- pred_results$SIGMA_Y
    SIGMA_T <- SIGMA_T_mis_ordered[order(mis_ordered$country,mis_ordered$year),order(mis_ordered$country,mis_ordered$year)]
    SIGMA_Y <- SIGMA_Y_mis_ordered[order(mis_ordered$country,mis_ordered$year),order(mis_ordered$country,mis_ordered$year)]
    ordered <- mis_ordered[order(mis_ordered$country,mis_ordered$year),]
    COV_beta_b <- diag(pred_results$COV_beta_b)
    # calculate sd's and %-tiles to get CI's and PI's
    sigma_T_est <- sqrt((ordered$sigma_all))
    sigma_Y_est <- sqrt((ordered$sigma_Y_all))
    pred_new   <- ordered$pred
    pred_fixed   <- ordered$pred_all_fixed
    pred_fixpen   <- ordered$pred_all_fixpen
    ordered$resid <- (ordered$Y - ordered$pred)/sqrt((ordered$sigma_Y_all))
    
    if(!TRANS){
      lower_CI <- pred_new - 1.96*sigma_T_est
      upper_CI <- pred_new + 1.96*sigma_T_est
      lower_CI2<- pred_new - 1.96*sigma_Y_est
      upper_CI2<- pred_new + 1.96*sigma_Y_est
    }
    if(TRANS){ 
      pred_new   <- exp(pred_new)
      ordered$Y  <- exp(ordered$Y)
      
      lower_CI <- exp(log(pred_new)-1.96*sigma_T_est)
      upper_CI <- exp(log(pred_new)+1.96*sigma_T_est)
      lower_CI2<- exp(log(pred_new)-1.96*sigma_Y_est)
      upper_CI2<- exp(log(pred_new)+1.96*sigma_Y_est)
    }
    plot_data <- data.frame(Y = ordered$Y,SE_var = ordered$SE_var,country=ordered$country,year=ordered$year,pred=pred_new,lower_CI=lower_CI,upper_CI=upper_CI,lower_CI2=lower_CI2,upper_CI2=upper_CI2,sigma_T_est=sigma_T_est,sigma_Y_est=sigma_Y_est,resid = ordered$resid,pred_fixed=pred_fixed,pred_fixpen=pred_fixpen, Sex = ordered[,"Sex"])
    full_pars <- c(pred_results$b)
    fin_res <- list(result=result, plot_data = plot_data, full_pars=full_pars,COV_beta_b=COV_beta_b,df=df) #,pred_results=pred_results
    if(full_sig){
      fin_res <- list(result=result,plot_data = plot_data,full_pars=full_pars,COV_beta_b=COV_beta_b,df=df,SIGMA_T=SIGMA_T,SIGMA_Y=SIGMA_Y)
    }
  }
  
  
  return(fin_res)
}  


not_in <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

prediction <- function(result,newdata,plots=TRUE,boot){
  ##This function will predict for the original and new data.  Also,
  # it will give plots, and save plots for the original and new data.
  ##ARGUMENTS:
  # result: a fitted model and data from the model_fit function.
  # newdata: a matrix containing the new data to predict.  This must 
  #         include all of the covariates included in the final model,          
  #         along with a variable called 'year' and 'country'.
  # knots: the knots for the PT-splines.
  # knots.r: additional random knots.
  # save: set to TRUE to save plots.
  # plots: set to TRUE to see plots.
  
  # Extract the fitted model and original data set.
  fitted_model <- result$model
  data_frame <- result$data_frame
  Xi.matrix <- result$Xi.matrix
  DF_P <- result$DF_P
  DF_R <- result$DF_R
  B.knots <- result$B.knots
  q.order <- result$q.order
  COV <- result$COV
  P_COV <- result$P_COV
  
  # Make the X and Z matrices for the original fitted model
  N <- length(unique(data_frame$country))
  X_orig <- result$Xi.matrix
  if(length(result$COV)>1){X_orig <- cbind(X_orig,COV)}
  Z_orig <- make_Z(data_frame,data_frame,P_COV,DF_P,DF_R,B.knots,q.order)
  
  # Estimate EBLUPs, predicted means, and variances of the predicted means.
  pred_mat <- pred_fun(X_orig,Z_orig,Y=data_frame$Y,fitted_model,N,SE_var=data_frame$SE_var,full_b = 1)
  S <- pred_mat$S
  b <- pred_mat$b
  row.names(data_frame) <- 1:length(data_frame[,1])
  
  data_frame <- data.frame(data_frame,pred_vals=pred_mat$mu_T)
  df <- pred_mat$df
  
  
  # Plotting for the new data.
  cov_w_preds <- 0
  SIGMA_T <- NULL
  SIGMA_Y <- NULL
  if(!is.null(dim(newdata))){
    
    # Extracting only those countries that were used in the original analysis
    cov_data_used_only <- newdata[newdata$country %in% data_frame$country,]
    cov_data_not_used <- newdata[not_in(newdata$country,data_frame$country),]
    # Ording the data so the random effects align.
    efes=c(cov_data_used_only$country)
    cov_data_used_only <- cov_data_used_only[order(efes,cov_data_used_only$year),]
    # Making the X and Z matrices for the new data.
    X <- cbind(1,cov_data_used_only$year)
    if(length(COV)>1){X <- cbind(X,cov_data_used_only[,names(cov_data_used_only) %in% colnames(COV)])}
    new_P_COV <- cov_data_used_only[,names(cov_data_used_only) %in% colnames(P_COV)]
    
    Z <- make_Z(cov_data_used_only,data_frame,new_P_COV,DF_P,DF_R,B.knots,q.order)
    # Getting rid of NA's then making X and Z again
    cov_data_used_only <- cov_data_used_only[!c(apply((apply(cbind(X,new_P_COV,Z),1,is.na)),2,any)),]
    # Making the X and Z matrices for the new data.
    X <- cbind(1,cov_data_used_only$year)
    if(length(COV)>1){X <- cbind(X,cov_data_used_only[,names(cov_data_used_only) %in% colnames(COV)])}
    
    new_P_COV <- cov_data_used_only[,names(cov_data_used_only) %in% colnames(P_COV)]
    Z <- make_Z(cov_data_used_only,data_frame,new_P_COV,DF_P,DF_R,B.knots,q.order)
    # Setting the SE_var to the median value
    SE_val <- median(data_frame$SE_var)
    # Predicting the mean and variance for each new observation.
    pred_mat_new <- pred_fun(X,Z,Y=0,fitted_model,N,SE_var=SE_val,full_b = b,S=S)
    
    cov_data_all <- rbind(cov_data_used_only)
    X_all <- as.matrix(X)
    Z_all <- as.matrix(Z) 
    D_new <- pred_mat_new$D
    row.names(cov_data_all) <- NULL
    pred_all  <- pred_mat_new$mu_T
    pred_all_fixed  <- pred_mat_new$mu_F
    pred_all_fixpen  <- pred_mat_new$mu_RF
    sigma_all <- diag(pred_mat_new$sigma_T)
    sigma_Y_all <- diag(pred_mat_new$sigma_Y)
    which_pred<- rep(1,length(pred_mat_new$mu_T))
    COV_beta_b <- pred_mat_new$S
    SIGMA_T <- pred_mat_new$sigma_T
    SIGMA_Y <- pred_mat_new$sigma_Y
    
    if(length(cov_data_not_used$year)>0 & (!boot)){
      # Predicting for countries that had no data in the initial estimation procedure.
      X_NU <- cbind(1,cov_data_not_used$year)
      if(length(COV)>1){X_NU <- cbind(X_NU,cov_data_not_used[,names(cov_data_not_used) %in% colnames(COV)])}
      new_P_COV_not_used <- cov_data_not_used[,names(cov_data_not_used) %in% colnames(P_COV)]
      
      Z_NU <- make_Z(cov_data_not_used,data_frame,new_P_COV_not_used,DF_P,DF_R,B.knots,q.order)
      
      # Getting rid of NA's then making X and Z again
      cov_data_not_used <- cov_data_not_used[!c(apply((apply(cbind(X_NU,Z_NU),1,is.na)),2,any)),]
      
      if(length(cov_data_not_used$year)>0){
        # Making the X and Z matrices for the new data.
        X_NU <- cbind(1,cov_data_not_used$year)
        if(length(COV)>1){X_NU <- cbind(X_NU,cov_data_not_used[,names(cov_data_not_used) %in% colnames(COV)])}
        new_P_COV_not_used <- cov_data_not_used[,names(cov_data_not_used) %in% colnames(P_COV)]
        
        Z_NU <- make_Z(cov_data_not_used,data_frame,new_P_COV_not_used,DF_P,DF_R,B.knots,q.order)
        
        # Predicting the mean and variance for each new observation.
        fixed_col <- 1:(dim(X_NU)[2]+(length(DF_P)+3 - q.order)*(dim(new_P_COV_not_used)[2])+dim(new_P_COV_not_used)[2])
        fixed_b <- b[fixed_col]
        fixed_S <- S[fixed_col,fixed_col]
        N2      <- length(unique(cov_data_not_used$country))
        pred_mat_nob <- pred_fun_nob(X_NU,Z_NU,Y=0,fitted_model,N2,SE_val=SE_val,fixed_b = fixed_b,fixed_S=fixed_S)
        
        cov_data_all <- rbind(cov_data_used_only,cov_data_not_used)
        row.names(cov_data_all) <- NULL
        X_all <- rbind(X_all, as.matrix(X_NU))
        pred_all  <- c(pred_mat_new$mu_T,pred_mat_nob$mu_T)
        pred_all_fixed   <- c(pred_mat_new$mu_F,pred_mat_nob$mu_F)
        pred_all_fixpen  <- c(pred_mat_new$mu_RF,pred_mat_nob$mu_T)
        sigma_all <- c(sigma_all,diag(pred_mat_nob$sigma_T))
        sigma_Y_all <- c(sigma_Y_all,diag(pred_mat_nob$sigma_Y))
        SIGMA_T=blockMatrixDiagonal(SIGMA_T,pred_mat_nob$sigma_T)
        SIGMA_Y=blockMatrixDiagonal(SIGMA_Y,pred_mat_nob$sigma_Y)
      }
    }
    
    
    cov_w_preds <- data.frame(cov_data_all,pred=pred_all,sigma_all=sigma_all,sigma_Y_all=sigma_Y_all,pred_all_fixed=pred_all_fixed,pred_all_fixpen=pred_all_fixpen)
    
    
    # Plotting the newdata with the original data overlayed.  
    if(plots){
      p <- ggplot(data=cov_w_preds, aes(x=year,y=pred)) + facet_wrap(~country)
      q <- p + geom_line(data=cov_w_preds,aes(x=year,y=pred),color="blue")  + geom_point(aes(x=year,y=pred),color="blue",alpha=0.5)
      p2 <- q + geom_point(data=data_frame, aes(x=year,y=Y))
      suppressMessages(print(p2))
    }
    
  }
  # Exporting all of the prediction results.
  Ret_list <- list(newdata = cov_w_preds, pred_orig=pred_mat$mu_T, sigma_orig = pred_mat$sigma_T,b=b,D=pred_mat$D,D_new=D_new,X_new=X_all,Z_new=Z_all,X_orig=X_orig,Z_orig=Z_orig,COV_beta_b=COV_beta_b,SIGMA_T=SIGMA_T,SIGMA_Y=SIGMA_Y,df=df)
  return(Ret_list)
}


pred_fun <- function(X,Z,Y,fitted_model,N,SE_var=1,full_b=0,S=0,EBLUP=TRUE){  
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
  sigma2 <- fitted_model$sigma^2
  gamma <- c(coef(fitted_model$modelStruct$varStruct))
  
  C <- as.matrix(cbind(X,Z))
  R <- sigma2*(1 + gamma*SE_var^2)*diag(dim(Z)[1])
  df <- NULL
  if(length(full_b)==1){
    # Construct the block diagonal covariance matrix of the random effects.  
    u_var <- as.matrix(fitted_model$modelStruct$reStruct$ind)*sigma2
    t1 <- try(u2_var<- as.matrix(fitted_model$modelStruct$reStruct$ind2)*sigma2,silent=TRUE)
    b_var <- as.matrix(fitted_model$modelStruct$reStruct$country)*sigma2
    l_var <- as.matrix(fitted_model$modelStruct$reStruct$country2)*sigma2
    Dpr1 <- do.call(blockMatrixDiagonal,replicate(N, b_var, simplify=FALSE))
    Dpr2 <- do.call(blockMatrixDiagonal,replicate(N, l_var, simplify=FALSE))
    if(!is.null(attr(t1,"class"))) {D <- blockMatrixDiagonal(u_var,Dpr1,Dpr2)}
    if(is.null(attr(t1,"class"))) {D <- blockMatrixDiagonal(u_var,u2_var,Dpr1,Dpr2)}
    
    t1 <- try(R_inv <- Matrix::solve(R),silent = TRUE)
    if(!is.null(attr(t1,"class"))){
      t1 <- try(R_inv <- solve(R),silent = TRUE)
      if(!is.null(attr(t1,"class"))){
        cat("Error: Inverse failed. Consider revising model.\n")
        break
      }
    }
    t1 <- try(D_inv <- Matrix::solve(D),silent = TRUE)
    
    if(!is.null(attr(t1,"class"))){
      t1 <- try(D_inv <- solve(D),silent = TRUE)
      if(!is.null(attr(t1,"class"))){
        if(any(diag(D)<(10e-13))){cat("Warning: Some variance components are zero (with machine precision), setting them to a non-zero small value to increase computational stability. See lme output for to see which variance components are essentially zero and consider revising model. \n")}
        diag(u_var)[diag(u_var)<(10e-13)] <- 10e-13
        diag(b_var)[diag(b_var)<(10e-13)] <- 10e-13
        diag(l_var)[diag(l_var)<(10e-13)] <- 10e-13 
        Dpr1 <- do.call(blockMatrixDiagonal,replicate(N, b_var, simplify=FALSE))
        Dpr2 <- do.call(blockMatrixDiagonal,replicate(N, l_var, simplify=FALSE))
        t1 <- try(u2_var<- as.matrix(fitted_model$modelStruct$reStruct$ind2)*sigma2,silent=TRUE)
        if(!is.null(attr(t1,"class"))) {D <- blockMatrixDiagonal(u_var,Dpr1,Dpr2)}
        if(is.null(attr(t1,"class"))) {
          diag(u2_var)[diag(u2_var)<(10e-13)] <- 10e-13
          D <- blockMatrixDiagonal(u_var,u2_var,Dpr1,Dpr2)
        }
        t1 <- try(D_inv <- Matrix::solve(D),silent = TRUE)
        if(!is.null(attr(t1,"class"))){
          t1 <- try(D_inv <- solve(D),silent = TRUE)
          if(!is.null(attr(t1,"class"))){
            cat("Error: Inverse failed. Consider revising model.\n")
            break
          }
        }
      }
    }
    
    G <- t(C)%*%R_inv%*%C
    B <- blockMatrixDiagonal(0*diag(dim(X)[2]),D_inv)
    
    t1 <- try(S <- Matrix::solve(G + B),silent = TRUE)
    if(!is.null(attr(t1,"class"))){
      t1 <- try(S <- solve(G + B),silent = TRUE)
      if(!is.null(attr(t1,"class"))){
        stop("S found to be singular.\n")
      }
    }
    df <- sum(diag(S%*%G))
    
    Sigma <- Z%*%D%*%t(Z) + R
    t1 <- try(Sigma_inv <- Matrix::solve(Sigma),silent = TRUE)
    if(!is.null(attr(t1,"class"))){
      t1 <- try(Sigma_inv <- solve(Sigma),silent = TRUE)
      if(!is.null(attr(t1,"class"))){
        cat("Error: Inverse failed. Consider revising model.\n")
        break
      }
    }
    if(!is.null(attr(t1,"class"))){stop("Inverse of covariance matrix failed.")}
    
    t1 <- try(beta_hat <- Matrix::solve(t(X)%*%Sigma_inv%*%X)%*%t(X)%*%Sigma_inv%*%Y,silent = TRUE)
    if(!is.null(attr(t1,"class"))){
      t1 <- try(beta_hat <- solve(t(X)%*%Sigma_inv%*%X)%*%t(X)%*%Sigma_inv%*%Y,silent = TRUE)
      if(!is.null(attr(t1,"class"))){
        cat("Error: Inverse failed. Consider revising model.\n")
        break
      }
    }
    full_b <- c(beta_hat,D%*%t(Z)%*%Sigma_inv%*%(Y - X%*%beta_hat))
  }
  
  # Calculate the mean and covariance matrix of the predictions.
  mu_F <- as.matrix(X)%*%full_b[1:(dim(X)[2])] #fixed mean
  
  u2_var <- NULL
  t1 <- try(u2_var<- as.matrix(fitted_model$modelStruct$reStruct$ind2)*sigma2,silent=TRUE)
  u2_dim <- 0
  if(!is.null(u2_var)){u2_dim <- dim(u2_var)[2]}
  num <- (dim(X)[2]) + dim(as.matrix(fitted_model$modelStruct$reStruct$ind))[2] + u2_dim
  mu_RF <- C[,1:num]%*%full_b[1:num] #fixed with penalized coefficients
  mu_T <- C%*%full_b
  sigma_T <- C%*%S%*%t(C)
  sigma_Y <- R+sigma_T 
  
  
  # Export the mean, covariance matrix, and random effects 
  preds <- list(mu_T=mu_T,mu_F=mu_F,mu_RF=mu_RF,sigma_T=sigma_T,sigma_Y=sigma_Y,b=full_b,D=D,S=S,df=df)
  return(preds)
}


pred_fun_nob <- function(X_NU,Z_NU,Y,fitted_model,N2,SE_val=1,fixed_b=0,fixed_S=0){  
  
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
  sigma2 <- fitted_model$sigma^2
  gamma <- c(coef(fitted_model$modelStruct$varStruct))
  
  # Construct the block diagonal covariance matrix of the random effects.  
  b_var <- as.matrix(fitted_model$modelStruct$reStruct$country)*sigma2
  l_var <- as.matrix(fitted_model$modelStruct$reStruct$country2)*sigma2
  Dpr1 <- do.call(blockMatrixDiagonal,replicate(N2, b_var, simplify=FALSE))
  Dpr2 <- do.call(blockMatrixDiagonal,replicate(N2, l_var, simplify=FALSE))
  D <- blockMatrixDiagonal(Dpr1,Dpr2)
  
  Z_pen <- Z_NU[,(1:(length(fixed_b)-dim(X_NU)[2]))]
  Z_ran <- Z_NU[,-(1:(length(fixed_b)-dim(X_NU)[2]))]
  C <- as.matrix(cbind(X_NU,Z_pen))
  R <- sigma2*(1 + gamma*SE_val^2)*diag(dim(Z_NU)[1])
  
  
  mu_F    <- as.matrix(X_NU)%*%fixed_b[1:(dim(X_NU)[2])]
  mu_T    <- C%*%fixed_b
  sigma_T <- Z_ran%*%D%*%t(Z_ran)+C%*%fixed_S%*%t(C)
  sigma_Y <- R+sigma_T 
  
  # Export the mean, covariance matrix, and random effects 
  preds <- list(mu_T=mu_T,mu_F=mu_F,sigma_T=sigma_T,sigma_Y=sigma_Y)
  return(preds)
}



make_Z <- function(newdata,data_frame,P_COV,DF_P,DF_R,B.knots,q.order){
  ## This function constructs the Z design matrix for a new dataset.
  
  ## Cleaning the inputted dataset to match the original.
  cov_data_used_only <- newdata
  new_country <- cov_data_used_only$country
  country2 <- new_country
  cov_data_used_only <- cbind(cov_data_used_only,country2)
  
  ## Calculating the fixed and random spline matrices. 
  bsplines <- bs(cov_data_used_only$year,knots=DF_P,Boundary.knots=B.knots,intercept= F)
  B.matrix <- as.matrix(bsplines)
  b.col <- ncol(B.matrix)
  b.row <- nrow(B.matrix)
  D.matrix <- diff(diag(b.col),differences = q.order)
  P.matrix <- crossprod(D.matrix)
  Pi.matrix.svd <- svd(P.matrix)
  
  Ui.matrix <-(Pi.matrix.svd$u)[,1:(b.col-q.order)]  		#matrix of eigenvectors
  eigen.vec <-(Pi.matrix.svd$d)[1:(b.col-q.order)]   		#vector of eigenvalues
  Sigmai.inv <- diag(1/sqrt(eigen.vec))  	  				#diagonal matrix of eigenvalues
  all_spline_pred <- B.matrix%*%Ui.matrix%*%Sigmai.inv           		#instead of Z=BU
  rand_spline_pred <- as.matrix(bs(cov_data_used_only$year,knots=DF_R,Boundary.knots = B.knots,intercept = F) )
  
  P_COV <- as.matrix(P_COV)
  ## Creating the Z matrix.
  Z_u <- model.matrix(~0 + all_spline_pred + P_COV[,-1]:all_spline_pred)
  Z_u2<- model.matrix(~0 + P_COV)
  Z_b <- NULL
  Z_l <- NULL
  unq_cnt <- unique(new_country)
  count <- 1
  for(k in 1:length(unq_cnt)){
    mat1 <- matrix(0,dim(rand_spline_pred)[1],1)
    mat2 <- matrix(0,dim(rand_spline_pred)[1],dim(rand_spline_pred)[2])
    ni <- length(new_country[new_country==unq_cnt[k]])
    count2 <- count + ni -1
    mat1[count:count2,] <- cbind(rep(1,ni))
    mat2[count:count2,] <- rand_spline_pred[count:count2,]
    Z_b <- cbind(Z_b,mat1)
    Z_l <- cbind(Z_l,mat2)
    count <- count2+1  
  }
  Z <- cbind(Z_u,Z_u2,Z_b,Z_l)
  
  return(Z)
}

























model_fit_noPcov <- function(data_w_out,all_data,DF_P,DF_R,B.knots,q.order,cov_data=0,symm=FALSE,diag=FALSE,plots=TRUE,boot=FALSE,B=2000,seed.val=NULL,TRANS=FALSE,EBLUP = TRUE,full_sig=FALSE){
  
  ##ARGUEMENTS:
  # data_w_covs: Data set that contains:
  #     -Y: Vector of the outcome of interest {length=sum(ni)}
  #     -COV: Matrix of covariate data not including year {dimension=sum(ni) x p}.
  #     -year: vector of year {length=sum(ni)}
  #     -SE_var: variable under consideration to effect the residual variability {length=sum(ni)}
  #     -country: country level ID variable {length=sum(ni)}
  # knots: the knots for the PT-splines.
  # knots.r: additional random knots.
  # threshold: cutoff to leave a variable in the model.
  ##This function will do backwards elimination of all the elements in "cov_data" with the proposed LMM
  COV = 0
  data_w_rm_NA <- data_w_out
  if(length(cov_data)>1){
    if(dim(data_w_out)[1] != dim(cov_data)[1]){stop("Outcome and Covariate data have different n.")}
    if(dim(cov_data)[2]>1){
      cnames = colnames(cov_data)
      data_w_rm_NA <- data_w_out[!c(apply((apply(cov_data,1,is.na)),2,any)),]
      COV <- as.matrix(cov_data[!c(apply((apply(cov_data,1,is.na)),2,any)),])
      colnames(COV) = cnames
    }
    if(dim(cov_data)[2]==1){
      cnames = colnames(cov_data)
      data_w_rm_NA <- data_w_out[!is.na(cov_data),]
      COV <- as.matrix(cov_data[!is.na(cov_data)])
      colnames(COV) = cnames
    }
  }
  
  Y <- data_w_rm_NA$Y
  SE_var <- (data_w_rm_NA$SE_var)
  year <- data_w_rm_NA$year
  country <- data_w_rm_NA$country
  country2 <- country
  
  bsplines <- bs(year,knots=DF_P,Boundary.knots=B.knots,intercept= F)
  B.matrix <- as.matrix(bsplines)
  b.col <- ncol(B.matrix)
  b.row <- nrow(B.matrix)
  
  D.matrix <- diff(diag(b.col),differences = q.order)
  P.matrix <- crossprod(D.matrix)
  Pi.matrix.svd <- svd(P.matrix)
  
  Ui1.matrix <-(Pi.matrix.svd$u)[,1:(b.col-q.order)]  	# matrix of eigenvectors corresponding to non-zero eigenvalues
  eigen.vec <-(Pi.matrix.svd$d)[1:(b.col-q.order)]   		# vector of eigenvalues
  Sigmai.inv <- diag(1/sqrt(eigen.vec))  	  				    # diagonal matrix of eigenvalues
  Zi.matrix <- B.matrix%*%Ui1.matrix%*%Sigmai.inv       # instead of Z=BU
  Xi.matrix <- cbind(1,year)
  rand_spline <- as.matrix(bs(year,knots=DF_R,Boundary.knots = B.knots,intercept = F) )
  
  
  # The "id" for the PT-splines
  ind <- rep(1,length(country))
  
  ctrl <- lmeControl(opt = c("optim"),maxIter = 1000, msMaxIter = 1000, niterEM = 500,  msMaxEval = 2000,tolerance = 1e-6)
  ## Fit the model
  if(length(cov_data)==1){
    if(!symm & !diag){fitted_model <- lme(Y ~ 0+Xi.matrix,random=(list(ind=pdIdent(~0+Zi.matrix),country=pdSymm(~1),country2=pdCompSymm(~0+rand_spline))),weights = varSum(form=~(SE_var^2)),control = ctrl)}
    if(symm){         fitted_model <- lme(Y ~ 0+Xi.matrix,random=(list(ind=pdIdent(~0+Zi.matrix),country=pdSymm(~1),country2=pdSymm(~0+rand_spline))) ,weights = varSum(form=~(SE_var^2)),control = ctrl)}
    if(!symm & diag){ fitted_model <- lme(Y ~ 0+Xi.matrix,random=(list(ind=pdIdent(~0+Zi.matrix),country=pdSymm(~1),country2=pdDiag(~0+rand_spline))) ,weights = varSum(form=~(SE_var^2)),control = ctrl)}
  }
  if(length(cov_data)>1){
    if(!symm & !diag){fitted_model <- lme(Y ~ 0+Xi.matrix+COV,random=(list(ind=pdIdent(~0+Zi.matrix),country=pdSymm(~1),country2=pdCompSymm(~0+rand_spline))),weights = varSum(form=~(SE_var^2)),control = ctrl)}
    if(symm){         fitted_model <- lme(Y ~ 0+Xi.matrix+COV,random=(list(ind=pdIdent(~0+Zi.matrix),country=pdSymm(~1),country2=pdSymm(~0+rand_spline))) ,weights = varSum(form=~(SE_var^2)),control = ctrl)}
    if(!symm & diag){ fitted_model <- lme(Y ~ 0+Xi.matrix+COV,random=(list(ind=pdIdent(~0+Zi.matrix),country=pdSymm(~1),country2=pdDiag(~0+rand_spline))) ,weights = varSum(form=~(SE_var^2)),control = ctrl)}
  }
  
  ## Exporting the results.
  data_frame <- data.frame(Y=Y,year=year,ind=ind,country=country,SE_var=SE_var,country2=country)
  result <- list(model=fitted_model,data_frame=data_frame,Xi.matrix=Xi.matrix,DF_P=DF_P,DF_R=DF_R,B.knots=B.knots,q.order=q.order,COV = COV)
  
  pred_results <- prediction_noPcov(result,all_data,plots,boot)
  df <- pred_results$df
  if(!boot){
    mis_ordered <- cbind(pred_results$newdata)
    SIGMA_T_mis_ordered <- pred_results$SIGMA_T
    SIGMA_Y_mis_ordered <- pred_results$SIGMA_Y
    SIGMA_T <- SIGMA_T_mis_ordered[order(mis_ordered$country,mis_ordered$year),order(mis_ordered$country,mis_ordered$year)]
    SIGMA_Y <- SIGMA_Y_mis_ordered[order(mis_ordered$country,mis_ordered$year),order(mis_ordered$country,mis_ordered$year)]
    ordered <- mis_ordered[order(mis_ordered$country,mis_ordered$year),]
    COV_beta_b <- diag(pred_results$COV_beta_b)
    # calculate sd's and %-tiles to get CI's and PI's
    sigma_T_est <- sqrt((ordered$sigma_all))
    sigma_Y_est <- sqrt((ordered$sigma_Y_all))
    pred_new   <- ordered$pred
    pred_fixed   <- ordered$pred_all_fixed
    pred_fixpen   <- ordered$pred_all_fixpen
    ordered$resid <- (ordered$Y - ordered$pred)/sqrt((ordered$sigma_Y_all))
    
    if(!TRANS){
      lower_CI <- pred_new - 1.96*sigma_T_est
      upper_CI <- pred_new + 1.96*sigma_T_est
      lower_CI2<- pred_new - 1.96*sigma_Y_est
      upper_CI2<- pred_new + 1.96*sigma_Y_est
    }
    if(TRANS){ 
      pred_new   <- exp(pred_new)
      ordered$Y  <- exp(ordered$Y)
      
      lower_CI <- exp(log(pred_new)-1.96*sigma_T_est)
      upper_CI <- exp(log(pred_new)+1.96*sigma_T_est)
      lower_CI2<- exp(log(pred_new)-1.96*sigma_Y_est)
      upper_CI2<- exp(log(pred_new)+1.96*sigma_Y_est)
      
    }
    plot_data <- data.frame(Y = ordered$Y,SE_var = ordered$SE_var,country=ordered$country,year=ordered$year,pred=pred_new,lower_CI=lower_CI,upper_CI=upper_CI,lower_CI2=lower_CI2,upper_CI2=upper_CI2,sigma_T_est=sigma_T_est,sigma_Y_est=sigma_Y_est,resid = ordered$resid,pred_fixed=pred_fixed,pred_fixpen=pred_fixpen, Sex = ordered[,"Sex"])
    full_pars <- c(pred_results$b)
    fin_res <- list(result=result,plot_data = plot_data,full_pars=full_pars,COV_beta_b=COV_beta_b,df=df) 
    if(full_sig){
      fin_res <- list(result=result,plot_data = plot_data,full_pars=full_pars,COV_beta_b=COV_beta_b,df=df,SIGMA_T=SIGMA_T,SIGMA_Y=SIGMA_Y)
    }
  }
  
  return(fin_res)
}  


prediction_noPcov <- function(result,all_data,plots=TRUE,boot = FALSE,EBLUP = TRUE){
  ##This function will predict for the original and new data.  Also,
  # it will give plots, and save plots for the original and new data.
  ##ARGUMENTS:
  # result: a fitted model and data from the model_fit function.
  # all_data: a matrix containing the new data to predict.  This must 
  #         include all of the covariates included in the final model,          
  #         along with a variable called 'year' and 'country'.
  # knots: the knots for the PT-splines.
  # knots.r: additional random knots.
  # save: set to TRUE to save plots.
  # plots: set to TRUE to see plots.
  
  # Extract the fitted model and original data set.
  fitted_model <- result$model
  data_frame <- result$data_frame
  Xi.matrix <- result$Xi.matrix
  DF_P <- result$DF_P
  DF_R <- result$DF_R
  B.knots <- result$B.knots
  q.order <- result$q.order
  
  # Make the X and Z matrices for the original fitted model
  N <- length(unique(data_frame$country))
  X_orig <- result$Xi.matrix
  if(length(result$COV)>1){X_orig <- cbind(X_orig,result$COV)}
  Z_orig <- make_Z_noPcov(data_frame,data_frame,DF_P,DF_R,B.knots,q.order)
  
  # Estimate EBLUPs, predicted means, and variances of the predicted means.
  pred_mat <- pred_fun_noPcov(X=X_orig,Z=Z_orig,Y=data_frame$Y,fitted_model,N,SE_var=data_frame$SE_var,full_b = 1)
  df<- pred_mat$df
  S <- pred_mat$S
  b <- pred_mat$b
  row.names(data_frame) <- 1:length(data_frame[,1])
  
  data_frame <- data.frame(data_frame,pred_vals=pred_mat$mu_T)
  
  # Plotting for the new data.
  cov_w_preds <- 0
  SIGMA_T=NULL
  SIGMA_Y=NULL
  if(!is.null(dim(all_data))){
    # Extracting only those countries that were used in the original analysis
    cov_data_used_only <- all_data[all_data$country %in% data_frame$country,]
    cov_data_not_used <- all_data[not_in(all_data$country,data_frame$country),]
    # Ording the data so the random effects align.
    div=c(cov_data_used_only$country)
    cov_data_used_only <- cov_data_used_only[order(div,cov_data_used_only$year),]
    # Making the X and Z matrices for the new data.
    X <- cbind(1,cov_data_used_only$year)
    if(length(result$COV)>1){X <- cbind(X,cov_data_used_only[,names(cov_data_used_only) %in% colnames(result$COV)])}
    Z <- make_Z_noPcov(cov_data_used_only,data_frame,DF_P,DF_R,B.knots,q.order)
    
    # Getting rid of NA's then making X and Z again
    cov_data_used_only <- cov_data_used_only[!c(apply((apply(cbind(X,Z),1,is.na)),2,any)),]
    # Making the X and Z matrices for the new data.
    X <- cbind(1,cov_data_used_only$year)
    if(length(result$COV)>1){X <- cbind(X,cov_data_used_only[,names(cov_data_used_only) %in% colnames(result$COV)])}
    Z <- make_Z_noPcov(cov_data_used_only,data_frame,DF_P,DF_R,B.knots,q.order)
    
    # Setting the SE_var to the median value
    SE_val <- median(data_frame$SE_var,na.rm = TRUE)
    
    # Predicting the mean and variance for each new observation.
    pred_mat_new <- pred_fun_noPcov(X,Z,Y=0,fitted_model,N,SE_var=SE_val,full_b = b,S=S)
    
    cov_data_all <- rbind(cov_data_used_only)
    X_all <- as.matrix(X)
    Z_all <- as.matrix(Z) 
    D_new <- pred_mat_new$D
    row.names(cov_data_all) <- NULL
    pred_all  <- pred_mat_new$mu_T
    pred_all_fixed  <- pred_mat_new$mu_F
    pred_all_fixpen  <- pred_mat_new$mu_RF
    sigma_all <- diag(pred_mat_new$sigma_T)
    sigma_Y_all <- diag(pred_mat_new$sigma_Y)
    SIGMA_T=pred_mat_new$sigma_T
    SIGMA_Y=pred_mat_new$sigma_Y
    which_pred<- rep(1,length(pred_mat_new$mu_T))
    COV_beta_b <- pred_mat_new$S
    
    if(length(cov_data_not_used$year)>0 & (!boot)){
      # Predicting for countries that had no data in the initial estimation procedure.
      X_NU <- cbind(1,cov_data_not_used$year)
      if(length(result$COV)>1){X_NU <- cbind(X_NU,cov_data_not_used[,names(cov_data_not_used) %in% colnames(result$COV)])}
      Z_NU <- make_Z_noPcov(cov_data_not_used,data_frame,DF_P,DF_R,B.knots,q.order)
      # Getting rid of NA's then making X and Z again
      cov_data_not_used <- cov_data_not_used[!c(apply((apply(cbind(X_NU,Z_NU),1,is.na)),2,any)),]
      
      if(length(cov_data_not_used$year)>0){
        # Making the X and Z matrices for the new data.
        X_NU <- cbind(1,cov_data_not_used$year)
        if(length(result$COV)>1){X_NU <- cbind(X_NU,cov_data_not_used[,names(cov_data_not_used) %in% colnames(result$COV)])}
        Z_NU <- make_Z_noPcov(cov_data_not_used,data_frame,DF_P,DF_R,B.knots,q.order)
        
        # Predicting the mean and variance for each new observation.
        fixed_b <- b[1:(dim(X_NU)[2]+length(DF_P)+3 - q.order)]
        fixed_S <- S[1:(dim(X_NU)[2]+length(DF_P)+3 - q.order),1:(dim(X_NU)[2]+length(DF_P)+3 - q.order)]
        N2      <- length(unique(cov_data_not_used$country))
        pred_mat_nob <- pred_fun_noPcov_nob(X_NU,Z_NU,Y=0,fitted_model,N=N2,SE_val=SE_val,fixed_b = fixed_b,fixed_S=fixed_S)
        
        cov_data_all <- rbind(cov_data_used_only,cov_data_not_used)
        row.names(cov_data_all) <- NULL
        X_all <- rbind(X_all, as.matrix(X_NU))
        pred_all  <- c(pred_mat_new$mu_T,pred_mat_nob$mu_T)
        pred_all_fixed   <- c(pred_mat_new$mu_F,pred_mat_nob$mu_F)
        pred_all_fixpen  <- c(pred_mat_new$mu_RF,pred_mat_nob$mu_T)
        sigma_all <- c(sigma_all,diag(pred_mat_nob$sigma_T))
        sigma_Y_all <- c(sigma_Y_all,diag(pred_mat_nob$sigma_Y))
        SIGMA_T=blockMatrixDiagonal(SIGMA_T,pred_mat_nob$sigma_T)
        SIGMA_Y=blockMatrixDiagonal(SIGMA_Y,pred_mat_nob$sigma_Y)
      }
    }
    
    
    cov_w_preds <- data.frame(cov_data_all,pred=pred_all,sigma_all=sigma_all,sigma_Y_all=sigma_Y_all,pred_all_fixed=pred_all_fixed,pred_all_fixpen=pred_all_fixpen)
    
    # Plotting the newdata with the original data overlayed.  
    if(plots){
      p <- ggplot(data=cov_w_preds, aes(x=year,y=pred)) + facet_wrap(~country)
      q <- p + geom_line(data=cov_w_preds,aes(x=year,y=pred),color="blue")  + geom_point(aes(x=year,y=pred),color="blue",alpha=0.5)
      p2 <- q + geom_point(data=data_frame, aes(x=year,y=Y))
      suppressMessages(print(p2))
    }
    
  }
  # Exporting all of the prediction results.
  
  Ret_list <- list(newdata = cov_w_preds, pred_orig=pred_mat$mu_T, sigma_orig = pred_mat$sigma_T,b=b,D=pred_mat$D,D_new=D_new,X_new=X_all,Z_new=Z_all,X_orig=X_orig,Z_orig=Z_orig,COV_beta_b=COV_beta_b,SIGMA_T=SIGMA_T,SIGMA_Y=SIGMA_Y,df=df)
  return(Ret_list)
}



pred_fun_noPcov <- function(X,Z,Y,fitted_model,N,SE_var=1,full_b=0,S=0,EBLUP=TRUE){  
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
  gamma <- c(coef(fitted_model$modelStruct$varStruct))
  sigma2 <- fitted_model$sigma^2
  
  
  
  C <- as.matrix(cbind(X,Z))
  R <- sigma2*(1 + gamma*SE_var^2)*diag(dim(Z)[1])
  df<- NULL
  if(length(full_b)==1){
    u_var <- as.matrix(fitted_model$modelStruct$reStruct$ind)*sigma2
    b_var <- as.matrix(fitted_model$modelStruct$reStruct$country)*sigma2
    l_var <- as.matrix(fitted_model$modelStruct$reStruct$country2)*sigma2
    
    # Construct the block diagonal covariance matrix of the random effects.  
    Dpr1 <- do.call(blockMatrixDiagonal,replicate(N, b_var, simplify=FALSE))
    Dpr2 <- do.call(blockMatrixDiagonal,replicate(N, l_var, simplify=FALSE))
    D <- blockMatrixDiagonal(u_var,Dpr1,Dpr2)
    
    t1 <- try(R_inv <- Matrix::solve(R),silent = TRUE)
    if(!is.null(attr(t1,"class"))){
      t1 <- try(R_inv <- solve(R),silent = TRUE)
      if(!is.null(attr(t1,"class"))){
        cat("Error: Inverse failed. Consider revising model.\n")
        break
      }
    }
    t1 <- try(D_inv <- Matrix::solve(D),silent = TRUE)
    if(!is.null(attr(t1,"class"))){
      t1 <- try(D_inv <- solve(D),silent = TRUE)
      if(!is.null(attr(t1,"class"))){
        cat("Error: Inverse failed. Consider revising model.\n")
        break
      }
    }
    
    if(!is.null(attr(t1,"class"))){
      diag(u_var)[diag(u_var)<(10e-13)] <- 10e-13
      D <- blockMatrixDiagonal(u_var,Dpr1,Dpr2)
      
      t1 <- try(D_inv <- Matrix::solve(D),silent = TRUE)
      if(!is.null(attr(t1,"class"))){
        t1 <- try(D_inv <- solve(D),silent = TRUE)
        if(!is.null(attr(t1,"class"))){
          cat("Error: Inverse failed. Consider revising model.\n")
          break
        }
      }
    }
    
    G <- t(C)%*%R_inv%*%C
    B <- blockMatrixDiagonal(0*diag(dim(X)[2]),D_inv)
    
    t1 <- try(S <- Matrix::solve(G + B),silent = TRUE)
    if(!is.null(attr(t1,"class"))){
      t1 <- try(S <- solve(G + B),silent = TRUE)
      if(!is.null(attr(t1,"class"))){
        stop("S found to be singular.\n")
      }
    }
    df <- sum(diag(S%*%G))
    
    
    
    Sigma <- Z%*%D%*%t(Z) + R
    t1 <- try(Sigma_inv <- Matrix::solve(Sigma),silent = TRUE)
    if(!is.null(attr(t1,"class"))){
      t1 <- try(Sigma_inv <- solve(Sigma),silent = TRUE)
      if(!is.null(attr(t1,"class"))){
        cat("Error: Inverse failed. Consider revising model.\n")
        break
      }
    }
    if(!is.null(attr(t1,"class"))){stop("Inverse of covariance matrix failed.")}
    
    t1 <- try(beta_hat <- Matrix::solve(t(X)%*%Sigma_inv%*%X)%*%t(X)%*%Sigma_inv%*%Y,silent = TRUE)
    if(!is.null(attr(t1,"class"))){
      t1 <- try(beta_hat <- solve(t(X)%*%Sigma_inv%*%X)%*%t(X)%*%Sigma_inv%*%Y,silent = TRUE)
      if(!is.null(attr(t1,"class"))){
        cat("Error: Inverse failed. Consider revising model.\n")
        break
      }
    }
    full_b <- c(beta_hat,D%*%t(Z)%*%Sigma_inv%*%(Y - X%*%beta_hat))
    
  }
  
  mu_F <- as.matrix(X)%*%full_b[1:(dim(X)[2])]
  num <- (dim(X)[2]) + dim(as.matrix(fitted_model$modelStruct$reStruct$ind))[2]
  mu_RF <- C[,1:num]%*%full_b[1:num]
  mu_T <- C%*%full_b
  sigma_T <- C%*%S%*%t(C)
  sigma_Y <- R+sigma_T 
  
  # Export the mean, covariance matrix, and random effects 
  preds <- list(mu_T=mu_T,mu_F=mu_F,mu_RF=mu_RF,sigma_T=sigma_T,sigma_Y=sigma_Y,b=full_b,D=D,S=S,df=df)
  return(preds)
}


pred_fun_noPcov_nob <- function(X_NU,Z_NU,Y,fitted_model,N2,SE_val=1,fixed_b=0,fixed_S=0){  
  
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
  
  # Export variance components from the fitted model.
  
  # Calculate the marginal covariance, fixed effect mean, and random effecs.
  gamma <- c(coef(fitted_model$modelStruct$varStruct))
  sigma2 <- fitted_model$sigma^2
  
  Z_pen <- Z_NU[,(1:(length(fixed_b)-dim(X_NU)[2]))]
  Z_ran <- Z_NU[,-(1:(length(fixed_b)-dim(X_NU)[2]))]
  C <- as.matrix(cbind(X_NU,Z_pen))
  R <- sigma2*(1 + gamma*SE_val^2)*diag(dim(Z_NU)[1])
  b_var <- as.matrix(fitted_model$modelStruct$reStruct$country)*sigma2
  l_var <- as.matrix(fitted_model$modelStruct$reStruct$country2)*sigma2
  
  
  # Construct the block diagonal covariance matrix of the random effects.  
  Dpr1 <- do.call(blockMatrixDiagonal,replicate(N2, b_var, simplify=FALSE))
  Dpr2 <- do.call(blockMatrixDiagonal,replicate(N2, l_var, simplify=FALSE))
  D <- blockMatrixDiagonal(Dpr1,Dpr2)  
  
  
  mu_F <- as.matrix(X_NU)%*%fixed_b[1:(dim(X_NU)[2])]
  mu_T <- C%*%fixed_b
  sigma_T <- Z_ran%*%D%*%t(Z_ran)+C%*%fixed_S%*%t(C)
  sigma_Y <- R+sigma_T 
  
  # Export the mean, covariance matrix, and random effects 
  preds <- list(mu_T=mu_T,mu_F=mu_F,sigma_T=sigma_T,sigma_Y=sigma_Y)
  return(preds)
}


make_Z_noPcov <- function(all_data,data_frame,DF_P,DF_R,B.knots,q.order){
  ## This function constructs the Z design matrix for a new dataset.
  
  ## Cleaning the inputted dataset to match the original.
  cov_data_used_only <- all_data
  new_country <- cov_data_used_only$country
  country2 <- new_country
  cov_data_used_only <- cbind(cov_data_used_only,country2)
  
  ## Calculating the fixed and random spline matrices. 
  bsplines <- bs(cov_data_used_only$year,knots=DF_P,Boundary.knots=B.knots,intercept= F)
  B.matrix <- as.matrix(bsplines)
  b.col <- ncol(B.matrix)
  b.row <- nrow(B.matrix)
  D.matrix <- diff(diag(b.col),differences = q.order)
  P.matrix <- crossprod(D.matrix)
  Pi.matrix.svd <- svd(P.matrix)
  
  Ui.matrix <-(Pi.matrix.svd$u)[,1:(b.col-q.order)]  		#matrix of eigenvectors
  eigen.vec <-(Pi.matrix.svd$d)[1:(b.col-q.order)]   		#vector of eigenvalues
  Sigmai.inv <- diag(1/sqrt(eigen.vec))  	  				#diagonal matrix of eigenvalues
  all_spline_pred <- B.matrix%*%Ui.matrix%*%Sigmai.inv           		#instead of Z=BU
  rand_spline_pred <- as.matrix(bs(cov_data_used_only$year,knots=DF_R,Boundary.knots = B.knots,intercept = F) )
  
  ## Creating the Z matrix.
  Z_u <-all_spline_pred
  Z_b <- NULL
  Z_l <- NULL
  unq_cnt <- unique(new_country)
  count <- 1
  for(k in 1:length(unq_cnt)){
    mat1 <- matrix(0,dim(rand_spline_pred)[1],1)
    mat2 <- matrix(0,dim(rand_spline_pred)[1],dim(rand_spline_pred)[2])
    ni <- length(new_country[new_country==unq_cnt[k]])
    count2 <- count + ni -1
    mat1[count:count2,] <- cbind(rep(1,ni))
    mat2[count:count2,] <- rand_spline_pred[count:count2,]
    Z_b <- cbind(Z_b,mat1)
    Z_l <- cbind(Z_l,mat2)
    count <- count2+1  
  }
  Z <- cbind(Z_u,Z_b,Z_l)
  
  return(Z)
}













######## This programs takes off some of the penalized and random splines from the end and adds a random slope.  We do this to encourage linearity towards the end of the predictions.  Also, it allows for different penalizations by P_COV.

model_fit_slope <- function(data_w_out,all_data,DF_P,DF_R,B.knots,q.order,cov_data=0,Pcov_data,symm=FALSE,diag=FALSE,plots=TRUE,TRANS=FALSE,full_sig=FALSE){
  
  ##ARGUEMENTS:
  # data_w_covs: Data set that contains:
  #     -Y: Vector of the outcome of interest {length=sum(ni)}
  #     -COV: Matrix of covariate data not including year {dimension=sum(ni) x p}.
  #     -year: vector of year {length=sum(ni)}
  #     -SE_var: variable under consideration to effect the residual variability {length=sum(ni)}
  #     -country: country level ID variable {length=sum(ni)}
  # knots: the knots for the PT-splines.
  # knots.r: additional random knots.
  # threshold: cutoff to leave a variable in the model.
  ##This function will do backwards elimination of all the elements in "cov_data" with the proposed LMM
  COV = 0
  data_w_rm_NA <- data_w_out
  P_COV <- as.matrix(Pcov_data)
  if(length(cov_data)>1){if(dim(data_w_out)[1] != dim(cov_data)[1]){stop("Outcome and Covariate data have different n.")}
    
    if(dim(cov_data)[2]==1){
      cnames = colnames(cov_data)
      COV <- as.matrix(cov_data[!is.na(cov_data)])
      P_COV <- as.matrix(Pcov_data[!is.na(cov_data),])
      colnames(COV) = cnames
    }
    if(dim(cov_data)[2]>1){
      data_w_rm_NA <- data_w_out[!c(apply((apply(cov_data,1,is.na)),2,any)),]
      P_COV <- as.matrix(Pcov_data[!c(apply((apply(cov_data,1,is.na)),2,any)),])
      COV <- as.matrix(cov_data[!c(apply((apply(cov_data,1,is.na)),2,any)),])
    }
  }
  
  Y <- data_w_rm_NA$Y
  SE_var <- (data_w_rm_NA$SE_var)
  year <- data_w_rm_NA$year
  country <- data_w_rm_NA$country
  country2 <- country
  
  bsplines <- bs(year,knots=DF_P,Boundary.knots=B.knots,intercept= F)
  B.matrix <- as.matrix(bsplines)
  b.col <- ncol(B.matrix)
  b.row <- nrow(B.matrix)
  D.matrix <- diff(diag(b.col),differences = q.order)
  P.matrix <- crossprod(D.matrix)
  Pi.matrix.svd <- svd(P.matrix)
  
  Ui.matrix <-(Pi.matrix.svd$u)[,1:(b.col-q.order)]  		# matrix of eigenvectors
  eigen.vec <-(Pi.matrix.svd$d)[1:(b.col-q.order)]   		# vector of eigenvalues
  Sigmai.inv <- diag(1/sqrt(eigen.vec))  	  				    # diagonal matrix of eigenvalues
  Zi.matrix <- B.matrix%*%Ui.matrix%*%Sigmai.inv        # instead of Z=BU
  Xi.matrix <- cbind(1,year)
  rand_spline <- as.matrix(bs(year, knots=DF_R, Boundary.knots = B.knots, intercept = F) )
  rand_spline <- as.matrix(rand_spline)
  
  # The "id" for the PT-splines
  ind <- rep(1,length(country))
  ind2<- ind
  
  ctrl <- lmeControl(opt = c("nlminb"), maxIter = 1000, msMaxIter = 1000, 
                     niterEM = 500,  msMaxEval = 2000,tolerance = 1e-6)
  c_year <- year - median(year)
  all_data$c_year <- all_data$year - median(year)
  
  
  
  #int <-  P_COV[,1]*Zi.matrix
  need1 <- list(pdIdent(~0+Zi.matrix))
  for(j in 2:(dim(P_COV)[2])){
    if(j==2){
      int2 <-  P_COV[,j]*Zi.matrix
      need1 <- list(need1[[1]],pdIdent(~0+int2))}
    if(j==3){
      int3 <-  P_COV[,j]*Zi.matrix
      need1 <- list(need1[[1]],need1[[2]],pdIdent(~0+int3))}
    if(j==4){
      int4 <-  P_COV[,j]*Zi.matrix
      need1 <- list(need1[[1]],need1[[2]],need1[[3]],pdIdent(~0+int4))}
    if(j==5){
      int5 <-  P_COV[,j]*Zi.matrix
      need1 <- list(need1[[1]],need1[[2]],need1[[3]],need1[[4]],pdIdent(~0+int5))}
    if(j==6){
      int6 <-  P_COV[,j]*Zi.matrix
      need1 <- list(need1[[1]],need1[[2]],need1[[3]],need1[[4]],need1[[5]],pdIdent(~0+int6))}
    if(j==7){
      int7 <-  P_COV[,j]*Zi.matrix
      need1 <- list(need1[[1]],need1[[2]],need1[[3]],need1[[4]],need1[[5]],need1[[6]],pdIdent(~0+int7))}
    if(j==8){
      int8 <-  P_COV[,j]*Zi.matrix
      need1 <- list(need1[[1]],need1[[2]],need1[[3]],need1[[4]],need1[[5]],need1[[6]],need1[[7]],pdIdent(~0+int8))}
    if(j==9){
      int9 <-  P_COV[,j]*Zi.matrix
      need1 <- list(need1[[1]],need1[[2]],need1[[3]],need1[[4]],need1[[5]],need1[[6]],need1[[7]],need1[[8]],pdIdent(~0+int9))}
    if(j==10){
      int10 <-  P_COV[,j]*Zi.matrix
      need1 <- list(need1[[1]],need1[[2]],need1[[3]],need1[[4]],need1[[5]],need1[[6]],need1[[7]],need1[[8]],need1[[9]],pdIdent(~0+int10))}
    if(j==11){
      int11 <-  P_COV[,j]*Zi.matrix
      need1 <- list(need1[[1]],need1[[2]],need1[[3]],need1[[4]],need1[[5]],need1[[6]],need1[[7]],need1[[8]],need1[[9]],need1[[10]],pdIdent(~0+int11))}
    if(j==12){
      int12 <-  P_COV[,j]*Zi.matrix
      need1 <- list(need1[[1]],need1[[2]],need1[[3]],need1[[4]],need1[[5]],need1[[6]],need1[[7]],need1[[8]],need1[[9]],need1[[10]],need1[[11]],pdIdent(~0+int12))}
    if(j==13){
      int13 <-  P_COV[,j]*Zi.matrix
      need1 <- list(need1[[1]],need1[[2]],need1[[3]],need1[[4]],need1[[5]],need1[[6]],need1[[7]],need1[[8]],need1[[9]],need1[[10]],need1[[11]],need1[[12]],pdIdent(~0+int13))}
  }
  test1 = pdBlocked(need1)
  
  
  ## Fit the model
  if(length(cov_data)==1){
    if(!symm & diag){ fitted_model <- lme(Y ~ 0+Xi.matrix,random=(list(ind=test1,country=pdDiag(~1),country2=pdDiag(~0+rand_spline))) ,weights = varSum(form=~(SE_var^2)),control = ctrl)}
    if(symm){         fitted_model <- lme(Y ~ 0+Xi.matrix,random=(list(ind=test1,country=pdDiag(~1),country2=pdSymm(~0+rand_spline))) ,weights = varSum(form=~(SE_var^2)),control = ctrl)}
    if(!symm & !diag){fitted_model <- lme(Y ~ 0+Xi.matrix,random=(list(ind=test1,country=pdDiag(~1),country2=pdCompSymm(~0+rand_spline))) ,weights = varSum(form=~(SE_var^2)),control = ctrl)}
  }
  if(length(cov_data)>1){
    if(!symm & diag){ fitted_model <- lme(Y ~ 0+Xi.matrix+COV,random=(list(ind=test1,country=pdDiag(~1),country2=pdDiag(~0+rand_spline))) ,weights = varSum(form=~(SE_var^2)),control = ctrl)}
    if(symm){         fitted_model <- lme(Y ~ 0+Xi.matrix+COV,random=(list(ind=test1,country=pdDiag(~1),country2=pdSymm(~0+rand_spline))) ,weights = varSum(form=~(SE_var^2)),control = ctrl)}
    if(!symm & !diag){fitted_model <- lme(Y ~ 0+Xi.matrix+COV,random=(list(ind=test1,country=pdDiag(~1),country2=pdCompSymm(~0+rand_spline))) ,weights = varSum(form=~(SE_var^2)),control = ctrl)
    }
  }
  
  
  
  ## Exporting the results.
  data_frame <- data.frame(Y=Y,year=year,ind=ind,country=country,SE_var=SE_var,country2=country,c_year=c_year)
  result <- list(model=fitted_model,data_frame=data_frame,Xi.matrix=Xi.matrix,DF_P=DF_P,DF_R=DF_R,B.knots=B.knots,q.order=q.order,COV = COV,P_COV = P_COV)
  pred_results <- prediction_slope(result,all_data,plots)
  
  df <- pred_results$df
  plot_data <- NULL
  mis_ordered <- cbind(pred_results$newdata)
  SIGMA_T_mis_ordered <- pred_results$SIGMA_T
  SIGMA_Y_mis_ordered <- pred_results$SIGMA_Y
  SIGMA_T <- SIGMA_T_mis_ordered[order(mis_ordered$country,mis_ordered$year),order(mis_ordered$country,mis_ordered$year)]
  SIGMA_Y <- SIGMA_Y_mis_ordered[order(mis_ordered$country,mis_ordered$year),order(mis_ordered$country,mis_ordered$year)]
  ordered <- mis_ordered[order(mis_ordered$country,mis_ordered$year),]
  COV_beta_b <- diag(pred_results$COV_beta_b)
  # calculate sd's and %-tiles to get CI's and PI's
  sigma_T_est <- sqrt((ordered$sigma_all))
  sigma_Y_est <- sqrt((ordered$sigma_Y_all))
  pred_new   <- ordered$pred
  pred_fixed   <- ordered$pred_all_fixed
  pred_fixpen   <- ordered$pred_all_fixpen
  ordered$resid <- (ordered$Y - ordered$pred)/sqrt((ordered$sigma_Y_all))
  
  if(!TRANS){
    lower_CI <- pred_new - 1.96*sigma_T_est
    upper_CI <- pred_new + 1.96*sigma_T_est
    lower_CI2<- pred_new - 1.96*sigma_Y_est
    upper_CI2<- pred_new + 1.96*sigma_Y_est
  }
  if(TRANS){ 
    pred_new   <- exp(pred_new)
    ordered$Y  <- exp(ordered$Y)
    
    lower_CI <- exp(log(pred_new)-1.96*sigma_T_est)
    upper_CI <- exp(log(pred_new)+1.96*sigma_T_est)
    lower_CI2<- exp(log(pred_new)-1.96*sigma_Y_est)
    upper_CI2<- exp(log(pred_new)+1.96*sigma_Y_est)
    
    #upper_CI[upper_CI>1] <- 1
    #upper_CI2[upper_CI2>1] <- 1
  }
  plot_data <- data.frame(Y = ordered$Y,SE_var = ordered$SE_var,country=ordered$country,year=ordered$year,pred=pred_new,lower_CI=lower_CI,upper_CI=upper_CI,lower_CI2=lower_CI2,upper_CI2=upper_CI2,sigma_T_est=sigma_T_est,sigma_Y_est=sigma_Y_est,resid = ordered$resid,pred_fixed=pred_fixed,pred_fixpen=pred_fixpen, Sex = ordered[,"Sex"])
  full_pars <- c(pred_results$b)
  fin_res <- list(result=result,plot_data = plot_data,full_pars=full_pars,COV_beta_b=COV_beta_b,df=df) #,pred_results=pred_results
  if(full_sig){
    fin_res <- list(result=result,plot_data = plot_data,full_pars=full_pars,COV_beta_b=COV_beta_b,df=df,SIGMA_T=SIGMA_T,SIGMA_Y=SIGMA_Y)
  }
  
  
  return(fin_res)
}  


prediction_slope <- function(result,newdata,plots=TRUE){
  ##This function will predict for the original and new data.  Also,
  # it will give plots, and save plots for the original and new data.
  ##ARGUMENTS:
  # result: a fitted model and data from the model_fit function.
  # newdata: a matrix containing the new data to predict.  This must 
  #         include all of the covariates included in the final model,          
  #         along with a variable called 'year' and 'country'.
  # knots: the knots for the PT-splines.
  # knots.r: additional random knots.
  # save: set to TRUE to save plots.
  # plots: set to TRUE to see plots.
  
  # Extract the fitted model and original data set.
  fitted_model <- result$model
  data_frame <- result$data_frame
  Xi.matrix <- result$Xi.matrix
  DF_P <- result$DF_P
  DF_R <- result$DF_R
  B.knots <- result$B.knots
  q.order <- result$q.order
  COV <- result$COV
  P_COV <- result$P_COV
  
  # Make the X and Z matrices for the original fitted model
  N <- length(unique(data_frame$country))
  X_orig <- result$Xi.matrix
  if(length(result$COV)>1){X_orig <- cbind(X_orig,COV)}
  Z_orig <- make_Z_slope(data_frame,data_frame,P_COV,DF_P,DF_R,B.knots,q.order)
  
  # Estimate EBLUPs, predicted means, and variances of the predicted means.
  pred_mat <- pred_fun(X_orig,Z_orig,Y=data_frame$Y,fitted_model,N,SE_var=data_frame$SE_var,full_b = 1)
  S <- pred_mat$S
  b <- pred_mat$b
  row.names(data_frame) <- 1:length(data_frame[,1])
  data_frame <- data.frame(data_frame, 
                           pred_vals=pred_mat$mu_T, 
                           Sex = X_orig[,"Sex"])
  df <- pred_mat$df
  
  
  # Plotting for the new data.
  cov_w_preds <- 0
  SIGMA_T <- NULL
  SIGMA_Y <- NULL
  if(!is.null(dim(newdata))){
    
    # Extracting only those countries that were used in the original analysis
    cov_data_used_only <- newdata[newdata$country %in% data_frame$country,]
    cov_data_not_used <- newdata[not_in(newdata$country,data_frame$country),]
    
    efes=c(cov_data_used_only$country)
    cov_data_used_only <- cov_data_used_only[order(efes,cov_data_used_only$year),]
    # Making the X and Z matrices for the new data.
    X <- cbind(1,cov_data_used_only$year)
    if(length(COV)>1){X <- cbind(X,cov_data_used_only[,names(cov_data_used_only) %in% colnames(COV)])}
    new_P_COV <- cov_data_used_only[,names(cov_data_used_only) %in% colnames(P_COV)]
    
    Z <- make_Z_slope(cov_data_used_only,data_frame,new_P_COV,DF_P,DF_R,B.knots,q.order)
    # Getting rid of NA's then making X and Z again
    cov_data_used_only <- cov_data_used_only[!c(apply((apply(cbind(X,new_P_COV,Z),1,is.na)),2,any)),]
    # Making the X and Z matrices for the new data.
    X <- cbind(1,cov_data_used_only$year)
    if(length(COV)>1){X <- cbind(X,cov_data_used_only[,names(cov_data_used_only) %in% colnames(COV)])}
    
    new_P_COV <- cov_data_used_only[,names(cov_data_used_only) %in% colnames(P_COV)]
    Z <- make_Z_slope(cov_data_used_only,data_frame,new_P_COV,DF_P,DF_R,B.knots,q.order)
    # Setting the SE_var to the median value
    SE_val <- median(data_frame$SE_var)
    # Predicting the mean and variance for each new observation.
    pred_mat_new <- pred_fun(X,Z,Y=0,fitted_model,N,SE_var=SE_val,full_b = b,S=S)
    
    cov_data_all <- rbind(cov_data_used_only)
    X_all <- as.matrix(X)
    Z_all <- as.matrix(Z) 
    D_new <- pred_mat_new$D
    row.names(cov_data_all) <- NULL
    pred_all  <- pred_mat_new$mu_T
    pred_all_fixed  <- pred_mat_new$mu_F
    pred_all_fixpen  <- pred_mat_new$mu_RF
    sigma_all <- diag(pred_mat_new$sigma_T)
    sigma_Y_all <- diag(pred_mat_new$sigma_Y)
    which_pred<- rep(1,length(pred_mat_new$mu_T))
    COV_beta_b <- pred_mat_new$S
    SIGMA_T <- pred_mat_new$sigma_T
    SIGMA_Y <- pred_mat_new$sigma_Y
    
    if(length(cov_data_not_used$year)>0){
      # Predicting for countries that had no data in the initial estimation procedure.
      X_NU <- cbind(1,cov_data_not_used$year)
      if(length(COV)>1){X_NU <- cbind(X_NU,cov_data_not_used[,names(cov_data_not_used) %in% colnames(COV)])}
      new_P_COV_not_used <- cov_data_not_used[,names(cov_data_not_used) %in% colnames(P_COV)]
      
      Z_NU <- make_Z_slope(cov_data_not_used,data_frame,new_P_COV_not_used,DF_P,DF_R,B.knots,q.order)
      
      # Getting rid of NA's then making X and Z again
      cov_data_not_used <- cov_data_not_used[!c(apply((apply(cbind(X_NU,Z_NU),1,is.na)),2,any)),]
      if(length(cov_data_not_used$year)>0){
        # Making the X and Z matrices for the new data.
        X_NU <- cbind(1,cov_data_not_used$year)
        if(length(COV)>1){X_NU <- cbind(X_NU,cov_data_not_used[,names(cov_data_not_used) %in% colnames(COV)])}
        new_P_COV_not_used <- cov_data_not_used[,names(cov_data_not_used) %in% colnames(P_COV)]
        
        Z_NU <- make_Z_slope(cov_data_not_used,data_frame,new_P_COV_not_used,DF_P,DF_R,B.knots,q.order)
        
        # Predicting the mean and variance for each new observation.
        fixed_col <- 1:(dim(X_NU)[2]+(length(DF_P)+3 - q.order)*(dim(new_P_COV_not_used)[2]))
        
        fixed_b <- b[fixed_col]
        fixed_S <- S[fixed_col,fixed_col]
        N2      <- length(unique(cov_data_not_used$country))
        pred_mat_nob <- pred_fun_nob(X_NU,Z_NU,Y=0,fitted_model,N2,SE_val=SE_val,fixed_b = fixed_b,fixed_S=fixed_S)
        
        cov_data_all <- rbind(cov_data_used_only,cov_data_not_used)
        row.names(cov_data_all) <- NULL
        X_all <- rbind(X_all, as.matrix(X_NU))
        pred_all  <- c(pred_mat_new$mu_T,pred_mat_nob$mu_T)
        pred_all_fixed   <- c(pred_mat_new$mu_F,pred_mat_nob$mu_F)
        pred_all_fixpen  <- c(pred_mat_new$mu_RF,pred_mat_nob$mu_T)
        sigma_all <- c(sigma_all,diag(pred_mat_nob$sigma_T))
        sigma_Y_all <- c(sigma_Y_all,diag(pred_mat_nob$sigma_Y))
        SIGMA_T=blockMatrixDiagonal(SIGMA_T,pred_mat_nob$sigma_T)
        SIGMA_Y=blockMatrixDiagonal(SIGMA_Y,pred_mat_nob$sigma_Y)
      }
    }
    
    
    
    cov_w_preds <- data.frame(cov_data_all,pred=pred_all,sigma_all=sigma_all,sigma_Y_all=sigma_Y_all, pred_all_fixed=pred_all_fixed,pred_all_fixpen=pred_all_fixpen)
    
    
    # Plotting the newdata with the original data overlayed.  
    if(plots){
      p <- ggplot(data=cov_w_preds, aes(x=year,y=pred)) + facet_wrap(~country)
      q <- p + geom_line(data=cov_w_preds,aes(x=year,y=pred),color="blue")  + geom_point(aes(x=year,y=pred),color="blue",alpha=0.5)
      p2 <- q + geom_point(data=data_frame, aes(x=year,y=Y))
      suppressMessages(print(p2))
    }
    
  }
  # Exporting all of the prediction results.
  Ret_list <- list(newdata = cov_w_preds, pred_orig=pred_mat$mu_T, sigma_orig = pred_mat$sigma_T,b=b,D=pred_mat$D,D_new=D_new,X_new=X_all,Z_new=Z_all,X_orig=X_orig,Z_orig=Z_orig,COV_beta_b=COV_beta_b,SIGMA_T=SIGMA_T,SIGMA_Y=SIGMA_Y,df=df)
  return(Ret_list)
}


make_Z_slope <- function(newdata,data_frame,P_COV,DF_P,DF_R,B.knots,q.order){
  ## This function constructs the Z design matrix for a new dataset.
  
  ## Cleaning the inputted dataset to match the original.
  cov_data_used_only <- newdata
  new_country <- cov_data_used_only$country
  country2 <- new_country
  cov_data_used_only <- cbind(cov_data_used_only,country2)
  
  ## Calculating the fixed and random spline matrices. 
  bsplines <- bs(cov_data_used_only$year,knots=DF_P,Boundary.knots=B.knots,intercept= F)
  B.matrix <- as.matrix(bsplines)
  b.col <- ncol(B.matrix)
  b.row <- nrow(B.matrix)
  D.matrix <- diff(diag(b.col),differences = q.order)
  P.matrix <- crossprod(D.matrix)
  Pi.matrix.svd <- svd(P.matrix)
  
  Ui.matrix <-(Pi.matrix.svd$u)[,1:(b.col-q.order)]  		#matrix of eigenvectors
  eigen.vec <-(Pi.matrix.svd$d)[1:(b.col-q.order)]   		#vector of eigenvalues
  Sigmai.inv <- diag(1/sqrt(eigen.vec))  	  				#diagonal matrix of eigenvalues
  all_spline_pred <- B.matrix%*%Ui.matrix%*%Sigmai.inv           		#instead of Z=BU
  rand_spline_pred <- as.matrix(bs(cov_data_used_only$year,knots=DF_R,Boundary.knots = B.knots,intercept = F) )
  rand_spline_pred <- as.matrix(rand_spline_pred)
  
  
  P_COV <- as.matrix(P_COV)
  ## Creating the Z matrix.
  Z_u <- model.matrix(~0 + all_spline_pred + P_COV[,-1]:all_spline_pred)
  Z_b <- NULL
  Z_l <- NULL
  unq_cnt <- unique(new_country)
  count <- 1
  for(k in 1:length(unq_cnt)){
    mat1 <- matrix(0,length(cov_data_used_only$year),1)
    mat2 <- matrix(0,length(cov_data_used_only$year),dim(rand_spline_pred)[2])
    ni <- length(new_country[new_country==unq_cnt[k]])
    count2 <- count + ni -1
    mat1[count:count2,] <- cbind(rep(1,ni))
    mat2[count:count2,] <- rand_spline_pred[count:count2,]
    Z_b <- cbind(Z_b,mat1)
    Z_l <- cbind(Z_l,mat2)
    count <- count2+1  
  }
  Z <- cbind(Z_u,Z_b,Z_l)
  
  return(Z)
}






















model_fit_noPcov_slope <- function(data_w_out,all_data,DF_P,DF_R,B.knots,q.order,cov_data=0,symm=FALSE,diag=FALSE,plots=TRUE,boot=FALSE,B=2000,seed.val=NULL,TRANS=FALSE,EBLUP = TRUE,full_sig=FALSE){
  
  ##ARGUEMENTS:
  # data_w_covs: Data set that contains:
  #     -Y: Vector of the outcome of interest {length=sum(ni)}
  #     -COV: Matrix of covariate data not including year {dimension=sum(ni) x p}.
  #     -year: vector of year {length=sum(ni)}
  #     -SE_var: variable under consideration to effect the residual variability {length=sum(ni)}
  #     -country: country level ID variable {length=sum(ni)}
  # knots: the knots for the PT-splines.
  # knots.r: additional random knots.
  # threshold: cutoff to leave a variable in the model.
  ##This function will do backwards elimination of all the elements in "cov_data" with the proposed LMM
  COV = 0
  data_w_rm_NA <- data_w_out
  if(length(cov_data)>1){
    if(dim(data_w_out)[1] != dim(cov_data)[1]){stop("Outcome and Covariate data have different n.")}
    if(dim(cov_data)[2]>1){
      cnames = colnames(cov_data)
      data_w_rm_NA <- data_w_out[!c(apply((apply(cov_data,1,is.na)),2,any)),]
      COV <- as.matrix(cov_data[!c(apply((apply(cov_data,1,is.na)),2,any)),])
      colnames(COV) = cnames
    }
    if(dim(cov_data)[2]==1){
      cnames = colnames(cov_data)
      data_w_rm_NA <- data_w_out[!is.na(cov_data),]
      COV <- as.matrix(cov_data[!is.na(cov_data)])
      colnames(COV) = cnames
    }
  }
  
  Y <- data_w_rm_NA$Y
  SE_var <- (data_w_rm_NA$SE_var)
  year <- data_w_rm_NA$year
  country <- data_w_rm_NA$country
  country2 <- country
  
  bsplines <- bs(year,knots=DF_P,Boundary.knots=B.knots,intercept= F)
  B.matrix <- as.matrix(bsplines)
  b.col <- ncol(B.matrix)
  b.row <- nrow(B.matrix)
  
  D.matrix <- diff(diag(b.col),differences = q.order)
  P.matrix <- crossprod(D.matrix)
  Pi.matrix.svd <- svd(P.matrix)
  
  Ui1.matrix <-(Pi.matrix.svd$u)[,1:(b.col-q.order)]  	# matrix of eigenvectors corresponding to non-zero eigenvalues
  eigen.vec <-(Pi.matrix.svd$d)[1:(b.col-q.order)]   		# vector of eigenvalues
  Sigmai.inv <- diag(1/sqrt(eigen.vec))  	  				    # diagonal matrix of eigenvalues
  Zi.matrix <- B.matrix%*%Ui1.matrix%*%Sigmai.inv       # instead of Z=BU
  Xi.matrix <- cbind(1,year)
  rand_spline <- as.matrix(bs(year,knots=DF_R,Boundary.knots = B.knots,intercept = F) )
  rand_spline <- rand_spline[,-c((length(DF_R)+2):(length(DF_R)+3))]
  
  
  # The "id" for the PT-splines
  ind <- rep(1,length(country))
  c_year <- year - median(year)
  all_data$c_year <- all_data$year - median(year)
  
  ctrl <- lmeControl(opt = c("optim"),maxIter = 1000, msMaxIter = 1000, niterEM = 500,  msMaxEval = 2000,tolerance = 1e-6)
  ## Fit the model
  if(length(cov_data)==1){
    if(!symm & !diag){fitted_model <- lme(Y ~ 0+Xi.matrix,random=(list(ind=pdIdent(~0+Zi.matrix),country=pdDiag(~1+c_year),country2=pdCompSymm(~0+rand_spline))),weights = varSum(form=~(SE_var^2)),control = ctrl)}
    if(symm){         fitted_model <- lme(Y ~ 0+Xi.matrix,random=(list(ind=pdIdent(~0+Zi.matrix),country=pdDiag(~1+c_year),country2=pdSymm(~0+rand_spline))) ,weights = varSum(form=~(SE_var^2)),control = ctrl)}
    if(!symm & diag){ fitted_model <- lme(Y ~ 0+Xi.matrix,random=(list(ind=pdIdent(~0+Zi.matrix),country=pdDiag(~1+c_year),country2=pdDiag(~0+rand_spline))) ,weights = varSum(form=~(SE_var^2)),control = ctrl)}
  }
  if(length(cov_data)>1){
    if(!symm & !diag){fitted_model <- lme(Y ~ 0+Xi.matrix+COV,random=(list(ind=pdIdent(~0+Zi.matrix),country=pdDiag(~1+c_year),country2=pdCompSymm(~0+rand_spline))),weights = varSum(form=~(SE_var^2)),control = ctrl)}
    if(symm){         fitted_model <- lme(Y ~ 0+Xi.matrix+COV,random=(list(ind=pdIdent(~0+Zi.matrix),country=pdDiag(~1+c_year),country2=pdSymm(~0+rand_spline))) ,weights = varSum(form=~(SE_var^2)),control = ctrl)}
    if(!symm & diag){ fitted_model <- lme(Y ~ 0+Xi.matrix+COV,random=(list(ind=pdIdent(~0+Zi.matrix),country=pdDiag(~1+c_year),country2=pdDiag(~0+rand_spline))) ,weights = varSum(form=~(SE_var^2)),control = ctrl)}
  }
  
  ## Exporting the results.
  data_frame <- data.frame(Y=Y,year=year,ind=ind,country=country,SE_var=SE_var,country2=country,c_year=c_year)
  result <- list(model=fitted_model,data_frame=data_frame,Xi.matrix=Xi.matrix,DF_P=DF_P,DF_R=DF_R,B.knots=B.knots,q.order=q.order,COV = COV)
  
  pred_results <- prediction_noPcov_slope(result,all_data,plots,boot)
  df <- pred_results$df
  if(!boot){
    mis_ordered <- cbind(pred_results$newdata)
    SIGMA_T_mis_ordered <- pred_results$SIGMA_T
    SIGMA_Y_mis_ordered <- pred_results$SIGMA_Y
    SIGMA_T <- SIGMA_T_mis_ordered[order(mis_ordered$country,mis_ordered$year),order(mis_ordered$country,mis_ordered$year)]
    SIGMA_Y <- SIGMA_Y_mis_ordered[order(mis_ordered$country,mis_ordered$year),order(mis_ordered$country,mis_ordered$year)]
    ordered <- mis_ordered[order(mis_ordered$country,mis_ordered$year),]
    COV_beta_b <- diag(pred_results$COV_beta_b)
    # calculate sd's and %-tiles to get CI's and PI's
    sigma_T_est <- sqrt((ordered$sigma_all))
    sigma_Y_est <- sqrt((ordered$sigma_Y_all))
    pred_new   <- ordered$pred
    pred_fixed   <- ordered$pred_all_fixed
    pred_fixpen   <- ordered$pred_all_fixpen
    ordered$resid <- (ordered$Y - ordered$pred)/sqrt((ordered$sigma_Y_all))
    
    if(!TRANS){
      lower_CI <- pred_new - 1.96*sigma_T_est
      upper_CI <- pred_new + 1.96*sigma_T_est
      lower_CI2<- pred_new - 1.96*sigma_Y_est
      upper_CI2<- pred_new + 1.96*sigma_Y_est
    }
    if(TRANS){ 
      pred_new   <- exp(pred_new)
      ordered$Y  <- exp(ordered$Y)
      
      lower_CI <- exp(log(pred_new)-1.96*sigma_T_est)
      upper_CI <- exp(log(pred_new)+1.96*sigma_T_est)
      lower_CI2<- exp(log(pred_new)-1.96*sigma_Y_est)
      upper_CI2<- exp(log(pred_new)+1.96*sigma_Y_est)
      
    }
    plot_data <- data.frame(Y = ordered$Y,SE_var = ordered$SE_var,country=ordered$country,year=ordered$year,pred=pred_new,lower_CI=lower_CI,upper_CI=upper_CI,lower_CI2=lower_CI2,upper_CI2=upper_CI2,sigma_T_est=sigma_T_est,sigma_Y_est=sigma_Y_est, resid = ordered$resid,pred_fixed=pred_fixed,pred_fixpen=pred_fixpen, Sex = ordered[,"Sex"])
    full_pars <- c(pred_results$b)
    fin_res <- list(result=result,plot_data = plot_data,full_pars=full_pars,COV_beta_b=COV_beta_b,df=df) 
    if(full_sig){
      fin_res <- list(result=result,plot_data = plot_data,full_pars=full_pars,COV_beta_b=COV_beta_b,df=df,SIGMA_T=SIGMA_T,SIGMA_Y=SIGMA_Y)
    }
  }
  
  return(fin_res)
}  


prediction_noPcov_slope <- function(result,all_data,plots=TRUE,boot = FALSE,EBLUP = TRUE){
  ##This function will predict for the original and new data.  Also,
  # it will give plots, and save plots for the original and new data.
  ##ARGUMENTS:
  # result: a fitted model and data from the model_fit function.
  # all_data: a matrix containing the new data to predict.  This must 
  #         include all of the covariates included in the final model,          
  #         along with a variable called 'year' and 'country'.
  # knots: the knots for the PT-splines.
  # knots.r: additional random knots.
  # save: set to TRUE to save plots.
  # plots: set to TRUE to see plots.
  
  # Extract the fitted model and original data set.
  fitted_model <- result$model
  data_frame <- result$data_frame
  Xi.matrix <- result$Xi.matrix
  DF_P <- result$DF_P
  DF_R <- result$DF_R
  B.knots <- result$B.knots
  q.order <- result$q.order
  
  # Make the X and Z matrices for the original fitted model
  N <- length(unique(data_frame$country))
  X_orig <- result$Xi.matrix
  if(length(result$COV)>1){X_orig <- cbind(X_orig,result$COV)}
  Z_orig <- make_Z_noPcov_slope(data_frame,data_frame,DF_P,DF_R,B.knots,q.order)
  
  # Estimate EBLUPs, predicted means, and variances of the predicted means.
  pred_mat <- pred_fun_noPcov(X=X_orig,Z=Z_orig,Y=data_frame$Y,fitted_model,N,SE_var=data_frame$SE_var,full_b = 1)
  df<- pred_mat$df
  S <- pred_mat$S
  b <- pred_mat$b
  row.names(data_frame) <- 1:length(data_frame[,1])
  data_frame <- data.frame(data_frame,pred_vals=pred_mat$mu_T)
  
  # Plotting for the new data.
  cov_w_preds <- 0
  SIGMA_T=NULL
  SIGMA_Y=NULL
  if(!is.null(dim(all_data))){
    # Extracting only those countries that were used in the original analysis
    cov_data_used_only <- all_data[all_data$country %in% data_frame$country,]
    cov_data_not_used <- all_data[not_in(all_data$country,data_frame$country),]
    # Ording the data so the random effects align.
    div=c(cov_data_used_only$country)
    cov_data_used_only <- cov_data_used_only[order(div,cov_data_used_only$year),]
    # Making the X and Z matrices for the new data.
    X <- cbind(1,cov_data_used_only$year)
    if(length(result$COV)>1){X <- cbind(X,cov_data_used_only[,names(cov_data_used_only) %in% colnames(result$COV)])}
    Z <- make_Z_noPcov_slope(cov_data_used_only,data_frame,DF_P,DF_R,B.knots,q.order)
    
    # Getting rid of NA's then making X and Z again
    cov_data_used_only <- cov_data_used_only[!c(apply((apply(cbind(X,Z),1,is.na)),2,any)),]
    # Making the X and Z matrices for the new data.
    X <- cbind(1,cov_data_used_only$year)
    if(length(result$COV)>1){X <- cbind(X,cov_data_used_only[,names(cov_data_used_only) %in% colnames(result$COV)])}
    Z <- make_Z_noPcov_slope(cov_data_used_only,data_frame,DF_P,DF_R,B.knots,q.order)
    
    # Setting the SE_var to the median value
    SE_val <- median(data_frame$SE_var,na.rm = TRUE)
    
    # Predicting the mean and variance for each new observation.
    pred_mat_new <- pred_fun_noPcov(X,Z,Y=0,fitted_model,N,SE_var=SE_val,full_b = b,S=S)
    
    cov_data_all <- rbind(cov_data_used_only)
    X_all <- as.matrix(X)
    Z_all <- as.matrix(Z) 
    D_new <- pred_mat_new$D
    row.names(cov_data_all) <- NULL
    pred_all  <- pred_mat_new$mu_T
    pred_all_fixed  <- pred_mat_new$mu_F
    pred_all_fixpen  <- pred_mat_new$mu_RF
    sigma_all <- diag(pred_mat_new$sigma_T)
    sigma_Y_all <- diag(pred_mat_new$sigma_Y)
    SIGMA_T=pred_mat_new$sigma_T
    SIGMA_Y=pred_mat_new$sigma_Y
    which_pred<- rep(1,length(pred_mat_new$mu_T))
    COV_beta_b <- pred_mat_new$S
    
    if(length(cov_data_not_used$year)>0){
      # Predicting for countries that had no data in the initial estimation procedure.
      X_NU <- cbind(1,cov_data_not_used$year)
      if(length(result$COV)>1){X_NU <- cbind(X_NU,cov_data_not_used[,names(cov_data_not_used) %in% colnames(result$COV)])}
      Z_NU <- make_Z_noPcov_slope(cov_data_not_used,data_frame,DF_P,DF_R,B.knots,q.order)
      # Getting rid of NA's then making X and Z again
      cov_data_not_used <- cov_data_not_used[!c(apply((apply(cbind(X_NU,Z_NU),1,is.na)),2,any)),]
      
      if(length(cov_data_not_used$year)>0){
        # Making the X and Z matrices for the new data.
        X_NU <- cbind(1,cov_data_not_used$year)
        if(length(result$COV)>1){X_NU <- cbind(X_NU,cov_data_not_used[,names(cov_data_not_used) %in% colnames(result$COV)])}
        Z_NU <- make_Z_noPcov_slope(cov_data_not_used,data_frame,DF_P,DF_R,B.knots,q.order)
        
        # Predicting the mean and variance for each new observation.
        fixed_col <- 1:(dim(X_NU)[2]+length(DF_P)+3 - q.order)
        
        fixed_b <- b[fixed_col]
        fixed_S <- S[fixed_col,fixed_col]
        N2      <- length(unique(cov_data_not_used$country))
        pred_mat_nob <- pred_fun_noPcov_nob(X_NU,Z_NU,Y=0,fitted_model,N=N2,SE_val=SE_val,fixed_b = fixed_b,fixed_S=fixed_S)
        
        cov_data_all <- rbind(cov_data_used_only,cov_data_not_used)
        row.names(cov_data_all) <- NULL
        X_all <- rbind(X_all, as.matrix(X_NU))
        pred_all  <- c(pred_mat_new$mu_T,pred_mat_nob$mu_T)
        pred_all_fixed   <- c(pred_mat_new$mu_F,pred_mat_nob$mu_F)
        pred_all_fixpen  <- c(pred_mat_new$mu_RF,pred_mat_nob$mu_T)
        sigma_all <- c(sigma_all,diag(pred_mat_nob$sigma_T))
        sigma_Y_all <- c(sigma_Y_all,diag(pred_mat_nob$sigma_Y))
        SIGMA_T=blockMatrixDiagonal(SIGMA_T,pred_mat_nob$sigma_T)
        SIGMA_Y=blockMatrixDiagonal(SIGMA_Y,pred_mat_nob$sigma_Y)
      }
    }
    
    cov_w_preds <- data.frame(cov_data_all,pred=pred_all,sigma_all=sigma_all,sigma_Y_all=sigma_Y_all,pred_all_fixed=pred_all_fixed,pred_all_fixpen=pred_all_fixpen)
    
    # Plotting the newdata with the original data overlayed.  
    if(plots){
      p <- ggplot(data=cov_w_preds, aes(x=year,y=pred)) + facet_wrap(~country)
      q <- p + geom_line(data=cov_w_preds,aes(x=year,y=pred),color="blue")  + geom_point(aes(x=year,y=pred),color="blue",alpha=0.5)
      p2 <- q + geom_point(data=data_frame, aes(x=year,y=Y))
      suppressMessages(print(p2))
    }
    
  }
  # Exporting all of the prediction results.
  Ret_list <- list(newdata = cov_w_preds, pred_orig=pred_mat$mu_T, sigma_orig = pred_mat$sigma_T,b=b,D=pred_mat$D,D_new=D_new,X_new=X_all,Z_new=Z_all,X_orig=X_orig,Z_orig=Z_orig,COV_beta_b=COV_beta_b,SIGMA_T=SIGMA_T,SIGMA_Y=SIGMA_Y,df=df)
  return(Ret_list)
}


make_Z_noPcov_slope <- function(all_data,data_frame,DF_P,DF_R,B.knots,q.order){
  ## This function constructs the Z design matrix for a new dataset.
  
  ## Cleaning the inputted dataset to match the original.
  cov_data_used_only <- all_data
  new_country <- cov_data_used_only$country
  country2 <- new_country
  cov_data_used_only <- cbind(cov_data_used_only,country2)
  
  ## Calculating the fixed and random spline matrices. 
  bsplines <- bs(cov_data_used_only$year,knots=DF_P,Boundary.knots=B.knots,intercept= F)
  B.matrix <- as.matrix(bsplines)
  b.col <- ncol(B.matrix)
  b.row <- nrow(B.matrix)
  D.matrix <- diff(diag(b.col),differences = q.order)
  P.matrix <- crossprod(D.matrix)
  Pi.matrix.svd <- svd(P.matrix)
  
  Ui.matrix <-(Pi.matrix.svd$u)[,1:(b.col-q.order)]  		#matrix of eigenvectors
  eigen.vec <-(Pi.matrix.svd$d)[1:(b.col-q.order)]   		#vector of eigenvalues
  Sigmai.inv <- diag(1/sqrt(eigen.vec))  	  				#diagonal matrix of eigenvalues
  all_spline_pred <- B.matrix%*%Ui.matrix%*%Sigmai.inv           		#instead of Z=BU
  rand_spline_pred <- as.matrix(bs(cov_data_used_only$year,knots=DF_R,Boundary.knots = B.knots,intercept = F) )
  rand_spline_dim <- dim(rand_spline_pred)
  rand_spline_dim[2] <- rand_spline_dim[2] - ncol(rand_spline_pred[,c((length(DF_R)+2):(length(DF_R)+3))])
  rand_spline_pred <- matrix(
    rand_spline_pred[,-c((length(DF_R)+2):(length(DF_R)+3))], 
    nrow = rand_spline_dim[1],
    ncol = rand_spline_dim[2])
  
  ## Creating the Z matrix.
  Z_u <-all_spline_pred
  Z_b <- NULL
  Z_l <- NULL
  unq_cnt <- unique(new_country)
  count <- 1
  for(k in 1:length(unq_cnt)){
    mat1 <- matrix(0,rand_spline_dim[1],2)
    mat2 <- matrix(0,rand_spline_dim[1],rand_spline_dim[2])
    ni <- length(new_country[new_country==unq_cnt[k]])
    count2 <- count + ni -1
    mat1[count:count2,] <- cbind(rep(1,ni),cov_data_used_only$c_year[count:count2])
    mat2[count:count2,] <- rand_spline_pred[count:count2,]
    Z_b <- cbind(Z_b,mat1)
    Z_l <- cbind(Z_l,mat2)
    count <- count2+1  
  }
  Z <- cbind(Z_u,Z_b,Z_l)
  
  return(Z)
}





varSum <-
  function(value = numeric(0), form = ~ fitted(.), fixed = NULL)
  {
    value <- unlist(value)		# may be given as a list
    fixed <- attr(value, "fixed") <- unlist(fixed)
    attr(value, "formula") <- form <- asOneSidedFormula(form)
    if (length(all.vars(getCovariateFormula(form))) == 0) {
      stop("\"form\" must have a covariate")
    }
    if (!is.null(getGroupsFormula(form))) {
      if (is.null(grpNames <- names(value)) && (length(value) > 1)) {
        stop("Initial values must have groups names in varPower")
      }
      if (!is.null(fixed)) {
        if (is.null(names(fixed))) {
          stop("Fixed parameters must have groups names in varPower")
        }
      }
      attr(value, "groupNames") <- c(grpNames, names(fixed))
    } else {
      attr(value, "whichFix") <- !is.null(fixed)
    }
    class(value) <- c("varSum", "varFunc")
    value
  }

###*# Methods for standard generics

coef.varSum <-
  function(object, unconstrained = TRUE, allCoef = FALSE, ...)
  {
    if (((length(object) == 0) &&
         (!allCoef || is.null(attr(object, "fixed")))) ||
        is.null( wPar <- attr(object, "whichFix"))) {
      return(numeric(0))
    }
    val <- double(length(wPar))
    if (any(wPar)) {
      val[wPar] <- attr(object, "fixed")
    }
    if (any(!wPar)) {
      val[!wPar] <- as.vector(object)
    }
    if (!is.null(getGroupsFormula(object))) {
      ##different values per group
      names(val) <- attr(object, "groupNames")
    } else {
      names(val) <- "covarCoef"
    }
    if (!allCoef) {
      val <- val[!wPar]
    }
    val
  }

"coef<-.varSum" <-
  function(object, ..., value)
  {
    if (length(object) > 0) {		# varying parameters
      value <- as.numeric(value)
      if (length(value) != length(object)) {
        stop(paste("Cannot change the length of the varStruct",
                   "parameter after initialization"))
      }
      object[] <- value
      aux <- coef(object, FALSE, allCoef = TRUE)
      if (!is.null(grps <- getGroups(object))) {
        aux <- aux[grps]
      }
      attr(object, "logLik") <-
        sum(log(attr(object, "weights") <- 1/sqrt(1+exp(aux)* getCovariate(object))))
    } else {
      stop(paste("Cannot change coefficients before initialization or",
                 "when all parameters are fixed"))
    }
    object
  }

Initialize.varSum <-
  function(object, data, ...)
  {
    form <- formula(object)
    if (all(!is.na(match(all.vars(getCovariateFormula(form)), names(data))))) {
      ## can evaluate covariate on data
      attr(object, "needUpdate") <- FALSE
      attr(object, "covariate") <- getCovariate(data, form)
    } else {
      attr(object, "needUpdate") <- TRUE
    }
    if (!is.null(grpForm <- getGroupsFormula(form))) {
      strat <- as.character(getGroups(data, form,
                                      level = length(splitFormula(grpForm, sep = "*")),
                                      sep = "*"))
      uStrat <- unique(strat)
      if (length(uStrat) > 1) {		# multi-groups
        attr(object, "groups") <- strat
        if (!is.null(attr(object, "fixed"))) {
          fixNames <- names(attr(object, "fixed"))
          if (is.null(fixNames)) {
            stop("Fixed parameters must have group names")
          }
          if (any(is.na(match(fixNames, uStrat)))) {
            stop("Mismatch between group names and fixed values names")
          }
        } else {
          fixNames <- NULL
        }
        uStratVar <- uStrat[is.na(match(uStrat, fixNames))]
        nStratVar <- length(uStratVar)
        attr(object, "whichFix") <- !is.na(match(uStrat, fixNames))
        if (nStratVar > 0) {
          if (length(object) <= 1) {
            ## repeat for all groups
            names(object) <- NULL
            oldAttr <- attributes(object)
            if (length(object) > 0) {
              object <- rep(as.vector(object), nStratVar)
            } else {
              object <- rep(0, nStratVar)
            }
            attributes(object) <- oldAttr
            attr(object, "groupNames") <- uStrat
            names(object) <- uStratVar
          } else {
            if (length(as.vector(object)) != nStratVar) {
              stop(paste("Initial value for \"varSum\" should be of length",
                         nStratVar))
            }
            stN <- attr(object, "groupNames") #must have names
            if ((length(stN) != length(uStrat)) ||
                any(sort(stN) != sort(uStrat))) {
              stop("Nonexistent groups names for initial values in varSum")
            }
          }
        } else {
          if (all(attr(object, "fixed") == 0)) {
            ## equal variances structure
            return(Initialize(varIdent(), data))
          } else {
            oldAttr <- attributes(object)
            object <- numeric(0)
            attributes(object) <- oldAttr
            attr(object, "groupNames") <- uStrat
          }
        }
      } else {                            # single stratum
        attr(object, "formula") <- getCovariateFormula(formula(object))
        attr(object, "whichFix") <- !is.null(attr(object, "fixed"))
      }
    }
    if (is.null(getGroupsFormula(object))) {
      ## single stratum
      if (attr(object, "whichFix")) {
        if (!attr(object, "fixed")) {
          ## equal variances structure
          return(Initialize(varIdent(), data))
        } else {
          oldAttr <- attributes(object)
          object <- numeric(0)
          attributes(object) <- oldAttr
        }
      } else {
        len <- length(as.vector(object))
        if (len == 0) {			# uninitialized
          oldAttr <- attributes(object)
          object <- 0
          attributes(object) <- oldAttr
        } else if (len > 1) {
          stop("Initial value for \"varSum\" should be of length 1.")
        }
      }
    }
    if (!is.null(covar <- getCovariate(object))) {
      natPar <- coef(object, allCoef = TRUE)
      if (!is.null(grps <- getGroups(object))) {
        natPar <- natPar[grps]
      }
      attr(object, "logLik") <-
        sum(log(attr(object, "weights") <-  1/sqrt(1+exp(natPar)* getCovariate(object))))
      object
    } else {
      NextMethod()
    }
  }


summary.varSum <-
  function(object, structName = "Sum of one and constant times variance covariate", ...)
  {
    if (!is.null(getGroupsFormula(object))) {
      structName <- paste(structName, " different strata", sep = ",")
    }
    summary.varFunc(object, structName)
  }

update.varSum <-
  function(object, data, ...)
  {
    val <- NextMethod()
    if (length(val) == 0) {		# chance to update weights
      aux <- coef(val, allCoef = TRUE)
      if (!is.null(grps <- getGroups(val))) {
        aux <- aux[grps]
      }
      attr(val, "logLik") <-
        sum(log(attr(val, "weights") <- 1/sqrt(1+exp(aux)* getCovariate(object))))
    }
    val
  }

summary.varFunc <-
  function(object, structName = class(object)[1], ...)
  {
    attr(object, "structName") <- structName
    attr(object, "oClass") <- class(object)
    class(object) <- "summary.varFunc"
    object
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



output_function <- function(plot_data, Estimation, data_one, marker, date,path){
  
  cov_data <- Estimation$cov_data
  zero_covs <- Estimation$zero_covs
  gamma <- Estimation$gamma
  sigma2 <- Estimation$sigma2
  smart_cov <- Estimation$smart_cov
  
  plot_data <- plot_data %>% 
    group_by(country, year, Sex) %>%
    filter(row_number()==1) %>% 
    ungroup() %>% 
    mutate(
      Sex = case_when(
        Sex ==  0 ~ "Both",
        Sex ==  1 ~ "Female",
        Sex == -1 ~ "Male"
      )
    ) %>% 
    arrange(country, year, Sex) %>% 
    select(Y:Sex, SE_mean_pred, SE_pred) 
  
  
  ########### Writing summaries to external csv files. ###########
  #### Make sure data_one from the program is loaded in ######
  
  ### Outputting data to .csv file.
  PP_plot_data <- data_one %>% 
    mutate(
      Region = factor(Region), 
      "SMART" = case_when(
        ShortSource=="SMART" ~ 1,
        TRUE ~ 0
      ),
      "Surveillance" = case_when(
        ShortSource=="Surveillance" ~ 1,
        TRUE ~ 0
      )
    ) %>% 
    select(-"SE_pred") %>% 
    arrange(country, year, Sex) %>% 
    left_join(plot_data, relationship = "many-to-one", by = c("country", "year", "Sex"))  %>% 
    mutate(
      all_SE_var = SE_val*((1/Point.Estimate.Imp) + 1/(1-Point.Estimate.Imp)), 
      all_Y = log(Point.Estimate.Imp/(1-Point.Estimate.Imp))
    ) %>% 
    mutate(
      SE_var = case_when(
        SE_var == 0 ~ all_SE_var,
        TRUE ~ SE_var
      ),
      Y_NS = Y,
      Y = case_when(
        is.na(Y) ~ all_Y,
        TRUE ~ Y
      )
    ) 
  if(!is.null(zero_covs)){ 
    PP_plot_data <- PP_plot_data %>%  
      mutate(
        adj_pred_ns = Y_NS - smart_cov*SMART,
        adj_pred = Y - smart_cov*SMART
      ) %>% 
      mutate(
        adj_pred = exp(adj_pred)/(1+exp(adj_pred)),
        adj_pred_ns = exp(adj_pred_ns)/(1+exp(adj_pred_ns))
      )
  }
  
  PP_plot_data <- PP_plot_data %>% 
    mutate(
      Point.Estimate = exp(Y)/(1+exp(Y)),
      Prediction = exp(pred)/(1+exp(pred)),
      pred_fixed = exp(pred_fixed)/(1+exp(pred_fixed)),
      pred_fixpen = exp(pred_fixpen)/(1+exp(pred_fixpen)),
      SE_var = sqrt(sigma2*(1 + gamma*SE_var^2))
    ) %>% 
    mutate(
      Std.Err = Point.Estimate*(1-Point.Estimate)*(SE_var), 
      SE_mean_pred = Prediction*(1-Prediction)*(sigma_T_est),
      SE_pred = Prediction*(1-Prediction)*(sigma_Y_est),
      lower_CI = exp(lower_CI)/(1+exp(lower_CI)),
      upper_CI = exp(upper_CI)/(1+exp(upper_CI)),
      lower_PI = exp(lower_CI2)/(1+exp(lower_CI2)),
      upper_PI = exp(upper_CI2)/(1+exp(upper_CI2))
    ) %>% 
    select(
      -c("Y","pred","lower_CI2","upper_CI2","sigma_T_est","sigma_Y_est","SE_var","Y_NS", "all_SE_var","all_Y")
    ) %>% 
    mutate(check = case_when(
      is.na(Point.Estimate) ~ 0,
      is.na(Point.Estimate.Imp) ~ 0,
      abs(Point.Estimate.Imp -Point.Estimate) < 0.000001 ~ 0,
      TRUE ~ 1)
    ) %>% 
    filter(check == 0) %>% 
    mutate(
      Outlier=1*I(abs(resid)>3),
      SE_Imput_Ind = 1*I(Standard.Error==SE_val), 
      Standard.Error = Standard.Error/100) %>% 
    rename( 
      SSE = Standard.Error, 
      SSE_imp = SE_val,
      Total_SSE = Std.Err
    )   %>% 
    select(c(ISO.code, country, year, Sex, Point.Estimate.Orig, Point.Estimate.Imp, 
             Prediction, lower_CI, upper_CI, lower_PI, upper_PI, SSE, SSE_imp, 
             Total_SSE, Region, N, LowerLimit, UpperLimit, ShortSource),everything()) %>% 
    select(-c(check)) %>% 
    mutate(ShortSource = as.character(ShortSource)) 
  
  both_data <- PP_plot_data %>% 
    select(country, year, Sex, Prediction) %>% 
    filter(Sex == "Both") %>% 
    rename(Prediction_overall = Prediction) %>% 
    group_by(country, year, Sex) %>%
    filter(row_number()==1) %>% 
    ungroup()%>% 
    select(-Sex) 
  male_data <- PP_plot_data %>% 
    select(country, year, Sex, Prediction) %>% 
    filter(Sex == "Male") %>% 
    rename(Prediction_male = Prediction) %>% 
    group_by(country, year, Sex) %>%
    filter(row_number()==1) %>% 
    ungroup()%>% 
    select(-Sex) 
  female_data <- PP_plot_data %>% 
    select(country, year, Sex, Prediction) %>% 
    filter(Sex == "Female") %>% 
    rename(Prediction_female = Prediction) %>% 
    group_by(country, year, Sex) %>%
    filter(row_number()==1) %>% 
    ungroup() %>% 
    select(-Sex) 
  
  # wide_data <- both_data %>% 
  #   left_join(male_data) %>% 
  #   left_join(female_data) %>% 
  #   mutate(diff = Prediction_overall - 
  #            (Prediction_male + Prediction_female)/2) %>% 
  #   select(country:year, diff)
  
  PP_plot_data <- PP_plot_data %>% 
    # left_join(wide_data) %>% 
    # mutate(
    #   diff = case_when(
    #     Sex ==  "Both" ~ 0,
    #     TRUE ~ diff
    #   )
    # ) %>% 
    mutate(
      Sex = case_when(
        Sex ==  "Both" ~ "Overall",
        TRUE ~ Sex
      )
    ) %>%
    # mutate(
    #   Point.Estimate = Point.Estimate + diff,
    #   Prediction = Prediction + diff,
    #   pred_fixed = pred_fixed + diff,
    #   pred_fixpen = pred_fixpen + diff,
    #   lower_CI = lower_CI + diff,
    #   upper_CI = upper_CI + diff,
    #   lower_PI = lower_PI + diff,
    #   upper_PI = upper_PI + diff
    # ) %>% 
    # mutate(
    #   resid = (Point.Estimate - Prediction)/SE_pred
    # ) %>% 
    select(-c("Point.Estimate","SMART","Surveillance"))
  # 
  # min_test <- PP_plot_data %>% select(Prediction,
  #                                     pred_fixed,
  #                                     pred_fixpen,
  #                                     lower_CI,
  #                                     upper_CI,
  #                                     lower_PI,
  #                                     upper_PI)
  # if(min(min_test,na.rm = TRUE) < 0){stop("Minimum test failed")}
  
  write.csv(PP_plot_data,paste0(path,marker," results ",date,".csv",sep=""),row.names = FALSE)
  
}




