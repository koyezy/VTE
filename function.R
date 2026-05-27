######### THIS IS THE LATEST VERSION OF FUNCTION R SCRIPT #########

# RNDR
# functions to add p-value to table 1

rndr <- function(x, ...) {
  y <- render.default(x, ...)
  if (is.factor(x) & length(unique(x))==2) y[3] else y
}


pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times = sapply(x, length)))
  
  if (is.numeric(y)) {
    # Nonparametric alternative to 2-sample t-test
    p <- wilcox.test(y ~ g)$p.value
  } else {
    # Categorical variables: chi-squared test
    p <- chisq.test(table(y, g))$p.value
  }
  
  # Format p-value
  c(sub("<", "&lt;", format.pval(p, digits = 3, eps = 0.001)))
}


rndr.strat <- function(label, n, ...) {
  ifelse(n==0, label, render.strat.default(label, n, ...))
}

# EXPIT, LOGIT

expit = function(x){
  1/(1+exp(-x))
}

logit = function(x){
  log(x/(1-x))
}

# BIGSUMMARY

BigSummary <- function (data, lev = NULL, model = NULL) {
  pr_auc <- try(MLmetrics::PRAUC(data[, lev[2]],
                                 ifelse(data$obs == lev[2], 1, 0)),
                silent = TRUE)
  brscore <- try(mean((data[, lev[2]] - ifelse(data$obs == lev[2], 1, 0)) ^ 2),
                 silent = TRUE)
  rocObject <- try(pROC::roc(ifelse(data$obs == lev[2], 1, 0), data[, lev[2]],
                             direction = "<", quiet = TRUE), silent = TRUE)
  if (inherits(pr_auc, "try-error")) pr_auc <- NA
  if (inherits(brscore, "try-error")) brscore <- NA
  rocAUC <- if (inherits(rocObject, "try-error")) {
    NA
  } else {
    rocObject$auc
  }
  tmp <- unlist(e1071::classAgreement(table(data$obs,
                                            data$pred)))[c("diag", "kappa")]
  out <- c(Acc = tmp[[1]],
           Kappa = tmp[[2]],
           AUCROC = rocAUC,
           AUCPR = pr_auc,
           Brier = brscore,
           Recall = caret:::recall.default(data = data$pred,
                                           reference = data$obs,
                                           relevant = lev[2]))
  out
}

# DATA AUGMENTATION
library(EMgaussian)

data_augment_em <- function(data, max_iter = 100, tol = 1e-6, verbose = TRUE) {
  # Load necessary library
  if (!requireNamespace("MASS", quietly = TRUE)) stop("Install 'MASS' package.")
  library(MASS)
  
  # Identify columns
  vars <- colnames(data)
  n <- nrow(data)
  
  # Identify variable with missing values
  miss_col <- vars[colSums(is.na(data)) > 0]
  if (length(miss_col) != 1) stop("Only one variable with missing values is allowed.")
  #y <- log(data[[miss_col]]) # log normal -> normal
  
  # Transform exponential distribution to normal distribution
  u = 1 - exp(-data[[miss_col]]/80)
  y=qnorm(u)
  
  # Split data into observed and missing parts
  obs_idx <- !is.na(y)
  mis_idx <- is.na(y)
  y_obs <- y[obs_idx]
  X_obs <- data[obs_idx, vars != miss_col]
  X_mis <- as.matrix(data[mis_idx, vars != miss_col])
  
  # Initialize using regression on observed data
  obs_data = data.frame(y_obs = y_obs, X_obs)
  fit <- lm(y_obs ~ ., data=obs_data)
  beta <- coef(fit)
  sigma2 <- var(resid(fit))
  
  iter <- 1
  diff <- Inf
  
  while (iter <= max_iter && diff > tol) {
    # -------------------
    # E-step: Impute missing Y
    # -------------------
    y_mis_hat <- cbind(1, X_mis) %*% beta
    
    # Optional stochastic DA step: add noise from posterior
    y_mis <- rnorm(length(y_mis_hat), mean = y_mis_hat, sd = sqrt(sigma2))
    
    # Replace missing values with imputed
    y_new <- y
    y_new[mis_idx] <- y_mis
    
    # -------------------
    # M-step: Update parameters
    # -------------------
    fit_new <- lm(y_new ~ as.matrix(data[, vars != miss_col]))
    beta_new <- coef(fit_new)
    sigma2_new <- var(resid(fit_new))
    
    # -------------------
    # Check convergence
    # -------------------
    diff <- max(abs(beta_new - beta))
    #if (verbose) cat(sprintf("Iteration %d | diff = %.6f\n", iter, diff))
    
    # Update parameters
    beta <- beta_new
    sigma2 <- sigma2_new
    iter <- iter + 1
  }
  
  # Final imputed dataset
  data_imp <- data
  data_imp[[miss_col]][mis_idx] <- y_new[mis_idx]
  u=pnorm(data_imp[[miss_col]])
  data_imp[[miss_col]] =-80*log(1-u)
  data_imp[[miss_col]] = ifelse(data_imp[[miss_col]]>1000, 1000, data_imp[[miss_col]])
  #data_imp[[miss_col]] <- exp(data_imp[[miss_col]]) # normal -> log normal
  
  # list(
  #   imputed_data = data_imp,
  #   beta = beta,
  #   sigma2 = sigma2,
  #   iterations = iter - 1
  # )
  return(data_imp)
}

data_augment_em_mult <- function(data, max_iter = 100, tol = 1e-6, verbose = TRUE) {
  
  data_work <- data
  # ----- Transform to normal -----
  data_work$BUN_FIRST_VAL <- log(data_work$BUN_FIRST_VAL) # log normal -> normal
  u = 1 - exp(-data_work$CRP_FIRST_VAL/80)
  data_work$CRP_FIRST_VAL=qnorm(u)
  
  vars <- colnames(data_work)
  n <- nrow(data_work)
  
  # Identify variables with missing
  miss_cols <- vars[colSums(is.na(data_work)) > 0]
  if (length(miss_cols) == 0) stop("No missing values found.")
  
  # Save missingness pattern
  mis_mask <- is.na(data_work)
  
  # Initialize missing values (mean imputation)
  for (v in miss_cols) {
    data_work[[v]][mis_mask[, v]] <- mean(data_work[[v]], na.rm = TRUE)
  }
  
  iter <- 1
  diff <- Inf
  
  while (iter <= max_iter && diff > tol) {
    old_data <- data_work
    
    for (v in miss_cols) {
      mis_idx <- mis_mask[, v]
      if (!any(mis_idx)) next
      
      # Regression without using variable v as predictor
      preds <- setdiff(vars, v)
      form <- as.formula(paste(v, "~", paste(preds, collapse = "+")))
      
      fit <- lm(form, data = data_work)
      
      # Predicted means for missing entries
      pred <- predict(fit, newdata = data_work[mis_idx, , drop = FALSE])
      
      # Residual variance
      sigma2 <- summary(fit)$sigma^2
      
      # Stochastic draw
      imp <- rnorm(sum(mis_idx), mean = pred, sd = sqrt(sigma2))
      
      data_work[[v]][mis_idx] <- imp
    }
    
    # Check convergence only on imputed cells
    diff <- max(abs(data_work[mis_mask] - old_data[mis_mask]))
    
    if (verbose) cat(sprintf("Iteration %d | max change = %.6f\n", iter, diff))
    
    iter <- iter + 1
  }
  
  # Back-transform to lognormal or exponential 
  data_work$BUN_FIRST_VAL <- exp(data_work$BUN_FIRST_VAL) 
  data_work$BUN_FIRST_VAL = ifelse(data_work$BUN_FIRST_VAL>500, 500, data_work$BUN_FIRST_VAL)
  u=pnorm(data_work$CRP_FIRST_VAL)
  data_work$CRP_FIRST_VAL =-80*log(1-u)
  data_work$CRP_FIRST_VAL = ifelse(data_work$CRP_FIRST_VAL>1000, 1000, data_work$CRP_FIRST_VAL)
  
  return(data_work)
}


# JOINT.MI

joint.MI <- function(data, mu, sigma, n.imp) {
  dat.imputed <- array(rep(as.matrix(data), n.imp), dim=c(nrow(data),ncol(data),n.imp))
  colnames(dat.imputed) <- colnames(data)  
  missing.col <- c()
  for (i in 1:nrow(data)) {
    if(any(is.na(data[i,]))){ # check down the rows
      x <- data[i,] # match names of x and mu. 
      dep   <- names(x[which(is.na(x))])
      given <- names(x[which(!is.na(x))])
      missing.col <- which(is.na(x))
      x.obs <- as.numeric(x[which(names(x) %in% given)])
      
      condMVN <- rcmvnorm(n=n.imp, mean=mu, sigma=sigma, dep=dep, given=given, X=x.obs)
    }
    
    for(l in missing.col){ 
      index <- which(missing.col == l)
      dat.imputed[i, l, ] <- ifelse(condMVN[,index]<0,0,condMVN[,index]) # put 0's to any negative values
    }
  }
  return(dat.imputed)
}

# CONDITIONAL.MI

micesym <- function(x) {
  (x + t(x))/2
}

conditional.estimation <- function(training_data) {
  model.estimation <- list()
  # fit logistic regressions for binary predictors
  binnames = c("THROMBOSIS_HX", "CENTRAL_LINE_DURING_ENC", "HYPOXEMIA", 
               "Emergency", "Trauma", "Urgent", "Elective", "Missing")
  for(i in 1:ncol(training_data)) {
      formula <- paste(names(training_data[-i]), collapse = " + ")
      formula <- paste(c(names(training_data[i]), formula), collapse = " ~ ")
      if(names(training_data)[i] %in% binnames) {
        model.estimation[[i]] <- glm(formula, training_data, family="binomial")
      } else {
        # We use ridge regression as adopted in mice to facilitate implementation of 
        # imputation processes. An additional advantage is that the ridge penalty 
        # accommodates for some estimation problems in sparse datasets
        X <- cbind(1, as.matrix(training_data[,-i]))
        fit <- tryCatch({
          estimice(x = X, y = training_data[,i])
        }, error = function(e) {
          estimice(x = X, y = training_data[,i], ls.meth = "ridge")
        }, finally = {
        })
        fit$df.residual <- fit$df # Facilitate compatibility with glm objects
        class(fit) <- "estimice" # Facilitate compatibility with glm objects
        model.estimation[[i]] <- fit
      }
  }
  return(model.estimation)
}

# Multiple Imputation of a single missing value for a single patient.
conditional.MI.single <- function(test_case, model_estimation, mu, sigma, data_classes, n.imp) {
  
  missing.col <- which(colSums(is.na(test_case)) > 0)
  if (length(missing.col) > 1) {
    stop("This function is only allowed for imputation of a single missing value.")
  }
  
  out <- rep(NA, n.imp)
  
  # observed variables (explanatory variables)
  predictors <- paste(names(test_case[,-missing.col]), collapse = " + ")
  
  # model matrix of missing values explained by explanatory variables
  pred.data <- model.matrix(formula(paste("~", predictors)), data = test_case)
  
  if(data_classes[missing.col]=="logical") {
    
    # Directly draw all 'n.imp' draws for the regression coefficients
    beta_star <- rmvt(n = n.imp,
                      delta = model_estimation[[missing.col]]$c, 
                      sigma = model_estimation[[missing.col]]$v,
                      df = model_estimation[[missing.col]]$df.residual)
    
    prob <- rep(NA, n.imp)
    for(j in 1:n.imp) {
      prob[j] <- 1/(1+exp(-as.numeric(pred.data[1,]) %*% beta_star[j,]))
    }
    out <- rbinom(n = n.imp, size = 1, prob = prob)
  } else {
    # In case we are dealing with the imputation of a continuous variable, we will use the same functionalities
    # as mice.impute.norm, which corresponds to Bayesian Linear Regression. The code below is adapted from the 
    # function mice::.norm.draw
    
    p <- model_estimation[[missing.col]] # This should be an object of class 'estimice'
    
    sigma_star <- sqrt(sum((p$r)^2)/rchisq(n.imp, p$df))
    for(j in 1:n.imp) {
      beta_star <- p$c + (t(chol(micesym(p$v))) %*% rnorm(length(p$c))) * sigma_star[j]
      out[j] <- (pred.data[1,] %*% beta_star) + rnorm(1) * sigma_star[j] 
    }
  }
  return(out)
}

conditional.MI <- function(data, # Data frame with the patient data
                           model_estimation, # List of conditional imputation models (multivariable models for each predictor)
                           mu = rep(0, ncol(data)), # Mean vector to initialize Gibbs sampler 
                           sigma = diag(rep(10, ncol(data))), # Covariance matrix to initialize Gibbs sampler 
                           data_classes, 
                           n.imp, # Number of required imputed datasets
                           maxit = 15) # set maximum iterations for convergence imputations
{
  
  dat.imputed <- array(rep(as.matrix(data), n.imp), dim=c(nrow(data),ncol(data),n.imp))
  colnames(dat.imputed) <- colnames(data)
  
  # get variable for which models were fitted with family "binomial"
  binnames <- colnames(dat.imputed)[unlist(lapply(model_estimation, function(x) x$family$family == "binomial"))]
  
  # each row/patient is taken separately
  for(i in 1:nrow(data)) {
    # identify missing variables
    missing.col <- which(is.na(data[i,]))
    
    # No Gibbs sampler needed if a patient has one missing value
    if(length(missing.col) == 1) {
      dat.imputed[i,missing.col,1:n.imp] <- conditional.MI.single(test_case = data[i,],  
                                                                  model_estimation, 
                                                                  mu, 
                                                                  sigma, 
                                                                  data_classes, 
                                                                  n.imp)
    } else if(length(missing.col) > 1) {
      
      # initialize vector for patient -> these will for be the concurrent imputations
      gibbsdata <- lapply(data_classes, vector, length = 1)
      
      # iterate over till all multiple imputations are done for this patient
      for(l in 1:n.imp) {  
        # Generate initial draws for imputation & iterate from 1 to 23 variables
        for(j in 1:length(data)) { 
          gibbsdata[[j]][1] <- data[i,j]
          # itial value = draw from multivariate normal distribution
          gibbsdata[[j]][2] <- ifelse(j %in% missing.col, mvtnorm::rmvnorm(1, mean=mu, sigma=sigma)[j], gibbsdata[[j]][1])
        }
        # convert vector to dataframe
        temp_data <- as.data.frame(gibbsdata)
        
        # first 2 iterations already done in intialization, iterate till max iterations for convergence (15)
        
        # Iterate over the imputation cycles
        for (iter in 3:maxit) {
          temp_data[iter,] <- temp_data[iter-1,]
          
          # Iterate ver the different variables with missing values
          # Each time, use the most recently imputed value from the remaining variables
          
          for(k in 1:length(missing.col)){
            # Extract last available data
            test_case <- temp_data[iter,]
            test_case[missing.col[k]] <- NA
            test_case[missing.col[k]] <- conditional.MI.single(test_case = test_case,  
                                                               model_estimation, 
                                                               mu, 
                                                               sigma, 
                                                               data_classes, 
                                                               n.imp = 1)
            # fill temporary data with last iteration
            temp_data[iter,] <- test_case
          }
        }
        
        # Save most recent imputation as imputed dataset 'l' for patient 'i',
        dat.imputed[i,missing.col,l] <- unlist(temp_data[maxit, missing.col])
      }
    } # end if-structure of the Gibbs sampler
  } # end iteration over the patients
  return(dat.imputed)
}


simple.conditional.mi = function(data, outcome){
  p <- model_estimation[[missing.col]] # This should be an object of class 'estimice'
  
  sigma_star <- sqrt(sum((p$r)^2)/rchisq(n.imp, p$df))
  for(j in 1:n.imp) {
    beta_star <- p$c + (t(chol(micesym(p$v))) %*% rnorm(length(p$c))) * sigma_star[j]
    out[j] <- (pred.data[1,] %*% beta_star) + rnorm(1) * sigma_star[j] 
  }
}

# create a function to find k to adjust absence rate and event rate
find_k = function(input, target.rate){
  k = -1000
  while(mean(expit(input+k)) <= target.rate){
    k = k + 0.01
  }
  return(k)
}

generate_ogm <- function(data, dependence, outcome.rate, interaction) {
  # coefficient sets
  coefs <- list(
    weak     = seq(0.05, 0.5, length.out = 16),
    moderate = seq(0.05, 1.0, length.out = 16),
    strong   = seq(0.05, 2.0, length.out = 16)
  )
  
  form_no <- ~ WBC_FIRST_VAL + THROMBOSIS_HX + HYPOXEMIA + BMI_FIRST_VALUE +
    BUN_FIRST_VAL + HEART_RATE_FIRST_VALUE + CRP_FIRST_VAL + FIRST_BRADEN_SCALE_SCORE +
    admit.type + CENTRAL_LINE_DURING_ENC
  
  form_int <- ~ WBC_FIRST_VAL + THROMBOSIS_HX + HYPOXEMIA + THROMBOSIS_HX:BUN_FIRST_VAL +
    THROMBOSIS_HX:WBC_FIRST_VAL + BMI_FIRST_VALUE + THROMBOSIS_HX:CRP_FIRST_VAL +
    BUN_FIRST_VAL + HEART_RATE_FIRST_VALUE + CRP_FIRST_VAL + FIRST_BRADEN_SCALE_SCORE +
    admit.type + CENTRAL_LINE_DURING_ENC
  
  scale_cols <- if(interaction==TRUE) -c(1,3,4,10:17) else -c(1,3,4,10:14) # exclude intercept & categorical
  coef_idx <- if(interaction==TRUE) 1:16 else 1:13
  
  # build matrix and apply robust scaling
  mat <- model.matrix(if(interaction==TRUE) form_int else form_no, data)
  mat[, scale_cols] <- RobScale(mat[, scale_cols], center = TRUE, scale = TRUE)
  
  # get dependence-specific coefficients
  coef_set <- coefs[[dependence]]
  if (is.null(coef_set)) stop("dependence must be 'weak', 'moderate', or 'strong'")
  
  pred <- mat %*% c(0, coef_set[coef_idx])
  
  # adjust intercept and simulate outcomes
  k <- find_k(pred, outcome.rate)
  out <- rbinom(nrow(data), 1, expit(pred + k))
  
  return(out)
}


# create a function to generate missingness for varying absence rate
generate_mgm = function(data, dependence, absence.rate){
  coef_weak = -seq(from=0.05, to=0.5, length.out=12) # weak
  coef_mod = -seq(from=0.05, to=1, length.out=12) # moderate
  coef_strong = -seq(from=0.05, to=2, length.out=12) # strong

  mat = model.matrix(~ FIRST_BRADEN_SCALE_SCORE + BMI_FIRST_VALUE + THROMBOSIS_HX + WBC_FIRST_VAL + 
                       HYPOXEMIA + BUN_FIRST_VAL + HEART_RATE_FIRST_VALUE + 
                       CENTRAL_LINE_DURING_ENC + admit.type + vte.hac2, data)
  mat[, -c(1,4,6,9:14)] <- RobScale(mat[, -c(1,4,6,9:14)], center=T, scale=T)
  
  if(dependence=="weak"){
    pred = mat %*% c(0,coef_weak, -0.5) # vte and missing have negative association
  }
  else if(dependence=="moderate"){
    pred = mat %*% c(0,coef_mod, -1)
  }
  else if(dependence=="strong"){
    pred = mat %*% c(0,coef_strong, -2)
  }
  k = find_k(pred, absence.rate)
  out = rbinom(nrow(data),1,expit(pred+k)) # rbinom
  return(out)
}

# create a function to generate plasmode simulated dataset
plasmode <- function(data, size=10000, dep.outcome, outcome.rate, outcome.interaction, dep.absence, absence.rate){
  # resample covariate information
  resampled = sample(nrow(data), size) # use complete data
  data_sample = data[resampled,]
  
  # assign outcome
  data_sample$vte.hac2 = generate_ogm(data_sample, dep.outcome, outcome.rate, outcome.interaction)
  data_sample$vte.hac2 = as.factor(data_sample$vte.hac2)
  
  # assign missingness
  data_sample$miss = generate_mgm(data_sample, dep.absence, absence.rate)
  
  # split data into derivation and implementation cohort
  smp_size <- floor(0.70 * nrow(data_sample)) # 70% of sample
  train_ind <- sample(seq_len(nrow(data_sample)), size = smp_size)
  sample_der <- data_sample[train_ind, ]
  sample_imp <- data_sample[-train_ind, ]
  
  sample_der$imp = 0
  sample_imp$imp = 1 # imp=1 indicates implementation cohort
  data_sample = data.frame(rbind(sample_der, sample_imp))

  return(data_sample)
}


# create a function to develop five prediction models all at once on a derivation cohort
develop_mod = function(data, mod){
  # drop "miss" column
  data = data[1:(ncol(data)-1)]
  
  if(mod=="logistic regression"){
    # 1. Logistic reg not allowing interaction
    fitControl <- trainControl(method = "cv", number = 5,
                               classProbs = T, summaryFunction = BigSummary,
                               seeds = set.seed(1234))
    data$vte.hac2 = ifelse(data$vte.hac2==1, "yes", "no")
    glm = caret::train(vte.hac2 ~ .,
                       data=data, method = "glm", trControl = fitControl,
                       metric = "AUCROC")
    out = glm
  }
  else if(mod=="random forest"){
    # 2. Random forest
    n <- nrow(data)
    s.size <- n / 5 # 5-fold
    swor <- TRUE # cross validation: sampling without replacement
    samp <- randomForestSRC:::make.sample(ntree=3000, n, s.size, swor)
    rf.src.mod = rfsrc(vte.hac2 ~ .,
                       data=as.data.frame(data), ntree = 3000,
                       rfq = TRUE, nodesize=1, mtry=10,
                       bootstrap = "by.user", samp = samp)
    out = rf.src.mod
  }
  else if(mod=="xg boost"){
    # 3. XG boost
    fitControl <- trainControl(method = "cv", number = 5,
                               classProbs = T, summaryFunction = BigSummary,
                               seeds = set.seed(1234))
    data$vte.hac2 = ifelse(data$vte.hac2==1, "yes", "no")
    xgb = caret::train(vte.hac2 ~ .,
                data=data, method = "xgbTree", trControl = fitControl,
                metric = "AUCROC", tuneLength = 3, verbosity = 0)
    out = xgb
  }
 else if(mod=="neural network"){
   # 4. Deep neural network
   x.mod = data %>% dplyr::select(-vte.hac2) %>% 
     fastDummies::dummy_cols(select_columns = "admit.type",
                             remove_selected_columns = TRUE,
                             remove_first_dummy = FALSE) %>%
     dplyr::rename_with(~ gsub("admit.type_", "", .x)) %>%
     mutate(Missing = ifelse("Missing" %in% names(.), Missing, 0)) %>%
     dplyr::rename(Trauma = `Trauma `) %>% 
     mutate(across(where(is.factor), as.numeric)) %>%
     mutate(across(where(is.character), as.numeric)) %>%
     as.matrix()
   # y.mod = ifelse(data$vte.hac2=="Yes",1,0)
   y.mod = data$vte.hac2
   tensorflow::set_random_seed(1234)
   
   dnn = keras_model_sequential()
   dnn %>%
     layer_dense(units = 64, activation = "relu",
                 kernel_initializer = "he_normal") %>%
     layer_batch_normalization() %>%
     layer_dropout(rate=0.4) %>%
     layer_dense(units = 32, activation = "relu", kernel_initializer = "he_normal") %>%
     layer_batch_normalization() %>%
     layer_dropout(rate=0.3) %>%
     layer_dense(units = 16, activation = "relu", kernel_initializer = "he_normal") %>%
     layer_batch_normalization() %>%
     layer_dropout(rate=0.2) %>%
     layer_dense(units = 8, activation = "relu", kernel_initializer = "he_normal") %>%
     layer_batch_normalization() %>%
     layer_dropout(rate=0.2) %>%
     layer_dense(units = 1, activation ="sigmoid")
   dnn %>% compile(
     loss = 'binary_crossentropy',
     optimizer = optimizer_adam(),
     metrics = c('accuracy'))
   history = dnn %>% fit(x = x.mod, y = y.mod, 
                         epochs = 20, batch_size = 10, validation_split = 0.4, verbose = 0) 
   out = dnn
 }
  return(out)
}

simulate.impute <- function(data){
  sample_der <- data[data$imp==0, ]
  sample_imp <- data[data$imp==1, ]
  sample_imp$vte.hac2 <- as.numeric(as.character(sample_imp$vte.hac2))
  out <- list()
  
  # True values
  out$true <- sample_imp$CRP_FIRST_VAL[sample_imp$miss==1]
  
  # 1. Median
  med <- median(sample_der$CRP_FIRST_VAL)
  # med <- median(sample_imp[sample_imp$miss==0,"CRP_FIRST_VAL"])
  sample_imp[sample_imp$miss==1,"CRP_FIRST_VAL"] <- med
  out$median <- sample_imp[sample_imp$miss==1,"CRP_FIRST_VAL"]
  
  # 2. MICE (save all 5 imputations)
  sample_imp[sample_imp$miss==1,"CRP_FIRST_VAL"] <- NA
  m_times = round(100*sum(sample_imp$miss)/nrow(sample_imp)) # number of imputation = percentage of absence 
  mice_obj <- mice(data = sample_imp[, 2:11], m = m_times, method = 'pmm', print = FALSE)
  
  out$mice_list <- lapply(1:mice_obj$m, function(m){
    tmp <- sample_imp
    tmp[sample_imp$miss==1,"CRP_FIRST_VAL"] <- complete(mice_obj, m)[sample_imp$miss==1,"CRP_FIRST_VAL"]
    tmp
  })
  out$mice <- complete(mice_obj, 1)[sample_imp$miss==1,"CRP_FIRST_VAL"]  # take the first set for structure consistency
  
  # 3. Data augmentation
  sample_imp[sample_imp$miss==1,"CRP_FIRST_VAL"] <- NA
  
  out$da_list <- lapply(1:5, function(m){
    tmp <- sample_imp
    tmp[sample_imp$miss==1,"CRP_FIRST_VAL"] <- data_augment_em(sample_imp[c(3,7:11)])[sample_imp$miss==1,"CRP_FIRST_VAL"]
    tmp
  })
  # out$da <- sample_imp[sample_imp$miss==1,"CRP_FIRST_VAL"]
  out$da <- out$da_list[[1]][sample_imp$miss==1,"CRP_FIRST_VAL"] # take the first set
  
  # 4. Conditional modeling
  sample_imp[sample_imp$miss==1,"CRP_FIRST_VAL"] <- NA
  sample_imp[sample_imp$miss==1,"CRP_FIRST_VAL"] <-
    predict(lm(CRP_FIRST_VAL ~ ., data = sample_der[2:11]), sample_imp[2:11])[sample_imp$miss==1]
  out$cm <- sample_imp[sample_imp$miss==1,"CRP_FIRST_VAL"]
  
  # 5. Joint modeling
  sample_imp[sample_imp$miss==1,"CRP_FIRST_VAL"] <- NA
  x <- sample_imp[, 2:11]
  x_wide <- x %>%
    fastDummies::dummy_cols(select_columns = "admit.type",
                            remove_selected_columns = TRUE,
                            remove_first_dummy = FALSE) %>%
    dplyr::rename_with(~ gsub("admit.type_", "", .x)) %>%
    mutate(
      Missing = ifelse("Missing" %in% names(.), Missing, 0),
      `Trauma ` = ifelse("Trauma " %in% names(.), `Trauma `, 0)
    ) %>%
    dplyr::rename(Trauma = `Trauma `)
  
  mean <- colMeans(x_wide[, 1:length(x_wide)], na.rm = TRUE)
  var <- cov(x_wide[, 1:length(x_wide)], use = "complete.obs")
  var[which(var != diag(var))] <- 0
  cov <- var + diag(ncol(var)) * 1e-5
  
  jm_imps <- joint.MI(x_wide[, 1:length(x_wide)], mean, cov, n.imp = m_times)
  out$jm_list <- lapply(1:m_times, function(m){
    tmp <- sample_imp
    tmp[sample_imp$miss==1,"CRP_FIRST_VAL"] <- jm_imps[sample_imp$miss==1,"CRP_FIRST_VAL", m]
    tmp
  })
  out$jm <- jm_imps[sample_imp$miss==1,"CRP_FIRST_VAL", 1] # take the first set
  
  return(out)
}

simulate.predict <- function(data, data.impute){
  varlist <- c("THROMBOSIS_HX", "BMI_FIRST_VALUE", "HYPOXEMIA",
                   "WBC_FIRST_VAL", "BUN_FIRST_VAL",
                   "CENTRAL_LINE_DURING_ENC", "admit.type", "HEART_RATE_FIRST_VALUE", 
                   "FIRST_BRADEN_SCALE_SCORE", "CRP_FIRST_VAL", "vte.hac2", "miss")
  
  sample_der <- data[data$imp==0, varlist]
  sample_imp <- data[data$imp==1, varlist]
  sample_imp$vte.hac2 <- as.numeric(as.character(sample_imp$vte.hac2))
  
  # Train models once
  glm <- develop_mod(sample_der, mod = "logistic regression")
  rf  <- develop_mod(sample_der, mod = "random forest")
  xgb <- develop_mod(sample_der, mod = "xg boost")
  dnn <- develop_mod(sample_der, mod = "neural network")
  
  # Helper: get predictions from any model
  predict_probs <- function(model, type, newdata, xmod=NULL){
    newdata <- newdata[, varlist]
    if(type %in% c("glm", "xgb")){
      predict(model, newdata=newdata[, 1:(ncol(newdata)-2)], type="prob")[newdata$miss==1, 2]
    } else if(type=="rf"){
      predict(model, newdata[newdata$miss==1, 1:(ncol(newdata)-2)])$predicted[,2]
    } else if(type=="dnn"){
      xmod <- newdata %>%
        dplyr::select(-miss,-vte.hac2) %>%
        fastDummies::dummy_cols(select_columns = "admit.type",
                                remove_selected_columns = TRUE,
                                remove_first_dummy = FALSE) %>%
        dplyr::rename_with(~ gsub("admit.type_", "", .x)) %>%
        mutate(Missing = ifelse("Missing" %in% names(.), Missing, 0),
               `Trauma ` = ifelse("Trauma " %in% names(.), `Trauma `, 0)) %>%
        dplyr::rename(Trauma = `Trauma `)
      predict(model, data.matrix(xmod))[newdata$miss==1]
    }
  }
  
  # Helper: average predictions across multiple imputations
  average_preds <- function(model, type, imp_list){
    preds <- sapply(1:length(imp_list), function(m){
      newdata <- imp_list[[m]]
      predict_probs(model, type, newdata)
    })
    rowMeans(preds)
  }
  
  # Construct simple imputations
  imputations <- list(
    true   = {tmp <- sample_imp; tmp},
    median = {tmp <- sample_imp; tmp[sample_imp$miss==1,"CRP_FIRST_VAL"] <- data.impute$median; tmp},
    cm     = {tmp <- sample_imp; tmp[sample_imp$miss==1,"CRP_FIRST_VAL"] <- data.impute$cm; tmp}
  )
  
  # Collect predictions
  out <- list()
  
  for(mod_name in c("glm", "rf", "xgb", "dnn")){
    model <- get(mod_name)
    out[[mod_name]] <- data.frame(
      true   = predict_probs(model, mod_name, imputations$true),
      median = predict_probs(model, mod_name, imputations$median),
      mice   = average_preds(model, mod_name, data.impute$mice_list),
      da     = average_preds(model, mod_name, data.impute$da_list),
      cm     = predict_probs(model, mod_name, imputations$cm),
      jm     = average_preds(model, mod_name, data.impute$jm_list)
    )
  }
  # out[["ens"]] = Reduce('+',list(out[["glm"]],out[["rf"]],out[["xgb"]],out[["dnn"]]))/4
  return(out)
}

simulate.all = function(data, dep.outcome, outcome.rate, outcome.interaction, dep.absence, absence.rate, nsim){
  
  sdat = lapply(seq_len(nsim), function(x) plasmode(data, size=10000, dep.outcome, outcome.rate, 
                                                    outcome.interaction, dep.absence, absence.rate))
  out = list()
  for (i in 1:nsim){
    idat = simulate.impute(sdat[[i]])
    out$impute[[i]] = idat
    pdat = simulate.predict(sdat[[i]], idat)
    out$predict[[i]] = pdat
    print(paste("simulation",i,"is complete!"))
  }

  return(out)
}
