# all Machine learning (ML) algorithms
RunML <- function(method, Train_expr, Train_surv, mode = "Model", timeVar = "OS.time", statusVar = "OS", ...){
  # for example: Enet [alpha=0.4]
  method = gsub(" ", "", method) 
  method_name = gsub("(\\w+)\\[(.+)\\]", "\\1", method)  # get name of ML algorithm, e.g., Enet
  method_param = gsub("(\\w+)\\[(.+)\\]", "\\2", method) # get parameter of ML algorithm, e.g., alpha=0.4
  
  method_param = switch(
    EXPR = method_name,
    "Enet" = list("alpha" = as.numeric(gsub("alpha=", "", method_param))),
    "StepCox" = list("direction" = method_param),
    NULL
  )
  message("Run ", method_name, " algorithm for ", mode, "; ",
          method_param, ";",
          " using ", ncol(Train_expr), " Variables")
  
  args = list("Train_expr" = Train_expr,
              "Train_surv" = Train_surv,
              "mode" = mode,
              "timeVar" = timeVar, "statusVar" = statusVar)
  args = c(args, method_param)
  
  obj <- do.call(what = paste0("Run", method_name),
                 args = args) 
  
  if(mode == "Variable"){
    message(length(obj), " Variables retained;\n")
  }else{message("\n")}
  return(obj)
}

RunEnet <- function(Train_expr, Train_surv, mode, timeVar, statusVar, alpha){
  cv.fit = cv.glmnet(x = Train_expr,
                     y = Surv(Train_surv[[timeVar]], Train_surv[[statusVar]]),
                     family = "cox", alpha = alpha, nfolds = 10)
  fit = glmnet(x = Train_expr,
               y = Surv(Train_surv[[timeVar]], Train_surv[[statusVar]]),
               family = "cox", alpha = alpha, lambda = cv.fit$lambda.min)
  fit$subFeature = colnames(Train_expr)
  if (mode == "Model") return(fit)
  co = coef(fit)[, 1]
  if (mode == "Variable") return(names(co)[co!=0])
}

RunLasso <- function(Train_expr, Train_surv, mode, timeVar, statusVar){
  RunEnet(Train_expr, Train_surv, mode, timeVar, statusVar, alpha = 1)
}

RunRidge <- function(Train_expr, Train_surv, mode, timeVar, statusVar){
  RunEnet(Train_expr, Train_surv, mode, timeVar, statusVar, alpha = 0)
}

RunStepCox <- function(Train_expr, Train_surv, mode, timeVar, statusVar, direction){
  fit <- step(coxph(formula = Surv(Train_surv[[timeVar]], Train_surv[[statusVar]]) ~ .,
                    data = as.data.frame(Train_expr)),
              direction = direction, trace = 0)
  fit$subFeature = colnames(Train_expr)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(names(coef(fit)))
}

RunsurvivalSVM <- function(Train_expr, Train_surv, mode, timeVar, statusVar){
  fit = survivalsvm(formula = Surv(Train_surv[[timeVar]], Train_surv[[statusVar]]) ~ .,
                    data= as.data.frame(Train_expr),
                    gamma.mu = 1, opt.meth = "ipop")
  fit$subFeature = colnames(Train_expr)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(fit["var.names"])
}

RunCoxBoost <- function(Train_expr, Train_surv, mode, timeVar, statusVar){
  pen <- optimCoxBoostPenalty(time = Train_surv[[timeVar]],
                              status = Train_surv[[statusVar]],
                              x = Train_expr,
                              trace = F, start.penalty=500, parallel = F)
  cv.res <- cv.CoxBoost(time = Train_surv[[timeVar]],
                        status = Train_surv[[statusVar]],
                        x = Train_expr,
                        maxstepno=500, K=10, type="verweij", penalty=pen$penalty)
  fit <- CoxBoost(time = Train_surv[[timeVar]],
                  status = Train_surv[[statusVar]],
                  x = Train_expr,
                  stepno = cv.res$optimal.step,
                  penalty = pen$penalty)
  fit$subFeature = colnames(Train_expr)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(names(coef(fit)[abs(coef(fit))>0]))
}

RunSuperPC <- function(Train_expr, Train_surv, mode, timeVar, statusVar){
  data <- list(x = t(Train_expr),
               y = Train_surv[[timeVar]],
               censoring.status = Train_surv[[statusVar]],
               featurenames = colnames(Train_expr))
  fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) # default
  cv.fit <- suppressWarnings(superpc.cv(fit, data,
                                        n.threshold = 20,# default
                                        n.fold = 10,
                                        n.components = 3,
                                        min.features = 5,
                                        max.features = nrow(data$x),
                                        compute.fullcv = TRUE,
                                        compute.preval = TRUE))
  fit$threshold <- cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])]
  fit$data <- data
  fit$subFeature = colnames(Train_expr)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(names(fit$feature.scores)[abs(fit$feature.scores)>0.5])
}

RunplsRcox <- function(Train_expr, Train_surv, mode, timeVar, statusVar){
  data <- list(x = Train_expr,
               time = Train_surv[[timeVar]],
               status = Train_surv[[statusVar]])
  cv.plsRcox.res = cv.plsRcox(data = data,
                              nt=10, verbose = FALSE)
  fit <- plsRcox(Xplan = data$x,
                 time = data$time,
                 event = data$status,
                 nt = as.numeric(cv.plsRcox.res[5]),
                 verbose = F, sparse = T)
  fit$subFeature = colnames(Train_expr)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(rownames(fit$Coeffs)[fit$Coeffs!=0])
}

RunRSF <- function(Train_expr, Train_surv, mode, timeVar, statusVar){
  rf_nodesize = 2 # may modify
  fit <- rfsrc(formula = formula(paste0("Surv(", timeVar, ", ", statusVar, ")", "~.")),
               data = cbind(Train_expr, Train_surv),
               ntree = 1000, nodesize = rf_nodesize,
               splitrule = 'logrank',
               importance = T,
               proximity = T,
               forest = T)
  fit$subFeature = colnames(Train_expr)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(var.select(fit, verbose = F)$topvars)
}

RunGBM <- function(Train_expr, Train_surv, mode, timeVar, statusVar){
  fit <- gbm(formula = Surv(Train_surv[[timeVar]], Train_surv[[statusVar]]) ~ .,
             data = as.data.frame(Train_expr),
             distribution = 'coxph',
             n.trees = 10000,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 10,n.cores = 6)
  best <- which.min(fit$cv.error)
  fit <- gbm(formula = Surv(Train_surv[[timeVar]], Train_surv[[statusVar]]) ~ .,
             data = as.data.frame(Train_expr),
             distribution = 'coxph',
             n.trees = best,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001, n.cores = 8)
  fit$subFeature = colnames(Train_expr)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(rownames(summary.gbm(fit))[summary.gbm(fit)$rel.inf>0])
}

quiet <- function(..., messages=FALSE, cat=FALSE){
  if(!cat){
    sink(tempfile())
    on.exit(sink())
  }
  out <- if(messages) eval(...) else suppressMessages(eval(...))
  out
}

RunEval <- function(fit, 
                    Test_expr = NULL, 
                    Test_surv = NULL, 
                    Train_expr = NULL, 
                    Train_surv = NULL, 
                    Train_name = NULL,
                    cohortVar = "Cohort",
                    timeVar = "OS.time", 
                    statusVar = "OS"){
  
  if(!is.element(cohortVar, colnames(Test_surv))) {
    stop(paste0("There is no [", cohortVar, "] indicator, please fill in one more column!"))
  } 
  
  if((!is.null(Train_expr)) & (!is.null(Train_surv))) {
    new_data <- rbind.data.frame(Train_expr[, fit$subFeature],
                                 Test_expr[, fit$subFeature])
    
    if(!is.null(Train_name)) {
      Train_surv$Cohort <- Train_name
    } else {
      Train_surv$Cohort <- "Training"
    }
    colnames(Train_surv)[ncol(Train_surv)] <- cohortVar
    Test_surv <- rbind.data.frame(Train_surv[,c(cohortVar, timeVar, statusVar)],
                                  Test_surv[,c(cohortVar, timeVar, statusVar)])
    Test_surv[,1] <- factor(Test_surv[,1], 
                            levels = c(unique(Train_surv[,cohortVar]), setdiff(unique(Test_surv[,cohortVar]),unique(Train_surv[,cohortVar]))))
  } else {
    new_data <- Test_expr[, fit$subFeature]
  }
  
  RS <- quiet(switch(
    EXPR = class(fit)[1],
    "coxnet"      = predict(fit, type = 'link', as.matrix(new_data)),
    "coxph"       = predict(fit, type = 'risk', as.data.frame(new_data)),
    "survivalsvm" = predict(fit, as.data.frame(new_data))$predicted,
    "CoxBoost"    = predict(fit, type = "lp", as.data.frame(new_data)),
    "superpc"     = superpc.predict(object = fit,
                                    data = fit$data,
                                    newdata = list(x = t(as.matrix(new_data))),
                                    threshold = fit$threshold,
                                    n.components = 1)$v.pred,
    "plsRcoxmodel" = predict(fit, type = "lp", as.data.frame(new_data)),
    "rfsrc"        = predict(fit, as.data.frame(new_data))$predicted,
    "gbm"          = predict(fit, type = 'link', as.data.frame(new_data))
  ))
  Predict.out <- Test_surv
  Predict.out$RS <- as.vector(RS)
  Predict.out <- split(x = Predict.out, f = Predict.out[,cohortVar])
  f <- as.formula(paste0("Surv(", timeVar,",",statusVar,")~RS"))
  unlist(lapply(Predict.out, function(data){
    unname(summary(coxph(formula = f,
                         data = data))$concordance["C"])
  }))
}

