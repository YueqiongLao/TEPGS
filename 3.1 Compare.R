# core functions
makeCox <- function(Features, 
                    coefs = NULL, 
                    SIGname,
                    Train_expr, 
                    Train_surv, 
                    unmatchR = 0.2,
                    statusVar = "OS", 
                    timeVar = "OS.time"){
  unmatch.Features <- setdiff(Features, colnames(Train_expr))
  
  if (all(is.na(coefs))) coefs = NULL
  if (!is.null(coefs) & is.null(names(coefs))) names(coefs) = Features
  
  if (length(unmatch.Features)/length(Features) > unmatchR){
    message(SIGname, ": Warnings, ", length(unmatch.Features), " of ", length(Features), 
            " (", round(length(unmatch.Features) / length(Features) * 100) , "%)", 
            " Features were not matched in train set. Skip this signature.")
    return(NULL)
  }else{
    if (length(unmatch.Features) > 0){
      message(SIGname, ": ", length(unmatch.Features), " of ", length(Features), 
              " (", round(length(unmatch.Features) / length(Features) * 100) , "%)", 
              " Features were not matched in train set. ")
    }
    subFeature <- intersect(Features, colnames(Train_expr))
    
    if (!is.null(coefs)){
      model = setNames(object = coefs[subFeature], nm = subFeature)
    }else{
      model <- coxph(formula = Surv(Train_surv[[timeVar]], Train_surv[[statusVar]]) ~ .,
                     data = as.data.frame(Train_expr[, subFeature]))$coefficients
    }
    return(model)
  }
  
 
}

calCindex <- function(model, name,
                      Test_expr, 
                      Test_surv, 
                      Train_expr = NULL,
                      Train_surv = NULL,
                      Train_name = NULL,
                      CohortVar = "Cohort", 
                      metaCohort = TRUE,
                      timeVar = "OS.time", 
                      statusVar = "OS"){
  
  CohortSet <- list(); CohortSurv <- list()
  CohortSet[["Test"]] <- Test_expr
  Test_surv$Sample <- rownames(Test_surv)
  CohortSurv[["Test"]] <- Test_surv[,c("Sample", CohortVar, timeVar, statusVar)]
  
  if((!is.null(Train_expr)) & (!is.null(Train_surv))){
    newSurv <- Train_surv
    
    if(!is.null(Train_name)) newSurv[[CohortVar]] <- Train_name 
    else newSurv[[CohortVar]] <- "Training"
    newSurv$Sample <- rownames(newSurv)
    
    CohortSet[["Train"]] <- Train_expr
    CohortSurv[["Train"]] <- newSurv[,c("Sample", CohortVar, timeVar, statusVar)]
  }
  
  if(metaCohort){
    CohortSet[["Meta"]] <- Test_expr
    CohortSurv[["Meta"]] <- Test_surv
    CohortSurv[["Meta"]][[CohortVar]] <- "MetaCohort"
  }
  
  Cohortlevel <- c(unique(CohortSurv[["Train"]][[CohortVar]]), unique(CohortSurv[["Test"]][[CohortVar]]))
  if(metaCohort) Cohortlevel <- c(Cohortlevel, unique(CohortSurv[["Meta"]][[CohortVar]]))
  CohortSet <- do.call(rbind, CohortSet)
  CohortSurv <- do.call(rbind, CohortSurv)
  CohortSurv[[CohortVar]] <- factor(CohortSurv[[CohortVar]], levels = Cohortlevel)

  if (is.null(model)){
    return(NULL)
  }else if(length(intersect(names(model), colnames(CohortSet))) == 0){
    message(name, ": estimate C-index using calculated RS score")
    if (length(setdiff(CohortSurv$Sample, names(model)))>0) 
      message("There are no match RS score for ", length(setdiff(CohortSurv$Sample, names(model))), " samples")
    Predict.out <- CohortSurv
    Predict.out$RS <- as.vector(model[Predict.out$Sample])
    Predict.out <- split(x = Predict.out, f = Predict.out[,CohortVar])
    f <- as.formula(paste0("Surv(", timeVar,",",statusVar,")~RS"))
    res <- do.call(rbind, lapply(Predict.out, function(data){
      s = summary(coxph(formula = f, data = data))
      c(s$concordance, n = s$n)
    }))
    res <- as.data.frame(res)
    res$Cohort <- factor(rownames(res), levels = Cohortlevel)
    colnames(res) <- c("C", "se", "n", "Cohort")
    res$method = name
    res$RS <- lapply(Predict.out, function(x) x$RS)
    return(res)
  }else{
    message(name, ": estimate C-index using Cox model fitted score")
    
    model <- model[names(model) %in% colnames(CohortSet)]
    model.mat <- as.matrix(model)
    RS <- as.matrix(CohortSet[, names(model)]) %*% model.mat
    
    Predict.out <- CohortSurv
    Predict.out$RS <- as.vector(RS)
    Predict.out <- split(x = Predict.out, f = Predict.out[,CohortVar])
    f <- as.formula(paste0("Surv(", timeVar,",",statusVar,")~RS"))
    res <- do.call(rbind, lapply(Predict.out, function(data){
      s = summary(coxph(formula = f, data = data))
      c(s$concordance, n = s$n)
    }))
    res <- as.data.frame(res)
    res$Cohort <- factor(rownames(res), levels = Cohortlevel)
    colnames(res) <- c("C", "se", "n", "Cohort")
    res$method = name
    res$RS <- lapply(Predict.out, function(x) x$RS)
    return(res)
  }
}

# other functions
standarize.fun <- function(indata, centerFlag, scaleFlag) {  
  scale(indata, center=centerFlag, scale=scaleFlag)
}

scaleData <- function(data, cohort = NULL, centerFlags = NULL, scaleFlags = NULL){
  samplename = rownames(data)
  if (is.null(cohort)){
    data <- list(data); names(data) = "training"
  }else{
    data <- split(as.data.frame(data), cohort)
  }
  
  if (is.null(centerFlags)){
    centerFlags = F; message("No centerFlags found, set as FALSE")
  }
  if (length(centerFlags)==1){
    centerFlags = rep(centerFlags, length(data)); message("set centerFlags for all cohort as ", unique(centerFlags))
  }
  if (is.null(names(centerFlags))){
    names(centerFlags) <- names(data); message("match centerFlags with cohort by order\n")
  }
  
  if (is.null(scaleFlags)){
    scaleFlags = F; message("No scaleFlags found, set as FALSE")
  }
  if (length(scaleFlags)==1){
    scaleFlags = rep(scaleFlags, length(data)); message("set scaleFlags for all cohort as ", unique(scaleFlags))
  }
  if (is.null(names(scaleFlags))){
    names(scaleFlags) <- names(data); message("match scaleFlags with cohort by order\n")
  }
  
  centerFlags <- centerFlags[names(data)]; scaleFlags <- scaleFlags[names(data)]
  outdata <- mapply(standarize.fun, indata = data, centerFlag = centerFlags, scaleFlag = scaleFlags, SIMPLIFY = F)
  # lapply(out.data, function(x) summary(apply(x, 2, var)))
  outdata <- do.call(rbind, outdata)
  outdata <- outdata[samplename, ]
  return(outdata)
}

quiet <- function(..., messages=FALSE, cat=FALSE){
  if(!cat){
    sink(tempfile())
    on.exit(sink())
  }
  out <- if(messages) eval(...) else suppressMessages(eval(...))
  out
}