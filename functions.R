library(tidyverse)
library(epitools)
library(abind)
library(DescTools)
library(logistf)
library(pROC)
library(psych)
library(lmtest)

make_table=function(dats, nrow=2, ncol=2, col_names=c("case", "control"), row_names=c("outcome", "no outcome"), include_totals=FALSE){
  
  suppressWarnings({ 
    data = matrix(dats, nrow=nrow, ncol=ncol, byrow=TRUE)
    table = as.table(data)
    colnames(table)=col_names
    rownames(table)=row_names
    
    if (include_totals == TRUE){
      rowsums=rowSums(table)
      table=cbind(table, rowsums)
      colsums=colSums(table)
      table=rbind(table, rowsums)
      colnames(table) = c(col_names, 'Totals')
      rownames(table) = c(row_names, 'Totals')
    }
    
    print(table)
    return(table)
  }) 
}

chisquared=function(t){
  return(chisq.test(t))
}

fisher_exact=function(t){
  return(fisher.test(t))
}

crude_odd_ratio=function(t){
  return(epitools::oddsratio(t))
}

crude_risk_ratio=function(t){
  return(epitools::riskratio(t, rev="both"))
}

mcnemars_test=function(t){
  return(mcnemar.test(t))
}

mh_test=function(...){
  # take in a list of 2x2 tables and calculate Cochran-Mantel-Haenszel Test
  dat=abind::abind(..., along=3)
  return(mantelhaen.test(dat))
}

breslow_day_test=function(...){
  # take in a list of 2x2 tables and calculate the Breslow-Day test of homogeneity
  dat=abind::abind(..., along=3)
  return(DescTools::BreslowDayTest(dat))
}


breslowday.test = function(x, OR=NA, printORi.s=TRUE){
  ## function to compute the Breslow Day test of homogeneity for a 2 by 2 by k table
  ## x is a three dim array, 2x2xk
  ## tests to see if all strata have the same OR
  ## if OR is not given, the Mantel-Haenszel estimate is used.
  if(is.na(OR)) {
    OR = mantelhaen.test(x)$estimate
    names(OR) = ""
  } 
  OR.i <- apply(x, 3,  function(x) x[1,1] * x[2,2] / x[1,2] /x[2,1])
  k = dim(x)[3]
  n11k = x[1,1,]
  n21k = x[2,1,]
  n12k = x[1,2,]
  n22k = x[2,2,]
  row1sums = n11k + n12k
  row2sums = n21k + n22k
  col1sums = n11k + n21k
  Amax = apply(cbind(row1sums,col1sums),1,min)
  ## Astar must be no more than col1sums and no more than row1sums
  bb = row2sums +row1sums * OR - col1sums*(1-OR)
  determ = sqrt(bb^2 + 4 * (1-OR) *  OR * row1sums * col1sums)
  Astar = (-bb + cbind( -determ, determ))/ (2 -2*OR)
  Astar = ifelse(Astar[,1] <= Amax & Astar[,1] >= 0, Astar[,1], Astar[,2])
  ## print(Astar)
  Bstar = row1sums - Astar
  Cstar = col1sums - Astar
  Dstar = row2sums - col1sums + Astar
  Var = apply(1 / (.5+cbind(Astar,Bstar,Cstar,Dstar)), 1, sum)^(-1)
  ## print(Var)
  X2 = sum( (x[1,1,] - Astar)^2/Var )
  pvalue = 1 - pchisq(X2,k-1)
  if(printORi.s) {
    out <- rbind(log(OR.i), 1/Var)
    dimnames(out)[[1]] <- c("log OR","Weight")
    print(out)
  }
  return(unlist(list(OR = OR, Stat = X2, df = k-1, pvalue = pvalue)))
}

ORfunction = function(confidence,a,b,c,d) {
  #enter confidence percent as a whole number, e.g. "95"
  
  OR = (a*d) / (b*c)
  
  lnOR = log(OR)
  error = sqrt(1/a + 1/b + 1/c + 1/d)
  Z = -qnorm((1 - confidence/100)/2)
    #gives left hand Z score, multiply by negative
  
  lower = exp(lnOR - Z*error)
  upper = exp(lnOR + Z*error)
  
  output = c(OR, lower, upper)
  names(output) = c("OR", "lower", "upper")
  return(output)
}

logit.plot = function(data,predictor,outcome,numIntervals=10,spline=TRUE,spar=0.7) {
  #Check inputs
  if (!is.character(predictor)) {stop("The predictor must be a string.")}
  if (!is.character(outcome)) {stop("The outcome must be a string.")}
  
  #Initialize vectors
  logits = rep(0,numIntervals)
  means = rep(0,numIntervals)
  
  #Define outcome and predictor columns (renames them to make the code below simpler)
  outcomeIndex = which(colnames(data)==outcome)
  predictorIndex = which(colnames(data)==predictor)
  names(data)[outcomeIndex] = "outcome"
  names(data)[predictorIndex] = "predictor"
  
  #Check for missing data
  if (any(is.na(data$predictor))) {stop("The predictor contains missing data.")}
  if (any(is.na(data$outcome))) {stop("The outcome contains missing data.")}
  
  #Sort data in ascending order of predictor
  data=data[order(data$predictor),]
  
  #Define the intervals
  intervalLength = floor(dim(data)[1] / numIntervals)
  #floor() enforces rounding down, in case numIntervals doesn't divide evenly into dim(data)[1]
  
  if ((dim(data)[1] / numIntervals) - floor(dim(data)[1] / numIntervals) != 0) {
    warning("The number of intervals does not divide evenly into the size of the dataset. Some data points will be omitted.")
  }
  
  #Define the starting index of each section
  intervalStarts = c()
  for (k in 1:numIntervals) {
    intervalStarts[k] = 1 + (k-1)*intervalLength
  }
  
  #Loop over each section
  for (j in intervalStarts) {
    positive=0
    negative=0
    sum=0
    
    #Loop over each data point in the section
    for (i in c(j:(j+intervalLength-1))) {
      if (data$outcome[i] == 1) {             #outcome variable is 1
        positive = positive + 1
      } else if (data$outcome[i] == 0) {      #outcome variable is 0
        negative = negative + 1
      } else {
        stop("The outcome column must be binary.")
      }
      
      sum = sum + data$predictor[i]   #adding the predictor values so we can calculate the mean
    }
    
    x = positive/(positive+negative)        #calculates the probability
    logits[1+round(j/intervalLength)] = log(x/(1-x))        #puts the logit in the index of j
    means[1+round(j/intervalLength)] = sum/intervalLength   #puts the mean in the index of j 
    
  }
  
  #Check for infinite logits
  if (any(is.infinite(logits))) {
    warning("Infinite logits generated. Not all points will be plotted.")
  }
  
  #Plot each logit at the corresponding predictor value
  plot(means,logits,xlab=paste(predictor))
  
  #Plot a smoothing spline
  if (spline==TRUE) {
    

    spline = smooth.spline(means[!is.infinite(logits)], logits[!is.infinite(logits)], spar=spar)
    lines(spline)  
  }
}

exact_logreg=function(formula,dataset){
  return(logistf(formula, family = "binomial", dataset)
  )
}

roc.curve = function(glm.fit,data,outcome,returnValues=FALSE,returnPred=FALSE,returnROC=FALSE,weights=NULL) {
  #Make sure pROC package is present
  if (!"package:pROC" %in% search()) {stop("This function requires the pROC package.")}
  
  #Check input
  if (!is.character(outcome)) {stop("Outcome must be a string.")}
  if (!is.null(weights) & !is.character(weights)) {stop("Weights argument must be a string.")}
  if (sum(returnValues,returnPred,returnROC)>1) {stop("Select one output at a time.")}
  
  #Define outcome column by renaming it
  outcomeIndex = which(colnames(data)==outcome)
  names(data)[outcomeIndex] = "outcome"
  
  #Make prediction based on the fit and add to data frame
  prediction = predict(glm.fit,type=c("response"))
  data = data.frame(data,prediction)
  
  #Multiply rows based on weights variable if present
  if (!is.null(weights)) {
    data = data[rep(row.names(data), unlist(data[weights])), -(names(data)==weights)]
  }
  
  #Calculate ROC curve
  curve = roc(outcome~prediction,data=data)
  plot(curve)    
  print(curve)   #Summarizes the result, including the area under the curve
  #By default, the plot function labels the x-axis as "Specificity" (from 1 to 0) instead of "1-specificity" from 0 to 1
  
  #Output, if requested
  if (returnValues==TRUE) {
    output = data.frame(curve$thresholds,curve$sensitivities,curve$specificities)
    names(output) = c("threshold","sensitivity","specificity")
    return(round(output,digits=10))   #rounding is to prevent display as scientific notation
  }
  
  if (returnPred==TRUE) {
    #For confidence intervals: calculate the predictions on the linear ("link") scale
    linkPred = predict(glm.fit, type=c("link"), se.fit=TRUE)
    
    #Calculate the 95% confidence intervals, then convert back to the logistic scale using the inverse logit function
    lower.limit = linkPred$fit - 1.96 * linkPred$se.fit
    lower.limit = exp(lower.limit)/(1+exp(lower.limit))
    upper.limit = linkPred$fit + 1.96 * linkPred$se.fit
    upper.limit = exp(upper.limit)/(1+exp(upper.limit))
    
    output = data.frame(prediction, lower.limit, upper.limit)
    return(round(output,digits=10))   #rounding is to prevent display as scientific notation
  }
  
  if (returnROC==TRUE) {return(curve)}
}


plot_ROC <- function(outcomes, pred_probs){
  # takes in a list of outcomes as 0/1 and a list of predicted probabilities from the model
  # e.g. outcomes=c(1,1,1,0,0,0,0,0), pred_probs=c(0.759, 0.759, 0.209, 0.154, 0.139, 0.685, 0.163, 0.132)
  # outputs an ROC curve and the AUC
  auc=pROC::roc(outcomes, pred_probs, plot=TRUE, auc=TRUE, legacy.axes = TRUE)$auc
  return(auc)
}


standardize = function(data){
  return(scale(data))
}

correlations = function(data){
  return(cor(data, use = "complete.obs"))
}

calculate_pca = function(data, no_pcs, adjustment = "varimax"){
  return(psych::principal(data, nfactors = no_pcs, rotate = adjustment))
}

scree_plot = function(pca_df){
  # The input is the really any dataframe, but it can take the output from the "calculate_pca" function
  return(plot(pca_df$values/sum(pca_df$values),
              xlab="Prinical component",
              ylab="Proportion of variance explained"))
}

pc_values = function(data, pca_df, chosen_pc){
  # The goal of this function is to take the original dataframe and previously calculated PCs
  # And then calculates individuals PC values
  factor = apply(data, chosen_pc, function(x) sum(x*pca_df$weights[,chosen_pc]))
  return(factor)
}

likelihood_ratio = function(glm_model){
  # This takes a glm model and runs the likelihood ratio on that model
  return(lrtest(glm_model))
}


propensity_scores = function(data, primary_predictor, covariates){
  # The function calculates the propensity score for for a primary predictor given a vector of covariate names

  # Build formula for propensity score 
  # First combine the vector of covariates
  # Then include the primary predictor as the outcome
  formula = str_c(covariates, collapse = "+") %>% 
    str_c(primary_predictor,"~",. )
  
  # now run the model then get the probabilities
  prop_model = glm(formula,
                   family = "binomial",
                   data)
  pscores = predict(prop_model, type = "response")
  
  return(pscores)
}