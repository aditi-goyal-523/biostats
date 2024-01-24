library(epitools)
library(abind)
library(DescTools)

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
