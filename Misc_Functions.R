# This script contains miscellaneous functions for interpreting 
# some of the more obscure categorical regression models


## Lecture 3 (Ordinal Logit Models)

polrLogitStdCoef <- function(mod, help = FALSE, digits = 4){
  # Like Stata's listcoef for ordinal logit in R
  # mod is a polr object
  # Replicates all beta's in listcoef except bStdY, and bStdXY
  # See the math here: https://www3.nd.edu/~rwilliam/stats3/L04.pdf
  
  ## z score
  zScore <- summary(mod)$coefficients[1:(ncol(mod$model)-1), 3]
  
  ## p_value (two tailed)
  mod.coef <-data.frame(coef(summary(mod)))[1:(ncol(mod$model))-1,]
  pVal <- round((pnorm(abs(mod.coef$t.value), lower.tail= FALSE) * 2),2)
  
  ## bStdX 
  bStdX <- mod$coef * apply(mod$model[, 2:(ncol(mod$model))], 2, sd)
  
  ## bStdY
  ## bStdY <- mod$coef/sd(as.numeric(mod$model[,1])) ## need latent dep var
  
  ## SDofX
  SDofX <- apply(mod$model[, 2:(ncol(mod$model))], 2, sd)
  
  
  outputDF <- data.frame(b.LogOdds = coef(mod), b.OddsRatio = exp(coef(mod)) ,z = zScore, p_value = pVal,
                         bStdX = bStdX, SDofX = SDofX)
  outputDF <- as.data.frame(lapply(outputDF, round, digits), row.names = row.names(outputDF)) # ROUND output
  
  helpText <- "b.LogOdds = raw coefficient in Log Odds
  b.OddsRatio = raw coefficient in Log Odds
  z = z-score for test of b=0
  P>|z| = p-value for z-test
  bStdX = x-standardized coefficient
  SDofX = standard deviation of X\n\n\n"
  
  if(help == "FALSE"){return(outputDF)} 
  else if(help == "TRUE"){cat(helpText) ; return(outputDF)}
  
}



## Lecture 4 (Multinomial Models)

mlogTestR <- function(unrestrictedModel, reflevel, type = "wald", digits = 4){
  # Inputs: unrestrictedmodel is a mlogit object
  #         reflevel is the reference level of the dependent variable
  #         type default is the wald test, It can be changed to "lr" i.e. likelihood ratio
  #         Rounded to 4 decimal places by default.
  # Output: Datafram containing the results of the wald/lr test similar to Stata's mlogtest, wald lr
  depVarList <- gsub(" \\+", ",", terms(unrestrictedModel)[[3]][3]) # get dep vars
  depVarList <- unlist(strsplit(depVarList, ", "))
  
  resultDF1 <- data.frame(chi2 = numeric(), df = numeric(), Pval = numeric())
  
  for (i in 1:length(depVarList)){
    restrictedFormula = gsub(" ", " + ", paste0(depVarList[-i], "",collapse = " "))
    restrictedFormula = paste("1 | (", restrictedFormula, ")", sep = "")
    restrictedFormula = formula(paste(terms(unrestrictedModel)[[2]], restrictedFormula, sep = " ~ "))
    restrictedModel = mlogit(restrictedFormula, data = unrestrictedModel$model, reflevel = reflevel)
    
    if(type == "wald"){
      waldDF <-  waldtest(unrestrictedModel, restrictedModel, test = "Chisq")
      resultDF1 <- rbind(resultDF1, c(waldDF[2,3], abs(waldDF[2,2]), waldDF[2,4]))
    }
    else if(type == "lr"){
      lrDF <- lrtest(unrestrictedModel, restrictedModel)
      resultDF1 <- rbind(resultDF1, c(lrDF[2,4], abs(lrDF[2,3]), lrDF[2,5]))
    }
  }
  
  if(type == "wald"){
    cat(paste("Wald tests for independent variables (N=", length(unrestrictedModel$fitted.values),")\n",sep = ""))
    cat("Ho: All coefficients associated with given variable(s) are 0.\n\n")
  }
  else if(type == "lr"){
    cat(paste("Likelihood-ratio tests for independent variables (N=", length(unrestrictedModel$fitted.values),")\n",sep = ""))
    cat("Ho: All coefficients associated with given variable(s) are 0.\n\n")
  }
  
  names(resultDF1) <- c("chi2", "df", "P>chi2")
  resultDF1 <- as.data.frame(resultDF1, row.names = depVarList)
  resultDF1 <- as.data.frame(lapply(resultDF1, round, digits), row.names = row.names(resultDF1))
  
  # Output
  return(resultDF1)
  
}
