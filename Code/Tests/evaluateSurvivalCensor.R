# Return survival metrics to evaluate the clustering results.
# clinicalData: data frame containing TCGA clinical data
# labels: numeric vector of cluster assignments for the patients

evaluateSurvivalCensor <- function(clinicalData, labels) {
  # Retrieve the patient survival times and death status
  survTime <- clinicalData$days_to_death
  deathStatus <- clinicalData$vital_status == "dead"
  
  # Replace missing survival times with days to last follow up
  missingSurvInd <- is.na(survTime)
  lastFollowup <- clinicalData$days_to_last_followup
  survTime[missingSurvInd] <- lastFollowup[missingSurvInd]
  
  # Calculate the p-value of the log-rank test and the concordance
  # index of the Cox proportional hazards model
  survObject <- Surv(survTime, deathStatus)
  survDiff <- survdiff(survObject~labels)
  pval <- 1 - pchisq(survDiff$chisq, length(survDiff$n)-1)
  coxphFit <- coxph(survObject~labels)
  ci <- summary(coxphFit)$concordance[1]
  results <- c(pval, ci)
  
  return(results)
}
