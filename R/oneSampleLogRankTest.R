#' @title Find Matched Cumulative Survival Probability
#' 
#' @param time follow up length
#' @param ageDiag age at diagnosis
#' @param sex sex
#' @param race race
#' @param dataPop Population level mortality data
#' @param maxFollowUp maximum follow-up, if max follow-up not provided then the 
#' time would be considered until death or censoring
#' 
#' @return matched survival probability
#' @examples 
#' # load data
#' data(dataSurv_small)
#' data(dataPop_2018_2021)
#' 
#' # Extract info for the first subject
#' time_vec <- dataSurv_small$time[1]
#' age_vec <- dataSurv_small$age[1]
#' sex_vec <- dataSurv_small$sex[1]
#' race_vec <- dataSurv_small$race[1]
#' 
#' # Generate cumulative survival probability
#' findMatchedCumuSurvProb(time = time_vec, ageDiag = age_vec, sex = sex_vec,
#' race = race_vec, dataPop = dataPop_2018_2021)
#'
#' #If maximum followup is determined to be 20 years
#' findMatchedCumuSurvProb(time = time_vec, ageDiag = age_vec, sex = sex_vec,
#' race = race_vec, dataPop = dataPop_2018_2021, maxFollowUp = 20)
#' 
#' @import dplyr
#' @export
#' @importFrom rlang .data 
#' @importFrom magrittr %>% 
#'
findMatchedCumuSurvProb <- function(time, ageDiag, sex, race, dataPop, maxFollowUp = NULL) {
  
  colSexRace <- paste0(race, '_', sex)
  if (!colSexRace %in% colnames(dataPop))
    stop('Selected sex-race column not found in population data!')
  
  if (is.null(maxFollowUp)) {
    
    time_floor <- floor(time)
    time_vec <- rep(1, time_floor)
    time_vec <- append(time_vec, time - time_floor)
    
    ageEnd <- ageDiag + time
    dataPop %>%
      filter(!!as.symbol("age") >= ageDiag, !!as.symbol("age") <= ageEnd) %>%
      select(.data[[colSexRace]]) %>%
      transmute(.data[[colSexRace]] * time_vec) %>%
      # divided by 100k because population deta rate is per 100 k
      sum() / 100000
  }
  else {
    ageEnd <- ageDiag + maxFollowUp
    dataPop %>%
      filter(!!as.symbol("age") >= ageDiag, !!as.symbol("age") < ageEnd) %>%
      pull(.data[[colSexRace]]) %>%
      # divided by 100k because population deta rate is per 100 k
      cumsum() / 100000
  }
}


#' @title Calculate One-Sample Log-Rank Test
#' 
#' @param dataSurv Survival data
#' @param dataPop Population data
#' @param type Type of test
#' 
#' @return p-value for one-sample log-rank test
#' 
#' @examples 
#' # load data
#' data(dataSurv_small)
#' data(dataPop_2018_2021)
#' 
#' # Since the dataset is small run an exact test
#' oneSampleLogRankTest(dataSurv_small, dataPop_2018_2021, type = "exact")
#' 
#' @export
#'
oneSampleLogRankTest <- function(dataSurv, dataPop, type = c("exact", "approximate")) {
  
  type <- match.arg(type)
  
  totalExpected <- 0
  for (i in 1:nrow(dataSurv)) {
    totalExpected <- totalExpected +  findMatchedCumuSurvProb(dataSurv$time[i],
                                                              dataSurv$age[i],
                                                              dataSurv$sex[i],
                                                              dataSurv$race[i],
                                                              dataPop)
  } 
  totalObserved <- sum(dataSurv$status)
  
  
  std_mort_ratio <- totalObserved / totalExpected
  
  
  if (type == "exact") {
    if (totalExpected > totalObserved){
      pvalue <- 2 * (1 - stats::pchisq(2 * totalExpected, df = (2 * totalObserved) + 2))
    }
    
    if (totalExpected < totalObserved){
      pvalue <- 2 * (stats::pchisq(2 * totalExpected, df = (2 * totalObserved)))
      
    }
    
    chisq_lower <- stats::qchisq(0.025, (2 * totalObserved))
    chisq_upper <- stats::qchisq(0.975, (2 * totalObserved) + 2)
    
    smr_lower <- (0.5 * chisq_lower) / totalExpected
    smr_upper <- (0.5 * chisq_upper) / totalExpected
    
    conf_table <- data.frame("std_mort_ratio_est" = std_mort_ratio,
                             "lwr" = smr_lower,
                             "upr" = smr_upper)
    
    
  }
  
  else {
    pvalue = 1 - stats::pchisq((totalObserved - totalExpected) ^ 2 / totalExpected, df=1)
    
    # assumption is alpha = 0.05, hence 95% CI; could be changed to user defined input
    crit_value <- stats::qchisq(0.975, 1)
    conf_band_term1 <- (crit_value/(2 * totalExpected)) 
    conf_band_term2 <- sqrt(crit_value * (4 * totalObserved + crit_value)) / (2 * totalExpected)
    
    conf_table <- data.frame("std_mort_ratio_est" = std_mort_ratio,
                             "lwr" = std_mort_ratio + conf_band_term1 - conf_band_term2,
                             "upr" = std_mort_ratio + conf_band_term1 + conf_band_term2)
  }
  
  
  
  
  test_obj <- list("p.value" = pvalue, "estimate"  = conf_table)
  return(test_obj)
}



#' @title Plot Kaplan-Meier Curve against Population
#'
#' @param dataSurv Survival data
#' @param dataPop Population data
#' @param type Type of test to conduct in order to display p-value
#' 
#' @return ggplot object
#' 
#' @importFrom survival survfit Surv
#' @importFrom ggplot2 geom_line aes
#' @importFrom survminer ggsurvplot
#' @importFrom stats median
#' 
#' @examples 
#' # load data
#' data(dataSurv_small)
#' data(dataPop_2018_2021)
#' 
#' plotKM(dataSurv_small, dataPop_2018_2021, type = "exact")
#' 
#' @export
#'
plotKM <- function(dataSurv, dataPop, type = c("exact", "approximate")) {
  
  type <- match.arg(type)
  
  n <- nrow(dataSurv)
  maxFollowUp <- max(dataSurv$time)
  dataCumuHaz <- NULL
  for (i in 1:n) {
    dataCumuHaz <- cbind(dataCumuHaz,
                         findMatchedCumuSurvProb(dataSurv$time[i],
                                                 dataSurv$age[i],
                                                 dataSurv$sex[i],
                                                 dataSurv$race[i],
                                                 dataPop,
                                                 maxFollowUp)
                         )
  }
  dataCumuSurvProb <- exp(-dataCumuHaz)
  curveDataPop <- data.frame(
    avgSurvProb = c(1, apply(dataCumuSurvProb, 1, mean)),
    time = 0:ceiling(maxFollowUp)
  )
  
  if (type == "exact"){
    test_output <- oneSampleLogRankTest(dataSurv, dataPop, "exact")
    
  } else{
    test_output <- oneSampleLogRankTest(dataSurv, dataPop, "approximate")
  }
  
  
  fit <- survival::survfit(survival::Surv(time, status) ~ 1, data=dataSurv)
  x_pos <- stats::median(fit$time)
  p_val <- round(test_output$p.value, 4)
  survminer::ggsurvplot(fit, legend.title = "", legend.labs = "Study Sample", conf.int = TRUE)$plot + 
    ggplot2::geom_line(ggplot2::aes(x = .data[["time"]], y = .data[["avgSurvProb"]], 
                       colour = "Demographically-matched U.S. population"), 
                       data = curveDataPop) + 
    ggplot2::annotate("text", x = x_pos, y = 0, label = paste0(type, " test p-value: ", p_val), size = 6) +
    ggplot2::theme(legend.text = ggplot2::element_text(size=16))
}









