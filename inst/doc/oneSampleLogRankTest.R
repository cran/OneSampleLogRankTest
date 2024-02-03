## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(OneSampleLogRankTest)

## -----------------------------------------------------------------------------
data(dataSurv, package = "OneSampleLogRankTest")
data(dataPop_2018_2021_race_sex_eth, package = "OneSampleLogRankTest")

## approximate one sample log rank test
oneSampleLogRankTest(dataSurv, dataPop_2018_2021_race_sex_eth,
                     type = "approximate")



## ---- fig.height=6, fig.width=8-----------------------------------------------

plotKM(dataSurv, dataPop_2018_2021_race_sex_eth, type = "approximate")


## -----------------------------------------------------------------------------

data(dataSurv_small, package = "OneSampleLogRankTest")

## exact test version of one sample log rank test for small sample size
oneSampleLogRankTest(dataSurv_small, dataPop_2018_2021_race_sex_eth)


## ---- fig.height=6, fig.width=8-----------------------------------------------

plotKM(dataSurv_small, dataPop_2018_2021_race_sex_eth, type = "exact")


## -----------------------------------------------------------------------------
data(simulated_clinical_data)

oneSampleLogRankTest(simulated_clinical_data, dataPop_2018_2021_race_sex_eth,
                     type = "approximate")


## ---- fig.height=6, fig.width=8-----------------------------------------------
plotKM(simulated_clinical_data, dataPop_2018_2021_race_sex_eth, type = "approximate")


## -----------------------------------------------------------------------------
sessionInfo()

