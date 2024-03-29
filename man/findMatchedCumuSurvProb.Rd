% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oneSampleLogRankTest.R
\name{findMatchedCumuSurvProb}
\alias{findMatchedCumuSurvProb}
\title{Find Matched Cumulative Survival Probability}
\usage{
findMatchedCumuSurvProb(time, ageDiag, sex, race, dataPop, maxFollowUp = NULL)
}
\arguments{
\item{time}{follow up length}

\item{ageDiag}{age at diagnosis}

\item{sex}{sex}

\item{race}{race}

\item{dataPop}{Population level mortality data}

\item{maxFollowUp}{maximum follow-up, if max follow-up not provided then the
time would be considered until death or censoring}
}
\value{
matched survival probability
}
\description{
Find Matched Cumulative Survival Probability
}
\examples{
# load data
data(dataSurv_small)
data(dataPop_2018_2021)

# Extract info for the first subject
time_vec <- dataSurv_small$time[1]
age_vec <- dataSurv_small$age[1]
sex_vec <- dataSurv_small$sex[1]
race_vec <- dataSurv_small$race[1]

# Generate cumulative survival probability
findMatchedCumuSurvProb(time = time_vec, ageDiag = age_vec, sex = sex_vec,
race = race_vec, dataPop = dataPop_2018_2021)

#If maximum followup is determined to be 20 years
findMatchedCumuSurvProb(time = time_vec, ageDiag = age_vec, sex = sex_vec,
race = race_vec, dataPop = dataPop_2018_2021, maxFollowUp = 20)

}
