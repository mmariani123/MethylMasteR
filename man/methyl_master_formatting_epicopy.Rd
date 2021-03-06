% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methyl_master_formatting_epicopy.R
\name{methyl_master_formatting_epicopy}
\alias{methyl_master_formatting_epicopy}
\title{methyl_master_formatting_epicopy}
\usage{
methyl_master_formatting_epicopy(
  epi.form.seg = NULL,
  epi.form.output.dir = getwd(),
  epi.form.save.seg = FALSE,
  epi.form.comparison = NULL,
  epi.form.thresholds = NULL
)
}
\arguments{
\item{epi.form.seg}{The input Epicopy seg results}

\item{epi.form.output.dir}{The output dir for the formatted seg object}

\item{epi.form.save.seg}{Whether to save the formatted Epicopy seg results}

\item{epi.form.comparison}{The comparison vector used earlier in the anlysis}

\item{epi.form.thresholds}{The thresholds used to determine the CNV
state, if NULL, the equation  seg.state <- round(2^seg.means * 2) is used}
}
\value{
Formatted epicopy result as a MethylMaster seg object
for downstream analysis
}
\description{
Formatting the Epicopy output into MethylMaster seg object
for downstream analyses
}
