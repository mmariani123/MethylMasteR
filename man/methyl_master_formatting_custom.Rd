% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methyl_master_formatting_custom.R
\name{methyl_master_formatting_custom}
\alias{methyl_master_formatting_custom}
\title{methyl_master_formatting_custom}
\usage{
methyl_master_formatting_custom(
  custom.form.seg = NULL,
  custom.form.output.dir = getwd(),
  custom.form.sample.sheet.path = NULL,
  custom.form.comparison = NULL,
  custom.form.save.seg = FALSE,
  custom.form.thresholds = NULL,
  ...
)
}
\arguments{
\item{custom.form.seg}{The input sesame segmentation list for formatting}

\item{custom.form.output.dir}{The output dir for sesame formatting}

\item{custom.form.sample.sheet.path}{The path to the MethylMaster
sample sheet}

\item{custom.form.comparison}{The MethylMaster comparison vector used in the
custom routine analysis}

\item{custom.form.save.seg}{Whether or not to save the formatted seg output}

\item{custom.form.thresholds}{The thresholds used to determine the CNV
state, if NULL, the equation  seg.state <- round(2^seg.means * 2) is used}

\item{...}{additional parameters to pass to methyl_master_formatting_sesame}
}
\value{
Formatted seg list object for visualizing etc.
}
\description{
Formatting the results of the custom functionality to prepare
for comparison and output
}
