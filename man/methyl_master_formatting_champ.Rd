% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methyl_master_formatting_champ.R
\name{methyl_master_formatting_champ}
\alias{methyl_master_formatting_champ}
\title{methyl_master_formatting_champ}
\usage{
methyl_master_formatting_champ(
  champ.form.seg = NULL,
  champ.form.output.dir = getwd(),
  champ.form.save.seg = FALSE,
  champ.form.comparison = NULL,
  champ.form.padj = NULL,
  champ.form.thresholds
)
}
\arguments{
\item{champ.form.seg}{The input ChAMP CNV segmentation results}

\item{champ.form.output.dir}{The output directory
used}

\item{champ.form.save.seg}{Whether to save the formatted ChAMP results}

\item{champ.form.comparison}{The MethylMaster comparison vector}

\item{champ.form.padj}{The padj that was used in the ChAMP routine}

\item{champ.form.thresholds}{The thresholds used to determine the CNV
state, if NULL, the equation  seg.state <- round(2^seg.means * 2) is used}

\item{...}{Additional parameters passed to methyl_master_champ_formatting}
}
\value{
Formatted champ result as MethylMaster seg object
for downstream analysis
}
\description{
formatting the champ seg results for
overlaps and visualization etc.a
Michael Mariani PhD Dartmouth College 2021
}
