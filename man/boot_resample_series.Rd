\name{boot_resample_series}
\alias{boot_resample_series}
\title{boot interface for resampling series}
\usage{
boot_resample_series(data, statistic, R, ..., include.other = TRUE,
  new.index.name = "resample")
}
\arguments{
  \item{data}{time.table to resample from}

  \item{statistic}{statistic to apply to each resampled
  time.table}

  \item{R}{number of bootstrap replicates}

  \item{...}{additonal parameters to \code{boot} (and
  therefore indirectly to \code{statistic})}

  \item{include.other}{whether to also pass (as second
  argument to \code{statistic}) a time.table containing
  those entities not in the bootstrapped one}

  \item{new.index.name}{column name to use for the new
  index (defaults to "resample"), old index columns are
  stored as auxiliary variables}
}
\description{
boot interface for resampling series
}

