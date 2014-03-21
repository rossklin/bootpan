\name{resample_series}
\alias{resample_series}
\title{Resample series from a \code{time.table}}
\usage{
resample_series(tt, num, resample = TRUE, new.index.name = "resample")
}
\arguments{
  \item{tt}{\code{time.table} to resample from}

  \item{num}{number of series to resample (defaults to the
  same number as in \code{tt})}

  \item{resample}{whether to sample with replacement
  (default)}

  \item{new.index.name}{name to use for additional index
  column (for disambigating separate resamples of the same
  data, defaults to "resample", set to NULL if new index is
  not to be included)}
}
\description{
Resample whole series from a \code{time.table},
disambiguating names if necessary (while retaining the
types of each column).
}

