\name{resample_by_cols}
\alias{resample_by_cols}
\title{Resample by a number of columns}
\usage{
resample_by_cols(dt, cols = colnames(dt), num = NULL, resample = TRUE,
  unique.name = NULL)
}
\arguments{
  \item{dt}{\code{data.table} to resample from}

  \item{cols}{columns to resample from (defaults to all
  columns of \code{dt})}

  \item{num}{number of *keys* to resample (defaults to
  number of unique values in \code{tt})}

  \item{resample}{whether to sample with replacement
  (default)}

  \item{unique.name}{column name to use for disambiguating
  resampled values (set to NULL, default, not to include
  one)}
}
\description{
Resample observations from a \code{data.frame} by
resampling from the unique values of a number of columns
}
\details{
Gives a \code{data.table} with \code{num} rows and
\code{length(cols)} columns, with resampled values from
\code{dt}.
}

