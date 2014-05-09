\name{resample_dynamics}
\alias{resample_dynamics}
\title{Bootstrap from dynamics}
\usage{
resample_dynamics(tt, num = nrow(unique(index(tt))), k,
  new.index.name = "resample", weight.name = NULL,
  weight.fun = resample_dynamics_default_weightfun, frequency = attr(tt,
  "frequency"), sample.from = frequency$from, resample.steps = TRUE, ...)
}
\arguments{
  \item{tt}{\code{time.table} to bootstrap from}

  \item{num}{number of series to resample (defaults to name
  number as \code{tt})}

  \item{k}{number of nearest neighbours to resample from}

  \item{new.index.name}{column name for the new unambiguous
  indices (defaults to "resample", set to NULL not to
  include one)}

  \item{weight.name}{Name of new (auxiliary) column
  containing the weight of the sampled point}

  \item{weight.fun}{weight function for local resampling}

  \item{frequency}{list containing starting (from) and
  stopping (to) times as well as the lenght of the timestep
  (delta), defaults to using information from \code{tt}.}

  \item{sample.from}{initial year of \code{tt} to sample
  from, defaults to the year specified in
  \code{frequency}.}

  \item{resample.steps}{whether to resample relative (to
  current position) transitions instead of actual
  transitions (defaults to TRUE), use FALSE if you want the
  resampled series to consisit only of data points in the
  original data set}

  \item{...}{additional arguments to pass to
  \code{weight.fun}}
}
\description{
Bootstrap samples by simulating a random walk with steps
sampled (locally) from a procided time.table.
}
\details{
Note that only complete cases from \code{tt} are used. If
no "new.index.name" is provided the procedure returns a
\code{data.table} rather than a \code{time.table}.
}

