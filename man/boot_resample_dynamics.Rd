\name{boot_resample_dynamics}
\alias{boot_resample_dynamics}
\title{boot interface to resample_dynamics}
\usage{
boot_resample_dynamics(data, statistic, R, k, ...,
  new.index.name = "resample", weight.name = NULL,
  weight.fun = resample_dynamics_default_weightfun)
}
\arguments{
  \item{data}{time.table the dynamics of which are reused}

  \item{statistic}{the statistic to compute for each new
  dataset}

  \item{R}{number of bootstrap replicates}

  \item{k}{the number of nearest neighbours to resample
  from}

  \item{...}{additonal parameters to \code{boot} (and
  therefore indirectly to \code{statistic})}

  \item{new.index.name}{column name for the new unambiguous
  indices (defaults to "resample", set to NULL not to
  include one)}

  \item{weight.name}{Name of new (auxiliary) column
  containing the weight of the sampled point}

  \item{weight.fun}{weight function for local resampling}
}
\description{
boot interface to resample_dynamics
}

