\name{sample_markov}
\alias{sample_markov}
\title{Generate panel/time series samples by "Markov" procedure}
\usage{
sample_markov(times, rinitial, rtransition, time.column.name = "time", ...)
}
\arguments{
  \item{times}{vector of timepoints at which to sample}

  \item{rinitial}{function to sample initial
  data.table/values}

  \item{rtransition}{function to sample transitions given a
  data.frame containing current values and the current and
  proceeding time point}

  \item{time.column.name}{column name to use store time
  values in (defaults to "time")}

  \item{...}{additional values to pass to rinitial and
  rtransition}
}
\description{
Generate bootstrap samples given a sampling procedure for
initial values and a procedure for sampling transitions
(conditional only on the current value of the process).
}
\details{
\code{rinitial} should take no paramters not specified in
\code{...} while \code{rtransition} takes three arguments:
the time to transition from, the time to transition to, and
a data.table containing the current set of values, in that
order.

\code{rtransition} *is* allowed to destructively update the
\code{data.table} it is passed, but should even then return
that destructively updated data table.
}

