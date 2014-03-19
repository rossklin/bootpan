## Software License Agreement (BSD License)
##
## Copyright (c) 2014, Tilo Wiklund (tilo@wiklund.co)
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
##
##     Redistributions of source code must retain the above copyright
##     notice, this list of conditions and the following disclaimer.
##
##     Redistributions in binary form must reproduce the above copyright
##     notice, this list of conditions and the following disclaimer in
##     the documentation and/or other materials provided with the
##     distribution.
##
##     The names of its contributors may not be used to endorse or promote products
##     derived from this software without specific prior written permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#' Resample by a number of columns
#'
#' Resample observations from a \code{data.frame} by resampling from the unique
#' values of a number of columns
#'
#' @param dt \code{data.table} to resample from
#' @param cols columns to resample from (defaults to all columns of \code{dt})
#' @param num number of *keys* to resample (defaults to number of unique values in \code{tt})
#' @param resample whether to sample with replacement (default)
#' @param unique.name column name to use for disambiguating resampled values (set to NULL, default, not to include one)
#'
#' @details Gives a \code{data.table} with \code{num} rows and
#' \code{length(cols)} columns, with resampled values from \code{dt}.
resample_by_cols <- function( dt, cols=colnames(dt), num=NULL
                            , resample=TRUE
                            , unique.name = NULL ) {
    available <- if(setequal(cols, key(dt))) {
        unique(dt)
    } else if(setequal(cols, colnames(cols))) {
        unique(dt)
    } else {
        unique(dt[,cols,with=FALSE])
    }
    #
    navailable <- nrow(available)
    num <- maybe(num, navailable)
    result.indices <- sample.int(navailable, num, resample)
    #
    result <- available[result.indices]
    if(!is.null(unique.name))
        result[,eval(unique.name):=.I]
    result
}

#' Generate panel/time series samples by "Markov" procedure
#'
#' Generate bootstrap samples given a sampling procedure for initial values and
#' a procedure for sampling transitions (conditional only on the current value
#' of the process).
#'
#' @param times vector of timepoints at which to sample
#' @param rinitial function to sample initial data.table/values
#' @param rtransition function to sample transitions given a data.frame containing current values and the current and proceeding time point
#' @param time.column.name column name to use store time values in (defaults to "time")
#' @param ... additional values to pass to rinitial and rtransition
#' 
#' @details \code{rinitial} should take no paramters not specified in \code{...}
#' while \code{rtransition} takes three arguments: the time to transition from,
#' the time to transition to, and a data.table containing the current set of
#' values, in that order.
#'
#' \code{rtransition} *is* allowed to destructively update the \code{data.table}
#' it is passed, but should even then return that destructively updated data
#' table.
#'
#' @export
sample_markov <- function( times, rinitial, rtransition
                         , time.column.name = "time", ... ) {
    first.frame <- rinitial(...)
    current.frame <- first.frame
    #
    frames <- list(copy(current.frame)[,eval(time.column.name):=times[1]])
    for(i in seq_along(times)[-1]) {
        current.frame <- rtransition(times[i-1], times[i], current.frame, ...)
        frames[[i]] <- copy(current.frame)[,eval(time.column.name):=times[i]]
    }
    #
    do.call(rbind, frames)
}
