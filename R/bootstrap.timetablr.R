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

#' non-boot interface for resampling series
#'
#' Resample whole series from a \code{time.table}, disambiguating names if
#' necessary (while retaining the types of each column).
#'
#' @param tt \code{time.table} to resample from
#' @param num number of series to resample (defaults to the same number as in \code{tt})
#' @param resample whether to sample with replacement (default)
#' @param new.index.name name to use for additional index column (for disambigating separate resamples of the same data, defaults to "resample", set to NULL if new index is not to be included)
#'
#' @export
resample_series <- function(tt, num=NULL, resample=TRUE, new.index.name="resample") {
    #TODO: Add option to disambiguate existing index columns instead
    resampled <-
        resample_by_cols( tt, index_names(tt), num=num, resample=resample
                        , unique.name = new.index.name )
    if(!is.null(new.index.name)) {
        as.time.table( resampled, new.index.name, time_name(tt), measurement_names(tt)
                     , c(index_names(tt), auxiliary_names(tt)) )
    } else {
        resampled
    }
}

#' boot interface for resampling series
#'
#' @param data time.table to resample from
#' @param statistic statistic to apply to each resampled time.table
#' @param R number of bootstrap replicates
#' @param ... additonal parameters to \code{boot} (and therefore indirectly to \code{statistic})
#' @param include.other whether to also pass (as second argument to \code{statistic}) a time.table containing those entities not in the bootstrapped one
#' @param new.index.name column name to use for the new index (defaults to "resample"), old index columns are stored as auxiliary variables
#' 
#' @export
boot_resample_series <- function(data, statistic, R, ..., include.other=TRUE, new.index.name="resample") {
    require(boot)
    boot.args <- list(...)
    boot.args$data <- unique(index(data))
    boot.args$R <- R
    #
    boot.args$statistic <- if(is.null(boot.args$m)) {
        function(indices, ents, ...) {
            ss <- as.data.table(indices[ents,])[,eval(new.index.name):=.I]
            statistic(promote( subset(data, index=ss)
                             , c("index", rep("auxiliary", length(index_names(data))))
                             , c(new.index.name, index_names(data)) ), ...)
        }
    } else {
        function(indices, ents, other, ...) {
            ss <- as.data.table(indices[ents,])[,eval(new.index.name):=.I]
            ssother <- as.data.table(indices[other,])
            statistic( promote( subset(data, index=ss)
                              , c("index", rep("auxiliary", length(index_names(data))))
                              , c(new.index.name, index_names(data)) )
                     , subset(data, index=ssother), ...)
        }
    }
    #
    result <- do.call(boot, boot.args)
    # print.boot inspects this...
    result$call <- match.call()
    # Otherwise print.boot thinks we're doing case resampling for censored data...
    result$call[[1L]] <- quote(boot)
    result
}

resample_dynamics_default_weightfun <- function(d, h=0.5) {
    expd <- exp(-d/h)
    result <-
        if(sum(expd) < .Machine$double.eps) as.numeric(d == min(d))
        else expd
    result/sum(result)
}

#' Bootstrap from dynamics
#'
#' Bootstrap samples by simulating a random walk with steps sampled (locally)
#' from a procided time.table.
#'
#' @param tt \code{time.table} to bootstrap from
#' @param num number of series to resample (defaults to name number as \code{tt})
#' @param k number of nearest neighbours to resample from
#' @param new.index.name column name for the new unambiguous indices (defaults to "resample", set to NULL not to include one)
#' @param weight.name Name of new (auxiliary) column containing the weight of the sampled point
#' @param weight.fun weight function for local resampling
#' @param frequency list containing starting (from) and stopping (to) times as well as the lenght of the timestep (delta), defaults to using information from \code{tt}. 
#' @param sample.from initial year of \code{tt} to sample from, defaults to the year specified in \code{frequency}.
#' @param resample.steps whether to resample relative (to current position) transitions instead of actual transitions (defaults to TRUE), use FALSE if you want the resampled series to consisit only of data points in the original data set
#' @param ... additional arguments to pass to \code{weight.fun}
#'
#' @details Note that only complete cases from \code{tt} are used. If no
#' "new.index.name" is provided the procedure returns a \code{data.table} rather
#' than a \code{time.table}.
#' 
#' @export
resample_dynamics <- function( tt, num = nrow(unique(index(tt))), k
                             , new.index.name = "resample"
                             , weight.name = NULL
                             , weight.fun = resample_dynamics_default_weightfun
                             , frequency=attr(tt, "frequency")
                             , sample.from=frequency$from
                             , resample.steps=TRUE
                             , ... ) {
    if(inherits(tt[[time_name(tt)]], "numeric"))
        warning("resample_dynamics does not play well with floating point times")
    times <- with(frequency, seq(from, to, delta))
    nms <- safe_name(tt, num=2)
    index.name <- maybe(new.index.name, nms[1])
    tmp.dist.name <- maybe(weight.name, nms[2])
    #
    tt.complete <- subset(tt, expr=complete.cases(measurement(tt)))
    # TODO: Make subset work with only times specified...
    tt.initial <- tt.complete[unlist(time(tt.complete))==sample.from]
    if(nrow(tt.initial) == 0) stop("No complete cases from this timepoint available")
    rinitial <- function() {
        result <- tt.initial[sample.int(nrow(tt.initial),num,TRUE)][,eval(index.name):=.I]
        if(!is.null(weight.name))
            result[,eval(weight.name):=NA]
        result
    }
    #
    rtransition <- if(resample.steps) {
        tt.diff <- diff(tt.complete)
        complete.diff <- complete.cases(measurement(tt.diff))
        tt.complete <- subset(tt.complete, expr=complete.diff)
        tt.diff <- subset(tt.diff, expr=complete.diff)
        stopifnot(all.equal(index(with.time=T, tt.complete), index(with.time=T, tt.diff)))
        ##
        tt.complete <- flanner_by_measurement(tt.complete)
        ##
        function(tfrom, tto, current) {
            # Find neighbouring point row numbers
            neighbours <- knn_lookup_rows(tt.complete, current, k)
            # Pick the step differences corresponding to the neighbouring points
            diffs <-
            # Adjoint distance and index information and key by index
                setkeyv( tt.diff[neighbours][,
                             eval(tmp.dist.name):=attr(neighbours,"distance")][,
                         eval(index.name):=rep(current[[index.name]], each=k)]
                       , eval(index.name) )[,
                       # Apply weight per index
                       eval(tmp.dist.name):=weight.fun(.SD[[tmp.dist.name]]),by=eval(index.name)][,
                       # Sample one row for each index
                      .SD[sample.int(nrow(.SD), 1, prob=.SD[[tmp.dist.name]])],by=eval(index.name)]
            setkeyv(current, index.name)
            # Step the current values
            for(col in measurement_names(tt))
                diffs[,eval(col):=.SD[[col]]+current[[col]]]
            diffs[,colnames(current),with=FALSE]
        }
    } else {
        tt.lag <- timetablr:::lag.time.table(tt.complete)
        complete.lag <- complete.cases(measurement(tt.lag))
        tt.lag <- subset(tt.lag, expr=complete.lag)
        tt.complete <- subset(tt.complete, expr=complete.lag)
        stopifnot(all.equal(index(with.time=T, tt.complete), index(with.time=T, tt.lag)))
        ##
        tt.lag <- flanner_by_measurement(tt.lag)
        ##
        function(tfrom, tto, current) {
            neighbours <- knn_lookup_rows(tt.lag, current, k)
            nxt <- setkeyv( tt[neighbours][,
                         eval(tmp.dist.name):=attr(neighbours,"distance")][,
                         eval(index.name):=rep(current[[index.name]], each=k)]
                       , eval(index.name) )[,
                       # Apply weight per index
                       eval(tmp.dist.name):=weight.fun(.SD[[tmp.dist.name]]),by=eval(index.name)][,
                       # Sample one row for each index
                      .SD[sample.int(nrow(.SD), 1, prob=.SD[[tmp.dist.name]])],by=eval(index.name)]
            nxt[,colnames(current),with=FALSE]
        }
    }
    #
    resampled <- sample_markov(times, rinitial, rtransition, time.column.name=time_name(tt))
    if(is.null(new.index.name)) {
        resampled[,eval(index.name):=NULL]
    } else {
        as.time.table( resampled, index.name, time_name(tt), measurement_names(tt)
                     , c(index_names(tt), auxiliary_names(tt), maybe(weight.name, c())) )
    }
}

#' boot interface to resample_dynamics
#'
#' @param data time.table the dynamics of which are reused
#' @param statistic the statistic to compute for each new dataset
#' @param R number of bootstrap replicates
#' @param k the number of nearest neighbours to resample from
#' @param ... additonal parameters to \code{boot} (and therefore indirectly to \code{statistic})
#' @param new.index.name column name for the new unambiguous indices (defaults to "resample", set to NULL not to include one)
#' @param weight.name Name of new (auxiliary) column containing the weight of the sampled point
#' @param weight.fun weight function for local resampling
#'
#' @export
boot_resample_dynamics <- function( data, statistic, R, k, ...
                                  , new.index.name = "resample"
                                  , weight.name = NULL
                                  , weight.fun = resample_dynamics_default_weightfun ) {
    require(boot)
    #
    boot.args <- list(...)
    boot.args$data <- data
    boot.args$statistic <- statistic
    boot.args$R <- R
    boot.args$sim <- "parametric"
    boot.args$ran.gen <- function(unused1, unused2) {
        resample_dynamics( data, k=k
                         , new.index.name=new.index.name
                         , weight.name=weight.name
                         , weight.fun=weight.fun )
    }
    #
    result <- do.call(boot, boot.args)
    # print.boot inspects this...
    result$call <- match.call()
    # Otherwise print.boot thinks we're doing case resampling for censored data...
    result$call[[1L]] <- quote(boot)
    result
}
