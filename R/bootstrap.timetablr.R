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

#' Resample series from a \code{time.table}
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
                             , ... ) {
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
    tt.diff <- diff(tt.complete)
    complete.diff <- complete.cases(measurement(tt.diff))
    tt.complete <- subset(tt.complete, expr=complete.diff)
    tt.diff <- subset(tt.diff, expr=complete.diff)
    stopifnot(all.equal(index(with.time=T, tt.complete), index(with.time=T, tt.diff)))
    #
    tt.complete <- flanner_by_measurement(tt.complete)
    #
    rtransition <- function(tfrom, tto, current) {
        # Find neighbouring point row numbers
        neighbours <- knn_lookup_rows(tt.complete, current, k)
        neighbours
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
    #
    resampled <- sample_markov(times, rinitial, rtransition, time.column.name=time_name(tt))
    if(is.null(new.index.name)) {
        resampled[,eval(index.name):=NULL]
    } else {
        as.time.table( resampled, index.name, time_name(tt), measurement_names(tt)
                     , c(index_names(tt), auxiliary_names(tt), maybe(weight.name, c())) )
    }
}
