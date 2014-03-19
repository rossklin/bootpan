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
#'
#' @details Gives a \code{data.table} with \code{num} rows and
#' \code{length(cols)} columns, with resampled values from \code{dt}.
resample_by_cols <- function(dt, cols=colnames(dt), num=NULL, resample=TRUE) {
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
    available[result.indices]
}
