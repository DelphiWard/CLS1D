##' Compute k-nearest neighbour predictions for the last dimension of
##' an embedding.
##'
##' This function is based on the \code{\link[RANN]{nn2}} function from \pkg{RANN}.
##' Given an array representing an embedding, it constructs inverse
##' distance weighted k-nearest neighbour predictions for the last
##' column of the embedded array.
##'
##' @title K-Nearest Neighbour Prediction
##' @param emb An embedding constructed with \code{\link{window50}}.
##' @param k The number of neighbours.
##' @return Returns a list with two components
##' \item{\code{obs}}{the observed responses.}
##' \item{\code{pred}}{the predicted responses.}
##' @seealso \code{\link{window50}},\code{\link{embed50}},\code{\link[RANN]{nn2}}
##' @importFrom RANN nn2
knn <- function(emb,k=7) {
  X <- emb[,-ncol(emb),drop=F]
  y <- emb[,ncol(emb)]
  ## Find the k+1 nearest neighbours
  nn <- nn2(as.matrix(X),k=k+1)
  ## Remove the self neighbour or the last.
  keep <- nn$nn.idx!=row(nn$nn.idx)
  keep[rowSums(keep)==k+1,k+1] <- FALSE
  nn$nn.idx <- t(matrix(t(nn$nn.idx)[t(keep)],k,length(y)))
  nn$nn.dists <- t(matrix(t(nn$nn.dists)[t(keep)],k,length(y)))
  ## Construct inverse distance weights
  w <- 1/pmax(nn$nn.dists,1.0E-12)
  w <- w/rowSums(w)
  ## Compute difference in actual and predicted
  pred <- rowSums(w*matrix(y[nn$nn.idx],nrow(nn$nn.idx),ncol(nn$nn.idx)))
  list(obs=y,pred=pred)
}



##' Compute averages in windows of length \code{2*lhalf} along a transect,
##' where each window overlaps the previous window by lhalf.
##'
##' Computes averages in sliding windows along transect by taking differences of
##' cumulative sums.
##' @title Compute averages in sliding window
##' @param x One-dimensional spatial presence absence data.
##' @param lhalf The window size increment.
##' @return A numeric vector containing mean density in each window of size \code{2*lhalf}.
window50 <- function(x,lhalf) {
  S <- c(0,cumsum(x))
  k <- seq(2*lhalf+1,length(x)+1,by=lhalf)
  (S[k] - S[k-2*lhalf])/(2*lhalf)
}



##' Form all embeddings of dimension \code{d} for window size \code{2*lhalf} where each window
##' along the transect overlaps the previous by \code{lhalf}.
##'
##' Takes one-dimensional spatial (transect) presence-absence or concentration data and creates a
##' delay embedding of density along the transect. The transect \code{x} must be formatted
##' so that each element in the vector x represents the occupancy of one species,
##' morphotype or habitat in one position along the transect. Returns a matrix where
##' each row is the trajectory of species density along \code{d+1} windows.
##' @title Create a Delay Embedding
##' @param x The spatial transect data.
##' @param lhalf The window size increment.
##' @param dim The embedding dimension + 1.
##' @return An embedding as a matrix.
##' @seealso \code{\link{window50}}
##' @export
##' @examples
##' x <- seq(1, 100, by=1)
##' embed50(x, 5, 4)
embed50 <- function(x,lhalf,dim) {
  w <- window50(x,lhalf)
  N <- length(w)
  k <- outer(seq_len(N-dim+1),seq_len(dim)-1,'+')
  matrix(w[k],nrow(k),ncol(k))
}


##' Compute the prediction error spectra for a range of window sizes using the sliding window approach,
##' for Characteristic Length Scale (CLS) estimation.
##'
##' The transect \code{x} must be numeric with each element corresponding to the presence \code{1}
##' (or concentration > 0) and absence \code{0} of the target species at one point along the transect.
##' The spatial resolution must be constant along the length of the transect and optimally at as fine
##' a resolution as practical.
##' The density of the target species in windows of size \code{2} to \code{2*Lhalf} is computed with
##' the \code{\link{window50}} function. Using the \code{\link{embed50}} function an embedding is created for each
##' window size into \code{d} dimensions. Using the \code{\link{knn}} function the \code{k} nearest neighbours
##' of each \code{d} dimensional point in the embedding are weighted according to their inverse distance,
##' and the trajectory to the \code{d+1} dimension predicted. Difference between the predicted and observed
##' trajectory is calculated either as Error X \code{"errX"}  or Prediction r-squared \code{"PRsq"}.
##'
##' \deqn{Error X = \sqrt{l\times E_t\times \bigg[\Big(X_l^t - \hat{X_l^t}\Big)^2\bigg]}}
##'
##' Where \eqn{l} is the window size,  \eqn{X_l^t} is the observed and \eqn{X_l^t} is the predicted
##' density of species \eqn{X}, and \eqn{E_t} is the expectation of their difference (mean difference).
##'
##' \deqn{Prediction  r^2 = 1 - \frac{E_t\times \bigg[\Big(X_l^t - \hat{X_l^t}\Big)^2\bigg]}{Var(X_l^t)}}
##'
##' Note that Error X is the preferred metric for this sliding window approach.
##' Choice of Lhalf will be constrained by the length of the transect. For example, if \code{d=3}, Lhalf should be
##' 0.49\% of the transect length to have 200 replicates or 0.96\% to obtain 100 replicates of the largest
##' window size.
##'
##' The CLS is defined as the window length at which the spectra reaches a plateau (when plotted).
##' @title Sliding window prediction error computation
##' @param x One-dimensional transect data.
##' @param metric Metric to calculate. Options are \code{"errX"} (Error X) or \code{"PRSq"} (prediction r-squared). Error X is the preferred option for this sliding window approach.
##' @param Lhalf The maximum half window size.
##' @param d The embedding dimension.
##' @param n.samples The number of resamples to draw from the embedding. Default is \code{10}.
##' @param k The number of neighbours to use in prediction. Default is \code{7}.
##' @param replace Sample from the embedding with replacement. Default is \code{TRUE}.
##' @return A matrix of Prediction r-squared or Errror X estimates in which each row is a different window
##' size and each column is a different subsample of the delay embedding.
##' @references Ward D, Wotherspoon S, Melbourne-Thomas J, Haapkyla J, Johnson CR (submitted).
##' Detecting ecological regime shifts from transect data.
##' @seealso \code{\link{window50}},\code{\link{embed50}}, \code{\link{knn}}
##' @export
##' @examples
##' #Calculate Error X from Blue Bowl reef in 2007
##' set.seed(2)
##' x <- CLS1D(ifelse(bb2007 == "CF", 1, 0), Lhalf=50, d=3, metric="errX", n.samples=20)
##' plotCLS(x, metric="errX")
##' #The CLS is 60-70 cm
##'
##' #Compare that to the Error X spectra from the same reef in 2010
##' set.seed(2)
##' x <- CLS1D(ifelse(bb2010 == "CF", 1, 0), Lhalf=50, d=3, metric="errX", n.samples=20)
##' plotCLS(x, metric="errX")
##' #The CLS has declined from 60-70cm to around 30-40 cm.
##'
##' #Now check the Prediction r-squared spectra
##' x <- CLS1D(ifelse(bb2007 == "CF", 1, 0), Lhalf=50, d=3, metric="PRSq", n.samples=20)
##' plotCLS(x, metric="PRSq")
CLS1D<- function(x, Lhalf, d, metric=c("PRSq","errX"), n.samples=10, k=7, replace=TRUE) {

  metric <- match.arg(metric)

  t(sapply(1:Lhalf,function(lhalf){
    ## Create an embedding for this window size
    emb <- embed50(x,lhalf,dim=d+1)
    ## Common sample size of the resamples.
    M <- nrow(embed50(x,Lhalf,d+1))
    ## Resample embedding
    sapply(seq_len(n.samples),
           function(.) {
             ## Subsample with replacement
             smp <- emb[sample(1:nrow(emb),M,replace=replace),,drop=F]
             #gather nearest neighbours and calculate predicted and observed distances
             nn <- knn(smp,k)
             switch(metric,
                    PRSq = {
                      1 - mean((nn$obs-nn$pred)^2)/mean((nn$obs-mean(nn$obs))^2)
                    },
                    errX = {
                      sqrt(2*lhalf)*sqrt(abs(mean((nn$obs-nn$pred)^2)))
                    })
           })
  })
  )
}

##' Compute the prediction error spectra for a range of window sizes using the short time-series
##' approach, for Characteristic Length Scale (CLS) estimation.
##'
##' Takes a matrix or dataframe \code{x} where each column contains transect data from a different time step,
##' and each row contains numeric data on the presence \code{1} (or concentration > 0) and absence \code{0}
##' of the target species at one point along the transect. The spatial resolution must be constant along the
##' length of the transect and optimally at as fine a resolution as practical.
##' The density of the target species in each window is computed with the \code{\link{window50sts}} function.
##' Then, using the \code{\link{knn}} function the \code{k} nearest neighbours of each point in the embedding are
##' weighted according to their inverse distance, and the trajectory to the final dimension (column) predicted.
##' Difference between the predicted and observed trajectory is calculated either as Error X \code{"errX"}  or
##' Prediction r-squared \code{"PRsq"}.
##'
##' \deqn{Error X = \sqrt{l\times E_t\times \bigg[\Big(X_l^t - \hat{X_l^t}\Big)^2\bigg]}}
##'
##' Where \eqn{l} is the window size,  \eqn{X_l^t} is the observed and \eqn{X_l^t} is the predicted
##' density of species \eqn{X}, and \eqn{E_t} is the expectation of their difference (mean difference).
##'
##' \deqn{Prediction  r^2 = 1 - \frac{E_t\times \bigg[\Big(X_l^t - \hat{X_l^t}\Big)^2\bigg]}{Var(X_l^t)}}
##'
##' @title Short time-series prediction error computation
##' @param x One-dimensional transect data.
##' @param metric Metric to calculate. Options are \code{"errX"} (Error X) or \code{"PRSq"} (prediction r-squared). Error X is the preferred option for this sliding window approach.
##' @param Lhalf The maximum half window size.
##' @param n.samples The number of resamples to draw from the embedding. Default is \code{10}.
##' @param k The number of neighbours to use in prediction. Default is \code{7}.
##' @param replace Sample from the embedding with replacement. Default is \code{TRUE}.
##' @return A matrix of Prediction r-squared or Errror X estimates in which each row is a different window
##' size and each column is a different subsample of the delay embedding.
##' @references Ward D, Wotherspoon S, Melbourne-Thomas J, Haapkyla J, Johnson CR (submitted).
##' Detecting ecological regime shifts from transect data.
##' @seealso \code{\link{window50sts}}, \code{\link{knn}}
##' @importFrom RANN nn2
##' @export
##' @examples
##' x <- CLS1Dsts(Sp2R1_sts, Lhalf=20, metric="errX", n.samples=20)
##' plotCLS(x, metric="errX")
##'
##' x <- CLS1Dsts(Sp2R2_sts, Lhalf=20, metric="errX", n.samples=20)
##' plotCLS(x, metric="errX")
##'
##' x <- CLS1Dsts(Sp2R1_sts, Lhalf=20, metric="PRSq", n.samples=20)
##' plotCLS(x, metric="PRSq")
##'
##' x <- CLS1Dsts(Sp2R2_sts, Lhalf=20, metric="PRSq", n.samples=20)
##' plotCLS(x, metric="PRSq")
CLS1Dsts<- function(x, Lhalf, metric=c("PRSq","errX"), n.samples=10, k=7, replace=TRUE){

  metric <- match.arg(metric)

  t(sapply(1:Lhalf,
           function(lhalf) {
             emb <- window50sts(x,lhalf)
             M <- nrow(window50sts(x,Lhalf))
             sapply(seq_len(n.samples), function(.) {
               ## Draw a random sample
               smp <- emb[sample(1:nrow(emb),M,replace=replace),,drop=F]
               ## gather nearest neighbours and calculate predicted and observed distances
               nn <- knn(smp,k)
               switch(metric,
                      PRSq = {
                        1 - mean((nn$obs-nn$pred)^2)/mean((nn$obs-mean(nn$obs))^2)
                      },
                      errX = {
                        sqrt(2*lhalf)*sqrt(abs(mean((nn$obs-nn$pred)^2)))
                      })
             })}))

}


##' Compute averages in windows of length \code{2*lhalf} along each transect in a set of time-delayed
##' transects. Each window overlaps the previous window by lhalf within the same transect. Between transects
##' the windows are in the same position spatially.
##'
##' Computes averages in windows along transects by taking differences of
##' cumulative sums.
##' @title Compute averages for short time-series
##' @param x Data frame or matrix of 1-dimensional spatial transects
##' @param lhalf the window size increment.
##' @return a matrix or dataframe containing mean density in each window of size \code{2*lhalf}
window50sts <- function(x,lhalf) {
  S <- rbind(0,x)
  for(k in 1:ncol(x)) S[,k] <- cumsum(S[,k])
  k <- seq(2*lhalf+1,nrow(x)+1,by=lhalf)
  (S[k,] - S[k-2*lhalf,])/(2*lhalf)
}


##' Plot the prediction error spectra for a range of window sizes, for estimation of the Characteristic
##' Length Scale (CLS).
##'
##' Takes matrix of Error X or Prediction r-squared estimates from \code{\link{CLS1D}} or \code{\link{CLS1Dsts}}.
##' Uses \code{\link[zoo]{rollapply}} from package \pkg{zoo} to estimate the rolling mean (rolling window width = 3)
##' of the estimates for each window size, and plots with pointwise 95 percent confidence intervals.
##' @title Plot prediction error spectra
##' @param X1 Matrix of prediction error estimates returned by \code{\link{CLS1D}} or \code{\link{CLS1Dsts}}
##' @param Lhalf The maximum half window size for which prediction error was estimated. Should be equal to \code{nrow(X1)}.
##' @param metric Optional. Specify metric for y-axis label. Options are \code{"errX"} (Error X) or \code{"PRSq"} (prediction r-squared).
##' @param ylim Optional. Specify y-axis range. Default is \code{range(X1)}.
##' @param xlab Optional. Specify x-axis label. Default is "Window length".
##' @param xticks Optional. Specify spacing of ticks on x-axis. Default is by 10.
##' @return Returns a plot of the prediction error spectra.
##' @importFrom graphics axis matlines mtext plot.new plot.window polygon
##' @importFrom grDevices adjustcolor
##' @importFrom stats var
##' @importFrom zoo rollapply
##' @export
##' @examples
##' x <- CLS1D(ifelse(bb2007 == "CF", 1, 0), Lhalf=50, d=3, metric="errX", n.samples=20)
##' plotCLS(x, metric="errX")
plotCLS <- function(X1, Lhalf=nrow(X1), metric, ylim, xlab, xticks=10){
  plot.new()

  if(missing(ylim)){plot.window(xlim=c(0, 2*Lhalf), ylim=c(range(X1)))
  } else {plot.window(xlim=c(0, 2*Lhalf), ylim=ylim)}

  axis(2, col="darkblue", lwd=2, cex.axis=0.8, las=1, tcl=-0.23, mgp=c(3, 0.28,0), col.axis="darkblue")

  se <- function(x) sqrt(var(x)/length(x))
  ciX1<- vector()
  for(i in 1:nrow(X1)){
    ciX1[i] <- as.vector(se(X1[i,]))
  }
  polygon(c(2*(1:Lhalf),rev(2*(1:Lhalf))),c(rollapply(rowMeans(X1)- 1.96*ciX1, mean, width=3, partial=TRUE),
                                    rev(rollapply(rowMeans(X1)+ 1.96*ciX1, mean, width=3, partial=TRUE))), col=adjustcolor("royalblue1", alpha.f = 0.35), border=NA)

  matlines(2*(1:Lhalf),rollapply(rowMeans(X1),mean, width=3, partial=TRUE),type = "l", col="darkblue", lwd = 2)
  axis(1, at=c(seq(0,2*Lhalf,by=xticks)), lwd=2, tcl=-0.23,las=1, cex.axis=1, mgp=c(3, 0.28,0))

  if(missing(xlab)){mtext("Window length",  side=1, line=2, cex=1.2)
  } else {mtext(xlab, side=1, line=2, cex=1.2)}
  if(missing(metric)){}
  else if(metric=="errX"){mtext("Error X", side=2, line=2, cex=1.1)}
  else if(metric=="PRSq"){mtext("Prediction r-sq", side=2, line=2, cex=1.1)}

}
