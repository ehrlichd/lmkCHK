% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/package-support.R
\name{LMK_swap}
\alias{LMK_swap}
\title{Swap Landmarks}
\usage{
LMK_swap(a, l1, l2)
}
\arguments{
\item{a}{a p x k x n matrix of landmark coordinates}

\item{l1}{a numeric vector of landmark index (sequence number) to be changed}

\item{l2}{a numeric vector of landmark index (Sequence number) to change to}
}
\description{
Function to swap landmark sequences within an individual or across an entire array
}
\examples{

lmks <- cbind(1:10, 1:10, 1:10)
plot(lmks, col = rainbow(10), pch = 16) ## plot sequence (in 2 dimensions)
lmks2 <- LMK_swap(lmks, l1 = c(2,8), l2 = c(8,2)) ## flip landmarks 2,8
plot(lmks2, col = rainbow(10), pch = 16) ## show flipped landmarks


}
\author{
Daniel Ehrlich
}
