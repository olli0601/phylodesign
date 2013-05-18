\name{acute.lkl.tree.xk.ik}
\alias{acute.lkl.tree.xk.ik}
\title{Compute the likelihood of a transmission chain with \code{nx} transmissions from the index case and \code{ni} transmissions from non-index cases under the \code{Acute} model
There are 1+nx+ni NODES in this tree.}
\usage{
  acute.lkl.tree.xk.ik(nx, ni, rx, ri, dT, log = 0)
}
\arguments{
  \item{nx}{number of transmissions from donor with risk
  group X, nx=0 is possible}

  \item{ni}{number of transmissions from donor with risk
  group I, ni=0 is possible}
}
\description{
  Compute the likelihood of a transmission chain with
  \code{nx} transmissions from the index case and \code{ni}
  transmissions from non-index cases under the \code{Acute}
  model There are 1+nx+ni NODES in this tree.
}

