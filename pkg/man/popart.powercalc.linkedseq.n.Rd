\name{popart.powercalc.linkedseq.n}
\alias{popart.powercalc.linkedseq.n}
\title{performs power calculations and computes binomial confidence intervals for the "linkedseq" approach}
\usage{
  popart.powercalc.linkedseq.n(inc, idx.A, idx.C,
    test.prop0, test.prop1, test.alpha, verbose = 0)
}
\arguments{
  \item{inc}{computed linked sequences of incident cases}

  \item{idx.A}{row index of arm A sites}

  \item{idx.C}{row index of arm C sites}

  \item{test.prop0}{null proportion of transmission during
  acute phase}

  \item{test.prop1}{alternate proportion of transmission
  during acute phase}

  \item{test.alpha}{level of test}
}
\value{
  list containing power and confidence intervals
}
\description{
  performs power calculations and computes binomial
  confidence intervals for the "linkedseq" approach
}

