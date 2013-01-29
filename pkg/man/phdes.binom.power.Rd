\name{phdes.binom.power}
\alias{phdes.binom.power}
\title{performs power calculations and computes binomial confidence intervals for a given number of acute 2 acute transmissions}
\usage{
  phdes.binom.power(x2i, i2i, idx.A, idx.C, test.prop0.A,
    test.prop1.A, test.prop0.C = NULL, test.prop1.C = NULL,
    test.alpha = 0.05, verbose = 0, method.pw = "homebrew",
    method.ci = "bayes")
}
\arguments{
  \item{x2i}{number of transmissions}

  \item{i2i}{number of acute 2 acute transmissions}

  \item{idx.A}{row index of arm A sites}

  \item{idx.C}{row index of arm C sites}

  \item{test.prop0.A}{null proportion of transmission
  during acute phase for arm A}

  \item{test.prop1.A}{alternate proportion of transmission
  during acute phase for arm A}

  \item{test.prop0.C}{null proportion of transmission
  during acute phase for arm C (if ommitted set to the one
  for arm A)}

  \item{test.prop1.C}{alternate proportion of transmission
  during acute phase for arm C (if ommitted set to the one
  for arm A)}

  \item{test.alpha}{level of test}

  \item{method.pw}{method used to calculate Binomial power}

  \item{method.ci}{method used to calculate Binomial
  confidence intervals}
}
\value{
  list containing power and confidence intervals
}
\description{
  performs power calculations and computes binomial
  confidence intervals for a given number of acute 2 acute
  transmissions
}

