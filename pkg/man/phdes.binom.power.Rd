\name{phdes.binom.power}
\alias{phdes.binom.power}
\title{performs power calculations and computes binomial confidence intervals for a given number of acute 2 acute transmissions}
\usage{
  phdes.binom.power(x2i, i2i, test.prop0, test.prop1,
    test.alpha = 0.05, verbose = 0, method.pw = "homebrew",
    method.ci = "bayes")
}
\arguments{
  \item{x2i}{number of transmissions, typically for a
  specific arm}

  \item{i2i}{number of acute 2 acute transmissions,
  typically for a specific arm}

  \item{test.prop0}{null proportion of transmission during
  acute phase}

  \item{test.prop1}{alternate proportion of transmission
  during acute phase}

  \item{test.alpha}{level of test}

  \item{method.pw}{method used to calculate Binomial power.
  Possible options are "homebrew","cloglog", "logit",
  "probit", "asymp", "lrt", "exact".}

  \item{method.ci}{method used to calculate Binomial
  confidence intervals. Possible options are "exact", "ac",
  "asymptotic", "wilson", "prop.test", "bayes", "logit",
  "cloglog", "probit".}
}
\value{
  \item{conf}{binomial confidence interval for each site in
  x2i and i2i} \item{is.conf}{boolean flag for each site,
  indicating if test.prop0 is below the confidence interval
  AND if the i2i exceed 5} \item{power}{binomial power for
  each site in x2i and i2i}
}
\description{
  performs power calculations and computes binomial
  confidence intervals for a given number of acute 2 acute
  transmissions
}
\examples{
sites<- popart.getdata.randomized.arm(1)
	sites<- popart.getdata.randomized.n(sites, 2500, 3)
	x2i<- popart.get.sampled.transmissions(sites)
	x2i<- phdes.get.linked.transmissions(x2i,0.8)
	tmp<- popart.pool(sites, x2i)
	x2i<- tmp[["transm"]]
	idx.A<- tmp[["idx.A"]]
	phdes.binom.power(	x2i[idx.A,,drop=0], round(x2i[idx.A,,drop=0]*test.prop0), test.prop0, test.prop1, test.alpha, verbose=0)
}

