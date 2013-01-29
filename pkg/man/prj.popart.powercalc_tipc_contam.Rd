\name{prj.popart.powercalc_tipc_contam}
\alias{prj.popart.powercalc_tipc_contam}
\title{compute the power in differentiating H0: acute to acute is <10\% vs H0: acute to acute is >40\% with the tipc cluster approach}
\usage{
  prj.popart.powercalc_tipc_contam(p.prev.instudy.clu.armC = 0.4,
    p.phylosignal = 0.7, opt.pooled = "no pooling",
    opt.sampling = "PC and HCC")
}
\arguments{
  \item{p.prev.instudy.clu.armC}{Probability that HIV+
  individuals visit an HCC in arm C (sensitivity
  analysis).}

  \item{p.phylosignal}{Probability that a true (potentially
  indirect) transmission is identified with a phylogenetic
  method. (Here not used)}

  \item{opt.pooled}{Pooling option for power analysis}

  \item{opt.sampling}{Sampling option for trial}
}
\description{
  This script varies the probability that transmission
  occurs outside the community.
}

