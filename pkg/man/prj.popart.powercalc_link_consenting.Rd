\name{prj.popart.powercalc_link_consenting}
\alias{prj.popart.powercalc_link_consenting}
\title{compute the power in differentiating H0: acute to acute is <10\% vs H0: acute to acute is >40\% with the linkage approach}
\usage{
  prj.popart.powercalc_link_consenting(p.phylosignal = 0.7,
    p.nocontam = 0.85, p.prev.instudy.clu.armC = 0.2,
    opt.pooled = "no pooling", opt.sampling = "PC and HCC")
}
\arguments{
  \item{p.phylosignal}{Probability that a true (potentially
  indirect) transmission is identified with a phylogenetic
  method.}

  \item{p.nocontam}{Frequency with which transmissions
  occur within the community}

  \item{p.prev.instudy.clu.armC}{Probability that HIV+
  individuals visit an HCC in arm C (sensitivity
  analysis).}

  \item{opt.pooled}{Pooling option for power analysis}

  \item{opt.sampling}{Sampling option for trial}
}
\description{
  This script varies the probability that HIV+ individuals
  consent to participation in the phylogenetic study at HCC
}

