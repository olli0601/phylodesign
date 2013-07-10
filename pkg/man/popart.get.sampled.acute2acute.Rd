\name{popart.get.sampled.acute2acute}
\alias{popart.get.sampled.acute2acute}
\title{compute the number of acute to acute transmissions from tip clusters up to order 3}
\usage{
  popart.get.sampled.acute2acute(x, tipc.p,
    method = "PC and HCC", rtn.int = 1,
    p.vhcc.prev.AB = 0.95, p.vhcc.prev.C = 0.5,
    p.vhcc.inc.AB = 0.8, p.vhcc.inc.C = 0.25,
    consent.PC = 0.9, consent.HCC = 0.5, p.lab = 0.8)
}
\arguments{
  \item{x}{data frame of trial sites with randomized arms
  and population numbers}

  \item{tipc.p}{expected frequencies of tip cluster counts
  under hypothesis}

  \item{method}{sampling scheme. Options are "PC after yr 1
  and HCC", "PC and HCC", "only HCC"}

  \item{rtn.int}{indicator if integer values are returned}

  \item{p.vhcc.prev.AB}{proportion of initially +ve
  individuals who visit HCC at some point in 3 year period
  in Arms A & B}

  \item{p.vhcc.prev.C}{proportion of initially +ve
  individuals who visit HCC at some point in 3 year period
  in Arm C}

  \item{p.vhcc.inc.AB}{proportion of incident cases who
  visit HCC at some point in 3 year period in Arms A and B}

  \item{p.vhcc.inc.C}{proportion of incident cases who
  visit HCC at some point in 3 year period in Arm C}

  \item{consent.PC}{consent rate to phylogenetics in
  population cohort}

  \item{consent.HCC}{consent rate to phylogenetics in
  health care centres}

  \item{p.lab}{proportion of isolates that make it through
  to full genomes}
}
\value{
  For each trial site, a vector containing
  \item{i2i.s}{number of sampled acute to acute
  transmissions} \item{i2i.c}{number of acute to acute
  transmissions} \item{x2i.s}{number of sampled
  transmissions} \item{x2i.c}{number of true transmissions}
}
\description{
  compute the number of acute to acute transmissions from
  tip clusters up to order 3
}
\examples{
sites<- popart.getdata.randomized.arm(1)
sites<- popart.getdata.randomized.n(sites, 2500, 3)
popart.get.sampled.transmissions(sites)
}

