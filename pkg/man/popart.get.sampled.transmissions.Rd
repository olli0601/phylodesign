\name{popart.get.sampled.transmissions}
\alias{popart.get.sampled.transmissions}
\title{compute the number of incident sequences that can be linked to a prevalent case with a perfect phylogenetic method}
\usage{
  popart.get.sampled.transmissions(x,
    method = "PC and HCC", rtn.int = 1,
    p.vhcc.prev.AB = 0.95, p.vhcc.prev.C = 0.5,
    p.vhcc.inc.AB = 0.8, p.vhcc.inc.C = 0.25,
    consent.PC = 0.9, consent.HCC = 0.5, p.lab = 0.8,
    p.community = 0.85)
}
\arguments{
  \item{x}{data frame of trial sites with randomized arms
  and population numbers}

  \item{method}{sampling scheme. Currently support "PC and
  HCC", "only HCC"}

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

  \item{p.community}{proportion of infections that come
  from within the community}
}
\value{
  x2i the number of sampled transmissions
}
\description{
  compute the number of incident sequences that can be
  linked to a prevalent case with a perfect phylogenetic
  method
}
\examples{
sites<- popart.getdata.randomized.arm(1)
				sites<- popart.getdata.randomized.n(sites, 2500, 3)
				popart.get.sampled.transmissions(sites)
}

