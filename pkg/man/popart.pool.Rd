\name{popart.pool}
\alias{popart.pool}
\title{pool the number of transmissions according to a pooling option}
\usage{
  popart.pool(opt.pooled, sites, inc)

  popart.pool(opt.pooled, sites, inc)
}
\arguments{
  \item{sites}{data frame of trial sites}

  \item{inc}{matrix of transmissions}

  \item{method}{poopling method. Supported methods are "no
  pooling", "pooled across country", "pooled across ZA",
  "pooled across SA", "pooled across trial".}

  \item{opt.pooled}{poopling option}

  \item{sites}{data frame of trial sites}

  \item{inc}{matrix of acute to acute transmissions}
}
\value{
  \item{inc}{pooled matrix of acute to acute transmissions}
  \item{idx.A}{row index of arm A sites} \item{idx.B}{row
  index of arm B sites} \item{idx.C}{row index of arm C
  sites}

  list with the components 'inc': pooled matrix of acute to
  acute transmissions, 'idx.A': row index of arm A sites,
  'idx.C': row index of arm C sites
}
\description{
  pool the number of transmissions according to a pooling
  option

  pool the number of acute to acute transmissions according
  to a pooling option
}
\examples{
sites<- popart.getdata.randomized.arm(1)
	sites<- popart.getdata.randomized.n(sites, 2500, 3)
	sampled.x2i<- popart.get.sampled.transmissions(sites)
	popart.pool(sites, sampled.x2i, "pooled across ZA")
}

