\name{phdes.get.linked.transmissions}
\alias{phdes.get.linked.transmissions}
\title{compute the number of incident sequences that can be linked to a prevalent case with an imperfect phylogenetic method}
\usage{
  phdes.get.linked.transmissions(x2i, p.linkage,
    rtn.exp = 0)
}
\arguments{
  \item{x2i}{number of incident sequences that could be
  linked to a prevalent case}

  \item{p.linkage}{proportion of incident sequences that
  are linked with a phylogenetic method}

  \item{rtn.exp}{indicator if expectations are returned
  instead of draws from binomial (for debugging)}
}
\value{
  linked.x2i number of incident sequences that are linked
  to a prevalent case under a Binomial link model
}
\description{
  compute the number of incident sequences that can be
  linked to a prevalent case with an imperfect phylogenetic
  method
}
\examples{
sites<- popart.getdata.randomized.arm(1)
	sites<- popart.getdata.randomized.n(sites, 2500, 3)
	x2i<- popart.get.sampled.transmissions(sites)
	phdes.get.linked.transmissions(x2i,0.8)
}

