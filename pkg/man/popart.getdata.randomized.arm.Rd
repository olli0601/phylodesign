\name{popart.getdata.randomized.arm}
\alias{popart.getdata.randomized.arm}
\title{randomize allocation of arms to triplets}
\usage{
  popart.getdata.randomized.arm(randomize.n, rtn.fixed = 0)
}
\arguments{
  \item{randomize.n}{number of times the randomization is
  repeated}

  \item{rtn.fixed}{indicator if a fixed allocation of arms
  is returned (for debugging)}
}
\value{
  data frame of trial sites with allocated arms
}
\description{
  randomize allocation of arms to triplets
}
\examples{
popart.getdata.randomized.arm(1)
}

