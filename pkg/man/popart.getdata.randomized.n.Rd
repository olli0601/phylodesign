\name{popart.getdata.randomized.n}
\alias{popart.getdata.randomized.n}
\title{randomize the number of individuals with HIV, on ART at each trial site and in the population cohort}
\usage{
  popart.getdata.randomized.n(x, pc.size, d.study,
    rtn.exp = 0)
}
\arguments{
  \item{x}{data frame of trial sites with allocated arms}

  \item{pc.size}{size of the population cohort}

  \item{d.study}{duration of the trial in years}

  \item{rtn.exp}{indicator if expectations are returned
  instead of draws from binomial (for debugging)}
}
\value{
  data frame of trial sites with allocated arms
}
\description{
  randomize the number of individuals with HIV, on ART at
  each trial site and in the population cohort
}
\examples{
sites<- popart.getdata.randomized.arm(1)
				popart.getdata.randomized.n(sites, 2500, 3)
}

