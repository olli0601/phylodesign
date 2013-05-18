\name{popart.getsampled.linkedseq.n}
\alias{popart.getsampled.linkedseq.n}
\title{compute the fraction of linked incident sequences}
\usage{
  popart.getsampled.linkedseq.n(x, opt.sampling,
    v_prev_AB = 0.95, v_prev_C = 0.5, v_inc_AB = 0.8,
    v_inc_C = 0.25, consent.PC = 0.9, consent.HCC = 0.5,
    p.lab = 0.8, p.community = 0.85)
}
\arguments{
  \item{x}{data frame}

  \item{opt.sampling}{sampling scheme}

  \item{v_prev_AB}{proportion of initially +ve individuals
  who visit HCC at some point in 3 year period in Arms A &
  B}

  \item{v_prev_C}{proportion of initially +ve individuals
  who visit HCC at some point in 3 year period in Arm C}

  \item{v_inc_AB}{proportion of incident cases who visit
  HCC at some point in 3 year period in Arms A and B}

  \item{v_inc_C}{proportion of incident cases who visit HCC
  at some point in 3 year period in Arm C}

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
  linked.incident
}
\description{
  compute the fraction of linked incident sequences
}

