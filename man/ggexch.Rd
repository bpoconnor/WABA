\name{ggexch}
\alias{ggexch}
\title{Correlational Analyses of Dyad-Level Data in the Exchangeable-Case}
\description{Conducts correlational analyses of dyad-level data in the 
	exchangeable-case, based on Griffin & Gonzalez (1995).}\arguments{
  \item{data}{
  \code{}An all-numeric dataframe where the rows are cases & the columns are the variables.
         Cases with missing values are not permitted in the data file.
         The first collumn of data must contain the dyad numbers (the dyad id# given to both
         dyad members). The second and subsequent collumns contain the variables to be analyzed.}
}
\details{The simplest nested data structures are based on dyads. The ggdist function 
	described elsewhere in this package can be used when the dyadic partners are 
	distinguishable (and not exchangeable or interchangeable), because they are drawn 
	from different classes or categories, such as men and women. The ggexch function 
	can be used for dyad-level data wherein the partners are exchangeable or 
	interchangeable, as in pairs of men, pairs of women, or pairs of same-sex children. 
	This function performs the statistical analyses described by Gonzalez and Griffin (1995).

    Produces the overall within-partner correlation, the overall cross-partner correlation, 
    intraclass correlations, and statistics for dyad-level and individual-level effects.
}
\references{Griffin, D., & Gonzalez, R. (1995). Correlational analysis of dyadlevel data 
	in the exchangeable case. Psychological Bulletin, 118, 430-439.
}
\author{Brian P. O'Connor }
\examples{ggexch(data_ggexch)}
