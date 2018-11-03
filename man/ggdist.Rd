\name{ggdist}
\alias{ggdist}
\title{Correlational Analysis of Dyad-Level Data in the Distinguishable-Case}
\description{Conducts correlational analyses of dyad-level data in the 
	distinguishable-case, based on Gonzalez & Griffin (1999).}
\usage{ggdist(data)}\arguments{
  \item{data}{
  \code{}An all-numeric dataframe where the rows are cases & the columns are the variables.
         Cases with missing values are not permitted in the data file.
         The first collumn of data must contain the dyad numbers (the dyad id# given to both
         dyad members). 
         The second collumn contains the person id #s (e.g., sex) coded as 1 or 2.
         The third and subsequent collumns contain the variables to be analyzed.}
}
\details{The simplest nested data structures are based on dyads. Imagine that neuroticism 
	and marital satisfaction scores are obtained from both husbands and wives from a 
	large number of couples. The dyadic partners in this case are distinguishable (and 
	not exchangeable or interchangeable), because they are drawn from different classes 
	or categories (men and women) that might well differ in their variable means, variances, 
	and covariances. A variety of potentially important questions can be examined using 
	such scores on two variables from each of the partners. The apparently simple question, 
	"Is neuroticism associated with marital satisfaction?" quickly becomes multifaceted in 
	dyadic datasets. This function performs the statistical analyses described by 
	Gonzalez and Griffin (1999).

    Produces the overall within-partner correlation, the overall cross-partner correlation, 
    intraclass correlations, and statistics for dyad-level and individual-level effects.}
\references{
 {Gonzalez, R., & Griffin, D. (1999). The correlational analysis of dyad-level data in the 
 	distinguishable case. Personal Relationships, 6, 449-469.}
}
\author{Brian P. O'Connor }
\examples{ggdist(data_ggdist)}
