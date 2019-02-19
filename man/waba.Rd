\name{waba}
\alias{waba}
\title{Within-And-Between-Analysis}
\description{Conducts Within-And-Between-Analyses and produces a wide range of WABA statistics}
\usage{waba(data)}
\arguments{
  \item{data}{
  \code{}An all-numeric dataframe where the rows are cases & the columns are the variables.

         Cases with missing values are not permitted in the data file. 

         The first value in each row (i.e., the first column of values
         in the data file) must be the individuals' group number/code,
         which must be an integer. The function sorts individuals into groups
         on the basis of these numbers/codes.

         Variable scores appear in subsequent columns. 

         Multiple Variable Analyses (MVA) are conducted when the number of
         variables = 3, which is a limit that is determined by
         Hotelling's t-test for dependent correlations. For MVAs involving
         more than three variables, simply run the function again using
         new/different combinations of three variables. }
}
\details{Within and between analysis (WABA) was developed by Dansereau, Alutto, and 
	Yammarino (1984). The procedure can process data from groups (not just dyads) 
	of varying sizes, and there are similarities with the Griffin and Gonzalez 
	approaches for dyad-level data described elsewhere in this package. However, 
	the computations in WABA are different, and there is some unique terminology 
	and statistical output.

    WABA involves three steps, as described by Yammarino and Markham (1992). First, 
    scores on each variable are examined to determine whether the variance occurs 
    primarily between groups, within groups, or both between and within groups. 
    Next, relationships between variables are examined to determine whether the 
    associations are primarily a function of between-groups covariance, within-groups 
    covariance, or both within- and between-groups covariances. Finally, the results 
    from these two steps are examined for consistency, and appropriate conclusions are drawn. 

    Produces multiple WABA statistics.
}
\references{
       {Dansereau, F., Alutto, J., & Yammarino, F. (1984). Theory testing in organizational 
       	behavior. Englewood Cliffs, NJ: Prentice-Hall.}
\cr\cr {Dansereau, F., Chandrasekaran, G., Dumas, M., Coleman, D., Ehrlich, S., & Bagchi, D. 
	(1986). DETECT: Data enquiry that tests entity and correlational/causal theories.Williamsville, 
	NY: Institute For Theory Testing.}
\cr\cr {Dansereau, F., & Yammarino, F. J. (2000). Within and between analysis: The varient 
	paradigm as an underlying approach to theory building. In K. J. Klein & S.W. J. Kozlowski 
	(Eds.), Multilevel theory, research, and methods in organizations: Foundations, extensions, 
	and new directions (pp. 425-466). San Francisco: Jossey-Bass.}
\cr\cr {Yammarino, F. J. (1998). Multivariate aspects of the varient /WABA approach. 
	Leadership Quarterly, 9, 203-227.}
\cr\cr {Yammarino, F. J., & Markham, S. (1992). On the application of within and between analysis: 
	Are absence and affect really group based. Journal of Applied Psychology, 77, 168-176.}
}
\author{Brian P. O'Connor}
\examples{
waba(data_Detect_Set_A)

waba(data_Bliese)

waba (data_jsp[c('school','english','maths','ravens' )])
}
