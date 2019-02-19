\name{wabamra}
\alias{wabamra}
\title{WABA - Multiple Relationship Analysis}
\description{Conducts WABA - Multiple Relationship Analyses, as described by 
	Dansereau et al. (1984) and Yammarino (1984).}
\usage{wabamra(data)}\arguments{
  \item{data}{
  \code{}An all-numeric dataframe where the rows are cases & the columns are the variables.

         Cases with missing values are not permitted in the data file. 

         This function conducts separate WABA analyses for each value
         of a specified "Condition" variable.  The function then conducts
         Multiple Relationship Analysis comparisons of the WABA
         coefficients for the different Condition values.

         The first value in each row (i.e., the first column of values
         in the data file) must be the Condition number.  Condition
         numbers must be integers. The lowest Condition number cannot
         be less than one. It is also best for there to be no missing
         values between the lowest and highest Condition numbers
         e.g., if the lowest value is 1 and the highest value is 5,
         then there should also be Condition values of 2, 3, and 4
         somewhere in the data file. Gaps in the integers may cause problems. 

         The second value in each row (i.e., the second column of values
         in the data file) must be the individual's group number/code.
         The function sorts individuals into groups on the basis of these
         numbers/codes. 

         Variable scores appear in subsequent columns. }
}
\details{Multiple relationship analyses (MRAs) can be conducted on the WABA results 
	for two or more conditions. The focus in this case is on the identification of 
	moderators or contingencies in WABA patterns. A condition is an additional group 
	variable across which WABA patterns may or may not be consistent. For example, 
	WABA patterns for a given set of variables may vary depending on sex, age group, 
	school type, or neighborhood type. MRA involves pairwise testing for possible 
	differences between WABA correlations for all possible combinations of values 
	of the condition variable of interest. Dansereau et al. (1984) and Yammarino (1998) 
	provided guidelines for interpreting the results from the many MRA contrasts.

    Produces multiple WABA MRA statistics.
}
\references{
       {Dansereau, F., Alutto, J., & Yammarino, F. (1984). Theory testing in organizational 
       	behavior. Englewood Cliffs, NJ: Prentice-Hall.}
\cr\cr {Yammarino, F. J. (1998). Multivariate aspects of the varient /WABA approach. 
	Leadership Quarterly, 9, 203-227.}
}
\author{Brian P. O'Connor }
\examples{
wabamra(data_Detect_Set_A_mra)

wabamra (data_jsp[c('sex','school','english','maths','ravens' )] )
}
