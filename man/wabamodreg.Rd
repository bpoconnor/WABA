\name{wabamodreg}
\alias{wabamodreg}
\title{WABA - Moderated Regression Analysis}
\description{This WABA function focuses on the interaction between two variables in 
	the prediction of a third. It is a WABA extension of regular moderated multiple regression.}
\usage{wabamodreg(data)}
\arguments{
  \item{data}{
  \code{}An all-numeric dataframe where the rows are cases & the columns are the variables.

         The first value in each row (i.e., the first column of values in
         the data file) must be the individuals' group number/code,
         which must be an integer. The function sorts individuals into groups
         on the basis of these numbers/codes.
         Variable scores appear in subsequent columns. 

         The second value in each row (i.e., the second column of values
         in the data file) are treated by the function as the
         Dependent/Outcome variable. 

         The third value in each row (i.e., the third column of values in
         the data file) are treated by the function as the
         Independent/Predictor variable. 

         The fourth value in each row (i.e., the fourth column of values
         in the data file) are treated by the function as the Moderator
         variable (technically also an IDV/predictor). 

         There is no need to enter product term variables or to enter
         composite variable scores. The function takes care of these
         aspects of the analyses internally. }
}
\details{In WABA analyses it is possibile to examine interactions between variables 
	(categorical or continuous) in the prediction of a designated dependent or outcome 
	variable. The analytic technique, described by Schriesheim (1995), is a direct 
	extension of familiar moderated multiple regression analysis to the WABA analyses 
	described above. Hierarchical regression analyses are conducted, focusing on whether 
	the product term for two variables (which carries the interaction) accounts for 
	significant variance in the outcome variable beyond the variation that is accounted 
	for by the two predictor variables (main effects). A significant increase in variance 
	accounted for by the product term indicates a significant interaction. These hierarchical 
	analyses are conducted using both the within-groups and the between-groups WABA correlations. 
	The wabamodreg function provides all of the statistics reported by Schriesheim for 
	his illustrative example.

    Produces multiple WABA statistics.}

\references{
       {Schriesheim, C. (1995). Multivariate and moderated within-and between-entity analysis 
       	(WABA) using hierarchical linear multiple regression. Leadership Quarterly, 6, 1-18.}
\cr\cr {Schriesheim, C. A., Cogliser, C. C., & Neider, L. L. (1995). Is it "trustworthy"? a 
	multiple levels-of-analysis reexamination of an ohio state leadership study, with implications 
	for future research. Leadership Quarterly, 6, 111-145.}
\cr\cr {Schriesheim C., Neider L.L., Scadura T. (1998), Delegation and leader-member exchange: 
	main effects, moderators, and measurement issues. Academy of Management Journal, 41(3), 298-318.}
}
\author{Brian P. O'Connor }
\examples{
wabamodreg(data_Detect_Set_A)

wabamodreg (data_jsp[c('school','english','maths','ravens' )])
}
