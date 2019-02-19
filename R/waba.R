

waba <- function (data) {

# WABA: Within And Between Groups Analysis  

# Sources:
# Dansereau, Alutto, & Yammarino (1984). Theory Testing in
# Organizational Behavior. Prentice-Hall
# Yammarino & Markham (1992). J. of Applied Psychology, 77, 168-176
# www.LevelsOfAnalysis.com 

# Prepare a raw data matrix for analysis,
# where rows = cases, & columns = variables
# Each row contains the data from a single individual
# Cases with missing values are not permitted in the data file. 

# The first value in each row (i.e., the first column of values
# in the data file) must be the individuals' group numbers/codes,
# which must be integers. The program sorts individuals into groups
# on the basis of these numbers/codes
# Variable scores appear in subsequent columns. 

# Multiple Variable Analyses are conducted when the number of
# variables = 3, which is a limit that is determined by
# Hotelling's t-test for dependent correlations.  For MVAs involving
# more than three variables, simply run this program again using
# new/different combinations of three variables. 


varnames <- colnames(data)[-1]

if ( all(complete.cases(data)) == 'FALSE' ) {
cat("\n\nERROR: There are missing values in the data matrix. Fix this problem and try again.\n\n") }

data <- as.matrix(data)

data <- data[order(data[,1]),] # sorting by 1st column (the group numbers/codes)

# removing cases for groups that have just one person/case
grpsize <- 1
dumped <- matrix(-9999,1,ncol(data))
for (luper1 in 2:nrow(data) ) {
	if ((data[luper1,1] == data[(luper1-1),1]))   grpsize <- grpsize + 1 
	if ((data[luper1,1] != data[(luper1-1),1] ) & grpsize == 1 ) { 
	    dumped <- rbind( dumped, data[(luper1-1),] )
		data[(luper1-1),2:ncol(data)] <- rep(NA, (ncol(data)-1)) 
	}
	if ((data[luper1,1] != data[(luper1-1),1] ) & grpsize > 1 )  grpsize <- 1
	if (luper1 == nrow(data) & data[luper1,1] != data[(luper1-1),1] ) 
		{ data[(luper1),2:ncol(data)] <- rep(NA, (ncol(data)-1))  }
}
if (nrow(dumped) > 1) {
	dumped <- dumped[2:nrow(dumped),]
	dimnames(dumped) <-list(rep("", dim(dumped)[1]))
	colnames(dumped) <- colnames(data)
	cat("\n\nThe following cases were removed from the data file because there were")
	cat("\nno other cases with the same group code:\n\n" )
	print(dumped)
}

data <- na.omit(data)

nss <- nrow(data)

totmn <- colSums(data[,2:ncol(data), drop = FALSE]) / nss

totdev  <- matrix( -9999, nss, ncol(data)-1 ) 
withdev <- matrix( -9999, nss, ncol(data)-1 )
betwdev <- matrix( -9999, nss, ncol(data)-1 )

grpns   <- matrix( -9999, 1, 1 )

first <- 1
ngrps <- 0

for (luper1 in 2:nss ) {
	if ((data[luper1,1] != data[(luper1-1),1] ) | luper1 == nss ) {
		if ( luper1 != nss ) { last <- luper1 - 1 }
		if ( luper1 == nss ) { last <- luper1 }

	ngrps <- ngrps + 1

	tempdat <- data[first:last,2:ncol(data), drop = FALSE]

	cmean <- colSums(tempdat) / nrow(tempdat)

	grpns <- rbind( grpns, ((last - first) + 1) )

	for (luper2 in first:last ) {
		totdev[luper2,]  <- data[luper2,2:ncol(data)] - totmn
		withdev[luper2,] <- data[luper2,2:ncol(data)] - cmean
		betwdev[luper2,] <- cmean - totmn
	}
first <- luper1
}}
grpns <- as.matrix(grpns[ 2:(nrow(grpns)),])

cormat <- cor( cbind( totdev, betwdev, withdev))

nvars <- nrow(cormat) / 3

# mean & sd 
tempdat <- data[,2:ncol(data), drop = FALSE]
mean <- colSums(tempdat) / nss
sd <- apply(tempdat, 2, FUN = sd)

cat("\n\nWABA")
cat("\n\nNumber of individual cases in the data file:",nss)
cat("\n\nNumber of groups in the data file:",ngrps)
cat("\n\nNumber of variables for the analyses:",nvars)
cat("\n\nAll significance tests are two-tailed\n")
cat("\n\nVariable means and standard deviations\n\n")
mnsd <- cbind( (mean), (sd) )
colnames(mnsd) <- cbind("   Mean","  Std Dev")
rownames(mnsd) <- varnames #colnames(data[,2:ncol(data)])
print(round(mnsd,3))


# Within- and Between-Groups Analysis of Variance 
if (nvars > 1) {
	etabetw <- as.matrix(diag( cormat[ (nvars+1):(nvars*2), 1:nvars]))
	etawith <- as.matrix(diag( cormat[ (nvars*2+1):(nvars*3),1:nvars]))
}
if (nvars == 1) {
	etabetw <- cormat[ (nvars+1):(nvars*2), 1:nvars, drop = FALSE]
	etawith <- cormat[ (nvars*2+1):(nvars*3),1:nvars, drop = FALSE]
}
etasq <- cbind( etabetw, etawith )**2
# E tests 
etests <- etabetw / etawith
# F tests 
ftrad <- (etests**2) * ( (nss - ngrps) / (ngrps - 1) )
fcrctd <- 1 / ftrad
# dfs for the F values (Dansereau et al 1984, p. 125, 128, 172) 
dfnumt <- matrix((ngrps - 1),nvars,1)
dfdemt <- matrix((nss - ngrps),nvars,1)
dfnumc <- matrix((nss - ngrps),nvars,1)
dfdemc <- matrix((ngrps - 1),nvars,1)
# sig levels for F 
pftrad  <- 1 - pf(ftrad,dfnumt[1,1],dfdemt[1,1])
pfcrctd <- 1 - pf(fcrctd,dfnumc[1,1],dfdemc[1,1])
# for displaying only the appropriate F results 
fetas <- cbind( as.matrix(1:nvars), as.matrix(ftrad), dfnumt, dfdemt, as.matrix(pftrad) )
for (lupev in 1:nvars ) {
	if ( etawith[lupev,1] > etabetw[lupev,1] ) {
	fetas[lupev,2] <- fcrctd[lupev,1]
	fetas[lupev,3] <- dfnumc[lupev,1]
	fetas[lupev,4] <- dfdemc[lupev,1]
	fetas[lupev,5] <- pfcrctd[lupev,1]
}
}

cat("\n\n\nWithin-Groups and Between-Groups Etas and Eta-squared values:\n\n")
labels <- cbind("   Eta-Betw","   Eta-With","   Eta-Bsq","   Eta-Wsq")
resebw <- cbind( etabetw, etawith, etasq )
dimnames(resebw) <-list(rep("", dim(resebw)[1]))
colnames(resebw) <- labels 
rownames(resebw) <- varnames
print (round(resebw,3))
cat("\n   The above Eta-Betw and Eta-With values for each variable are")
cat("\n   tested relative to one another with E tests of practical")
cat("\n   significance and F tests of statistical significance:\n\n")
labels <- cbind("   E-test","         F","    df num.","    df dem.","     p")
resfetas <- cbind( etests, fetas[, 2:5, drop = FALSE] )
dimnames(resfetas) <-list(rep("", dim(resfetas)[1]))
colnames(resfetas) <- labels
rownames(resfetas) <- varnames
print (round(resfetas,3))
cat("\n   E-test practical significance criteria & inductions:")
cat("\n   E >= 1.30323 = Wholes - 15      E >= 1.73205 = Wholes - 30")
cat("\n   E <= 0.76733 = Parts  - 15      E <= 0.57735 = Parts  - 30\n")

# Other indices of within-group agreement: Intraclass Correlations 
ssb <-  colSums(betwdev**2)
dfb <- ngrps - 1
msb <- ssb / dfb
ssw <-  colSums(withdev**2)
dfw <- colSums(grpns - 1)
msw <- ssw / dfw
icc1 <- as.matrix((msb - msw) / (msb+( (colSums(grpns)/nrow(grpns)) -1)*msw))
icc2 <- as.matrix((msb - msw) / msb)

cat("\n\n\nOther Indices of Within-Group Agreement: Intraclass Correlations\n\n")
labels <- cbind("     ICC(1)","     ICC(2)")
resicc <- cbind( icc1, icc2 )
dimnames(resicc) <-list(rep("", dim(resicc)[1]))
colnames(resicc) <- labels
rownames(resicc) <- varnames
print (round(resicc,3))

cat('\n\n\nThe following descriptions of ICC(1) & ICC(2) were provided by Bliese (2000).')
cat('\n\n\nICC(1):')
cat('\n\nICC(1) is the notation used by Bartko (1976), Bliese (2000),')
cat('\nJames, 1982, & McGraw & Wong, 1996). It is referred to as')
cat('\nthe ICC by Bryk & Raudenbush (1982) and as ICC(1,1) by Shrout &')
cat('\nFleiss (1979). Bryk & Raudenbush (1982) interpret it as the')
cat('\nproportion of the total variance that can be explained by group')
cat('\nmembership. For example, an ICC(1) of .10 on a leadership')
cat('\nvariable indicates that 10 percent of the variability in an')
cat('\nindividual,s ratings of leadership could be explained by or is')
cat('\nrelated to group membership. James (1982) interprets the ICC(1)')
cat('\nas the degree of reliability associated with a single assessment')
cat('\nof the group mean. That is, he interprets the ICC(1) as an index')
cat('\nof interrater reliability (the extent to which raters are')
cat('\nsubstitutable). ICC(1) is a way of estimating and contrasting the')
cat('\nbetween-group and within-group variance components from the ANOVA model.')

cat('\n\nThe range of the ICC(1) in the ANOVA model is from -l to +1.')
cat('\nNegative values occur when the within-group variance is greater')
cat('\nthan the between group variance. Negative ICC(1) values are of')
cat('\ntheoretical interest because they provide evidence of')
cat('\nwithin-group or frog-pond situations (Dansereau, Alutto,')
cat('\nYammarino, 1984). They suggest that individual variability,')
cat('\nrelative to a group mean, is an important source of variability.')
cat('\nAn example of where one might expect to find a frog-pond effect')
cat('\nwould be with pay satisfaction. One might expect the main')
cat('\npredictor of pay satisiaction to be an individual,s pay relative')
cat('\nto the average pay of his or her work group.')

cat('\n\nThe term non-independence refers to the degree to which responses')
cat('\nfrom individuals in the same group are influenced by, or depend')
cat('\non, or cluster by, group (Kenny & Judd, 1996). It may be')
cat('\nsomewhat confusing that ICC(1) can be considered both a measure')
cat('\nof reliability and a measure of non-independence. One way to help')
cat('\nclarify this distinction is to note that when ICC(1) is')
cat('\ncalculated on the dependent variable, it is generally considered')
cat('\na measure of non-independence (Bryk & Raudenbush, 1992; Kreft &')
cat('\nDeLeeuw, 1998). One is estimating ICC(1) to answer the question')
cat('\n"Is my outcome variable affected by or related to group')
cat('\nmembership? In contrast, when ICC(1) is calculated on an')
cat('\nindependent variable, it is generally being used as a measure of')
cat('\nreliability. In this case, one is typically attempting to answer')
cat('\nthe question "Can I aggregate this variable and analyze it as a')
cat('\ngroup mean?"')

cat('\n\nICC(2):')
cat('\n\nICC(2) is the notation used by Bartko (1976), Bliese (2000),')
cat('\n& James (1982). It is referred to as ICC(1,k) by Shrout &')
cat('\nFleiss (1979) and as ICC(k) by McGraw & Wong (1996). ICC(2) is an')
cat('\nestimate of the reliability of the group means.')
cat('\nThe formula is (MSB - MSW) / MSB')
cat('\nNegative values occur when MSW is greater than MSB.')
cat('\nICC(2) is often set to zero when this happens. "Negative')
cat('\nestimates are possible and can be interpreted as indicating that')
cat('\nthe true intraclass correlation is low, that is, two members')
cat('\nchosen randomly from any (group, dyad) vary almost as much as any')
cat('\ntwo randomly chosen members of the whole population" (Taylor, 2016).')

cat('\n\nICC(1) and ICC(2) are related to each other as a function of')
cat('\ngroup size (Bliese, 1998; Shrout & Fleiss, 1979).')
cat('\nThe relationship among ICC(1), ICC(2), and group size can be')
cat('\ndescribed as follows: ICC(l) may be considered a measure of the')
cat('\nreliability associated with a single assessment of the group mean')
cat('\n(James, 1982). When ICC(1) is large, a single rating from an')
cat('\nindividual is likely to provide a relatively reliable rating of')
cat('\nthe group mean; when ICC(1) is small, multiple ratings are')
cat('\nnecessary to provide reliable estimates of the group mean.\n')

cat("\nSources:\n")
cat("\n   Bartko, J. J. (1976). On various intraclass correlation reliability.")
cat("\n   Psychological Bulletin, 83, 762-765.")
cat("\n\n   Bliese, P. D. (2000). Within-group agreement, non-independence, and reliability:")
cat("\n   Implications for data aggregation and analysis. In K. J. Klein &")
cat("\n   S. W. J. Kozlowski (Eds.), Multilevel theory, research, and methods in organizations:")
cat("\n   Foundations, extensions, and new directions (pp. 349-381). San Francisco: Jossey-Bass.")
cat("\n\n   Bryk, A. S., & Raudenbush, S. W. (1982). Hiearchical linear models. Thousand Oaks, CA: Sage.")
cat("\n\n   Dansereau, F., Alutto,J. A., & Yammarino, F.J. (1984). Theory testing in organizational ")
cat("\n   behavior: The varient approach. Englewood Cliffs, NJ: Prentice Hall.")
cat("\n\n   James, L. R. (1982). Aggregation bias in estimates of perceptual agreement. ")
cat("\n   Psychological Bulletin, 67, 219-239.")
cat("\n\n   Kenny, D. A., & Judd, C. M. (1996). A general procedure for the estimation of ")
cat("\n   interdependence. Psychological Bulletin, 119, 138-148.")
cat("\n\n   Kreft, I., & Deleeuw,.J. (1993). Introducing multilevel modeling, Thousand Oaks, CA: Sage.")
cat("\n\n   McGraw, K. O., & Wong, S. P. (1996). Forming inferences about some intraclass ")
cat("\n   correlation coefficients. Psychological Methods, 1, 30-46.")
cat("\n\n   Shrout, P. E., & Fleiss, J. L. (1979). Intraclass correlations: Uses in assessing rater ")
cat("\n   reliability. Psychological Bulletin, 86, 420-428.")
cat("\n\n   Taylor (2016). www.faculty.umb.edu/peter_taylor/09b.pdf\n")

if (nvars > 1) {
	# Within- and Between-Groups Analysis of Covariance Analyses 
	rbetw <- cormat[ (nvars+1):(nvars*2),(nvars+1):(nvars*2)]
	rwith <- cormat[ (nvars*2+1):(nvars*3),(nvars*2+1):(nvars*3)]
		
	rbw <- matrix(-9999, 1, 4)
		for (luper in 1:(nvars-1) ) {
			for (lupec in (luper+1):nvars) {
			rbw <- rbind(rbw, cbind(luper,lupec,rbetw[luper,lupec],rwith[luper,lupec]) )
			}
		}
	#rbw <- t(as.matrix(rbw[2:nrow(rbw),]))
	if( nrow(rbw) > 2) rbw <- rbw[2:nrow(rbw),] else rbw <-  t(as.matrix(rbw[2:nrow(rbw),]))
	# A test 
	atests <- asin(sqrt(1-(rbw[,4]*rbw[,4])))-asin(sqrt(1-(rbw[,3]*rbw[,3])))
	# Z test 
	zbw <- matrix(-9999, nrow(rbw), 1)
	if ( ngrps > 3 ) {
		for (luper in 1:nrow(rbw)) {
			if ( abs(rbw[luper,3]) < 1  &  abs(rbw[luper,4]) < 1 ) {
			zbxy <-abs(0.5*((log(1+abs(rbw[luper,3])))-(log(1-abs(rbw[luper,3])))))
			zwxy <-abs(0.5*((log(1+abs(rbw[luper,4])))-(log(1-abs(rbw[luper,4])))))
			zbw[luper,1] <- (zbxy - zwxy) / sqrt( (1/(nss-ngrps-2))+(1/(ngrps-3)))
			}
		}
	}
	pzbw <- 1 - pnorm(abs(zbw)) 
	# R test of practical sigificance 
	rbig <- matrix(-9999,nrow(rbw), 2)
	for (luper in 1:nrow(rbw)) {
		for (lupec in 3:4) {
			if ( abs(rbw[luper,lupec]) < 1 ) {
			# rbig formula from Dansereau 1984 p 131 
			rbig[luper,(lupec-2)] <- abs( rbw[luper,lupec] / sqrt ( 1 - rbw[luper,lupec]*rbw[luper,lupec]))
			# rbig formula from Yamarino & Markham, 1992 p 170  = slightly diff 
			# rbig[luper,lupec-2] <- rbw[luper,lupec] / (((1 - (rbw[luper,lupec]*rbw[luper,lupec]) )**.6667) )
			}
		}
	}
	# t-tests of statistical significance (Yamarino & Markham, 1992 p 170) 
	tb <- as.matrix(rbig[,1] * sqrt(ngrps-2))
	dftb <- ngrps - 2
	tw <- as.matrix(rbig[,2] * sqrt(nss-ngrps-1))
	dftw <- nss - ngrps - 1
	ptb <- matrix(-9999, nrow(tb), 1 )
	ptw <- matrix(-9999, nrow(tw), 1 )
	for (luper in 1:nrow(tb)) {
		ptb[luper,1] <- (1 - pt(abs(tb[luper,1]),dftb)) * 2
		ptw[luper,1] <- (1 - pt(abs(tw[luper,1]),dftw)) * 2
	}
	
	cat("\n\n\nWithin-Groups and Between-Groups Correlations:\n\n")
	labels <- cbind("  Variable"," Variable","    r-Betw","    r-With","    A-test","    Z-test","       p")
	reswbc <- data.frame(round(cbind(rbw, atests, zbw, pzbw),3))
	reswbc[,1] <- varnames[reswbc[,1]]
	reswbc[,2] <- varnames[reswbc[,2]]
	colnames(reswbc)<-labels 
	print(reswbc, row.names = FALSE)
	cat("\n   A-tests represent practical, and Z-tests represent statistical, assessments of")
	cat("\n   the differences between the within-groups and between-groups correlations.\n")
	cat("\n   A-test practical significance criteria & inductions:")
	cat("\n   A >=  0.2618 = Wholes - 15      A >=  0.5236 = Wholes - 30")
	cat("\n   A <= -0.2618 = Parts  - 15      A <= -0.5236 = Parts  - 30\n")
	cat("\n   A value of -9999 indicates the Z test could not be computed.")


	cat("\n\n\n\nBetween-Groups Correlation Tests of Practical (R) and Statistical (t) Significance\n\n")
	labels <- cbind("  Variable"," Variable","       R","     t-test","       p")
	if( nrow(rbw) > 1) {resbcp <- data.frame(round(cbind( (rbw[,1:2]), rbig[,1], tb, ptb),3))
		} else {resbcp <- data.frame(round(cbind( t(rbw[,1:2]), rbig[,1], tb, ptb),3)) }
	resbcp[,1] <- varnames[resbcp[,1]]
	resbcp[,2] <- varnames[resbcp[,2]]
	colnames(resbcp)<-labels 
	print(resbcp, row.names = FALSE)
	cat("\n   R-test practical significance criteria & inductions:")
	cat("\n   R >=  0.26795 = S - 15      R >=  0.57735 = S - 30\n")
	cat("\n   Degrees of freedom for t-tests:", cbind(1, dftb),'\n')
	cat("\n   A value of -9999 indicates the t-test could not be computed.\n")


	cat("\n\n\nWithin-Groups Correlation Tests of Practical (R) and Statistical (t) Significance\n\n")
	labels <- cbind("  Variable"," Variable","       R","       t-test","       p")
	if( nrow(rbw) > 1) {reswcp <- data.frame(round(cbind( (rbw[,1:2]), rbig[,2], tw, ptw),3))
		} else {reswcp <- data.frame(round(cbind( t(rbw[,1:2]), rbig[,2], tw, ptw),3)) }
	reswcp[,1] <- varnames[reswcp[,1]]
	reswcp[,2] <- varnames[reswcp[,2]]
	colnames(reswcp)<-labels 
	print(reswcp, row.names = FALSE)
	cat("\n   R-test practical significance criteria & inductions:")
	cat("\n   R >=  0.26795 = S - 15      R >=  0.57735 = S - 30\n")
	cat("\n   Degrees of freedom for t-tests:", cbind( 1, dftw),'\n')
	cat("\n   A value of -9999 indicates the t-test could not be computed.\n")
	
	
	# Within- and Between-Groups Components and Raw (Total) Correlations 
	compons <- rbw
	for (luper in 1:nrow(rbw)) {
		compons[luper,3] <- etabetw[rbw[luper,1],1] * etabetw[rbw[luper,2],1] * rbw[luper,3]
		compons[luper,4] <- etawith[rbw[luper,1],1] * etawith[rbw[luper,2],1] * rbw[luper,4]
	}
	compons <- cbind( compons, (compons[,3] + compons[,4]) )
	
	# R test of practical sigificance of the total correlations 
	rbig2 <- matrix( -9999, nrow(compons), 1)
	for (luper in 1:nrow(compons)) {
		if ( abs(compons[luper,5]) < 1 ) {
			# rbig formula from Dansereau 1984 p 131 
			rbig2[luper,1] <- abs( compons[luper,5] / sqrt ( 1 - compons[luper,5]*compons[luper,5]))
			# rbig formula from Yamarino & Markham, 1992 p 170  = slightly diff 
			# rbig2(luper,lupec-2) <- rbw[luper,lupec] / (((1 - (rbw[luper,lupec]*rbw[luper,lupec]) )**.6667) ); 
		}
	}
	
	# t-tests for the total correlations Dansereau et al 1984 p 119 
	ttot <- rbig2 * sqrt(nss-2)
	dftot <- nss - 2
	ptot <- matrix( -9999,  nrow(ttot), 1 )
	for (luper in 1:nrow(ttot)) { ptot[luper,1] <- (1 - pt(abs(ttot[luper,1]),dftot)) * 2 }
	
	# A test of the between vs within components difference 
	aradians <- matrix( -9999, nrow(compons), 1)
	for (luper in 1:nrow(compons)) {
			aradians[luper,1] <- asin(sqrt(1-(compons[luper,4]*compons[luper,4]))) - 
								 asin(sqrt(1-(compons[luper,3]*compons[luper,3])))
	}
	
	# Z test of the between vs within components difference 
	zbwc <- matrix( -9999, nrow(compons), 1)
	if ( ngrps > 3 ) {
		for (luper in 1:nrow(compons)) {
			if (abs(compons[luper,3]) < 1 & abs(compons[luper,4]) < 1 ) {
				zbxyc <-abs( 0.5* ((log(1+abs(compons[luper,3]))) - (log(1-abs(compons[luper,3])))) )
				zwxyc <-abs( 0.5* ((log(1+abs(compons[luper,4]))) - (log(1-abs(compons[luper,4])))) )
				zbwc[luper,1] <- (zbxyc - zwxyc) / sqrt( (1/(nss-ngrps-2))+(1/(ngrps-3)) )
	}}}
	pzbwc <- (1 - pnorm(abs(zbwc))) * 2
	
	cat("\n\n\nWithin-Groups and Between-Groups WABA Components and Raw/Total Correlation(s):\n\n")
	labels <- cbind("  Variable"," Variable","      Betw","      With","    Raw/Tot.")
	componsDF <- data.frame(round(compons,3))
	componsDF[,1] <- varnames[componsDF[,1]]
	componsDF[,2] <- varnames[componsDF[,2]]
	colnames(componsDF) <- labels 
	print(componsDF, row.names = FALSE)
	
	
	cat("\nPractical (R) and Statistical (t) Significance of the Raw/Total correlations:\n\n")
	labels <- cbind("  Variable"," Variable","       R","      t-test","      p")
	if( nrow(compons) > 1) {resrtc <- cbind( (compons[,1:2]), rbig2, ttot, ptot )  
		} else {resrtc <- cbind( t(compons[,1:2]), rbig2, ttot, ptot ) }
	resrtc <- data.frame(round(resrtc,3))
	resrtc[,1] <- varnames[resrtc[,1]]
	resrtc[,2] <- varnames[resrtc[,2]]
	colnames(resrtc) <- labels 
	print(resrtc, row.names = FALSE)
	cat("\n   R-test practical significance criteria & inductions:")
	cat("\n   R >=  0.26795 = S - 15      R >=  0.57735 = S - 30\n")
	cat("\n   Degrees of freedom for t-tests:", cbind(1, dftot))
	cat("\n\n   A value of -9999 indicates the t-test could not be computed.\n")


	cat("\n\nThe Within-Groups and Between-Groups WABA Components are tested relative to")
	cat("\none another with A tests of practical Z tests of statistical significance.\n\n")	
	labels <- cbind("  Variable"," Variable","       A-test","       Z-test","       p")
	if( nrow(compons) > 1) {resaz <- cbind( (compons[,1:2]), aradians, zbwc, pzbwc)
		} else {resaz <- cbind( t(compons[,1:2]), aradians, zbwc, pzbwc) }
	resaz <- data.frame(round(resaz,3))
	resaz[,1] <- varnames[resaz[,1]]
	resaz[,2] <- varnames[resaz[,2]]
	colnames(resaz) <- labels 
	print(resaz, row.names = FALSE)
	cat("\n   A-test practical significance criteria & inductions:")
	cat("\n   A >=  0.2618 = Wholes - 15      A >=  0.5236 = Wholes - 30")
	cat("\n   A <= -0.2618 = Parts  - 15      A <= -0.5236 = Parts  - 30\n")
	cat("\n   A value of -9999 indicates the Z test could not be computed.\n\n\n")

}

# Multiple Variable Analysis 
#if (nvars != 3 ) {cat("\n\n\n\nMultiple Variable Analysis was not conducted because the number of variables was not = 3\n")}

if ( nvars == 3 ) {
cat("\n\nMultiple Variable Analysis Using Within-Groups Correlations\n")

cat("\nSources:\n")
cat("\n   Dansereau, Alutto, & Yammarino (1984). Theory Testing in Organizational Behavior. Prentice-Hall.\n")
cat("\n   Yammarino (1998). Multivariate aspects of the varient/WABA approach. Leadership Quarterly, 9, 203-227.\n")
cat("\n   www.LevelsOfAnalysis.com\n")
cat("\n\nCorrelations are compared using both A tests of practical")
cat("\nsignificance & Hotelling t-tests of statistical significance.\n")
cat("\nA-test values > 15 degrees are considered practically significant.\n")

dum <- rbw

for (luper in 1:3) {

	rxy <-  dum[1,4]
	rxz <-  dum[2,4]
	ryz <-  dum[3,4]

	cat("\n\nWithin correlations:\n\n")
	labels <- cbind("  Variable"," Variable","        r")
	reswc <- rbind( cbind( (dum[1,1]),(dum[1,2]),(rxy) ), cbind( (dum[2,1]),(dum[2,2]),(rxz) )  )
	reswc <- data.frame(round(reswc,3))
	reswc[,1] <- varnames[reswc[,1]]
	reswc[,2] <- varnames[reswc[,2]]
	colnames(reswc) <- labels 
	print(reswc, row.names = FALSE)

	# "A" test of the difference between the correlations: Dansereau et al, 1984, p 140 
	atest <- asin(abs(sqrt(1-(rxy*rxy)))) - asin(abs(sqrt(1-(rxz*rxz))))
	adegs <- atest * 57.29578
	# Hotelling's t-test of the difference between the correlations: Dansereau et al, 1984, p 141 
	hotest <- cbind( abs(rxy) - abs(rxz) ) * sqrt(  ( (ngrps-2)*(1+abs(ryz)) ) /
         ( 2 * ( 1 - ryz*ryz - rxy*rxy - rxz*rxz + 2 * abs(ryz) *
         abs(rxy) * abs(rxz))) )
	dfhotest <- ngrps-2
	photest <- (1 - pt(abs(hotest),dfhotest) ) * 2
	cat("\nA test (in radians & degrees)  &  Hotelling's t-test:\n\n")
	labels <- cbind( "   A-rads","   A-degs","        t","     df","     p")
	resah <- cbind(abs(atest), abs(adegs), hotest, dfhotest, photest)
	dimnames(resah) <-list(rep("", dim(resah)[1]))
	colnames(resah) <- labels
	print (round(resah,3))

	# rotating the rows of the correlation data matrix 
	dum <- rbind( dum[3,], dum )
	dum <- dum[1:3,]
}

for (luper in 1:3) {

	rxy <-  dum[1,3]
	rxz <-  dum[2,3]
	ryz <-  dum[3,3]

	cat("\n\n\nBetween correlations:\n\n")
	labels <- cbind("  Variable"," Variable","     r")
	resbc <- rbind( cbind((dum[1,1]), (dum[1,2]), (rxy) ), cbind( (dum[2,1]), (dum[2,2]), (rxz)))
	resbc <- data.frame(round(resbc,3))
	resbc[,1] <- varnames[resbc[,1]]
	resbc[,2] <- varnames[resbc[,2]]
	colnames(resbc) <- labels 
	print(resbc, row.names = FALSE)

	# "A" test of the difference between the correlations: Dansereau et al, 1984, p 140 
	atest <- asin(abs(sqrt(1-(rxy*rxy)))) - asin(abs(sqrt(1-(rxz*rxz))))
	adegs <- atest * 57.29578
	# Hotelling's t-test of the difference between the correlations: Dansereau et al, 1984, p 141 
	hotest <- ( abs(rxy) - abs(rxz) ) * sqrt(  ( (ngrps-3)*(1+abs(ryz)) ) /
        ( 2 * ( 1 - ryz*ryz - rxy*rxy - rxz*rxz + 2 *
        abs(ryz) * abs(rxy) * abs(rxz))) )
	dfhotest <- ngrps-3
	photest <- (1 - pt(abs(hotest),dfhotest) ) * 2
	cat("\nA test (in radians & degrees)  &  Hotelling's t-test:\n\n")
	labels <- cbind( "    A-rads","    A-degs","        t","     df","     p")
	resdc <- cbind( abs(atest), abs(adegs), hotest, dfhotest, photest )
	dimnames(resdc) <-list(rep("", dim(resdc)[1]))
	colnames(resdc) <- labels
	print (round(resdc,3))

	# rotating the rows of the correlation data matrix 
	dum <- rbind( dum[3,], dum )
	dum <- dum[1:3,]
}
}

} 

