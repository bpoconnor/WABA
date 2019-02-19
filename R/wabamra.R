

wabamra <- function (data) {


# WABA - Multiple Relationship Analysis 

# Sources (for interpretation of the results): 

# Dansereau, Alutto, & Yammarino (1984). Theory
# Testing in Organizational Behavior. Prentice-Hall
   
# Yammarino (1998). Multivariate aspects of the varient/WABA 
# approach. Leadership Quarterly, 9, 203-227

# www.LevelsOfAnalysis.com 

# This program conducts separate WABA analyses for each value
# of a specified "Condition" variable.  The program then conducts
# Multiple Relationship Analysis comparisons of the WABA
# coefficients for the different Condition values. 

# Prepare a raw data matrix for analysis, where rows = cases, &
# columns = variables.  Each row contains the data from a single
# individual. Cases with missing values are not permitted in the data file. 

# The first value in each row (i.e., the first column of values
# in the data file) must be the Condition number.  Condition
# numbers must be integers; the lowest Condition number cannot
# be less than one; it is also best for there to be no missing
# values between the lowest and highest Condition numbers
# e.g., if the lowest value is 1 and the highest value is 5,
# then there should also be Condition values of 2, 3, and 4
# somewhere in the data file. Gaps in the integers may cause problems. 

# The second value in each row (i.e., the second column of values
# in the data file) must be the individual's group number/code
# The program sorts individuals into groups on the basis of these
# numbers/codes. 

# Variable scores appear in subsequent columns. 


if ( all(complete.cases(data)) == 'FALSE' ) {
	cat("\n\nERROR: There are missing values in the data matrix. Fix this problem and try again") }

cat("\n\nWABA - Multiple Relationship Analysis\n")


# removing cases for groups that have just one person/case
grpsize <- 1
dumped <- matrix(-9999,1,ncol(data))
for (luper1 in 2:nrow(data) ) {
	if ( ( data[luper1,2] == data[(luper1-1),2] ) )   grpsize = grpsize + 1 
	if ( ( data[luper1,2] != data[(luper1-1),2] ) & grpsize == 1 )    { 
	    dumped <- rbind( dumped, data[(luper1-1),] )
		data[(luper1-1),3:ncol(data)] <- rep(NA, (ncol(data)-2)) }
	if ( ( data[luper1,2] != data[(luper1-1),2] ) & grpsize > 1 )  grpsize = 1 
	if (luper1 == nrow(data) & data[luper1,2] != data[(luper1-1),2] )  
		{data[(luper1),3:ncol(data)] <- rep(NA, (ncol(data)-2))}  
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

datacond <- as.matrix(data)

datacond <- datacond[order(datacond[,1],datacond[,2]),] # sorting by 1st & 2nd columns 

if (min(datacond[,1]) == 0) { 
	datacond[,1] <- datacond[,1] + 1
	cat("\n\nA value of zero for Condition was detected.")
	cat("\nAll Condition values were therefore incremented by 1.\n\n")
}

mincond <- min(datacond[,1])
maxcond <- max(datacond[,1])

rbwdf <- matrix( -9999, maxcond, 5)

# grand loop, for doing one condition at a time 
for (grand in mincond:maxcond) {

cat("\n\nResults for Condition #", grand, "\n")

data <- matrix( -9999, 1,3)
for (luper in 1:nrow(datacond)) {
	if ( datacond[luper,1] == grand ) data <- rbind( data, datacond[luper,2:4] )
}
data <- data[2:nrow(data),]

# regular WABA goes here 

nss <- nrow(data)

totmn <- colSums(data[,2:ncol(data)]) / nss

totdev  <- matrix( -9999, nss, ncol(data)-1 ) 
withdev <- matrix( -9999, nss, ncol(data)-1 )
betwdev <- matrix( -9999, nss, ncol(data)-1 )

grpns   <- matrix( -9999, 1, 1 )

first <- 1
ngrps <- 0

for (luper1 in 2:nss ) {
	if ((data[luper1,1] != data[(luper1-1),1] ) | luper1 == nss ) {
		if ( luper1 != nss ) last <- luper1 - 1 
		if ( luper1 == nss ) last <- luper1 
		ngrps <- ngrps + 1

		tempdat <- data[first:last,2:ncol(data)]
		cmean <- colSums(tempdat) / nrow(tempdat)

		grpns <- rbind( grpns, ((last - first) + 1))

		for (luper2 in first:last ) {
			totdev[luper2,]  <- data[luper2,2:ncol(data)] - totmn
			withdev[luper2,] <- data[luper2,2:ncol(data)] - cmean
			betwdev[luper2,] <- cmean - totmn
		}
	first <- luper1
	}
}

cormat <- cor( cbind( totdev, betwdev, withdev))

nvars <- nrow(cormat) / 3

# mean & sd 
tempdat <- data[,2:ncol(data)]
mean <- colSums(tempdat) / nss
sd <- apply(tempdat, 2, FUN = sd)

cat("\nWABA")
cat("\n\nNumber of individual cases in the data file:",nss)
cat("\n\nNumber of groups in the data file:",ngrps)
cat("\n\nNumber of variables for the analyses:",nvars)
cat("\n\nAll significance tests are two-tailed\n")
cat("\n\nVariable means and standard deviations\n\n")
mnsd <- cbind((mean), (sd))
varnames <- rownames(mnsd) # for later use
rownames(mnsd) <- varnames 
colnames(mnsd) <- cbind("Mean","Std Dev")
print(round(mnsd,3))

# Within- and Between-Groups Analysis of Variance 
etabetw <- as.matrix(diag( cormat[ (nvars+1):(nvars*2), 1:nvars]))
etawith <- as.matrix(diag( cormat[ (nvars*2+1):(nvars*3),1:nvars]))
etasq <- cbind( etabetw, etawith )**2
# E tests 
etests <- etabetw / etawith
# F tests 
ftrad <- (etests**2) * ((nss - ngrps) / (ngrps - 1))
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
fetas <- cbind( as.matrix(1:nvars), as.matrix(ftrad), dfnumt, dfdemt, as.matrix(pftrad))
for (lupev in 1:nvars ) {
	if ( etawith[lupev,1] > etabetw[lupev,1] ) {
	fetas[lupev,2] <- fcrctd[lupev,1]
	fetas[lupev,3] <- dfnumc[lupev,1]
	fetas[lupev,4] <- dfdemc[lupev,1]
	fetas[lupev,5] <- pfcrctd[lupev,1]
	}
}

cat("\n\nWithin-Groups and Between-Groups Etas and Eta-squared values:\n\n")
#labels <- cbind("Variable","Eta-Betw","Eta-With","Eta-Bsq","Eta-Wsq")
labels <- cbind("Eta-Betw","Eta-With","Eta-Bsq","Eta-Wsq")
#resebw <- cbind((1:nvars), etabetw, etawith, etasq )
resebw <- cbind( etabetw, etawith, etasq )
dimnames(resebw) <-list(rep("", dim(resebw)[1]))
colnames(resebw) <- labels 
rownames(resebw) <- varnames
print (round(resebw,3))
cat("\nThe above Eta-Betw and Eta-With values for each variable are\n")
cat("\ntested relative to one another with E tests of practical\n")
cat("\nsignificance and F tests of statistical significance:\n\n")
labels <- cbind("E-test","F","df num.","df dem.","p")
#resfetas <- cbind( fetas[,1], etests, fetas[, 2:5] )
resfetas <- cbind( etests, fetas[, 2:5] )
dimnames(resfetas) <-list(rep("", dim(resfetas)[1]))
colnames(resfetas) <- labels
rownames(resfetas) <- varnames
print (round(resfetas,3))
cat("\n   E-test practical significance criteria & inductions:")
cat("\n   E >= 1.30323 = Wholes - 15     E >= 1.73205 = Wholes - 30")
cat("\n   E <= 0.76733 = Parts  - 15     E <= 0.57735 = Parts  - 30\n")

# Other indices of within-group agreement: Intraclass Correlations 
grpns <- as.matrix(grpns[ 2:(nrow(grpns)),])
ssb <-  colSums(betwdev**2)
dfb <- ngrps - 1
msb <- ssb / dfb
ssw <-  colSums(withdev**2)
dfw <- colSums(grpns - 1)
msw <- ssw / dfw
icc1 <- as.matrix((msb - msw) / (msb+((colSums(grpns)/nrow(grpns)) -1)*msw))
icc2 <- as.matrix((msb - msw) / msb)

cat("\n\nOther Indices of Within-Group Agreement: Intraclass Correlations\n\n")
labels <- cbind("ICC(1)","ICC(2)")
resicc <- cbind( icc1, icc2 )
dimnames(resicc) <-list(rep("", dim(resicc)[1]))
colnames(resicc) <- labels
rownames(resicc) <- varnames
print (round(resicc,3))
cat("\n   Sources:\n")
cat("\n   Bliese (2000), in Klein et al, Multilevel theory, research, and methods in organizations.\n")
cat("\n   Shrout & Fleiss, 1979, Psychological Bulletin, 86, 420-428.\n")
cat("\n   McGraw & Wong, 1996, Psychological Methods, 1, 30-46.\n\n")

# Within- and Between-Groups Analysis of Covariance Analyses 
rbetw <- cormat[ (nvars+1):(nvars*2),(nvars+1):(nvars*2)]
rwith <- cormat[ (nvars*2+1):(nvars*3),(nvars*2+1):(nvars*3)]
rbw <- matrix(-9999, 1, 4)
for (luper in 1:(nvars-1)) {
	for (lupec in (luper+1):nvars) {rbw <- rbind(rbw, cbind(luper,lupec,rbetw[luper,lupec],rwith[luper,lupec]))
	}
}
rbw <- t(rbw[2:nrow(rbw),])
# A test 
atests <- asin(sqrt(1-(rbw[,4]*rbw[,4])))-asin(sqrt(1-(rbw[,3]*rbw[,3])))
# Z test 
zbw <- matrix(-9999, nrow(rbw), 1)
if ( ngrps > 3 ) {
	for (luper in 1:nrow(rbw)) {
		if ( abs(rbw[luper,3]) < 1  &  abs(rbw[luper,4]) < 1 ) {
			zbxy  <- abs(0.5*((log(1+abs(rbw[luper,3])))-(log(1-abs(rbw[luper,3])))))
			zwxy  <- abs(0.5*((log(1+abs(rbw[luper,4])))-(log(1-abs(rbw[luper,4])))))
			zbw[luper,1] <- (zbxy - zwxy) / sqrt((1/(nss-ngrps-2))+(1/(ngrps-3)))
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
if( nrow(rbw) > 1) {resbcp <- data.frame(round(cbind((rbw[,1:2]), rbig[,1], tb, ptb),3))
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
if( nrow(rbw) > 1) {reswcp <- data.frame(round(cbind((rbw[,1:2]), rbig[,2], tw, ptw),3))
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
compons <- cbind( compons, (compons[,3] + compons[,4]))
# R test of practical sigificance of the total correlations 
rbig2 <- matrix( -9999, nrow(compons), 1)
for (luper in 1:nrow(compons)) {
	if ( abs(compons[luper,5]) < 1 ) {
		# rbig formula from Dansereau 1984 p 131 
		rbig2[luper,1] <- abs( compons[luper,5] / sqrt ( 1 - compons[luper,5]*compons[luper,5]))
	}
}
# t-tests for the total correlations Dansereau et al 1984 p 119 
ttot <- rbig2 * sqrt(nss-2)
dftot <- nss - 2
ptot <- matrix( -9999,  nrow(ttot), 1 )
for (luper in 1:nrow(ttot)) { ptot[luper,1] <- (1 - pt(abs(ttot[luper,1]),dftot)) * 2
}
# A test of the between vs within components difference 
aradians <- asin(sqrt(1-(compons[,4]*compons[,4]))) - asin(sqrt(1-(compons[,3]*compons[,3])))
# Z test of the between vs within components difference 
zbwc <- matrix( -9999, nrow(compons), 1)
if ( ngrps > 3 ) {
	for (luper in 1:nrow(compons)) {
		if ( abs(compons[luper,3]) < 1 & abs(compons[luper,4]) < 1 ) {
			zbxyc  <- abs( 0.5* ((log(1+abs(compons[luper,3]))) - (log(1-abs(compons[luper,3])))))
			zwxyc  <- abs( 0.5* ((log(1+abs(compons[luper,4]))) - (log(1-abs(compons[luper,4])))))
			zbwc[luper,1] <- (zbxyc - zwxyc) / sqrt((1/(nss-ngrps-2))+(1/(ngrps-3)))
		}
	}
}
pzbwc <- (1 - pnorm(abs(zbwc))) * 2

	labels <- cbind("  Variable"," Variable","      Betw","      With","    Raw/Tot.")
	componsDF <- data.frame(round(compons,3))
	componsDF[,1] <- varnames[componsDF[,1]]
	componsDF[,2] <- varnames[componsDF[,2]]
	colnames(componsDF) <- labels 
	print(componsDF, row.names = FALSE)


	cat("\nPractical (R) and Statistical (t) Significance of the Raw/Total correlations:\n\n")
	labels <- cbind("  Variable"," Variable","       R","      t-test","      p")
	if( nrow(compons) > 1) {resrtc <- cbind((compons[,1:2]), rbig2, ttot, ptot )  
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
	if( nrow(compons) > 1) {resaz <- cbind((compons[,1:2]), aradians, zbwc, pzbwc)
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

rbwdf[grand,] <- cbind( t(rbw[1,3:4]), dftb, dftw, grand )

}



# Multiple Relationship Analysis 
cat("\n\n\nMultiple Relationship Analysis\n")
cat("\nSources:")
cat("\nDansereau, Alutto, & Yammarino (1984). Theory Testing in Organizational Behavior. Prentice-Hall.")
cat("\nYammarino (1998). Multivariate aspects of the varient/WABA approach. Leadership Quarterly, 9, 203-227.")
cat("\nwww.LevelsOfAnalysis.com")

# loops for comparing each condition with every other condition 
for (luper1 in 1:(nrow(rbwdf)-1)) {
	for (luper2 in (luper1+1):nrow(rbwdf)) {

		dumm <- rbind( rbwdf[luper1,], rbwdf[luper2,] )

		cat("\n\nComparisons of Coefficients From The Following Two Conditions:", dumm[1,5], "&", dumm[2,5] )

		cat("\n\n\nComparison of the two Between correlations:\n\n") 
		labels <- cbind("Cond.#","    r-Betw.","    Cond.#","    r-Betw.")
		res2b <- cbind(dumm[1,5],dumm[1,1],dumm[2,5],dumm[2,1])
		dimnames(res2b) <-list(rep("", dim(res2b)[1]))
		colnames(res2b) <- labels
		print (round(res2b,3))
			if ( abs(dumm[1,1]) < 1  &  abs(dumm[2,1]) < 1 ) {

			# "A" test of the difference between the correlations; Dansereau et al, 1984, p 142 
			atest <- asin(abs(sqrt(1-(dumm[1,1]*dumm[1,1])))) - asin(abs(sqrt(1-(dumm[2,1]*dumm[2,1]))))
			adegs <- atest * 57.29578

			# Z test of the difference between the correlations; Dansereau et al, 1984, p 142 
			z1  <- abs(0.5* ((log(1+abs(dumm[1,1]))) - (log(1-abs(dumm[1,1])))))
			z2  <- abs(0.5* ((log(1+abs(dumm[2,1]))) - (log(1-abs(dumm[2,1])))))
			zbb <- (z1 - z2) / sqrt((1/(dumm[1,3]-1))+(1/(dumm[2,3]-1)))
			pzbb <- (1 - pnorm(abs(zbb))) * 2 
			cat("\n     A test (in radians & degrees)  &  z-test:\n\n")
			labels <- cbind( "A-rads","    A-degs","       z","       p")
			resat <- cbind( abs(atest), abs(adegs), zbb, pzbb )
			dimnames(resat) <-list(rep("", dim(resat)[1]))
			colnames(resat) <- labels
			print (round(resat,3))
			}
		if (abs(dumm[1,1]) == 1 | abs(dumm[2,1]) == 1) cat("\nA and Z tests cannot be computed") 

		cat("\n\nComparison of the two Within correlations:\n\n") 
		labels <- cbind("Cond.#","    r-With.","    Cond.#","    r-With.")
		res2w <- cbind(dumm[1,5],dumm[1,2],dumm[2,5],dumm[2,2])
		dimnames(res2w) <-list(rep("", dim(res2w)[1]))
		colnames(res2w) <- labels
		print (round(res2w,3))
		if ( abs(dumm[1,2]) < 1  &  abs(dumm[2,2]) < 1 ) {
			# "A" test of the difference between the correlations; Dansereau et al, 1984, p 142 
			atest <- asin(abs(sqrt(1-(dumm[1,2]*dumm[1,2])))) -
        asin(abs(sqrt(1-(dumm[2,2]*dumm[2,2]))))
			adegs <- atest * 57.29578
			# Z test of the difference between the correlations; Dansereau et al, 1984, p 142 
			z1  <- abs(0.5* ((log(1+abs(dumm[1,2]))) - (log(1-abs(dumm[1,2])))))
			z2  <- abs(0.5* ((log(1+abs(dumm[2,2]))) - (log(1-abs(dumm[2,2])))))
			zww <- (z1 - z2) / sqrt((1/(dumm[1,4]-1))+(1/(dumm[2,4]-1)))
			pzww <- (1 - pnorm(abs(zww))) * 2 
			cat("\n     A test (in radians & degrees)  &  z-test:\n\n")
			labels <- cbind( "A-rads","    A-degs","       z","        p")
			res2aw <- cbind( abs(atest), abs(adegs), zww, pzww )
			dimnames(res2aw) <-list(rep("", dim(res2aw)[1]))
			colnames(res2aw) <- labels
			print (round(res2aw,3))
}
		if ((abs(dumm[1,2]) == 1 | abs(dumm[2,2]) == 1)) { cat("\nA and Z tests cannot be computed") }

		cat("\n\nComparison of the Between & Within correlations:\n\n") 
		labels <- cbind("Cond.#","    r-Betw.","    Cond.#","    r-With.")
		resbw <- cbind(dumm[1,5],dumm[1,1],dumm[2,5],dumm[2,2])
		dimnames(resbw) <-list(rep("", dim(resbw)[1]))
		colnames(resbw) <- labels
		print (round(resbw,3))
		if ( abs(dumm[1,1]) < 1  &  abs(dumm[2,2]) < 1 ) {
			# "A" test of the difference between the correlations; Dansereau et al, 1984, p 142 
			atest <- asin(abs(sqrt(1-(dumm[1,1]*dumm[1,1])))) - asin(abs(sqrt(1-(dumm[2,2]*dumm[2,2]))))
			adegs <- atest * 57.29578
			# Z test of the difference between the correlations; Dansereau et al, 1984, p 142 
			z1  <- abs(0.5* ((log(1+abs(dumm[1,1]))) - (log(1-abs(dumm[1,1])))))
			z2  <- abs(0.5* ((log(1+abs(dumm[2,2]))) - (log(1-abs(dumm[2,2])))))
			zbw <- (z1 - z2) / sqrt((1/(dumm[1,3]-1))+(1/(dumm[2,4]-1)))
			pzbw <- (1 - pnorm(abs(zbw))) * 2 
			cat("\n     A test (in radians & degrees)  &  z-test:\n\n")
			labels <- cbind( "A-rads","    A-degs","       z","        p")
			resbwc <- cbind( abs(atest), abs(adegs), zbw, pzbw )
			dimnames(resbwc) <-list(rep("", dim(resbwc)[1]))
			colnames(resbwc) <- labels
			print (round(resbwc,3))
		}
		if ((abs(dumm[1,1]) == 1 | abs(dumm[2,2]) == 1)) cat("\nA and Z tests cannot be computed") 

		cat("\n\nComparison of the other Between & Within correlations:\n\n") 
		labels <- cbind("Cond.#","    r-Betw.","    Cond.#","    r-With.")
		resobw <- cbind(dumm[2,5],dumm[2,1],dumm[1,5],dumm[1,2])
		dimnames(resobw) <-list(rep("", dim(resobw)[1]))
		colnames(resobw) <- labels
		print (round(resobw,3))
		if ( abs(dumm[2,1]) < 1  &  abs(dumm[1,2]) < 1 ) {
			# "A" test of the difference between the correlations; Dansereau et al, 1984, p 142 
			atest <- asin(abs(sqrt(1-(dumm[2,1]*dumm[2,1])))) - asin(abs(sqrt(1-(dumm[1,2]*dumm[1,2]))))
			adegs <- atest * 57.29578
			# Z test of the difference between the correlations; Dansereau et al, 1984, p 142 
			z1  <- abs(0.5* ((log(1+abs(dumm[2,1]))) - (log(1-abs(dumm[2,1])))))
			z2  <- abs(0.5* ((log(1+abs(dumm[1,2]))) - (log(1-abs(dumm[1,2])))))
			zbw <- (z1 - z2) / sqrt((1/(dumm[2,3]-1))+(1/(dumm[1,4]-1)))
			pzbw <- (1 - pnorm(abs(zbw))) * 2 
			cat("\n     A test (in radians & degrees)  &  z-test:\n\n")
			labels <- cbind( "A-rads","    A-degs","       z","        p")
			resobwa <- cbind( abs(atest), abs(adegs), zbw, pzbw )
			dimnames(resobwa) <-list(rep("", dim(resobwa)[1]))
			colnames(resobwa) <- labels
			print (round(resobwa,3))
		}
		if ((abs(dumm[2,1]) == 1 | abs(dumm[1,2]) == 1)) cat("\nA and Z tests cannot be computed") 

}
}

} # end of function wabamra







