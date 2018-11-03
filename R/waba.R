

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

if ( all(complete.cases(data)) == 'FALSE' ) {
cat("\n\nERROR: There are missing values in the data matrix. Fix this problem and try again.\n\n") }

data <- as.matrix(data)

data <- data[order(data[,1]),] # sorting by 1st column (the group numbers/codes)

# removing cases for groups that have just one person/case
grpsize <- 1
dumped <- matrix(-9999,1,ncol(data))
for (luper1 in 2:nrow(data) ) {
	if ( ( data[luper1,1] == data[(luper1-1),1] ) )                   { grpsize = grpsize + 1 }
	if ( ( data[luper1,1] != data[(luper1-1),1] ) & grpsize == 1 )    { 
	    dumped <- rbind( dumped, data[(luper1-1),] )
		data[(luper1-1),2:ncol(data)] <- rep(NA, (ncol(data)-1)) }
	if ( ( data[luper1,1] != data[(luper1-1),1] ) & grpsize > 1 )     { grpsize = 1 }
	if (luper1 == nrow(data) & data[luper1,1] != data[(luper1-1),1] ) { data[(luper1),2:ncol(data)] <- rep(NA, (ncol(data)-1))  }
}
if (nrow(dumped) > 1) {
	dumped <- dumped[2:nrow(dumped),]
	dimnames(dumped) <-list(rep("", dim(dumped)[1]))
	colnames(dumped) <- colnames(data)
	cat("\n\nThe following cases were removed from the data file because there were no other cases with the same group code:\n\n" )
	print(dumped)
}

data <- na.omit(data)

nss <- nrow(data)

totmn <- colSums(data[,2:ncol(data)]) / nss

totdev  <- matrix( -9999, nss, ncol(data)-1 ) 
withdev <- matrix( -9999, nss, ncol(data)-1 )
betwdev <- matrix( -9999, nss, ncol(data)-1 )

grpns   <- matrix( -9999, 1, 1 )

first <- 1
ngrps <- 0

for (luper1 in 2:nss ) {
	if ( ( data[luper1,1] != data[(luper1-1),1] ) | luper1 == nss ) {
		if ( luper1 != nss ) { last <- luper1 - 1 }
		if ( luper1 == nss ) { last <- luper1 }

	ngrps <- ngrps + 1

	tempdat <- data[first:last,2:ncol(data)]

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

cormat <- cor( cbind( totdev, betwdev, withdev ) )

nvars <- nrow(cormat) / 3

# mean & sd 
tempdat <- data[,2:ncol(data)]
mean <- colSums(tempdat) / nss
sd <- apply(tempdat, 2, FUN = sd)

cat("\n\nWABA")
cat("\n\nNumber of Individual Cases in the Data File:",nss)
cat("\n\nNumber of Groups in the Data File:",ngrps)
cat("\n\nNumber of Variables for the Analyses:",nvars)
cat("\n\nAll Significance Tests Are Two-Tailed\n")
cat("\n\nVariable Means and Standard Deviations\n\n")
mnsd <- cbind( (mean), (sd) )
varnames <- rownames(mnsd) # for later use
colnames(mnsd) <- cbind("Mean","Std Dev")
rownames(mnsd) <- colnames(data[,2:ncol(data)])
print(round(mnsd,3))


# Within- and Between-Groups Analysis of Variance 
etabetw <- as.matrix(diag( cormat[ (nvars+1):(nvars*2), 1:nvars]))
etawith <- as.matrix(diag( cormat[ (nvars*2+1):(nvars*3),1:nvars]))
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
labels <- cbind("Variable","Eta-Betw","Eta-With","Eta-Bsq","Eta-Wsq")
resebw <- cbind( (1:nvars), etabetw, etawith, etasq )
dimnames(resebw) <-list(rep("", dim(resebw)[1]))
colnames(resebw) <- labels 
rownames(resebw) <- varnames
print (round(resebw,3))
cat("\nThe above Eta-Betw and Eta-With values for each variable are\n")
cat("\ntested relative to one another with E tests of practical\n")
cat("\nsignificance and F tests of statistical significance:\n\n")
labels <- cbind("Variable","E-test","F","df num.","df dem.","p")
resfetas <- cbind( fetas[,1], etests, fetas[, 2:5] )
dimnames(resfetas) <-list(rep("", dim(resfetas)[1]))
colnames(resfetas) <- labels
rownames(resfetas) <- varnames
print (round(resfetas,3))
cat("\n   E-test practical significance criteria & inductions:\n")
cat("\n   E >= 1.30323 = Wholes - 15;     E >= 1.73205 = Wholes - 30;\n")
cat("\n   E <= 0.76733 = Parts  - 15;     E <= 0.57735 = Parts  - 30;\n")

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
labels <- cbind("Variable","ICC(1)","ICC(2)")
resicc <- cbind( (1:nvars), icc1, icc2 )
dimnames(resicc) <-list(rep("", dim(resicc)[1]))
colnames(resicc) <- labels
rownames(resicc) <- varnames
print (round(resicc,3))
cat("\n   Sources:\n")
cat("\n   Bliese (2000), in Klein et al, Multilevel theory, research, and methods in organizations.\n")
cat("\n   Shrout & Fleiss, 1979, Psychological Bulletin, 86, 420-428.\n")
cat("\n   McGraw & Wong, 1996, Psychological Methods, 1, 30-46.\n")

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
if( nrow(rbw) > 2) rbw = rbw[2:nrow(rbw),] else rbw =  t(as.matrix(rbw[2:nrow(rbw),]))
# A test 
atests <- asin(sqrt(1-(rbw[,4]*rbw[,4])))-asin(sqrt(1-(rbw[,3]*rbw[,3])))
# Z test 
zbw <- matrix(-9999, nrow(rbw), 1)
if ( ngrps > 3 ) {
	for (luper in 1:nrow(rbw)) {
		if ( abs(rbw[luper,3]) < 1  &  abs(rbw[luper,4]) < 1 ) {
		zbxy =abs(0.5*((log(1+abs(rbw[luper,3])))-(log(1-abs(rbw[luper,3])))))
		zwxy =abs(0.5*((log(1+abs(rbw[luper,4])))-(log(1-abs(rbw[luper,4])))))
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
		rbig[luper,(lupec-2)] <- abs( rbw[luper,lupec] / sqrt ( 1 - rbw[luper,lupec]*rbw[luper,lupec] ) )
		# compute rbig formula from Yamarino & Markham, 1992 p 170  = slightly diff 
		# compute rbig[luper,lupec-2] <- rbw[luper,lupec] / ( ( (1 - (rbw[luper,lupec]*rbw[luper,lupec]) )**.6667) ); 
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
labels <- cbind("Var.#","Var.#","r-Betw","r-With","A-test","Z-test","p")
reswbc <- cbind(rbw, atests, zbw, pzbw) 
dimnames(reswbc) <-list(rep("", dim(reswbc)[1]))
colnames(reswbc)=labels 
print (round(reswbc,3))
cat("\n   A-tests represent practical, and Z-tests represent statistical, assessments of")
cat("\n   the differences between the within-groups and between-groups correlations.\n")
cat("\n   A-test practical significance criteria & inductions:\n")
cat("\n   A >=  0.2618 = Wholes - 15;     A >=  0.5236 = Wholes - 30;\n")
cat("\n   A <= -0.2618 = Parts  - 15;     A <= -0.5236 = Parts  - 30;\n")
cat("\n   -9999 indicates the Z test could not be computed.")
cat("\n\n\n\nBetween-Groups Correlation Tests of Practical (R) and Statistical (t) Significance\n\n")
labels <- cbind("Var.#","Var.#","R","t-test","p")
#resbcp <- cbind( t(rbw[,1:2]), rbig[,1], tb, ptb )
#resbcp <- cbind( (rbw[,1:2]), rbig[,1], tb, ptb )
if( nrow(rbw) > 1) resbcp = cbind( (rbw[,1:2]), rbig[,1], tb, ptb )  else resbcp = cbind( t(rbw[,1:2]), rbig[,1], tb, ptb )
dimnames(resbcp) <-list(rep("", dim(resbcp)[1]))
colnames(resbcp) <- labels
print (round(resbcp,3))
cat("\n   R-test practical significance criteria & inductions:\n")
cat("\n   R >=  0.26795 = S - 15;     R >=  0.57735 = S - 30;\n")
cat("\n   Degrees of Freedom for t-tests:", cbind(1, dftb))
cat("\n   A value of -9999 indicates the t-test could not be computed.\n")
cat("\n\n\nWithin-Groups Correlation Tests of Practical (R) and Statistical (t) Significance\n\n")
labels <- cbind("Var.#","Var.#","R","t-test","p")
#reswcp <- cbind( t(rbw[,1:2]), rbig[,2], tw, ptw )
#reswcp <- cbind( (rbw[,1:2]), rbig[,2], tw, ptw )
if( nrow(rbw) > 1) reswcp = cbind( (rbw[,1:2]), rbig[,2], tw, ptw )  else reswcp = cbind( t(rbw[,1:2]), rbig[,2], tw, ptw )
dimnames(reswcp) <- list(rep("", dim(reswcp)[1]))
colnames(reswcp) <- labels
print (round(reswcp,3))
cat("\n   R-test practical significance criteria & inductions:\n")
cat("\n   R >=  0.26795 = S - 15;     R >=  0.57735 = S - 30;\n")
cat("\n   Degrees of Freedom for t-tests:", cbind( 1, dftw) )
cat("\n\n   A value of -9999 indicates the t-test could not be computed.\n")


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
		rbig2[luper,1] <- abs( compons[luper,5] / sqrt ( 1 - compons[luper,5]*compons[luper,5] ) )
		# rbig formula from Yamarino & Markham, 1992 p 170  = slightly diff 
		# rbig2(luper,lupec-2) <- rbw[luper,lupec] / ( ( (1 - (rbw[luper,lupec]*rbw[luper,lupec]) )**.6667) ); 
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
		aradians[luper,1] <- asin(sqrt(1-(compons[luper,4]*compons[luper,4]))) - asin(sqrt(1-(compons[luper,3]*compons[luper,3])))
}

# Z test of the between vs within components difference 
zbwc <- matrix( -9999, nrow(compons), 1)
if ( ngrps > 3 ) {
	for (luper in 1:nrow(compons)) {
		if (abs(compons[luper,3]) < 1 & abs(compons[luper,4]) < 1 ) {
			zbxyc =abs( 0.5* ((log(1+abs(compons[luper,3]))) - (log(1-abs(compons[luper,3])))) )
			zwxyc =abs( 0.5* ((log(1+abs(compons[luper,4]))) - (log(1-abs(compons[luper,4])))) )
			zbwc[luper,1] <- (zbxyc - zwxyc) / sqrt( (1/(nss-ngrps-2))+(1/(ngrps-3)) )
}}}
pzbwc <- (1 - pnorm(abs(zbwc))) * 2

cat("\n\n\nWithin-Groups and Between-Groups WABA Components and Raw/Total Correlation(s):\n\n")
labels <- cbind("Var.#","Var.#","Betw","With","Raw/Tot.")
dimnames(compons) <-list(rep("", dim(compons)[1]))
colnames(compons) <- labels
print (round(compons,3))

cat("\nPractical (R) and Statistical (t) Significance of the Raw/Total correlations:\n\n")
labels <- cbind("Var.#","Var.#","R","t-test","p")
#resrtc <- cbind( (compons[,1:2]), rbig2, ttot, ptot )
if( nrow(compons) > 1) resrtc <- cbind( (compons[,1:2]), rbig2, ttot, ptot )  else resrtc <- cbind( t(compons[,1:2]), rbig2, ttot, ptot )
dimnames(resrtc) <-list(rep("", dim(resrtc)[1]))
colnames(resrtc) <- labels
print (round(resrtc,3))

cat("\n   R-test practical significance criteria & inductions:\n")
cat("\n   R >=  0.26795 = S - 15;     R >=  0.57735 = S - 30;\n")
cat("\n   Degrees of Freedom for t-tests:", cbind(1, dftot))
cat("\n\n   A value of -9999 indicates the t-test could not be computed.\n")
cat("\n\nThe Within-Groups and Between-Groups WABA Components are tested relative to")
cat("\none another with A tests of practical Z tests of statistical significance.\n\n")

labels <- cbind("Var.#","Var.#","A-test","Z-test","p")
#resaz <- cbind( (compons[,1:2]), aradians, zbwc, pzbwc)
if( nrow(compons) > 1) resaz <- cbind( (compons[,1:2]), aradians, zbwc, pzbwc)  else resaz <- cbind( t(compons[,1:2]), aradians, zbwc, pzbwc)
dimnames(resaz) <-list(rep("", dim(resaz)[1]))
colnames(resaz) <- labels
print (round(resaz,3))
cat("\n   A-test practical significance criteria & inductions:\n")
cat("\n   A >=  0.2618 = Wholes - 15;     A >=  0.5236 = Wholes - 30;\n")
cat("\n   A <= -0.2618 = Parts  - 15;     A <= -0.5236 = Parts  - 30;\n")
cat("\n   -9999 indicates the Z test could not be computed.\n\n\n")



# Multiple Variable Analysis 
if (nvars != 3 ) {cat("\n\n\n\nMultiple Variable Analysis was not conducted because the number of variables was not = 3\n")}

if ( nvars == 3 ) {
cat("\n\nMultiple Variable Analysis Using Within-Groups Correlations\n")

cat("\nSources:\n")
cat("\nDansereau, Alutto, & Yammarino (1984). Theory Testing in Organizational Behavior. Prentice-Hall.\n")
cat("\nYammarino (1998). Multivariate aspects of the varient/WABA approach. Leadership Quarterly, 9, 203-227.\n")
cat("\nwww.LevelsOfAnalysis.com\n")
cat("\nCorrelations are compared using both A tests of practical")
cat("\nsignificance & Hotelling t-tests of statistical significance.\n")
cat("\nA-test values > 15 degrees are considered practically significant.\n")

dum <- rbw

for (luper in 1:3) {

	rxy <-  dum[1,4]
	rxz <-  dum[2,4]
	ryz <-  dum[3,4]

	cat("\nWithin correlations:\n")
	labels <- cbind("Var.#","Var.#","r")
	reswc <- rbind( cbind( (dum[1,1]),(dum[1,2]),(rxy) ), cbind( (dum[2,1]),(dum[2,2]),(rxz) )  )
	dimnames(reswc) <-list(rep("", dim(reswc)[1]))
	colnames(reswc) <- labels
	print (round(reswc,3))
	# "A" test of the difference between the correlations: Dansereau et al, 1984, p 140 
	atest <- asin(abs(sqrt(1-(rxy*rxy)))) - asin(abs(sqrt(1-(rxz*rxz))))
	adegs <- atest * 57.29578
	# Hotelling's t-test of the difference between the correlations: Dansereau et al, 1984, p 141 
	hotest <- cbind( abs(rxy) - abs(rxz) ) * sqrt(  ( (ngrps-2)*(1+abs(ryz)) ) /
         ( 2 * ( 1 - ryz*ryz - rxy*rxy - rxz*rxz + 2 * abs(ryz) *
         abs(rxy) * abs(rxz) ) ) )
	dfhotest <- ngrps-2
	photest <- (1 - pt(abs(hotest),dfhotest) ) * 2
	cat("\nA test (in radians & degrees)  &  Hotelling's t-test:\n\n")
	labels <- cbind( "A-rads","A-degs","t","df","p")
	resah <- cbind(abs(atest), abs(adegs), hotest, dfhotest, photest)
	dimnames(resah) <-list(rep("", dim(resah)[1]))
	colnames(resah) <- labels
	print (round(resah,3))

	# rotating the rows of the correlation data matrix 
	dum <- rbind( dum[3,], dum )
	dum <- dum[1:3,]
}

cat("\n\n\n\nMultiple Variable Analysis Using Between-Groups Correlations\n")

cat("\nSources:\n")
cat("\nDansereau, Alutto, & Yammarino (1984). Theory Testing in Organizational Behavior. Prentice-Hall.\n")
cat("\nYammarino (1998). Multivariate aspects of the varient/WABA approach. Leadership Quarterly, 9, 203-227.\n")
cat("\nwww.LevelsOfAnalysis.com\n")
cat("\nCorrelations are compared using both A tests of practical")
cat("\nsignificance & Hotelling t-tests of statistical significance.\n")
cat("\nA-test values > 15 degrees are considered practically significant.\n")
dum <- rbw

for (luper in 1:3) {

	rxy <-  dum[1,3]
	rxz <-  dum[2,3]
	ryz <-  dum[3,3]

	cat("\nBetween correlations:\n")
	labels <- cbind("Var.#","Var.#","r")
	resbc <- rbind( cbind((dum[1,1]), (dum[1,2]), (rxy) ), cbind( (dum[2,1]), (dum[2,2]), (rxz) ) )
	dimnames(resbc) <-list(rep("", dim(resbc)[1]))
	colnames(resbc) <- labels
	print (round(resbc,3))
	# "A" test of the difference between the correlations: Dansereau et al, 1984, p 140 
	atest <- asin(abs(sqrt(1-(rxy*rxy)))) - asin(abs(sqrt(1-(rxz*rxz))))
	adegs <- atest * 57.29578

	# Hotelling's t-test of the difference between the correlations: Dansereau et al, 1984, p 141 
	hotest <- ( abs(rxy) - abs(rxz) ) * sqrt(  ( (ngrps-3)*(1+abs(ryz)) ) /
        ( 2 * ( 1 - ryz*ryz - rxy*rxy - rxz*rxz + 2 *
        abs(ryz) * abs(rxy) * abs(rxz) ) ) )
	dfhotest <- ngrps-3
	photest <- (1 - pt(abs(hotest),dfhotest) ) * 2
	cat("\nA test (in radians & degrees)  &  Hotelling's t-test:\n\n")
	labels <- cbind( "A-rads","A-degs","t","df","p")
	resdc <- cbind( abs(atest), abs(adegs), hotest, dfhotest, photest )
	dimnames(resdc) <-list(rep("", dim(resdc)[1]))
	colnames(resdc) <- labels
	print (round(resdc,3))

	# rotating the rows of the correlation data matrix 
	dum <- rbind( dum[3,], dum )
	dum <- dum[1:3,]
}
}

} # end of function WABA



