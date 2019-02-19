

wabamodreg <- function (data) {

#  WABA - Moderated Regression Analysis 

   # Sources:

   # Schriesheim, C. (1995). Multivariate and moderated within-and between-entity 
   # analysis (WABA) using hierarchical linear multiple regression. Leadership Quarterly, 6, 1-18.

   # Schriesheim, C. A., Cogliser, C. C., & Neider, L. L. (1995). Is it "trustworthy"? 
   # a multiple levels-of-analysis reexamination of an ohio state leadership study, with
   # implications for future research. Leadership Quarterly, 6, 111-145.

   # Schriesheim C., Neider L.L., Scadura T. (1998), Delegation and leader-member
   # exchange: main effects, moderators, and measurement issues. Academy of Management
   # Journal, 41(3), 298-318.

   # This WABA program focuses on the interaction between two variables
   # in the prediction of a third.  It is a WABA extension of regular
   # moderated multiple regression.  All variables must be continuous. 

   # Prepare a raw data matrix for analysis,
   # where rows = cases, & columns = variables
   # Each row contains the data from a single individual.
   # Cases with missing values are not permitted in the data file. 

   # The first value in each row (i.e., the first column of values in
   # the data file) must be the individuals' group numbers/codes,
   # which must be integers. The program sorts individuals into groups
   # on the basis of these numbers/codes.
   # Variable scores appear in subsequent columns. 

   # The second value in each row (i.e., the second column of values
   # in the data file) are treated by the program as the
   # Dependent/Outcome variable. 

   # The third value in each row (i.e., the third column of values in
   # the data file) are treated by the program as the
   # Independent/Predictor variable. 

   # The fourth value in each row (i.e., the fourth column of values
   # in the data file) are treated by the program as the Moderator
   # variable (technically also an IDV/predictor). 

   # There is no need to enter product term variables or to enter
   # composite variable scores. The program takes care of these
   # aspects of the analyses internally. 

if ( all(complete.cases(data)) == 'FALSE' ) {
	cat("\n\nERROR: There are missing values in the data matrix. Fix this problem and try again") }

data <- as.matrix(data)

data <- data[order(data[,1]),] # sorting by 1st column (the group numbers/codes)

# removing cases for groups that have just one person/case
grpsize <- 1
dumped <- matrix(-9999,1,ncol(data))
for (luper1 in 2:nrow(data) ) {
	if ( ( data[luper1,1] == data[(luper1-1),1] ) )                   grpsize = grpsize + 1 
	if ( ( data[luper1,1] != data[(luper1-1),1] ) & grpsize == 1 )    { 
	    dumped <- rbind( dumped, data[(luper1-1),] )
		data[(luper1-1),2:ncol(data)] <- rep(NA, (ncol(data)-1)) }
	if ( ( data[luper1,1] != data[(luper1-1),1] ) & grpsize > 1 )     grpsize = 1 
	if (luper1 == nrow(data) & data[luper1,1] != data[(luper1-1),1] ) data[(luper1),2:ncol(data)] <- rep(NA, (ncol(data)-1))  
}
if (nrow(dumped) > 1) {
	dumped <- dumped[2:nrow(dumped),]
	dimnames(dumped) <-list(rep("", dim(dumped)[1]))
	colnames(dumped) <- colnames(data)
	cat("\n\nThe following cases were removed from the data file because there were no other cases with the same group code:\n\n" )
	print(dumped)
}

data <- na.omit(data)


if (ncol(data) > 4) {
	cat("\nThe are more than four columns in the data file.")
	cat("\nThe analyses will be conducted using the first four columns, as per the instructions.\n")
	data <- data[,1:4] 
}

# specifying/creating the data for each set of analyses 
# creating the IDV/Moderator Composite Variable 
x <- as.matrix(cbind( matrix(1,nrow(data),1), data[,3:4] ))
y <- data[,2]
b <- solve(t(x) %*% x) %*% t(x) %*% y
comp1 <- data[,3] * b[2,1] + data[,4] * b[3,1]
data <- cbind( data, comp1 )
# creating the IDV/Moderator Product Term Composite Variable 
product <- data[,3] * data[,4]
x <- as.matrix(cbind( matrix(1,nrow(data),1), data[,3:4], product ))
y <- as.matrix(data[,2])
b <- solve(t(x) %*% x) %*% t(x) %*% y
comp2 <- data[,3] * b[2,1] + data[,4] * b[3,1] + product * b[4,1]
data <- cbind( data, comp2 )

# start of regular WABA analyses 

nss <- nrow(data)

totmn <- colSums(data[,2:ncol(data)]) / nss

totdev  <- matrix( -9999,  nss, (ncol(data)-1))
withdev <- matrix( -9999,  nss, ncol(data)-1)
betwdev <- matrix( -9999,  nss, ncol(data)-1)

grpns    <- matrix( -9999,  1, 1)

first <- 1
ngrps <- 0

for (luper1 in 2:nss) {
	if ( data[(luper1),1] != data[(luper1-1),1]  |  luper1 == nss ) {
		if ( luper1 != nss ) last <- luper1 - 1
		if ( luper1 == nss ) last <- luper1

		ngrps <- ngrps + 1

		tempdat <- data[first:last,2:ncol(data)]
		cmean <- colSums(tempdat) / nrow(tempdat)

		grpns <- cbind( grpns, (last - first) + 1 )

		for (luper2 in first:last) {
			totdev[luper2,]  <- as.matrix(data[luper2,2:ncol(data)] - totmn)
			withdev[luper2,] <- as.matrix(data[luper2,2:ncol(data)] - cmean)
			betwdev[luper2,] <- cmean - totmn
		}
		first <- luper1
	}
}

cormat <- cor( cbind( totdev, betwdev, withdev ) )

nvars <- nrow(cormat) / 3

# mean & sd 
tempdat <- data[,2:4]
mean <- colSums(tempdat) / nss
sd <- apply(tempdat, 2, FUN = sd)

cat("\nWABA -- Moderated Regression\n")
cat("\n\nNumber of individual cases in the data file:",nss)
cat("\n\nNumber of groups in the data file:",ngrps)
cat("\n\nNumber of variables for the analyses:",nvars)
cat("\n\nAll significance tests are two-tailed\n")
cat("\n\nVariable means and standard deviations\n\n")
labels <- cbind("    Mean","   Std Dev")
resdes <- cbind( mean, sd )
colnames(resdes) <- labels
print (round(resdes,3))
rns <- as.matrix(row.names(resdes))
cat("\n\nThe dependent variable is:  ",rns[1,1])
cat("\n\nThe independent variable is:",rns[2,1])
cat("\n\nThe moderator variable is:  ",rns[3,1])

# Within- and Between-Groups Analysis of Variance 
etabetw <- as.matrix(diag( cormat[ (nvars+1):(nvars*2), 1:nvars]))
etawith <- as.matrix(diag( cormat[ (nvars*2+1):(nvars*3), 1:nvars]))
etasq <- cbind( etabetw, etawith )**2
# E tests 
etests <- etabetw / etawith
# dfs for the F values (Dansereau et al 1984, p. 125, 128, 172) 
dfnumt <- matrix((ngrps - 1),5,1)
dfdemt <- matrix((nss - ngrps),5,1)
dfnumc <- matrix((nss - ngrps),5,1)
dfdemc <- matrix((ngrps - 1),5,1)
# adjusting the dfs for the composite variables (Shriesheim 1995 p 11) 
dfdemt[4,1] <- dfdemt[4,1] - 1
dfdemt[5,1] <- dfdemt[5,1] - 2
dfdemc[4,1] <- dfdemc[4,1] - 1
dfdemc[5,1] <- dfdemc[5,1] - 2
# F tests 
ftrad <- (etests**2) * (dfdemt / dfnumt)
fcrctd <- ((etawith /etabetw)**2) * (dfdemc / dfnumc)
# sig levels for F 
pftrad  <- 1 - pf(ftrad,dfnumt[1,1],dfdemt[1,1])
pfcrctd <- 1 - pf(fcrctd,dfnumc[1,1],dfdemc[1,1])
# for displaying only the appropriate F results 
fetas <- cbind( (1:nvars), ftrad, dfnumt, dfdemt, pftrad )
for (lupev in 1:nvars) {
	if ( etawith[lupev,1] > etabetw[lupev,1] ) {
	fetas[lupev,2] <- fcrctd[lupev,1]
	fetas[lupev,3] <- dfnumc[lupev,1]
	fetas[lupev,4] <- dfdemc[lupev,1]
	fetas[lupev,5] <- pfcrctd[lupev,1]
	}
}

# Within- and Between-Groups Analysis of Covariance Analyses 
rbetw <- cormat[ (nvars+1):(nvars*2)  , (nvars+1):(nvars*2)]
rwith <- cormat[ (nvars*2+1):(nvars*3), (nvars*2+1):(nvars*3)]
rbw <- matrix( -9999, 1, 4)
for (luper in 2:5) {rbw <- rbind(rbw,  cbind(1, luper, rbetw[1,luper], rwith[1,luper]) ) }
rbw <- rbw[2:nrow(rbw),]
# A test 
atests <- as.matrix(asin(sqrt(1-(rbw[,4]*rbw[,4])))-asin(sqrt(1-(rbw[,3]*rbw[,3]))))
# Z test 
zbw <- matrix( -9999, 4, 1)
if ( ngrps > 3 ) {
	for (luper in 1:4) {
		if ( abs(rbw[luper,3]) < 1  &  abs(rbw[luper,4]) < 1 ) {
			zbxy  <- abs( 0.5* ((log(1+abs(rbw[luper,3]))) - (log(1-abs(rbw[luper,3])))) )
			zwxy  <- abs( 0.5* ((log(1+abs(rbw[luper,4]))) - (log(1-abs(rbw[luper,4])))) )
			zbw[luper,1] <- (zbxy - zwxy) / sqrt( (1/(nss-ngrps-2))+(1/(ngrps-3)) )
				if ( luper == 3 ) zbw[luper,1] <- (zbxy - zwxy) / sqrt( (1/(nss-ngrps-2))+(1/(ngrps-3-1)) )
				if ( luper == 4 ) zbw[luper,1] <- (zbxy - zwxy) / sqrt( (1/(nss-ngrps-2))+(1/(ngrps-3-2)) )
		}
	}
}
pzbw <- 1 - pnorm(abs(zbw)) 
# R test of practical sigificance 
rbig <- matrix( -9999, nrow(rbw), 2)
for (luper in 1:nrow(rbw)) {
	for (lupec in 3:4) {
		if ( abs(rbw[luper,lupec]) < 1 ) {
			# rbig formula from Dansereau 1984 p 131 
			rbig[luper,lupec-2] <- abs( rbw[luper,lupec] / sqrt ( 1 - rbw[luper,lupec]*rbw[luper,lupec] ) )
			# rbig formula from Yamarino & Markham, 1992 p 170  = slightly diff 
			# rbig[luper,lupec-2]  <-  rbw[luper,lupec] / ( ( (1 - (rbw[luper,lupec]*rbw[luper,lupec]] )**.6667) ); 
		}
	}
}
# t-tests of statistical significance (Yamarino & Markham, 1992 p 170) 
dftb <- matrix((ngrps - 2),4,1)
dftw <- matrix((nss - ngrps - 1),4,1)
# adjusting the dfs for the composite variables (Shriesheim 1995 p 11) 
dftb[3,1] <- dftb[3,1] - 1
dftb[4,1] <- dftb[4,1] - 2
tb <- rbig[,1] * sqrt(dftb)
tw <- rbig[,2] * sqrt(dftw)
ptb <- matrix( -9999,  nrow(tb), 1 )
ptw <- matrix( -9999,  nrow(tw), 1 )
for (luper in 1:nrow(tb)) {
	ptb[luper,1] <- (1 - pt(abs(tb[luper,1]),dftb[luper,1])) * 2
	ptw[luper,1] <- (1 - pt(abs(tw[luper,1]),dftw[luper,1])) * 2
}

# Within- and Between-Groups Components and Raw (Total) Correlations 
compons <- rbw
for (luper in 1:nrow(rbw)) {
	compons[luper,3] <- etabetw[rbw[luper,1],1] * etabetw[rbw[luper,2],1] * rbw[luper,3]
	compons[luper,4] <- etawith[rbw[luper,1],1] * etawith[rbw[luper,2],1] * rbw[luper,4]
}
compons <- cbind( compons, (compons[,3] + compons[,4]) )
# A test of the between vs within components difference 
aradians <- asin(sqrt(1-(compons[,4]*compons[,4]))) - asin(sqrt(1-(compons[,3]*compons[,3])))
# Z test of the between vs within components difference 
zbwc <- matrix( -9999, nrow(compons), 1)
if ( ngrps > 3 ) {
	for (luper in 1:nrow(compons)) {
		if ( abs(compons[luper,3]) < 1  &  abs(compons[luper,4]) < 1 ) {
			zbxyc  <- abs( 0.5* ((log(1+abs(compons[luper,3]))) - (log(1-abs(compons[luper,3])))) )
			zwxyc  <- abs( 0.5* ((log(1+abs(compons[luper,4]))) - (log(1-abs(compons[luper,4])))) )
			zbwc[luper,1] <- (zbxyc - zwxyc) / sqrt( (1/(nss-ngrps-2))+(1/(ngrps-3)) )
				if ( luper == 3 ) zbw[luper,1] <- (zbxy - zwxy) / sqrt( (1/(nss-ngrps-2))+(1/(ngrps-3-1)) )
				if ( luper == 4 ) zbw[luper,1] <- (zbxy - zwxy) / sqrt( (1/(nss-ngrps-2))+(1/(ngrps-3-2)) )
		}
	}
}
pzbwc <- (1 - pnorm(abs(zbwc))) * 2

cat("\n\n\nWABA 1 Results -- Dependent Variable:\n\n")
labels <- cbind("Eta-Betw","   Eta-With","   E Test","      F","   df num","   dfdem","      p")
resw1dv <- cbind( cbind( etabetw[1,1], etawith[1,1], etests[1,1]), rbind(fetas[1,2:5]) )
dimnames(resw1dv) <-list(rep("", dim(resw1dv)[1]))
colnames(resw1dv) <- labels 
print (round(resw1dv,3))

cat("\nWABA 1 Results -- Independent Variable:\n\n")
labels <- cbind("Eta-Betw","   Eta-With","   E Test","      F","   df num","   dfdem","       p")
resw1iv <- cbind( cbind( etabetw[2,1], etawith[2,1], etests[2,1]), rbind(fetas[2,2:5]) )
dimnames(resw1iv) <-list(rep("", dim(resw1iv)[1]))
colnames(resw1iv) <- labels
print (round(resw1iv,3))

cat("\nWABA 1 Results -- Moderator Variable:\n\n")
labels <- cbind("Eta-Betw","   Eta-With","   E Test","      F","   df num","   dfdem","       p")
resw1m <- cbind( cbind( etabetw[3,1], etawith[3,1], etests[3,1]), rbind(fetas[3,2:5]) )
dimnames(resw1m) <-list(rep("", dim(resw1m)[1]))
colnames(resw1m) <- labels
print (round(resw1m,3))

cat("\nWABA 1 Results -- IDV+Moderator Composite Variable:\n\n")
labels <- cbind("Eta-Betw","   Eta-With","   E Test","      F","   df num","   dfdem","       p")
resw1ivm <- cbind( cbind(etabetw[4,1], etawith[4,1], etests[4,1]), rbind(fetas[4,2:5]) )
dimnames(resw1ivm) <-list(rep("", dim(resw1ivm)[1]))
colnames(resw1ivm) <- labels
print (round(resw1ivm,3))

cat("\nWABA 1 Results -- IDV+Moderator+Product Term Composite Variable:\n\n")
labels <- cbind("Eta-Betw","   Eta-With","   E Test","      F","   df num","   dfdem","       p")
resw1ivmp <- cbind( cbind(etabetw[5,1], etawith[5,1], etests[5,1]), rbind(fetas[5,2:5]) )
dimnames(resw1ivmp) <-list(rep("", dim(resw1ivmp)[1]))
colnames(resw1ivmp) <- labels
print (round(resw1ivmp,3))
cat("\n     E-test practical significance criteria & inductions:")
cat("\n     E >= 1.30323 = Wholes - 15     E >= 1.73205 = Wholes - 30")
cat("\n     E <= 0.76733 = Parts  - 15     E <= 0.57735 = Parts  - 30\n")

cat("\n\nWABA 2 Results: Independent Variable Predicting the Dependent Variable:\n")
labels <- cbind("r-Betw","   r-With","   A-test","   Z-test","       p")
resw2iv <- cbind(  rbind(rbw[1,3:4]), cbind(atests[1,1], zbw[1,1],  pzbw[1,1]) )
dimnames(resw2iv) <-list(rep("", dim(resw2iv)[1]))
colnames(resw2iv) <- labels
cat("\n")
print (round(resw2iv,3))
labels <- cbind("r-Betw","      R","   t-test","   df","       p")
resw2b <- cbind( rbw[1,3], rbig[1,1], tb[1,1], dftb[1,1], ptb[1,1] )
colnames(resw2b) <- labels
cat("\n")
print (round(resw2b,3))
labels <- cbind("r-With","      R","   t-test","   df","       p")
resw2w <- cbind( rbw[1,4], rbig[1,2], tw[1,1], dftw[1,1], ptw[1,1] )
colnames(resw2w) <- labels
cat("\n")
print (round(resw2w,3))

cat("\nWABA 2 Results: Moderator Variable Predicting the Dependent Variable:\n")
labels <- cbind("r-Betw","   r-With","   A-test","   Z-test","       p")
resw2m <- cbind(  rbind(rbw[2,3:4]), cbind(atests[2,1], zbw[2,1], pzbw[2,1]) )
dimnames(resw2m) <-list(rep("", dim(resw2m)[1]))
colnames(resw2m) <- labels
cat("\n")
print (round(resw2m,3))
labels <- cbind("r-Betw","      R","   t-test","   df","       p")
resw2b <- cbind( rbw[2,3], rbig[2,1], tb[2,1], dftb[2,1], ptb[2,1] )
colnames(resw2b) <- labels
cat("\n")
print (round(resw2b,3))
labels <- cbind("r-With","      R","   t-test","   df","       p")
resw2w <- cbind( rbw[2,4], rbig[2,2], tw[2,1], dftw[2,1], ptw[2,1] )
colnames(resw2w) <- labels
cat("\n")
print (round(resw2w,3))

cat("\nWABA 2 Results: IDV+Moderator Composite Variable Predicting the Dependent Variable:\n")
labels <- cbind("r-Betw","   r-With","   A-test","   Z-test","       p")
resw2im <- cbind(  rbind(rbw[3,3:4]), cbind(atests[3,1], zbw[3,1], pzbw[3,1]) )
dimnames(resw2im) <-list(rep("", dim(resw2im)[1]))
colnames(resw2im) <- labels
cat("\n")
print (round(resw2im,3))
labels <- cbind("r-Betw","      R","   t-test","   df","       p")
resw2b <- cbind ( rbw[3,3], rbig[3,1], tb[3,1], dftb[3,1], ptb[3,1] )
colnames(resw2b) <- labels
cat("\n")
print (round(resw2b,3))
labels <- cbind("r-With","      R","   t-test","   df","       p")
resw2w <- cbind( rbw[3,4], rbig[3,2], tw[3,1], dftw[3,1], ptw[3,1] )
colnames(resw2w) <- labels
cat("\n")
print (round(resw2w,3))

cat("\nWABA 2 Results: IDV+Moderator+Product Term Composite Variable Predicting the Dependent Variable:\n")
labels <- cbind("r-Betw","   r-With","   A-test","   Z-test","       p")
resw2imp <- cbind(  rbind(rbw[4,3:4]), cbind(atests[4,1], zbw[4,1], pzbw[4,1]) )
dimnames(resw2imp) <-list(rep("", dim(resw2imp)[1]))
colnames(resw2imp) <- labels
cat("\n")
print (round(resw2imp,3))
labels <- cbind("r-Betw","      R","   t-test","   df","       p")
resw2b <- cbind( rbw[4,3], rbig[4,1], tb[4,1], dftb[4,1], ptb[4,1] )
colnames(resw2b) <- labels
cat("\n")
print (round(resw2b,3))
labels <- cbind("r-With","      R","   t-test","   df","       p")
resw2w <- cbind( rbw[4,4], rbig[4,2], tw[4,1], dftw[4,1], ptw[4,1] )
colnames(resw2w) <- labels
cat("\n")
print (round(resw2w,3))

cat("\n\nA-test practical significance criteria & inductions:")
cat("\n    A >=  0.2618 = Wholes - 15     A >=  0.5236 = Wholes - 30")
cat("\n    A <= -0.2618 = Parts  - 15     A <= -0.5236 = Parts  - 30\n")
cat("\n-9999 indicates the Z test could not be computed.\n")
cat("\nR-test practical significance criteria & inductions:")
cat("\n    R >=  0.26795 = S - 15     R >=  0.57735 = S - 30")
cat("\nA value of -9999 indicates the t-test could not be computed.\n")

cat("\n\nComponents: Independent Variable Predicting the Dependent Variable:\n\n")
labels <- cbind("Betw","   With","   A-test","   Z-test","       p")
rescom <- cbind(  rbind(compons[1,3:4]), cbind(atests[1,1], zbw[1,1],  pzbw[1,1]) )
dimnames(rescom) <-list(rep("", dim(rescom)[1]))
colnames(rescom) <- labels
print (round(rescom,3))

cat("\nComponents: Moderator Variable Predicting the Dependent Variable:\n\n")
labels <- cbind("Betw","   With","   A-test","   Z-test","       p")
rescomm <- cbind( rbind(compons[2,3:4]), cbind(atests[2,1], zbw[2,1], pzbw[2,1]) )
dimnames(rescomm) <-list(rep("", dim(rescomm)[1]))
colnames(rescomm) <- labels
print (round(rescomm,3))

cat("\nComponents: IDV+Moderator Composite Variable Predicting the Dependent Variable:\n\n")
labels <- cbind("Betw","   With","   A-test","   Z-test","       p")
rescomb <- cbind(  rbind(compons[3,3:4]), cbind(atests[3,1], zbw[3,1],  pzbw[3,1]) )
dimnames(rescomb) <-list(rep("", dim(rescomb)[1]))
colnames(rescomb) <- labels
print (round(rescomb,3))

cat("\nComponents: IDV+Moderator+Product Term Composite Variable Predicting the Dependent Variable:\n\n")
labels <- cbind("Betw","   With","   A-test","   Z-test","       p")
rescomim <- cbind(  rbind(compons[4,3:4]), cbind(atests[4,1], zbw[4,1],  pzbw[4,1]) )
dimnames(rescomim) <-list(rep("", dim(rescomim)[1]))
colnames(rescomim) <- labels
print (round(rescomim,3))

cat("\nA-test practical significance criteria & inductions:")
cat("\n    A >=  0.2618 = Wholes - 15     A >=  0.5236 = Wholes - 30")
cat("\n    A <= -0.2618 = Parts  - 15     A <= -0.5236 = Parts  - 30")
cat("\n\n-9999 indicates the Z test could not be computed.\n")

cat("\n\n\nWABA - Moderated Hierarchical Regression Results for the Within Correlations:\n")
# Step 2 
r2idv  <- rbw[1,4] * rbw[1,4]
r2idvm <- rbw[3,4] * rbw[3,4]
Fidvm <- (r2idvm-r2idv) / ((1-r2idvm)/(nss-2-1))
dfnmidvm <- 2
dfdmidvm <- nss - 2 - 1
pidvm <- 1 - pf(abs(Fidvm),dfnmidvm,dfdmidvm)
cat("\nStep 2: Increment for the moderator variable, as a main effect, over the IDV")
cat("\n        i.e., for the IDV+Moderator Composite Variable over the IDV:\n")
cat("\n Step 1  Step 2 Incr.\n")
labels <- cbind("r-With","   r-With","   rsq ch","      F","   df num","   df dem","       p")
resmhr <- cbind( rbw[1,4], rbw[3,4], (r2idvm-r2idv), Fidvm, dfnmidvm, dfdmidvm, pidvm )
colnames(resmhr) <- labels
print (round(resmhr,3))

# Step 3 
r2xn  <- rbw[4,4] * rbw[4,4]
Fxn <- (r2xn-r2idvm) / ((1-r2xn)/(nss-3-1))
dfnmxn <- 3
dfdmxn <- nss - 3 - 1
pxn <- 1 - pf(abs(Fidvm),dfnmxn,dfdmxn)
cat("\n\nStep 3: Increment for the IDV+Moderator+Product Term Composite Variable over")
cat("\n        the two main effects (i.e., over the IDV+Moderator Composite Variable):\n")
cat("\n Step 1  Step 2 Incr.\n")
labels <- cbind("r-With","   r-With","   rsq ch","      F","   df num","   df dem","       p")
resinc <- cbind( rbw[3,4], rbw[4,4], (r2xn-r2idvm), Fxn, dfnmxn, dfdmxn, pxn )
colnames(resinc) <- labels
print (round(resinc,3))

cat("\n\n\nWABA - Moderated Hierarchical Regression Results for the Between Correlations:\n")
# Step 2 
r2idv  <- rbw[1,3] * rbw[1,3]
r2idvm <- rbw[3,3] * rbw[3,3]
Fidvm <- (r2idvm-r2idv) / ((1-r2idvm)/(ngrps-2-1))
dfnmidvm <- 2
dfdmidvm <- ngrps - 2 - 1
pidvm <- 1 - pf(abs(Fidvm),dfnmidvm,dfdmidvm)
cat("\nStep 2: Increment for the moderator variable, as a main effect, over the IDV")
cat("\n        i.e., for the IDV+Moderator Composite Variable over the IDV:\n")
labels <- cbind("r-Betw","   Betw","   rsq ch","      F","   df num","   df dem","       p")
cat("\n Step 1  Step 2 Incr.\n")
resincm <- cbind( rbw[1,3], rbw[3,3], (r2idvm-r2idv), Fidvm, dfnmidvm, dfdmidvm, pidvm )
colnames(resincm) <- labels
print (round(resincm,3))

# Step 3 
r2xn  <- rbw[4,3] * rbw[4,3]
Fxn <- (r2xn-r2idvm) / ((1-r2xn)/(ngrps-3-1))
dfnmxn <- 3
dfdmxn <- ngrps - 3 - 1
pxn <- 1 - pf(abs(Fidvm),dfnmxn,dfdmxn)
cat("\n\nStep 3: Increment for the IDV+Moderator+Product Term Composite Variable over")
cat("\n        the two main effects (i.e., over the IDV+Moderator Composite Variable):\n")
cat("\n Step 1  Step 2 Incr.\n")
labels <- cbind("r-Betw","   Betw","   rsq ch","      F","   df num","   df dem","       p")
resincim <- cbind( rbw[3,3], rbw[4,3], (r2xn-r2idvm), Fxn, dfnmxn, dfdmxn, pxn)
colnames(resincim) <- labels
print (round(resincim,3))

} # end of function wabamodreg 


