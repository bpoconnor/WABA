

ggdist <- function (data) {


# The Correlational Analysis of Dyad-Level Data in the Distinguishable-Case

# Gonzalez & Griffin, 1999, Personal Relationships, 6, 449-469 

# Create a data matrix where rows = cases, & columns = variables
# the first collumn of data must contain the dyad numbers
# (the dyad id# given to both dyad members); the second collumn
# contains the person id #s (e.g., sex) coded as 1 or 2
# the third and subsequent collumns contain the variables to be analyzed.

# Cases with missing values are not permitted in the data file

varnames <- colnames(data)[3:ncol(data)]

if ( all(complete.cases(data)) == 'FALSE' ) {
	cat("\n\nERROR: There are missing values in the data matrix. Fix this problem and try again") 
}

data <- as.matrix(data)

data <- data[order(data[,1]),] # sorting by 1st column (the group numbers/codes)


# removing cases for dyads that have just one person/case
grpsize <- 1
dumped <- matrix(-9999,1,ncol(data))
for (luper1 in 2:nrow(data) ) {
	if ( ( data[luper1,1] == data[(luper1-1),1] ) )                   { grpsize = grpsize + 1 }
	if ( ( data[luper1,1] != data[(luper1-1),1] ) & grpsize == 1 )    { 
 	   dumped <- rbind( dumped, data[(luper1-1),] )
		data[(luper1-1),2:ncol(data)] <- rep(NA, (ncol(data)-1)) 
	}
	if ( ( data[luper1,1] != data[(luper1-1),1] ) & grpsize > 1 )     { grpsize = 1 }
	if (luper1 == nrow(data) & data[luper1,1] != data[(luper1-1),1] ) { data[(luper1),2:ncol(data)] <- rep(NA, (ncol(data)-1))  }
}
if (nrow(dumped) > 1) {
	dumped <- dumped[2:nrow(dumped),]
	dimnames(dumped) <-list(rep("", dim(dumped)[1]))
	colnames(dumped) <- colnames(data)
	cat("\n\nThe following cases were removed from the data file because there were no other cases with the same dyad code:\n\n" )
	print(dumped)
}

data <- na.omit(data)



# mean & sd 
tempdat <- data[,3:ncol(data)]
mean <- colSums(tempdat) / nrow(data)
sd <- apply(tempdat, 2, FUN = sd)

norig <- nrow(data)
tailed <- 2

# selecting pairs that have the same dyad #s, and where sex=1 & sex=2  
datarev <- matrix(-9999,1,ncol(data))
for (aa in 1:(nrow(data)-1) ) {
	if ( (data[aa,1] == data[(aa+1),1]) & data[aa,2] == 1 & data[(aa+1),2] == 2 ) datarev <- rbind( datarev, data[aa:(aa+1),] )
}

data <- datarev[2:nrow(datarev),2:ncol(data)]

n <- nrow(data) / 2

# setting up pairwise ("double-counted") collumns of data  
xp <- matrix( -9999, nrow(data), (ncol(data)-1) )
for (aa in 1:(nrow(data)/2) ) {
	xp[aa*2-1,] <- data[(aa*2),2:ncol(data)]
	xp[aa*2,]   <- data[(aa*2-1),2:ncol(data)]
}
m <- cbind( data, xp )

cr <- cor(m)

# correlations between residuals -- Velicer, 1976 
c22 <- cr[2:nrow(cr),2:ncol(cr)]
c11 <- cr[1,1]
c12 <- cr[1,2:ncol(cr)]
c22part <- c22 - (c12) %*% solve(c11) %*% c12 
d <- diag(1 / (sqrt(diag(c22part))))
resd <- d %*% c22part %*% d 

overall <- resd[1:(nrow(resd)/2), 1:(ncol(resd)/2)]
iccross <- resd[((nrow(resd)/2)+1):nrow(resd),1:(nrow(resd)/2)]

# p-Levels of the intraclass correlations  
intra   <- as.matrix(diag(iccross))
zintra  <-  (sqrt(n) * intra) / (1 - intra*intra)
pzintra <- (1 - pnorm(abs(zintra))) * tailed

# tests of assumptions involving variances & covariances   

# matrix of pairwise ("double-counted") collumns of data for just one sex 
half <- matrix( -9999, 1, (ncol(m)-1))
for (aa in 1:nrow(m) ) {
	if ( (m[aa,1] == 1) ) {
		half <- rbind( half, m[aa,(2:ncol(m))]  )
	}
}
half <- half[2:nrow(half),]

# variance-covariance matrix & correlation matrix   
vcv2 <- cov(half)
cr2 <- cor(half)

first  <- vcv2[1:(nrow(vcv2)/2), 1:(ncol(vcv2)/2)]
second <- vcv2[(nrow(vcv2)/2+1):nrow(vcv2),(ncol(vcv2)/2+1):ncol(vcv2)]
iccross2 <- cr2[((nrow(cr2)/2)+1):nrow(cr2),1:(nrow(cr2)/2)]
iccrossv <- vcv2[((nrow(vcv2)/2)+1):nrow(vcv2),1:(nrow(vcv2)/2)]

# tests for equality of variances   
vartests <- matrix( -9999, nrow(first), 6)
for (bb in 1:nrow(first) ) {
	vartests[bb,1] <- ( bb )
	vartests[bb,2:3] <- cbind( first[bb,bb], second[bb,bb]  )
	vartests[bb,4] <- ((first[bb,bb]-second[bb,bb])*sqrt(nrow(half)-2)) / 
 					(2*sqrt(first[bb,bb]*second[bb,bb]*(1-(iccross2[bb,bb]*iccross2[bb,bb]))))
	vartests[bb,5] <- nrow(half)-2
	vartests[bb,6] <- (1-pt(abs(vartests[bb,4]),vartests[bb,5]))*tailed
}

# tests for equality of the within-category population covariances 
covarswi <- matrix( -9999, 1, 6)
for (cc in 1:(ncol(first)-1) ) {
	for (rr in (cc+1):nrow(first) ) {
	line     <- matrix( -9999, 1, 6)
	line[1,1:2] <- cbind( cc, rr )
	line[1,3:4] <- cbind( first[cc,rr], second[cc,rr] )
	line[1,5] <- ((first[cc,rr]-second[cc,rr])*sqrt(nrow(half))) / sqrt(first[cc,cc]*first[rr,rr]+second[cc,cc]*second[rr,rr]+
  			2*(((first[cc,rr]+second[cc,rr])/2)  **2) - 2*(iccrossv[cc,cc]*iccrossv[rr,rr]+iccrossv[cc,rr]*iccrossv[rr,cc]))
	line[1,6] <- (1 - pnorm(abs(line[1,5])))*tailed
	covarswi <- rbind( covarswi, line )
	}
}
covarswi <- covarswi[2:nrow(covarswi),]

# tests for equality of the across-category population covariances 
covarscr <- matrix( -9999, 1, 6)
for (cc in 1:(ncol(first)-1) ) {
	for (rr in (cc+1):nrow(first)) {
		line     <- matrix( -9999, 1, 6)
		line[1,1:2] <- cbind( cc, rr )
		line[1,3:4] <- cbind( iccrossv[cc,rr], iccrossv[rr,cc] )
		line[1,5] <- ((iccrossv[cc,rr]-iccrossv[rr,cc])*sqrt(nrow(half)))/ sqrt(first[cc,cc]*second[rr,rr]+second[cc,cc]*first[rr,rr]+
			  2*(((iccrossv[cc,rr]+iccrossv[rr,cc])/2)  **2) - 2*(first[cc,rr]*second[cc,rr]+iccrossv[cc,cc]*iccrossv[rr,rr]))
		line[1,6] <- (1 - pnorm(abs(line[1,5])))*tailed
		covarscr <- rbind( covarscr, line )
	}
}
covarscr <- covarscr[2:nrow(covarscr),]


# initializing 
nstarovr   <- matrix( -9999, nrow(overall),ncol(overall) )
zoverwi    <- nstarovr
pzoverwi   <- nstarovr
zovercr    <- nstarovr
pzovercr   <- nstarovr
indcorrs   <- nstarovr
tind       <- nstarovr
ptind      <- nstarovr
dyadcors   <- nstarovr
dstar      <- nstarovr
zdyad      <- nstarovr
pzdyad     <- nstarovr

for (ii in 1:nrow(overall)) {
	for (jj in 1:ncol(overall)) {
		if ( ii != jj ) {

		# p-Levels of the overall within-partner correlations  
		nstarovr[ii,jj] <- (2*n) / (1+iccross[ii,ii]*iccross[jj,jj]+iccross[ii,jj]*iccross[ii,jj])
		zoverwi[ii,jj] <- overall[ii,jj]*sqrt(nstarovr[ii,jj])
		pzoverwi[ii,jj] <- (1-pnorm(abs(zoverwi[ii,jj])))*tailed

		# p-Levels of the overall cross-partner correlations  
		zovercr[ii,jj] <- iccross[ii,jj]*sqrt(nstarovr[ii,jj])
		pzovercr[ii,jj] <- (1-pnorm(abs(zovercr[ii,jj])))*tailed

		# individual-level correlations  
		indcorrs[ii,jj] <- (overall[ii,jj]-iccross[ii,jj]) / (sqrt(1-iccross[ii,ii])*sqrt(1-iccross[jj,jj]))
		if ( abs(indcorrs[ii,jj]) < 1 ) {
		        sqrt(1-(indcorrs[ii,jj]*indcorrs[ii,jj]))
				ptind[ii,jj] <- (1 - pt(abs(tind[ii,jj]),n-2))*tailed
		}

# dyad-level correlations  
if ( (iccross[ii,ii] > 0 & iccross[jj,jj] > 0) ) {
	dyadcors[ii,jj] <- iccross[ii,jj]/(sqrt(iccross[ii,ii])*sqrt(iccross[jj,jj]))
	dstar[ii,jj]  <- (2*n / (1+iccross[ii,ii]*iccross[jj,jj]+overall[ii,jj]*overall[ii,jj]))*iccross[ii,ii]*iccross[jj,jj]
	zdyad[ii,jj]  <- dyadcors[ii,jj]*sqrt(dstar[ii,jj])
	pzdyad[ii,jj] <- (1 - pnorm(abs(zdyad[ii,jj])))*tailed
}

}
}
}

# preparing the output matrices 
owpart  <- matrix( -9999, 1,5)
ocpart  <- owpart
dyadrzp <- owpart
withrzp <- matrix(-9999, 1,6)
for (cc in 1:(ncol(overall)-1)) {
	for (rr in (cc+1):nrow(overall)) {
		owpart <-  rbind( owpart, cbind(cc, rr, overall[cc,rr], zoverwi[cc,rr], pzoverwi[cc,rr] ) )
		ocpart <- rbind( ocpart, cbind(cc, rr, iccross[cc,rr], zovercr[cc,rr], pzovercr[cc,rr] ))
		dyadrzp <- rbind( dyadrzp, cbind(cc, rr, dyadcors[cc,rr] ,zdyad[cc,rr], pzdyad[cc,rr] ))
		withrzp <- rbind( withrzp, cbind(cc, rr, indcorrs[cc,rr], tind[cc,rr], (n-2), ptind[cc,rr] ))
	}
}
owpart  <- owpart[2:nrow(owpart),]
ocpart  <- ocpart[2:nrow(ocpart),]
dyadrzp <- dyadrzp[2:nrow(dyadrzp),]
withrzp <- withrzp[2:nrow(withrzp),]


cat("\n\nThe Correlational Analysis of Dyad-Level Data in the Distinguishable-Case")
cat("\nGonzalez & Griffin, 1999, Personal Relationships, 6, 449-469.\n")

cat("\n\nAll significance tests are two-tailed\n")
cat("\n\nNumber of individual cases in the data file:", norig)
cat("\n\nNumber of individual cases that will be processed:", (nrow(data)))

cat("\n\n\nVariable Means and Standard Deviations\n\n")
mnsd <- cbind( (mean), (sd) )
#varnames <- rownames(mnsd) # for later use
colnames(mnsd) <- cbind("    Mean","   Std Dev")
print(round(mnsd,3))

cat("\n\nVariance-Covariance Matrix for the Primary Variables:\n")
dimnames(vcv2) <-list(rep("", dim(vcv2)[1]), rep("", dim(vcv2)[2]))
print(round(vcv2,3))

cat("\n\nTests for Equality of Variances:\n\n")
vartests <- vartests[,-1]
labels <- cbind( "   Variance","   Variance","      t","     df","       p")
dimnames(vartests) <-list(rep("", dim(vartests)[1]))
colnames(vartests) <- labels
rownames(vartests) <- varnames
print(round(vartests,3))

cat("\n\nTests for Equality of the Within-Person Covariances:\n\n")
labels <- cbind( "Variable","   Variable","   Covariance","   Covariance","       z","       p")
covarswi <- data.frame(round(covarswi,3))
covarswi[,1] <- varnames[covarswi[,1]]
covarswi[,2] <- varnames[covarswi[,2]]
colnames(covarswi) <- labels 
print(covarswi, row.names = FALSE)

labels <- cbind("Variable","   Variable","       r","       z","       p")
cat("\n\nOverall Within-Partner Correlations, z-values & p-Levels:\n\n")
owpart <- data.frame(round(owpart,3))
owpart[,1] <- varnames[owpart[,1]]
owpart[,2] <- varnames[owpart[,2]]
colnames(owpart) <- labels 
print(owpart, row.names = FALSE)

labels <- cbind("Variable","   Variable","   Covariance","   Covariance","       z","       p")
cat("\n\nTests for Equality of the Across-Person Covariances:\n\n")
covarscr <- data.frame(round(covarscr,3))
covarscr[,1] <- varnames[covarscr[,1]]
covarscr[,2] <- varnames[covarscr[,2]]
colnames(covarscr) <- labels 
print (covarscr, row.names = FALSE)

labels <- cbind("Variable","   Variable","       r","       z","       p")
cat("\n\nOverall Cross-Partner Correlations, z-values & p-Levels:\n\n")
ocpart <- data.frame(round(ocpart,3))
ocpart[,1] <- varnames[ocpart[,1]]
ocpart[,2] <- varnames[ocpart[,2]]
colnames(ocpart) <- labels 
print (ocpart, row.names = FALSE)

labels <- cbind("       r","       z","       p")
cat("\n\nIntraclass Correlations, z-values & p-Levels:\n\n")
resicc <- cbind( (1:nrow(intra)),intra,zintra,pzintra )
resicc <- data.frame(round(resicc,3))
resicc <- resicc[,-1]
colnames(resicc) <- labels 
rownames(resicc) <- varnames
print (resicc)

labels <- cbind("Variable","   Variable","          r","       z","       p")
cat("\n\nDyad-Level Correlations, z-values & p-Levels:\n\n")
dyadrzp <- data.frame(round(dyadrzp,3))
dyadrzp[,1] <- varnames[dyadrzp[,1]]
dyadrzp[,2] <- varnames[dyadrzp[,2]]
colnames(dyadrzp) <- labels 
print (dyadrzp, row.names = FALSE)

labels <- cbind("Variable","   Variable","          r","        t","      df","       p")
cat("\n\nIndividual-Level Correlations, t-values & p-Levels:\n\n")
# dimnames(withrzp) <-list(rep("", dim(withrzp)[1]))
# colnames(withrzp) <- labels
withrzp <- data.frame(round(withrzp,3))
withrzp[,1] <- varnames[withrzp[,1]]
withrzp[,2] <- varnames[withrzp[,2]]
colnames(withrzp) <- labels 
print (withrzp, row.names = FALSE)

cat("\n\n-9999 indicates that a value could not be computed.")
cat("\nThe computations sometimes produce dyad-level")
cat("\ncorrelations > 1.00 or < -1.00.  This is most likely to")
cat("\noccur when an intraclass correlations is weak.")
cat("\nSee Griffin & Gonzalez, 1995, p. 434-437 for details.")


} # end of function ggdist 

