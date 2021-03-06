


ggexch <- function (data) {


# Correlational Analyses of Exchangeable-Case Dyadic Data
# Griffin & Gonzalez, 1995, Psychological Bulletin, 118, 430-439. 

# Create a data matrix where rows = cases, & columns = variables
# the first collumn of data must contain the dyad numbers
# (the dyad id# given to both dyad members); the second and subsequent
# collumns contain the variables to be analyzed.

# Cases with missing values are not permitted in the data file

varnames <- colnames(data)[-1]

if ( all(complete.cases(data)) == 'FALSE' ) {
	cat("\n\nERROR: There are missing values in the data matrix. Fix this problem and try again") }

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
tempdat <- data[,2:ncol(data)]
mean <- colSums(tempdat) / nrow(data)
sd <- apply(tempdat, 2, FUN = sd)

norig <- nrow(data)
tailed <- 2

# selecting pairs that have the same dyad #s 
datarev <- matrix(-9999,1,ncol(data))
for (aa in 1:(nrow(data)-1)) {
	if ( (data[aa,1] == data[(aa+1),1]) )  datarev <- rbind( datarev, data[aa:(aa+1) ,] )
}
data <- datarev[2:nrow(datarev),2:ncol(data)]

n <- nrow(data) / 2

# setting up pairwise ("double-counted") collumns of data 
xp <- matrix( -9999, nrow(data),ncol(data))
for (aa in 1:(nrow(data)/2)) {
	xp[aa*2-1,] <- data[aa*2  ,]
	xp[aa*2  ,] <- data[aa*2-1,]
}
m <- cbind( data, xp )

cr <- cor(m)

overall <- cr[1:(nrow(cr)/2), 1:(ncol(cr)/2)]
iccross <- cr[((nrow(cr)/2)+1):nrow(cr), 1:(nrow(cr)/2)]

# significance levels of the intraclass correlations 
intra   <- as.matrix(diag(iccross))
zintra  <- intra * sqrt(n)
pzintra <- (1 - pnorm(abs(zintra))) * tailed

# initializing 
nstarovr   <- matrix( -9999, nrow(overall),ncol(overall))
zover      <- nstarovr
pzover     <- nstarovr
nstarcrs   <- nstarovr
zcrs       <- nstarovr
pzcrs      <- nstarovr
indcorrs   <- nstarovr
tind       <- nstarovr
ptind      <- nstarovr
dyadcors   <- nstarovr
dstar      <- nstarovr
zdyad      <- nstarovr
pzdyad     <- nstarovr
totcors    <- nstarovr

for (ii in 1:nrow(overall)) {
	for (jj in 1:ncol(overall)) {
		if ( ii != jj ) {

		# significance levels of the overall correlations 
		nstarovr[ii,jj] <- (2 * n) / (1+iccross[ii,ii]*iccross[jj,jj]+iccross[ii,jj]*iccross[ii,jj])
		zover[ii,jj] <- overall[ii,jj]*sqrt(nstarovr[ii,jj])
		pzover[ii,jj] <- (1 - pnorm(abs(zover[ii,jj])))*tailed

		# significance levels of the cross-intraclass correlations 
		nstarcrs[ii,jj] <- 2*n / (1+iccross[ii,ii]*iccross[jj,jj]+overall[ii,jj]*overall[ii,jj])
		zcrs[ii,jj] <- iccross[ii,jj]*sqrt(nstarcrs[ii,jj])
		pzcrs[ii,jj] <- (1 - pnorm(abs(zcrs[ii,jj])))*tailed

		# individual-level correlations 
		indcorrs[ii,jj] <- (overall[ii,jj]-iccross[ii,jj]) / (sqrt(1-iccross[ii,ii])*sqrt(1-iccross[jj,jj]))
			if ( (abs(indcorrs[ii,jj]) < 1) ) {
				tind[ii,jj] <- (indcorrs[ii,jj]*sqrt(n-1)) / sqrt(1 - (indcorrs[ii,jj]*indcorrs[ii,jj]))
				ptind[ii,jj] <- (1 - pt(abs(tind[ii,jj]),n-1))*tailed
			}

		# dyad-level correlations 
		if ( (iccross[ii,ii] > 0 & iccross[jj,jj] > 0) ) {
			dyadcors[ii,jj] <- iccross[ii,jj] / (sqrt(iccross[ii,ii])*sqrt(iccross[jj,jj]))
			dstar[ii,jj]  <- (2*n / (1+iccross[ii,ii]*iccross[jj,jj]+overall[ii,jj]*overall[ii,jj]))*iccross[ii,ii]*iccross[jj,jj]
			zdyad[ii,jj]  <- dyadcors[ii,jj]*sqrt(dstar[ii,jj])
			pzdyad[ii,jj] <- (1 - pnorm(abs(zdyad[ii,jj])))*tailed
		}

		# total or mean-level correlations 
		totcors[ii,jj] <- -9999
		if ( (iccross[ii,ii] != -1 & iccross[jj,jj] != -1) ) {
			totcors[ii,jj] <- (overall[ii,jj]*iccross[ii,jj]) / (sqrt(1+iccross[ii,ii])*sqrt(1+iccross[jj,jj]))
		}

	}
}
}

# preparing the output matrices 
overzp   <- matrix( -9999, 1,5)
iccroszp <- overzp
indrzp   <- matrix( -9999, 1,6)
dyadrzp  <- overzp
totr     <- matrix( -9999, 1,3)

for (cc in 1:(ncol(overall)-1)) {
for (rr in (cc+1):nrow(overall)) {
		overzp <-  rbind( overzp, cbind(cc, rr,overall[cc,rr],zover[cc,rr],pzover[cc,rr] ))
		iccroszp <-  rbind( iccroszp, cbind(cc,rr,iccross[cc,rr],zcrs[cc,rr],pzcrs[cc,rr] ))
		indrzp <-  rbind( indrzp, cbind(cc,rr,indcorrs[cc,rr],tind[cc,rr],(n-1),ptind[cc,rr] ))
		dyadrzp <-  rbind( dyadrzp, cbind(cc,rr,dyadcors[cc,rr],zdyad[cc,rr],pzdyad[cc,rr] ))
		totr <-  rbind( totr, cbind(cc,rr,totcors[cc,rr] ))
	}
}
overzp   <- overzp[2:nrow(overzp),]
iccroszp <- iccroszp[2:nrow(iccroszp),]
indrzp   <- indrzp[2:nrow(indrzp),]
dyadrzp  <- dyadrzp[2:nrow(dyadrzp),]
totr     <- totr[2:nrow(totr),]

cat("\n\n\nCorrelational Analyses of Exchangeable-Case Dyadic Data")
cat("\nGriffin & Gonzalez, 1995, Psychological Bulletin, 118, 430-439.")

cat("\n\nAll significance tests are two-tailed\n")
cat("\n\nNumber of individual cases in the data file:", norig)
cat("\n\nNumber of individual cases that will be processed:", (nrow(data)))

cat("\n\n\nVariable Means and Standard Deviations\n\n")
mnsd <- cbind( (mean), (sd) )
#varnames <- rownames(mnsd) # for later use
colnames(mnsd) <- cbind("    Mean","   Std Dev")
print(round(mnsd,3))

cat("\n\nIntraclass Correlations, z-values & p-Levels:\n\n")
labels <- cbind("       r","       z","       p")
resicc <- cbind( (1:nrow(intra)),intra,zintra,pzintra )
resicc <- data.frame(round(resicc,3))
resicc <- resicc[,-1]
colnames(resicc) <- labels 
rownames(resicc) <- varnames
print (resicc)

cat("\n\nOverall Correlations, z-values & p-Levels:\n\n")
labels <- cbind("Variable","   Variable","       r","       z","       p")
overzp <- data.frame(round(overzp,3))
overzp[,1] <- varnames[overzp[,1]]
overzp[,2] <- varnames[overzp[,2]]
colnames(overzp) <- labels 
rownames(overzp) <- varnames
print(overzp, row.names = FALSE)

cat("\n\nCross-Intraclass Correlations, z-values & p-Levels:\n\n")
labels <- cbind("Variable","   Variable","       r","       z","       p")
iccroszp <- data.frame(round(iccroszp,3))
iccroszp[,1] <- varnames[iccroszp[,1]]
iccroszp[,2] <- varnames[iccroszp[,2]]
colnames(iccroszp) <- labels 
rownames(iccroszp) <- varnames
print(iccroszp, row.names = FALSE)

cat("\n\nIndividual-Level Correlations, t-values & p-Levels:\n\n")
labels <- cbind("Variable","   Variable","       r","       t","      df","       p")
indrzp <- data.frame(round(indrzp,3))
indrzp[,1] <- varnames[indrzp[,1]]
indrzp[,2] <- varnames[indrzp[,2]]
colnames(indrzp) <- labels 
rownames(indrzp) <- varnames
print(indrzp, row.names = FALSE)

cat("\n\nDyad-Level Correlations, z-values & p-Levels:\n\n")
labels <- cbind("Variable","   Variable","       r","       z","       p")
dyadrzp <- data.frame(round(dyadrzp,3))
dyadrzp[,1] <- varnames[dyadrzp[,1]]
dyadrzp[,2] <- varnames[dyadrzp[,2]]
colnames(dyadrzp) <- labels 
rownames(dyadrzp) <- varnames
print(dyadrzp, row.names = FALSE)

cat("\n\nTotal or Mean-Level Correlations:\n\n")
labels <- cbind("Variable","   Variable","       r")
totr <- data.frame(round(totr,3))
totr[,1] <- varnames[totr[,1]]
totr[,2] <- varnames[totr[,2]]
colnames(totr) <- labels 
rownames(totr) <- varnames
print(totr, row.names = FALSE)

cat("\n\n-9999 indicates that a value could not be computed.")
cat("\nThe computations sometimes produce dyad-level")
cat("\ncorrelations > 1.00 or < -1.00.  This is most likely to")
cat("\noccur when an intraclass correlations is weak.")
cat("\nSee Griffin & Gonzalez, 1995, p. 434-437 for details.")

} # end of function ggexch


