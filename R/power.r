# Project: RVCLUST
# File     power.r
#
# Author: Eli Stahl
#         Harvard Med./Partners Healthcare
#       
# Author: R. Michael Sivley
#         Center for Human Genetics Research
#         Vanderbilt University Medical Center
# Email:  mike.sivley@vanderbilt.edu
#
# Description:
#	This is an R port for Eli Stahl's SeqVarSimul-T1_VT_tests-Power.c
# power calculation function.
#
#
## ------------------------------------------------------------------------- ##

# Verify the correct distribution methods (Q(x) normal, Qinv(x) normal, etc...)

nReps       <- 10
nSNPs       <- 20
n_assocSNPs <- 10
maxPerms    <- 400000

counts          <- rep(list(rep(0,nSNPs)),nReps)
cases_counts    <- rep(list(rep(0,nSNPs)),nReps)
controls_counts <- rep(list(rep(0,nSNPs)),nReps)
nHomozygs       <- rep(list(rep(0,nSNPs)),nReps)
Ts              <- rep(list(rep(0,nSNPs)),nReps)

mafs      <- rep(list(rep(0,nSNPs)),nReps)
VTs       <- rep(list(rep(0,nSNPs)),nReps)
VT_denoms <- rep(list(rep(1,nSNPs)),nReps)

nTs       <- rep(0,nReps)
obsT1     <- rep(0,nReps)
obsVT     <- rep(0,nReps)

T1_nPerms_gt_obs <- rep(0,nReps)
T1_nPerms_eq_obs <- rep(0,nReps)
VT_nPerms_gt_obs <- rep(0,nReps)

power <- function(effectsize,ncases,nctrls,maf=NA) {

	t0 <- proc.time()[[3]]
	print(paste('Start time:',t0))
	set.seed(as.numeric(Sys.time()))

	v   <- maf
	Nca <- ncases
	Nco <- nctrls
	N   <- Nca + Nco
	K   <- effectsize

	T = qnorm(K)

	phenotypes <- c(rep(1,Nca),rep(0,Nco))
	all_possible_Ts <- rep(0,2*N)

	meanF             <- 0
	meanGRR           <- 0
	meanMAF           <- 0
	meanCaseCounts    <- 0
	meanControlCounts <- 0
	T1        <- 0
	donePerms <- 0
	perm      <- 0

	case_counts <- rep(list(rep(list(rep(0,Nca)),nSNPs)),nReps)
	control_counts <- rep(list(rep(list(rep(0,Nco)),nSNPs)),nReps)

	for (repi in 1:nReps) {

	for (i in 1:nSNPs) {
		f <- -log(runif(1,0,1)) * 0.0015

		while ( f>0.05 | f < 0.00001) {
			f <- -log(runif(1,0,1)) * 0.0015
		}
		if ( i < n_assocSNPs ) {
			l  <- sqrt(v/(2*f*(1-f)))
			u0 <- (-2 * f) * l
			u1 <- (1 - 2 * f ) * l
			u2 <- (2 - 2 * f) * l
			R0 <- pnorm((T - u0) / sqrt(1 - v))
			R1 <- pnorm((T - u1) / sqrt(1 - v))
			R2 <- pnorm((T - u2) / sqrt(1 - v))
			meanF <- meanF + f
			meanGRR <- meanGRR + R1/R0

			fca <- f * (f * R2 + (1-f) * R1) / K
			fca <- min(fca,1)
			fca <- max(fca,0)

			fco <- (f - fca*K) / (1 - K)
			fco <- min(fco,1)
			fco <- max(fco,0)
		}
		else {
			fca <- f
			fco <- f
		}

		while ((cases_counts[[repi]][i] + controls_counts[[repi]][i]) == 0) {
			for (ii in 1:Nca) {
				case_counts[[repi]][[i]][ii] <- as.numeric((runif(1,0,1) < fca))
				if (runif(1,0,1) < fca) { 
					case_counts[[repi]][[i]][ii] <- case_counts[[repi]][[i]][ii] + 1}
				if (case_counts[[repi]][[i]][ii] == 1) {
					nHomozygs[[repi]][i] <- nHomozygs[[repi]][i] + 1}
				cases_counts[[repi]][i] <- case_counts[[repi]][[i]][ii]
			}
			for (ii in 1:Nca) {
				control_counts[[repi]][[i]][ii] <- as.numeric((runif(1,0,1) < fco))
				if (runif(1,0,1) < fco) { 
					control_counts[[repi]][[i]][ii] <- control_counts[[repi]][[i]][ii] + 1}
				if (control_counts[[repi]][[i]][ii] == 1) {
					nHomozygs[[repi]][i] <- nHomozygs[[repi]][i] + 1}
				controls_counts[[repi]][i] <- control_counts[[repi]][[i]][ii]
			}
		}

		meanCaseCounts    <- meanCaseCounts + cases_counts[[repi]][i]
		meanControlCounts <- meanControlCounts + controls_counts[[repi]][i]
		counts[[repi]][i] <- cases_counts[[repi]][i] + controls_counts[[repi]][i]
		all_possible_Ts[counts[[repi]][i]] <- all_possible_Ts[counts[[repi]][i]] + 1
		mafs[[repi]][i] <- (counts[[repi]][i] + 1) / (2 * N + 2)
		meanMAF <- meanMAF + mafs[[repi]][i]
	} # Finished nSNPs loop

	meanF             <- meanF / n_assocSNPs
	meanGRR           <- meanGRR / n_assocSNPs
	meanMAF           <- meanMAF / n_assocSNPs
	meanCaseCounts    <- meanCaseCounts / n_assocSNPs
	meanControlCounts <- meanControlCounts / n_assocSNPs

	nTs[repi] <- 1
	iii       <- 1
	for (i in 1:2*N) {
		if (all_possible_Ts[i] > 0) {
			Ts[[repi]][nTs[repi]] <- i
			nTs[repi] <- nTs[repi] + 1
			iii <- iii + all_possible_Ts[i]
			if (iii==nSNPs) {break}
		}
	}
	obsT1[repi] <- 0
	for (i in 1:nSNPs) {
		if (mafs[[repi]][i] < 0.01) { 
			obsT1[repi] <- obsT1[repi] + cases_counts[[repi]][i] }
		for (iii in nTs[repi]:1) {
			if (counts[[repi]][i] > Ts[[repi]][iii]) {break}
			VTs[[repi]][iii] <- VTs[[repi]][iii] + (cases_counts[[repi]][i] - controls_counts[[repi]][i]) * 0.5
			VT_denoms[[repi]][iii] <- VT_denoms[[repi]][iii] + 4*nHomozygs[[repi]][i] + (counts[[repi]][i] - 2 * nHomozygs[[repi]][i])
		}
	}
	obsVT[repi] <- -999
	for (i in 1:nTs[repi]) {
		VTs[[repi]][i] <- VTs[[repi]][i] / sqrt(VT_denoms[[repi]][i])
		if (VTs[[repi]][i] > obsVT[repi]) {
			obsVT[repi] <- VTs[[repi]][i]
		}
	}

	} # finished nReps loop

	print(paste('Simulated',n_assocSNPs,'SNPs with per-variant liability-scale effect',v,',',nSNPs-n_assocSNPs,'null SNPs, N=',N,', rep',repi,' meanF=',meanF,'meanGRR=',meanGRR,'meanMAF=',meanMAF,'meanCaseCounts=',meanCaseCounts,'meanControlCounts=',meanControlCounts,', obsT1=',obsT1[repi],'obsVT=',obsVT[repi],sep=' '))
	print(paste('Mid time:',proc.time()[[3]]))
	print(paste('Processing time:',((proc.time()[[3]]-t0) / 60),'minutes'))


	donePerms <- 0
	while (donePerms < nReps) {
		phenotypes <- sample(phenotypes,length(phenotypes))

		for (repi in 1:nReps) {
			T1 <- 0
			VTs[[repi]] <- rep(0,nTs[repi])

			for (i in 1:nSNPs) {
				cases_counts[[repi]][i] <- 0
				for (ii in 1:Nca) {
					if (phenotypes[ii] == 1) {
						cases_counts[[repi]][i] <- cases_counts[[repi]][i] + case_counts[[repi]][[i]][i]}
				}
				#FIXME: I THINK THERE IS A LOGICAL ERROR IN THIS BLOCK,
				# BUT I HAVE COPIED IT PER ELI'S CODE.
				# WHY WOULD I BE INCREMENTING CASES_COUNTS INSTEAD
				# OF CONTROLS_COUNTS?
				for (ii in 1:Nco) {
					if (phenotypes[ii] == 1) {
						cases_counts[[repi]][i] <- cases_counts[[repi]][i] + control_counts[[repi]][[i]][ii]}
				}
				controls_counts[[repi]][i] <- counts[[repi]][i] - cases_counts[[repi]][i]
				if (mafs[[repi]][i] < 0.01) {
					T1 <- T1 + cases_counts[[repi]][i]
				}
				for (iii in nTs[repi]:1) {
					if (counts[[repi]][i] > Ts[[repi]][iii]) {break}
					VTs[[repi]][iii] <- VTs[[repi]][iii] + (cases_counts[[repi]][i] - controls_counts[[repi]][i]) * 0.5
				}
			} # finished nSNPs loop

			VT = -999
			for (i in 1:nTs[repi]) {
				VTs[[repi]][i] <- VTs[[repi]][i] / sqrt(VT_denoms[[repi]][i])
				if (VTs[[repi]][i] > VT) {
					VT <- VTs[[repi]][i]}
			}
		
		if (T1 > obsT1[repi]) {
			T1_nPerms_gt_obs[repi] <- T1_nPerms_gt_obs[repi] + 1}
		else if (T1 == obsT1[repi]) {
			T1_nPerms_gt_obs[repi] <- T1_nPerms_eq_obs[repi] + 1}
		if (VT > obsVT[repi]) {
			VT_nPerms_gt_obs[repi] <- VT_nPerms_gt_obs[repi] + 1}

		} # finished nReps loop

		perm <- perm + 1
		permtest <- 0
		if (perm < 10000) {permtest <- 1000}
		else if (perm < 100000) {permtest <- 10000}
		else {permtest <- 100000}
		if (perm %% permtest == 0) {
			donePerms <- 0
			for (repi in 1:nReps) {
				if (T1_nPerms_gt_obs[repi]+0.5*T1_nPerms_eq_obs[repi] > 1 & VT_nPerms_gt_obs[repi] > 1) {donePerms <- donePerms + 1}
			}
			print(paste('donePerms ->',donePerms))
		}
		if (perm == maxPerms) {break}
	} # finished nReps loop

	for (repi in 1:nReps) {
		print(paste('Params(nAssocSNPs,nSNPs,v,N;rep;obsT1,obsVT)= ',n_assocSNPs,' ',nSNPs,' ',v,' ',N,' ',repi,' ',obsT1[repi],' ',obsVT[repi],'=>Pvals(T1,Vt;perms=',perm,')= ',(T1_nPerms_gt_obs[repi]+runif(1,0,1)*T1_nPerms_eq_obs[repi])/perm,' ',VT_nPerms_gt_obs[repi]/perm,sep=''))}
	print('Finished!')
	print(paste('Stop time:',proc.time()[[3]]))
	print(paste('Processing time:',((proc.time()[[3]]-t0) / 60),'minutes'))
	return(0)
}