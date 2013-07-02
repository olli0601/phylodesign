#' this file contains the likelihood functions for each model

acute.MAX.TIPC.SIZE<<- 15

popart.CLUSTERP.ACHG<<- matrix(c(0.2,0.05,0.05,0.4/3,0.1/3,0.1/3,0.2/3,0.05/3,0.05/3),3,3,dimnames=list(c("U","T","O"),c("1","2","3")))

popart.CLUSTERP.ACLW<<- matrix(c(825/1800,495/1800,165/1800,50/1800,30/1800,10/1800,25/1800,15/1800,5/1800),3,3,dimnames=list(c("U","T","O"),c("1","2","3")))

###############################################################################
acute.get.rates<- function(ibm.beta, ibm.pop= NULL, ibm.initpop= NULL, pop.n=nrow(ibm.pop), state.n= as.matrix(table(subset(ibm.pop,select=status))), per.capita.i= 0, debug=0)
{
	#in model 'Acute', rates are ' base * rel. infectiousness * S / N '
	if(!setequal( names(ibm.beta[['i']][[1]]),rownames(state.n) ))	stop("expected same covariates in init.pop and beta")
	state.n				<- state.n[names(ibm.beta[['i']][[1]]),]		#re-order covariates
	infecteds			<- state.n
	if(per.capita.i)	
		infecteds[]		<- 1
	#ibm.beta[['i']][[1]] is relative transmission rates per covariate, ie s i t u   --> compute beta per covariate
	propens				<- ibm.beta[['i']][[1]]*infecteds*ibm.beta[["base"]]	
	#ibm.beta[['s']][[1]] is relative susceptibility per covariate; ONLY s is susceptible and we ASSUME no re-infection
	if(!debug)
		propens			<- propens %*% t(ibm.beta[['s']][[1]] * state.n / pop.n)
	else
	{		
		state.n			<- as.matrix(table(subset(ibm.initpop,select=status)))[names(ibm.beta[['i']][[1]]),]
		propens			<- propens %*% t(ibm.beta[['s']][[1]] * state.n / pop.n)
	}		 
	rownames(propens)	<- colnames(propens)
	propens
}	
###############################################################################
acute.subtree.lkl.PartNT<- function(ntrs,rIdx,rE,cl.timewindow, log=0)
{
	if(log)
		ans		<- -rIdx*cl.timewindow  + ntrs*log(1-exp( -rE*cl.timewindow ))
	else
		ans		<- exp( -rIdx*cl.timewindow ) * (1-exp( -rE*cl.timewindow ))^ntrs	
	ans			<- matrix(ans,length(ntrs),length(ntrs),byrow=1,dimnames=list(paste('i',ntrs,sep=''),paste('n',ntrs,sep='')))
	ans
}
###############################################################################
acutesampled.subtree.lkl.PartNT<- function(ntrs,sample.prob,log=0)
{
	if(log)
		ans	<- ntrs * log(sample.prob)
	else
		ans	<- sample.prob^ntrs	
	ans		<- matrix(ans, nrow=length(ntrs), ncol=length(ntrs), byrow=1, dimnames= list(paste('i',ntrs,sep=''), paste('n',ntrs,sep='')))
	ans
}
###############################################################################
acutesampled.subtree.lkl.PartMiss<- function(ntrs,sample.prob,rIdx,rE,cl.timewindow,miss.max, log=0)
{
	part.powermissing	<- ( (1-sample.prob) * (1-exp(-rE*cl.timewindow)) )^seq.int(0,miss.max)
	#print(part.powermissing)
	ans					<- sapply(seq.int(0,length(ntrs)-1), function(n)
							{
								tmp<- sapply(seq.int(0,n), function(i)
										{
											#for debug, return i
											part.summissing<- sapply(seq.int(0,miss.max), function(m)
																{
																	sum(	sapply(seq.int(0,m), function(j)	dhyper(i, i+j, n+m-(i+j), n, log = FALSE) * (rIdx/rE)^j		)	)																												
																})
											#print(part.summissing); print(part.summissing * part.powermissing)
											tmp2	<- sum(part.summissing * part.powermissing)
											ifelse(log, log(tmp2), tmp2)											
										})
								c(tmp, rep(ifelse(log,-Inf,0), length(ntrs)-1 - n))
							})
	colnames(ans)			<- paste('n',ntrs,sep='')
	rownames(ans)			<- paste('i',ntrs,sep='')
	ans
}
###############################################################################
acutesampled.subtree.lkl<- function(lkl.complete, sample.prob, rIdx, rE, cl.timewindow, log=0)
{
	if(sum(lkl.complete)>1+EPS)	stop("acutesampled.subtree.lkl: complete lkl sums to > 1")
	
	closure				<- ncol(lkl.complete)-1
	ntrs				<- seq.int(0,closure)
	part.powermissing	<- (1-exp(-rE*cl.timewindow))^seq.int(0,closure)
	ans					<- sapply(ntrs, function(n)
			{
				tmp<- sapply(seq.int(0,n), function(i)
						{
							#for debug, return i							
							part.summissing<- sapply(seq.int(0,closure-n), function(m)
									{										
										sum(	sapply(seq.int(0,m), function(j)	dhyper(i, i+j, n+m-(i+j), n, log = FALSE) * (rIdx/rE)^j	* lkl.complete[i+j+1,n+m+1]	)	)																												
									})
							#print(part.summissing); print(part.summissing * part.powermissing[seq_along(part.summissing)]); print(dbinom(n, n+seq.int(0,closure-n), sample.prob))
							tmp2	<- sum( part.summissing * part.powermissing[seq_along(part.summissing)] * dbinom(n, n+seq.int(0,closure-n), sample.prob) )
							#tmp2	<- sum(part.summissing * part.powermissing)
							ifelse(log, log(tmp2), tmp2)							
						})
				c(tmp, rep(ifelse(log,-Inf,0), closure - n))
			})
	colnames(ans)			<- paste('n',ntrs,sep='')
	rownames(ans)			<- paste('i',ntrs,sep='')
	if(!log && sum(ans)>1+EPS)		stop("acutesampled.subtree.lkl: sampled lkl sums to > 1")
	if(log && sum(exp(ans))>1+EPS)	stop("acutesampled.subtree.lkl: sampled lkl sums to > 1")
	ans	
}
###############################################################################
acute.subtree.lkl.PartIdx<- function(ntr,rIdx,rE, log=0)
{
	if(log)
		ans			<- round(log(upper.tri(matrix(1, ntr, ntr),diag=T))) + log(rIdx/rE)*seq.int(1,ntr)
	else
		ans			<- round(upper.tri(matrix(1, ntr, ntr),diag=T)) * (rIdx/rE)^seq.int(1,ntr)
	ans				<- cbind(rep(ifelse(log,-Inf,0),nrow(ans)),ans)
	if(log)
		ans			<- rbind( c(0,rep(-Inf,ncol(ans)-1)), ans)
	else	
		ans			<- rbind( c(1,rep(0,ncol(ans)-1)), ans)	
	dimnames(ans)	<- list(paste('idx',seq.int(0,nrow(ans)-1),sep=''),paste('n',seq.int(0,ncol(ans)-1),sep=''))
	ans
}
###############################################################################
#' Compute the likelihood of a transmission chain with \code{nx} transmissions from the index case and \code{ni} transmissions from non-index cases under the \code{Acute} model
#' There are 1+nx+ni NODES in this tree.
#' @param nx	number of transmissions from donor with risk group X, nx=0 is possible
#' @param ni	number of transmissions from donor with risk group I, ni=0 is possible
#' @export 
acute.lkl.tree.xk.ik<- function(nx,ni,rx,ri,dT, log=0)
{	
	if(length(rx)!=1 || length(ri)!=1 || length(nx)!=1 || length(ni)!=1)	stop("acute.lkl.tree.xk.ik: not vectorized")
	k<- seq.int(0,nx+ni)
	if(log)
		ans<- nx * log(rx/ri) -rx*dT - lfactorial(nx+ni) + log(sum(choose(nx+ni,k) * (-1)^k * exp(-k*ri*dT)))
	else		
		ans<- ((rx/ri)^nx * exp(-rx*dT) / factorial(nx+ni)) * sum(choose(nx+ni,k) * (-1)^k * exp(-k*ri*dT))
	ans
}
###############################################################################
#' Compute the log likelihood of a tip cluster table under the \code{Acute} model
acute.loglkl<- function(tpc, rate.m, dT, clu.n=NULL)
{
	#fix mle for initial frequencies	
	init.freq.mle		<- apply(tpc,1,sum) / sum(tpc[c('u','t'),])
#print(init.freq.mle)	
	#compute partial mle for rate.m
	tpc.n.mx			<- ncol(tpc)-1														#max number of transmissions in tip cluster
	tpc.ncol			<- tpc.n.mx		#without sampling, need to compute probabilities only up to the largest number of transmissions + 1 in the any tip cluster
#print(tpc.ncol)	
	if(tpc.n.mx>tpc.ncol)	
		stop(paste("\nfound tip cluster sizes",tpc.n.mx,"while max supported is",tpc.ncol))
	if(is.null(clu.n))	
	{
		clu.n			<- clu.tipc.n(tpc.ncol)
		clu.n.sum		<- apply(clu.n,2,sum)	#cheaper than 'seq_len(ncol(clu.n))^seq.int(-1,ncol(clu.n)-2)'
		clu.n			<- clu.n / matrix(clu.n.sum, nrow=nrow(clu.n), ncol=length(clu.n.sum), byrow=1)										
	}
	else
		clu.n			<- clu.n[seq_len(tpc.ncol+1),seq_len(tpc.ncol+1)]
#print(clu.n); print(tpc.ncol)
	#transmission chain likelihoods for those that start with 'E' - always -Inf	
	tipc.lkl.E			<-	rep(-Inf,tpc.ncol)
	#get transmission chain likelihoods for those that start with 'U'
	part1				<- acute.subtree.lkl.PartNT(seq.int(0,ncol(clu.n)-1),rate.m['u','s'],rate.m['i','s'],dT, log=1)
	part1				<- matrix(part1,nrow(clu.n),length(part1),byrow=1,dimnames=list(paste('idx',seq.int(0,nrow(clu.n)-1),sep=''),paste('n',seq.int(0,length(part1)-1),sep='')))
	part2				<- acute.subtree.lkl.PartIdx(ncol(clu.n)-1,rate.m['u','s'],rate.m['i','s'], log=1)
	chain.lkl			<- part1 + part2 + log(clu.n)		
	dimnames(chain.lkl)	<- dimnames(part2)
	#check that the whole lot is <1 and adjust if necessary -- still trying to find where this error comes from!
	checklsum			<- log1p( sum( exp(chain.lkl) )-1 )
	if(checklsum>0)
		chain.lkl		<- chain.lkl - checklsum
	#integrate transmission chains out
	tipc.lkl.U			<- apply(chain.lkl,2,function(x) log(sum(exp(x[x!=-Inf]), na.rm=1)))
	#
	#get transmission chain likelihoods for those that start with 'T'
	#
	part1				<- acute.subtree.lkl.PartNT(seq.int(0,ncol(clu.n)-1),rate.m['t','s'],rate.m['i','s'],dT, log=1)
	part1				<- matrix(part1,nrow(clu.n),length(part1),byrow=1,dimnames=list(paste('idx',seq.int(0,nrow(clu.n)-1),sep=''),paste('n',seq.int(0,length(part1)-1),sep='')))
	part2				<- acute.subtree.lkl.PartIdx(ncol(clu.n)-1,rate.m['t','s'],rate.m['i','s'], log=1)
	chain.lkl			<- part1 + part2 + log(clu.n)		
	dimnames(chain.lkl)	<- dimnames(part2)
	#check that the whole lot is <1 and adjust if necessary -- still trying to find where this error comes from!
	checklsum			<- log1p( sum( exp(chain.lkl) )-1 )
	if(checklsum>0)
		chain.lkl		<- chain.lkl - checklsum
	
	#integrate transmission chains out
	tipc.lkl.T			<- apply(chain.lkl,2,function(x) log(sum(exp(x[x!=-Inf]), na.rm=1)))
#print(tipc.lkl.T)		
	tipc.lkl.U			<- tipc.lkl.U + log(init.freq.mle)['u']
	tipc.lkl.T			<- tipc.lkl.T + log(init.freq.mle)['t']

	tipc.lkl			<- rbind(tipc.lkl.E,rbind(tipc.lkl.U,tipc.lkl.T))	
	rownames(tipc.lkl)	<- c('i','u','t')	
	tpc					<- tpc[rownames(tipc.lkl),]	
#print(tipc.lkl); print(tpc)
	if(ncol(tipc.lkl)!=ncol(tpc))	stop("columns of tipc.lkl and tpc do not match")
	if(nrow(tipc.lkl)!=nrow(tpc))	stop("rows of tipc.lkk and tpc do not match")
#print(tipc.lkl); print(tpc); stop()
	table.lkl			<- tpc * tipc.lkl
#print(table.lkl)	
	table.lkl			<- sum( table.lkl[!is.nan(table.lkl) & tipc.lkl!=-Inf] )	
	ans					<- list(table.lkl= table.lkl, tipc.lkl=tipc.lkl)
	ans
}
###############################################################################
#' Compute the log likelihood of a tip cluster table under the \code{Acute} model
acutesampled.loglkl<- function(tpc, rate.m, sample.prob, dT, lclu.n=NULL)
{
	#fix mle for initial frequencies	
	init.freq.mle		<- apply(tpc,1,sum) / sum(tpc[c('u','t'),])					#TODO still the MLE for the initial values ? 
#print(init.freq.mle)	
	#compute partial mle for rate.m
	tpc.n.mx	<- ncol(tpc)-1														#max number of transmissions in tip cluster	
	tpc.ncol	<- max(acute.MAX.TIPC.SIZE,tpc.n.mx+2)									#without sampling, need to compute probabilities only up to the largest number of transmissions + 1 in the any tip cluster
#print(tpc.ncol)
	if(	is.null(lclu.n)  || 
		(!is.null(lclu.n) && ncol(lclu.n)<tpc.ncol)	)
	{
		clu.n			<- clu.tipc.n(tpc.ncol)	
#print(clu.n)		
		clu.n.sum		<- apply(clu.n,2,sum)																	#faster than seq_len(ncol(clu.n))^seq.int(-1,ncol(clu.n)-2)
		lclu.n			<- log( clu.n / matrix(clu.n.sum, nrow=nrow(clu.n), ncol=length(clu.n.sum), byrow=1) ) 
		#lclu.n.sum		<- log(seq_len(ncol(clu.n)))*seq.int(-1,ncol(clu.n)-2)									#this is less accurate
		#lclu.n			<- log(clu.n) - matrix(lclu.n.sum, nrow=nrow(clu.n), ncol=length(lclu.n.sum), byrow=1)	#this is less accurate													
	}
	else
		lclu.n			<- lclu.n[seq_len(tpc.ncol),seq_len(tpc.ncol)]
#print(exp(lclu.n)); print(tpc.ncol); stop()
	#get transmission chain likelihoods for those that start with 'U'
	part1				<- acute.subtree.lkl.PartNT(seq.int(0,ncol(lclu.n)-1),rate.m['u','s'],rate.m['i','s'], dT, log=1)
	part2				<- acute.subtree.lkl.PartIdx(ncol(lclu.n)-1,rate.m['u','s'],rate.m['i','s'], log=1)	
	chain.lkl			<- part1 + part2 + lclu.n		
	dimnames(chain.lkl)	<- dimnames(part2)
	checklsum			<- log1p( sum( exp(chain.lkl) )-1 )
	if(checklsum>0)
		chain.lkl		<- chain.lkl - checklsum	
	chain.lkl			<- acutesampled.subtree.lkl(exp(chain.lkl), sample.prob, rate.m['u','s'],rate.m['i','s'], dT, log=1)
	#the [idx0, seq.int(2,ncol(chain.lkl))] are tip clusters starting with an 'E' index case; account for initial frequency of the U -- currently -Inf because of clu.n
	tipc.lkl.E			<- exp( chain.lkl[1,seq.int(0,tpc.n.mx)+2] ) * init.freq.mle['u']
	chain.lkl[1,seq.int(2,ncol(chain.lkl))]	<- -Inf	
	chain.lkl			<- chain.lkl[,seq.int(0,tpc.n.mx)+1]	
	#integrate transmission chains out
	tipc.lkl.U			<- apply(chain.lkl,2,function(x) log(sum(exp(x))))
	#
	#get transmission chain likelihoods for those that start with 'T'
	#
	part1				<- acute.subtree.lkl.PartNT(seq.int(0,ncol(lclu.n)-1),rate.m['t','s'],rate.m['i','s'],dT, log=1)	
	part2				<- acute.subtree.lkl.PartIdx(ncol(lclu.n)-1,rate.m['t','s'],rate.m['i','s'], log=1)
	chain.lkl			<- part1 + part2 + lclu.n		
	dimnames(chain.lkl)	<- dimnames(part2)
	checklsum			<- log1p( sum( exp(chain.lkl) )-1 )
	if(checklsum>0)
		chain.lkl		<- chain.lkl - checklsum	
	chain.lkl			<- acutesampled.subtree.lkl(exp(chain.lkl), sample.prob, rate.m['t','s'],rate.m['i','s'], dT, log=1)	
	#the [idx0, seq.int(2,ncol(chain.lkl))] are tip clusters starting with an 'E' index case
	tipc.lkl.E			<- log(      tipc.lkl.E    +    exp(chain.lkl[1,seq.int(0,tpc.n.mx)+2]) * init.freq.mle['t']      )
	names(tipc.lkl.E)	<- paste("ns",seq.int(0,length(tipc.lkl.E)-1),sep='')
	chain.lkl[1,seq.int(2,ncol(chain.lkl))]	<- -Inf	
	chain.lkl			<- chain.lkl[,seq.int(0,tpc.n.mx)+1]
#print(chain.lkl)	
	#integrate transmission chains out
	tipc.lkl.T			<- apply(chain.lkl,2,function(x) log(sum(exp(x))))
#print(tipc.lkl.T)		
	tipc.lkl.U			<- tipc.lkl.U + log(init.freq.mle)['u']
	tipc.lkl.T			<- tipc.lkl.T + log(init.freq.mle)['t']
	
	tipc.lkl			<- rbind(tipc.lkl.E,rbind(tipc.lkl.U,tipc.lkl.T))	
	rownames(tipc.lkl)	<- c('i','u','t')	
	tpc					<- tpc[rownames(tipc.lkl),]	
#print(tipc.lkl); print(tpc); #stop()
	if(ncol(tipc.lkl)!=ncol(tpc))	stop("columns of tipc.lkl and tpc do not match")
	if(nrow(tipc.lkl)!=nrow(tpc))	stop("rows of tipc.lkk and tpc do not match")
#print("tipc.lkl"); print(tipc.lkl); print("tipc table"); print(tpc); #stop()
	table.lkl			<- tpc * tipc.lkl
#print(table.lkl)
	ans					<- list(table.lkl= sum( table.lkl[!is.nan(table.lkl) & is.finite(tipc.lkl)] ), tipc.lkl=tipc.lkl)
	ans
}
###############################################################################
#' randomize allocation of arms to triplets 
#' @param randomize.n number of times the randomization is repeated
#' @return data frame of trial sites with allocated arms
popart.getdata.random.arm<- function( randomize.n )
{
	sites	<- read.csv("triplet_summary_table_121212.csv")		
	# select studies for phylogenetics
	#in ZA chose small high prevalence triplets outside of Lusaka
	#5 is Khayelitsha triplet in SA 
	#6 is High prevalence triplet in SA			
	simul	<- sites[rep(which(sites$triplet %in% c(1,4,5,6)),randomize.n),]
	# randomise communities arms
	simul$random	<- runif(nrow(simul))
	simul$arm.code	<- unlist(tapply(simul$random,rep(1:(nrow(simul)/3),each=3),rank))
	simul$arm		<- ifelse(simul$arm.code==1,"A",simul$arm.code)
	simul$arm		<- ifelse(simul$arm.code==2,"B",simul$arm)
	simul$arm		<- ifelse(simul$arm.code==3,"C",simul$arm)
	drops 			<- c("random","arm.code")
	simul			<- simul[,!names(simul) %in% drops]
	simul$inc.rate	<- sapply(simul$arm,function(x){(x %in% "A")*0.0056+(x %in% "B")*0.01+(x %in% "C")*0.013})
	simul$p.adults	<- ifelse(simul$country==1,0.5,0.6)
	simul$rid		<- rep(1:randomize.n, each=3*4)
	simul
}
###############################################################################
popart.getdata.random.n<- function(x, pc.size, d.study)
{
	n			<- nrow(x)
	x$n.adults	<- rbinom(n,x$popsize,x$p.adults)
	x$n.prev	<- rbinom(n,x$n.adults,x$hivcomb/100)
	x$n.on.art	<- rbinom(n,x$n.prev,x$artadjust/100)
	x$n.not.art	<- x$n.prev-x$n.on.art	
	x$n.inc		<- rbinom(n,x$n.adults-x$n.prev,d.study*x$inc.rate)
	x$n.neg		<- with(x,n.adults-n.prev-n.inc)	
	x$p.in.PC	<- pc.size/x$n.adults	
	x$PC.on.art	<- rbinom(n,x$n.on.art,x$p.in.PC)
	x$PC.not.art<- rbinom(n,x$n.not.art,x$p.in.PC)
	x$PC.inc	<- rbinom(n,x$n.inc,x$p.in.PC)
	x$PC.neg	<- pc.size-x$PC.on.art-x$PC.not.art-x$PC.inc
	x
}
###############################################################################
popart.confint.plot2<- function(p0,p1,c0,c1,f.name=NA,xlab='s.consent',ylab='conf.int',cols=c("red","blue"),legend.loc="topright",legend.txt=c('',''))
{
	#print(p0); print(p1); print(c0); print(c1)
	def.par <- par(no.readonly = TRUE)
	if(!is.na(f.name))
		pdf(paste(f.name),version="1.4",width=6,height=6)
	par(mar=c(5,5.5,0.5,2))
	x<- as.numeric(rownames(c0))
	plot(1,1,type='n',xlim=range(x),ylim=range(rbind(c0,c1)),xlab=xlab,ylab=ylab)
	
	lines(x,p0,col=cols[1], lty=4)		
	lines(x,p1,col=cols[2], lty=4)		
	polygon( c(x,rev(x)),c(c0[,1],rev(c0[,2])), border=NA, col=my.fade.col(cols[1],0.3) )
	polygon( c(x,rev(x)),c(c1[,1],rev(c1[,2])), border=NA, col=my.fade.col(cols[2],0.3) )
	legend(legend.loc,fill=c(cols,"transparent","transparent"),legend=c(legend.txt,"","mean and 95% CI"),bty='n',border=NA)	
	if(!is.na(f.name))
		dev.off()
}
###############################################################################
popart.confint.plot<- function(p0,p1,c0,c1,f.name=NA,xlab='s.consent',ylab='conf.int',cols=c("red","blue"),legend.loc="topright",legend.txt=c('',''))
{
	#print(p0); print(p1); print(c0); print(c1)
	def.par <- par(no.readonly = TRUE)
	if(!is.na(f.name))
		pdf(paste(f.name),version="1.4",width=12,height=6)
	par(mar=c(5,5.5,0.5,2))
	layout( matrix(c(1,2),1,2) )
	
	x<- as.numeric(colnames(p0))	 
	plot(1,1,type='n',xlim=range(x),ylim=range(rbind(c0,c1)),xlab=xlab,ylab=ylab)
	lines(x,p0["cohort",],col=cols[1], lty=4)		
	lines(x,p1["cohort",],col=cols[2], lty=4)		
	polygon( c(x,rev(x)),c(c0["cohort.l",],rev(c0["cohort.u",])), border=NA, col=my.fade.col(cols[1],0.3) )
	polygon( c(x,rev(x)),c(c1["cohort.l",],rev(c1["cohort.u",])), border=NA, col=my.fade.col(cols[2],0.3) )
	legend(legend.loc,fill=c("transparent",cols),legend=c("cohort only",legend.txt),bty='n',border=NA)
	
	plot(1,1,type='n',xlim=range(x),ylim=range(rbind(c0,c1)),xlab=xlab,ylab='')
	lines(x,p0["pop",],col=cols[1],lty=2)
	lines(x,p1["pop",],col=cols[2],lty=2)
	polygon( c(x,rev(x)),c(c0["pop.l",],rev(c0["pop.u",])), border=NA, col=my.fade.col(cols[1],0.3) )	
	polygon( c(x,rev(x)),c(c1["pop.l",],rev(c1["pop.u",])), border=NA, col=my.fade.col(cols[2],0.3) )	
	legend(legend.loc,fill=c("transparent",cols),legend=c("cohort and HCC",legend.txt),bty='n',border=NA)	
	par(def.par)
	if(!is.na(f.name))
		dev.off()
}
###############################################################################
popart.powercalc.plot2<- function(m, m2, c, c2, f.name, xlab, ylab, legend.txt, legend.loc, cols=c("red","blue"), verbose= 0)
{	
	x<- as.numeric(names(m))
	x2<- as.numeric(names(m2))
	legend.txt<- c(legend.txt,"lines w/o dots cannot be trusted")
	pdf(paste(f.name),version="1.4",width=6,height=6)
	par(mar=c(4.5,5,1,1))
	plot(1,1,type='n',xlim=range(c(x,x2)),ylim=c(0,1),xlab=xlab,ylab=ylab)
	lines(x,m,col=cols[1], lty=4)	
	points(x[c],m[c],col=cols[1],pch=19)
	
	lines(x2,m2,col=cols[2],lty=2)
	points(x2[c2],m2[c2],col=cols[2],pch=22)
	
	if(!is.character(legend.loc))
		legend(x=legend.loc$x,y=legend.loc$y,fill=c(cols[1],cols[2],"transparent","transparent"),legend=legend.txt,bty='n', border= NA)
	else
		legend(legend.loc,fill=c(cols[1],cols[2],"transparent","transparent"),legend=legend.txt,bty='n', border= NA)	
	dev.off()
}
###############################################################################
popart.powercalc.plot<- function(m, m2, f.name, xlab, ylab, legend.txt, legend.loc, cols=c("red","blue"), no.conf.m=matrix(TRUE,ncol=ncol(m),nrow=nrow(m),dimnames=list(rownames(m),colnames(m))), no.conf.m2=matrix(TRUE,ncol=ncol(m2),nrow=nrow(m2),dimnames=list(rownames(m2),colnames(m2))), verbose= 0)
{	
	x<- as.numeric(colnames(m))
	legend.txt<- c(legend.txt,"lines w/o dots cannot be trusted")
	pdf(paste(f.name),version="1.4",width=8,height=8)
	par(mar=c(4.5,5,1,1))
	plot(1,1,type='n',xlim=range(x),ylim=c(0,1),xlab=xlab,ylab=ylab)
	lines(x,m["cohort",],col=cols[1], lty=4)	
	points(x[which(!no.conf.m["cohort",])],m["cohort",which(!no.conf.m["cohort",])],col=cols[1],pch=19)
	lines(x,m2["cohort",],col=cols[1],lty=2)
	points(x[which(!no.conf.m2["cohort",])],m2["cohort",which(!no.conf.m2["cohort",])],col=cols[1],pch=22)
	lines(x,m["pop",],col=cols[2], lty=4)	
	
	points(x[which(!no.conf.m["pop",])],m["pop",which(!no.conf.m["pop",])],col=cols[2],pch=19)
	lines(x,m2["pop",],col=cols[2], lty=2)
	points(x[which(!no.conf.m2["pop",])],m2["pop",which(!no.conf.m2["pop",])],col=cols[2],pch=22)	
	if(!is.character(legend.loc))
		legend(x=legend.loc$x,y=legend.loc$y,fill=c(cols[1],cols[2],"transparent","transparent"),legend=legend.txt,bty='n', border= NA)
	else
		legend(legend.loc,fill=c(cols[1],cols[2],"transparent","transparent"),legend=legend.txt,bty='n', border= NA)
		
	dev.off()
}
###############################################################################
popart.get.tipcprops<- function(tipc.p,f.contam)
{
	tipc.p<- tipc.p / sum(tipc.p)								#tipc.p may not be standardized
	if(rownames(tipc.p)[nrow(tipc.p)]!='O')	stop("popart.get.tipcprops: error at 1a")
	
	tmp<- apply(tipc.p,1,sum)
	tmp<- c(tmp,sum(tmp[-length(tmp)]))
	ans<- sapply(seq_len(nrow(tipc.p)-1),function(i)
				{
					tmp[i]/tmp[length(tmp)] * (1-f.contam) / sum(tipc.p[i,]) * tipc.p[i,]
				})
	ans<- rbind(t(ans),	f.contam / tmp[length(tmp)-1] * tipc.p[nrow(tipc.p),] )
	rownames(ans)<- rownames(tipc.p)	
	#print(apply(ans,1,sum))
	ans
}
###############################################################################
popart.getsampled.tipcseq<- function(tipc.p,pop.size,f.coh,f.susc,f.untreated,cohort.dur,perc.inc,p.lab,p.nocontam,p.consent.coh,p.consent.clu,p.prev.instudy.clu,p.inc.instudy.clu, opt.power="All")
{	
	if(!opt.power%in%c("All","PonlyPC","IonlyPC","PonlyPCandIonlyPC"))
		stop("popart.getsampled.tipcseq: unknown power method")
	verbose<- 0
	if(length(pop.size)!=length(f.susc))		stop("popart.getsampled.tipcseq: error at 1a")
	if(length(pop.size)!=length(f.untreated))	stop("popart.getsampled.tipcseq: error at 1b")
	if(length(pop.size)!=length(f.coh))			stop("popart.getsampled.tipcseq: error at 1c")
	
	inc.e<- perc.inc*cohort.dur*f.susc*pop.size					#total expected incidence in community		
	tipc.p<- tipc.p / sum(tipc.p)								#tipc.p may not be standardized
	#print(tipc.p); print(apply(tipc.p,1,sum))
	inc.p<- apply(tipc.p,2,sum)*seq_len(ncol(tipc.p))			#fraction of incidence in tipc column		
	tipc.e<- inc.e / sum(inc.p) 								#total expected tip clusters in community with at least 1 infection	
	if(verbose){	print(sum(inc.p)); print(inc.e); print(tipc.e)	}
	
	#get tip cluster sampling matrix for cohort and community, and transpose so we re ready for matrix multiplication
	#since there are higher order terms of si.XX, need to do cohort and community separately
	#cannot extract sequence from treated individuals
	sp.clu<-	(p.consent.coh*f.coh	+  p.consent.clu*(1-f.coh)*p.prev.instudy.clu)	*	p.lab*f.untreated/(1-f.susc)	
	si.clu<-	p.lab*p.consent.clu*p.inc.instudy.clu
	si.coh<-	p.lab*p.consent.coh		
	if(length(si.clu)!=1)	stop("popart.getsampled.tipcseq: cannot vectorize si.clu")
	if(length(si.coh)!=1)	stop("popart.getsampled.tipcseq: cannot vectorize si.coh")
	
	s.clu<- t(matrix(	c(	si.clu,		2*si.clu*(1-si.clu),		3*si.clu*(1-si.clu)^2,	
							0,			si.clu*si.clu,				3*si.clu*si.clu*(1-si.clu),
							0,			0,							si.clu*si.clu*si.clu 	),
						nrow=ncol(tipc.p),ncol=ncol(tipc.p),byrow=1))
	s.coh<- t(matrix(	c(	si.coh,		2*si.coh*(1-si.coh),		3*si.coh*(1-si.coh)^2,	
							0,			si.coh*si.coh,				3*si.coh*si.coh*(1-si.coh),
							0,			0,							si.coh*si.coh*si.coh 	),
						nrow=ncol(tipc.p),ncol=ncol(tipc.p),byrow=1))
	if(verbose){ print(sp.clu*s.clu); print(sp.clu*s.coh) }
	
	tipc.p<- tipc.p[-3,]				#do not observe transitions from outside community	
	ans<- sapply(seq_along(tipc.e),function(i)
			{
				out<- numeric(4)
				names(out)<- c("i2i.s","i2i.c","x2i.s","x2i.c")
				tipc<- tipc.e[i]*tipc.p				#all tip cluster counts before sampling			
				if(verbose) print(tipc)
				tmp<- c(apply(tipc,1,sum), sum(tipc[,2])+2*sum(tipc[,3]))		#all U->I, T->I, I->I before sampling
				if(verbose)	print(tmp)
				out[c("i2i.c","x2i.c")]<- c(tmp[3],sum(tmp))					#I->I and U->I + T->I + I->I	before sampling
				
				tipc<- tipc.e[i]*sp.clu[i]*( f.coh[i]*tipc.p%*%s.coh	+ (1-f.coh[i])*tipc.p%*%s.clu)	#tip clusters sampled in PC + tip clusters sampled in HCC
				if(verbose) print(tipc)	
				tmp<- c(apply(tipc,1,sum), sum(tipc[,2])+2*sum(tipc[,3]))		#all U->I, T->I, I->I after sampling
				if(verbose) print(tmp)	
				out[c("i2i.s","x2i.s")]<- c(tmp[3],sum(tmp))					#I->I and U->I + T->I + I->I	after sampling
				out				
			})
	if(verbose) print(ans)
	ans
}
###############################################################################
#' compute the fraction of linked incident sequences
#' @param x				data frame
#' @param opt.sampling	sampling scheme
#' @param v_prev_AB 	proportion of initially +ve individuals who visit HCC at some point in 3 year period in Arms A & B
#' @param v_prev_C  	proportion of initially +ve individuals who visit HCC at some point in 3 year period in Arm C
#' @param v_inc_AB		proportion of incident cases who visit HCC at some point in 3 year period in Arms A and B
#' @param v_inc_C		proportion of incident cases who visit HCC at some point in 3 year period in Arm C
#' @param consent.PC	consent rate to phylogenetics in population cohort
#' @param consent.HCC	consent rate to phylogenetics in health care centres
#' @param p.lab			proportion of isolates that make it through to full genomes
#' @param p.community	proportion of infections that come from within the community
#' @return linked.incident
popart.getsampled.linkedseq.n<- function(x, opt.sampling, v_prev_AB=0.95, v_prev_C=0.5, v_inc_AB=0.8, v_inc_C=0.25, consent.PC=0.9, consent.HCC=0.5, p.lab=0.8, p.community=0.85)
{
	visit.prev	<- numeric(nrow(x))
	visit.inc	<- numeric(nrow(x))
	visit.prev[ x$arm%in% c("A","B") ]	<- v_prev_AB
	visit.prev[ x$arm=="C" ]			<- v_prev_C
	visit.inc[ x$arm%in% c("A","B") ]	<- v_inc_AB
	visit.inc[ x$arm=="C" ]				<- v_inc_C
			
	if(opt.sampling=="PC and HCC")
	{
		seq.PC.prev		<- floor(x$PC.not.art*consent.PC*p.lab)
		seq.PC.inc		<- floor(x$PC.inc*consent.PC*p.lab)
	}
	else if(opt.sampling=="only HCC")
	{
		seq.PC.prev		<- floor(x$PC.not.art*visit.prev*consent.HCC*p.lab)
		seq.PC.inc		<- floor(x$PC.inc*visit.inc*consent.HCC*p.lab)
	}
	else
		stop("popart.getsampled.linkedseq.n: unknown method to sample")
	
	seq.nonPC.prev	<- floor((x$n.not.art-x$PC.not.art)*visit.prev*consent.HCC*p.lab)		
	seq.nonPC.inc	<- floor((x$n.inc-x$PC.inc)*visit.inc*consent.HCC*p.lab)
	
	seq.prev.coverage	<- (seq.PC.prev + seq.nonPC.prev)/x$n.prev	
	linked.incident		<- floor((seq.PC.inc + seq.nonPC.inc)*seq.prev.coverage*p.community)
	linked.incident
}
###############################################################################
#' compute the fraction of linked incident sequences based on overall proportions
#' this is replaced by "popart.getsampled.linkedseq.n"
popart.getsampled.linkedseq<- function(pop.size,f.coh,f.susc,f.untreated,cohort.dur,perc.inc,p.lab,p.nocontam,p.consent.coh,p.consent.clu,p.prev.instudy.clu,p.inc.instudy.clu, opt.power="All")
{	
	if(!opt.power%in%c("All","PonlyPC","IonlyPC","PonlyPCandIonlyPC"))
		stop("popart.getsampled.linkedseq: unknown power method")
	tmp<-			p.lab*p.consent.coh*f.coh	+	p.lab*p.consent.clu*(1-f.coh)*p.prev.instudy.clu 	
	seq.prev<- 		tmp*f.untreated*pop.size		#cannot extract sequence from treated
	tmp<- 			p.lab*p.consent.coh*f.coh	+	p.lab*p.consent.clu*(1-f.coh)*p.inc.instudy.clu
	seq.inc.clu<-	tmp*perc.inc*cohort.dur*f.susc*pop.size
	tmp<- 			seq.prev/((1-f.susc)*pop.size)	#fraction of prevalent cases for which sequences can be extracted
	linked<-		p.nocontam*tmp*seq.inc.clu		#treated and untreated equally likely to transmit
	linked
}	
###############################################################################
#'	pool the number of acute to acute transmissions according to a pooling option
#' 	@param opt.pooled poopling option
#' 	@param sites	data frame of trial sites
#' 	@param inc		matrix of acute to acute transmissions
#' 	@return list with the components 'inc': pooled matrix of acute to acute transmissions, 'idx.A': row index of arm A sites, 'idx.C': row index of arm C sites
popart.pool<- function(opt.pooled, sites, inc)
{
	if(opt.pooled=="no pooling")
	{
		idx.A<-	which(sites$arm=="A")
		idx.C<-	which(sites$arm=="C")				
	}
	else if(opt.pooled=="pooled across country")
	{
		inc<- do.call(rbind,tapply(seq_len(nrow(sites)),sites$rid, function(ris)
						{
							rbind(	apply(inc[ris[which(sites$arm[ris]=="A" & sites$country[ris]==1)],],2,sum),
									apply(inc[ris[which(sites$arm[ris]=="B" & sites$country[ris]==1)],],2,sum),
									apply(inc[ris[which(sites$arm[ris]=="C" & sites$country[ris]==1)],],2,sum),
									apply(inc[ris[which(sites$arm[ris]=="A" & sites$country[ris]==2)],],2,sum),
									apply(inc[ris[which(sites$arm[ris]=="B" & sites$country[ris]==2)],],2,sum),
									apply(inc[ris[which(sites$arm[ris]=="C" & sites$country[ris]==2)],],2,sum))					
						}))
		idx.A<-	seq(1,nrow(inc),by=3)
		idx.C<-	seq(3,nrow(inc),by=3)
	}
	else if(opt.pooled=="pooled across ZA")
	{
		inc<- do.call(rbind,tapply(seq_len(nrow(sites)),sites$rid, function(ris)
						{
							rbind(	apply(inc[ris[which(sites$arm[ris]=="A" & sites$country[ris]==1)],],2,sum),
									apply(inc[ris[which(sites$arm[ris]=="B" & sites$country[ris]==1)],],2,sum),
									apply(inc[ris[which(sites$arm[ris]=="C" & sites$country[ris]==1)],],2,sum))					
						}))
		idx.A<-	seq(1,nrow(inc),by=3)
		idx.C<-	seq(3,nrow(inc),by=3)
	}
	else if(opt.pooled=="pooled across SA")
	{
		inc<- do.call(rbind,tapply(seq_len(nrow(sites)),sites$rid, function(ris)
						{
							rbind(	apply(inc[ris[which(sites$arm[ris]=="A" & sites$country[ris]==2)],],2,sum),
									apply(inc[ris[which(sites$arm[ris]=="B" & sites$country[ris]==2)],],2,sum),
									apply(inc[ris[which(sites$arm[ris]=="C" & sites$country[ris]==2)],],2,sum))					
						}))
		idx.A<-	seq(1,nrow(inc),by=3)
		idx.C<-	seq(3,nrow(inc),by=3)
	}
	else if(opt.pooled=="pooled across trial")
	{
		inc<- do.call(rbind,tapply(seq_len(nrow(sites)),sites$rid, function(ris)
						{
							rbind(	apply(inc[ris[which(sites$arm[ris]=="A")],],2,sum),
									apply(inc[ris[which(sites$arm[ris]=="B")],],2,sum),
									apply(inc[ris[which(sites$arm[ris]=="C")],],2,sum)	)					
						}))
		idx.A<-	seq(1,nrow(inc),by=3)
		idx.C<-	seq(3,nrow(inc),by=3)				
	}
	else stop("popart.pool: pooling option not recognized")
	
	ans<- list(inc=inc, idx.A=idx.A, idx.C=idx.C)
	ans
}	
###############################################################################
popart.getpower.exactbinomtest<- function(x,test.prop0, test.prop1, test.alpha, method="homebrew")
{
	if(!method%in%c("homebrew","cloglog", "logit", "probit", "asymp", "lrt", "exact"))
		stop("popart.getpower.exactbinomtest: cannot find method")
	if(length(test.prop0)!=length(x))	test.prop0<- rep(test.prop0,length(x))
	if(length(test.prop1)!=length(x))	test.prop1<- rep(test.prop1,length(x))
	x<- round(x)
	if(method=="homebrew")
		ans<- sapply(seq_along(x),function(i)		#power based on exact binomial distribution
				{
					if(test.prop0[i]<test.prop1[i])
					{	
						tmp<- seq( round(test.prop0[i]*x[i]), x[i] )
						cu<- tmp[ which( pbinom( tmp, x[i], prob= test.prop0[i], lower.tail=0 )<test.alpha ) ][1]		#lower.tail=0 means P(X>x) so cu+1 is the lowest value at which T is rejected
	#print(c(cu+1,n))			
						tmp<- pbinom(cu,x[i],prob= test.prop1[i], lower.tail=0)	#evaluates P(X>cu)=P(X>=cu+1)=P(rej)
					}
					else
					{
						tmp<- seq( 0, round(x[i]) )
						cl<- tmp[ which(pbinom( tmp, x[i], prob= test.prop0[i] )>=test.alpha) ][1]		#cl is the lowest value at which T is still accepted
						tmp<- pbinom(cl-1,x[i],prob= test.prop1[i])		#by default, pbinom(-1,..) is 0						
					}
					tmp
				})
	else
	{
		require(binom)		
		ans<- sapply(seq_along(x),function(i)		#various approximations
			{
				alternative<- ifelse(test.prop0[i]<test.prop1[i],"greater","less")					
				binom.power(test.prop1[i], x[i], test.prop0[i], test.alpha, phi=1, alternative=alternative, method)
			})
	}
	ans
}	
###############################################################################
popart.getpower.exactbinomtest.randomarm<- function(m.type, loc.type, cohort.size, cohort.dur, s.prev.coh.armA, s.inc.coh.armA, s.prev.clu.armA, s.inc.clu.armA, s.prev.coh.armC, s.inc.coh.armC, s.prev.clu.armC, s.inc.clu.armC, test.prop0, test.prop1, test.alpha, verbose= 0)
{		
	if(loc.type=="SA")
		tmp<- sample(13:21,9)	#first 3 arm 'A', next 3 arm 'B', then arm 'C'
	else
		tmp<- sample(seq_len(nrow(POPART.SITES)),nrow(POPART.SITES))	#first 7 arm 'A', next 7 arm 'B', then arm 'C'
	loc.arm<- rep('C',nrow(POPART.SITES))
	loc.arm[ tmp[seq_len(nrow(POPART.SITES)/3)] ]<- 'A'
	loc.arm[ tmp[seq(nrow(POPART.SITES)/3+1,nrow(POPART.SITES)/3*2)] ]<- 'B'
	#compute incidence in cohort, enriched cohort, community that can be phylogenetically linked to index case
	#for each community and a random allocation of arm A,B,C
	good.inc<- sapply(seq_len(nrow(POPART.SITES)),function(i)
			{
				loc.type<- rownames(POPART.SITES)[i]
				pop.size<- POPART.SITES[loc.type,"pop"]*0.6 								
				pop.distr<- ibm.init.popinitdistributions(ibm.init.attributes("Acute"), loc.type)[["status"]]
				perc.inc<- ibm.init.incidence(paste(loc.type,loc.arm[i],sep='-'))*pop.distr['s']	#percent incidence in susceptible population per year					
				susc.clu<- pop.size * pop.distr['s']
				susc.coh<-  cohort.size * pop.distr['s']				
				if(loc.arm[i]!='C')
					out<- popart.getsampled.incidence(cohort.dur, susc.clu, susc.coh, perc.inc, s.prev.coh.armA, s.inc.coh.armA, s.prev.clu.armA, s.inc.clu.armA)
				else
					out<- popart.getsampled.incidence(cohort.dur, susc.clu, susc.coh, perc.inc, s.prev.coh.armC, s.inc.coh.armC, s.prev.clu.armC, s.inc.clu.armC)			
				out
			})
	colnames(good.inc)<- paste(rownames(POPART.SITES),loc.arm,sep='-') 
	#compute total number of detected links observed in cohort, enriched cohort, community, grouped by arm	
	good.inc.byarm<- sapply(c('A','B','C'),function(arm)
			{
				apply(good.inc[,which(loc.arm==arm)],1,sum)
			})
	good.inc.pw.byarm<- popart.getpower.exactbinomtest(good.inc.byarm,test.prop0, test.prop1, test.alpha)
	good.inc.pw.byarm<- matrix(good.inc.pw.byarm, nrow=nrow(good.inc.byarm), ncol=ncol(good.inc.byarm),dimnames=dimnames(good.inc.byarm))
	list(inc=good.inc.byarm,power=good.inc.pw.byarm)		
}
###############################################################################
popart.powercalc.pop<- function(pop.size, perc.inc, cohort.size, pc24.size, cohort.dur, sample.initcohort, sample.linktocluster, sample.pop, test.prop0, test.prop1, test.alpha, verbose= 0)
{
	ans<- numeric(3)
	names(ans)<- c("cohort","enriched","pop")
	
	ans[1]<- round( cohort.dur * cohort.size * perc.inc * sample.linktocluster * sample.initcohort) 
	
	sample.pc24<- cohort.size / pop.size * sample.initcohort
	ans[2]<- ans[1] + round( pc24.size * perc.inc / 2 * sample.linktocluster * sample.pc24 )
	
	tmp<- (sample.initcohort*cohort.size + sample.pop*(pop.size-cohort.size))/pop.size	#weighted average of finding an index case
	
	ans[3]<- round( cohort.dur * cohort.size * perc.inc * sample.linktocluster * tmp 
					+ 	pc24.size * perc.inc / 2 * sample.linktocluster * tmp
					+	(cohort.dur * (pop.size - cohort.size) - pc24.size/2 ) * perc.inc * tmp * sample.linktocluster ) 
	ans
}
###############################################################################
popart.powercalc.randomarm<- function(m.type, cohort.size, pc24.size, cohort.dur, sample.initcohort, sample.linktocluster, sample.pop, test.prop0, test.prop1, test.alpha, verbose= 0)
{		
	tmp<- sample(seq_len(nrow(POPART.SITES)),nrow(POPART.SITES))	#first 7 arm 'A', next 7 arm 'B', then arm 'C'
	loc.arm<- rep('C',nrow(POPART.SITES))
	loc.arm[ tmp[seq_len(nrow(POPART.SITES)/3)] ]<- 'A'
	loc.arm[ tmp[seq(nrow(POPART.SITES)/3+1,nrow(POPART.SITES)/3*2)] ]<- 'B'
	#compute incidence in cohort, enriched cohort, community that can be phylogenetically linked to index case
	#for each community and a random allocation of arm A,B,C
	good.inc<- sapply(seq_len(nrow(POPART.SITES)),function(i)
			{
				loc.type<- rownames(POPART.SITES)[i]
				pop.size<- POPART.SITES[loc.type,"pop"]
				perc.inc<- ibm.init.incidence(paste(loc.type,loc.arm[i],sep='-'))	#percent incidence in population (N) per year				
				popart.powercalc.pop(pop.size, perc.inc, cohort.size, pc24.size, cohort.dur, sample.initcohort, sample.linktocluster, sample.pop, test.prop0, test.prop1, test.alpha, verbose= 0)
			})
	colnames(good.inc)<- paste(rownames(POPART.SITES),loc.arm,sep='-') 
	#compute total number of detected links observed in cohort, enriched cohort, community, grouped by arm	
	good.inc.byarm<- sapply(c('A','B','C'),function(arm)
			{
				apply(good.inc[,which(loc.arm==arm)],1,sum)
			})
	if(verbose)
	{
		print(good.inc)
		print(good.inc.byarm)
	}		
	good.inc.pw.byarm<- sapply(as.vector(good.inc.byarm),function(n)
			{
				tmp<- seq( round(test.prop0*n), n )
				cu<- tmp[ which( pbinom( tmp, n, prob= test.prop0, lower.tail=0 )<test.alpha ) ][1]		#lower.tail=0 means P(X>x) so cu is the upper limit at which T is not rejected
				pbinom(cu+1,n,prob= test.prop1, lower.tail=0)
			})
	good.inc.pw.byarm<- matrix(good.inc.pw.byarm, nrow=nrow(good.inc.byarm), ncol=ncol(good.inc.byarm),dimnames=dimnames(good.inc.byarm))
	good.inc.pw.byarm	
}
###############################################################################
#' performs power calculations and computes binomial confidence intervals for the "linkedseq" approach
#' @param inc 			computed linked sequences of incident cases
#' @param idx.A			row index of arm A sites
#' @param idx.C			row index of arm C sites
#' @param test.prop0	null proportion of transmission during acute phase
#' @param test.prop1 	alternate proportion of transmission during acute phase
#' @param test.alpha	level of test
#' @return list containing power and confidence intervals
popart.powercalc.linkedseq.n<- function(inc, idx.A, idx.C, test.prop0, test.prop1, test.alpha, verbose=0)
{
	x<- inc[idx.A,]	
	x<- apply(x,2,mean)	
	if(verbose){	cat("\nphylog linked incident cases\n");	print(x)	}	
	conf.hg.med.armA<- binom.confint(as.vector(round(x*test.prop1)), as.vector(x), conf.level = 0.95, methods="bayes")[,c("lower","upper")]
	conf.lw.med.armA<- binom.confint(as.vector(round(x*test.prop0)), as.vector(x), conf.level = 0.95, methods="bayes")[,c("lower","upper")]
	rownames(conf.lw.med.armA)<- rownames(conf.hg.med.armA)<- colnames(inc)
	#print(conf.lw.med.armA); print(conf.hg.med.armA)
	is.conf.hg.med.armA<- test.prop0<conf.hg.med.armA[,"lower"] & x*test.prop1>=5
	power.hg.med.armA<- popart.getpower.exactbinomtest(x, test.prop1, test.prop0, test.alpha)
	power.lw.med.armA<- popart.getpower.exactbinomtest(x, test.prop0, test.prop1, test.alpha)
	names(power.hg.med.armA)<- names(power.lw.med.armA)<- colnames(inc)
	if(verbose){ cat(paste("\npower to distinguish",test.prop0, "from", test.prop1,"\n"));	print(power.lw.med.armA); print(power.hg.med.armA)	}
	
	x<- inc[idx.C,]
	x<- apply(x,2,mean)
	conf.hg.med.armC<- binom.confint(as.vector(round(x*test.prop1)), as.vector(x), conf.level = 0.95, methods="bayes")[,c("lower","upper")]
	conf.lw.med.armC<- binom.confint(as.vector(round(x*test.prop0)), as.vector(x), conf.level = 0.95, methods="bayes")[,c("lower","upper")]
	rownames(conf.hg.med.armC)<- rownames(conf.lw.med.armC)<- colnames(inc)	
	is.conf.hg.med.armC<- test.prop0<conf.hg.med.armC[,"lower"] & x*test.prop1>=5		
	power.hg.med.armC<- popart.getpower.exactbinomtest(x, test.prop1, test.prop0, test.alpha)
	power.lw.med.armC<- popart.getpower.exactbinomtest(x, test.prop0, test.prop1, test.alpha)
	names(power.hg.med.armC)<- names(power.lw.med.armC)<- colnames(inc)
	if(verbose){ cat(paste("\npower to distinguish",test.prop0, "from", test.prop1,"\n"));	print(power.lw.med.armC)	}
	
	ans<- list(power.lw.med.armA, power.lw.med.armC, is.conf.hg.med.armA, is.conf.hg.med.armC, conf.lw.med.armA, conf.hg.med.armA,conf.lw.med.armC, conf.hg.med.armC)
	ans		
}
