###############################################################################
prj.test<- function()
{	
	m.type<- "Acute"	
	my.mkdir(DATA,"test")
	
	if(1)		#test acute likelihood
	{
		loc.type<- "ZA-C"
		m.popsize<- 3e4
		theta<- c(8, 0.09)
		
		dir.name<- paste(DATA,"simudata",sep='/')
		f.name<- paste(dir.name,paste("tpc_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta[1],"_b",theta[2],sep=''),sep='/')
		print(f.name)
		if(1)
		{
			load(paste(f.name,"_tpc.R",sep=''))
			tpc.t<- tpc.tabulate(tpc[[1]][["tpc"]])
			ibm<- tpc[[1]][["ibm"]]
			cat(paste("\nsave",paste(f.name,"_example_contingencytable_upto3.R",sep='') ))
			save(tpc.t,ibm,file=paste(f.name,"_example_contingencytable_upto3.R",sep=''))
			print(tpc.t)
			stop()
		}
		load(paste(f.name,"_example_contingencytable_upto3.R",sep=''))
		print(tpc.t)
		#construct rates
		ibm<- ibm.as.data.table(ibm)
		ibm[["beta"]][['i']][["status"]]['i']<- theta[1]
		ibm[["beta"]][["base"]]<- theta[2]		
		
		
		cat(paste("\ntheta.acute",ibm[["beta"]][['i']][["status"]]['i'],"theta.base",ibm[["beta"]][["base"]],"\n"))
		lkl.rates<- acute.get.rates(ibm, per.capita.i=1)
		lkl<- acute.loglkl(tpc.t, lkl.rates, ibm[["beta"]][["dT"]])
		print(lkl.rates)
		print(lkl)
		
		ibm[["beta"]][['i']][["status"]]['i']<- 1
		ibm[["beta"]][["base"]]<- theta[2]		
		cat(paste("\ntheta.acute",ibm[["beta"]][['i']][["status"]]['i'],"theta.base",ibm[["beta"]][["base"]],"\n"))
		lkl.rates<- acute.get.rates(ibm, per.capita.i=1)
		lkl<- acute.loglkl(tpc.t, lkl.rates, ibm[["beta"]][["dT"]])
		print(lkl.rates)
		print(lkl)
		
		ibm[["beta"]][['i']][["status"]]['i']<- 4
		ibm[["beta"]][["base"]]<- theta[2]		
		cat(paste("\ntheta.acute",ibm[["beta"]][['i']][["status"]]['i'],"theta.base",ibm[["beta"]][["base"]],"\n"))
		lkl.rates<- acute.get.rates(ibm, per.capita.i=1)
		lkl<- acute.loglkl(tpc.t, lkl.rates, ibm[["beta"]][["dT"]])
		print(lkl.rates)
		print(lkl)
		
		ibm[["beta"]][['i']][["status"]]['i']<- 8
		ibm[["beta"]][["base"]]<- theta[2]		
		cat(paste("\ntheta.acute",ibm[["beta"]][['i']][["status"]]['i'],"theta.base",ibm[["beta"]][["base"]],"\n"))
		lkl.rates<- acute.get.rates(ibm, per.capita.i=1)
		lkl<- acute.loglkl(tpc.t, lkl.rates, ibm[["beta"]][["dT"]])
		print(lkl.rates)
		print(lkl)
		stop()
		
		ibm[["beta"]][['i']][["status"]]['i']<- theta[1]
		ibm[["beta"]][["base"]]<- 0.04
		cat(paste("\ntheta.acute",ibm[["beta"]][['i']][["status"]]['i'],"theta.base",ibm[["beta"]][["base"]],"\n"))
		lkl.rates<- acute.get.rates(ibm, per.capita.i=1)
		lkl<- acute.loglkl(tpc.t, lkl.rates, ibm[["beta"]][["dT"]])
		print(lkl.rates)
		print(lkl)
		print(sum(lkl))
		
		ibm[["beta"]][['i']][["status"]]['i']<- theta[1]
		ibm[["beta"]][["base"]]<- 0.11		
		cat(paste("\ntheta.acute",ibm[["beta"]][['i']][["status"]]['i'],"theta.base",ibm[["beta"]][["base"]],"\n"))
		lkl.rates<- acute.get.rates(ibm, per.capita.i=1)
		lkl<- acute.loglkl(tpc.t, lkl.rates, ibm[["beta"]][["dT"]])
		print(lkl.rates)
		print(lkl)
		print(sum(lkl))
		#print(ibm[["beta"]])
		#print(lkl.rates)
		#print(lkl)
		#compute likelihood
	}
	if(0)		#test acute likelihood
	{
		loc.type<- "ZA-A"
		m.popsize<- 6e4
		theta<- c(2, 0.098)
		
		dir.name<- paste(DATA,"simudata",sep='/')
		f.name<- paste(dir.name,paste("tpc_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta[1],"_b",theta[2],sep=''),sep='/')
		#print(f.name)
		#load(paste(f.name,"_tpc.R",sep=''))
		#tpc.t<- tpc.tabulate(tpc[[1]][["tpc"]])
		#ibm<- tpc[[1]][["ibm"]]
		#cat(paste("\nsave",paste(f.name,"_example_contingencytable_upto3.R",sep='') ))
		#save(tpc.t,ibm,file=paste(f.name,"_example_contingencytable_upto3.R",sep=''))
		#stop()
		load(paste(f.name,"_example_contingencytable_upto3.R",sep=''))
		print(tpc.t)
		#construct rates
		ibm<- ibm.as.data.table(ibm)
		ibm[["beta"]][['i']][["status"]]['i']<- theta[1]
		ibm[["beta"]][["base"]]<- theta[2]		
		
		
		cat(paste("\ntheta.acute",ibm[["beta"]][['i']][["status"]]['i'],"theta.base",ibm[["beta"]][["base"]],"\n"))
		lkl.rates<- acute.get.rates(ibm, per.capita.i=1)
		lkl<- acute.loglkl(tpc.t, lkl.rates, ibm[["beta"]][["dT"]])
		print(lkl.rates)
		print(lkl)
		print(sum(lkl))
		
		ibm[["beta"]][['i']][["status"]]['i']<- 1.5
		ibm[["beta"]][["base"]]<- theta[2]		
		cat(paste("\ntheta.acute",ibm[["beta"]][['i']][["status"]]['i'],"theta.base",ibm[["beta"]][["base"]],"\n"))
		lkl.rates<- acute.get.rates(ibm, per.capita.i=1)
		lkl<- acute.loglkl(tpc.t, lkl.rates, ibm[["beta"]][["dT"]])
		print(lkl.rates)
		print(lkl)
		print(sum(lkl))
		
		ibm[["beta"]][['i']][["status"]]['i']<- 6
		ibm[["beta"]][["base"]]<- theta[2]		
		cat(paste("\ntheta.acute",ibm[["beta"]][['i']][["status"]]['i'],"theta.base",ibm[["beta"]][["base"]],"\n"))
		lkl.rates<- acute.get.rates(ibm, per.capita.i=1)
		lkl<- acute.loglkl(tpc.t, lkl.rates, ibm[["beta"]][["dT"]])
		print(lkl.rates)
		print(lkl)
		print(sum(lkl))
		
		ibm[["beta"]][['i']][["status"]]['i']<- theta[1]
		ibm[["beta"]][["base"]]<- 0.04
		cat(paste("\ntheta.acute",ibm[["beta"]][['i']][["status"]]['i'],"theta.base",ibm[["beta"]][["base"]],"\n"))
		lkl.rates<- acute.get.rates(ibm, per.capita.i=1)
		lkl<- acute.loglkl(tpc.t, lkl.rates, ibm[["beta"]][["dT"]])
		print(lkl.rates)
		print(lkl)
		print(sum(lkl))
		
		ibm[["beta"]][['i']][["status"]]['i']<- theta[1]
		ibm[["beta"]][["base"]]<- 0.11		
		cat(paste("\ntheta.acute",ibm[["beta"]][['i']][["status"]]['i'],"theta.base",ibm[["beta"]][["base"]],"\n"))
		lkl.rates<- acute.get.rates(ibm, per.capita.i=1)
		lkl<- acute.loglkl(tpc.t, lkl.rates, ibm[["beta"]][["dT"]])
		print(lkl.rates)
		print(lkl)
		print(sum(lkl))
		#print(ibm[["beta"]])
		#print(lkl.rates)
		#print(lkl)
		#compute likelihood
	}
	if(0)		#test ibm attack rate
	{	
		m.type<- "Acute"
		loc.type<- "ZA-C"
		m.popsize<- 1e4
		resume<- 0
		verbose<- 1
		record.tpc<- 0
		theta<- c(16, 0.073)
		dir.name<- paste(DATA,"test",sep='/')
		f.name<- paste(dir.name,paste("tpc_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta[1],"_b",theta[2],sep=''),sep='/')
		print(f.name)
		
		
		ibm<- ibm.init.model(m.type,loc.type,m.popsize,theta,save=paste(f.name,"_ibm.R",sep=''), resume= resume)
		ssa.ibm<- ibm.collapse(ibm)
		n<- replicate(100, ssa.run(ssa.ibm, ssa.ctime= 0, ssa.etime= 1,record.tpc= record.tpc, verbose= 0, save=paste(f.name,"_tpc.R",sep=''), resume= resume))
		print("\n")
		print(summary(n))
	}
	if(1)		#test ibm sampling
	{	
		m.type<- "Acute"
		loc.type<- "ZA-C"
		m.popsize<- 6e4
		resume<- 1
		verbose<- 1
		record.tpc<- 1
		theta<- c(2, 0.112, 0.5, 0.5)
		names(theta)<- c("acute","base","m.st1","m.st2")
		
		dir.name<- paste(DATA,"test",sep='/')
		f.name<- paste(dir.name,paste("tpc_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta[1],"_b",theta[2],sep=''),sep='/')
		print(f.name)
		
		
		ibm<- ibm.init.model(m.type,loc.type,m.popsize,theta,save=paste(f.name,"_ibm.R",sep=''), resume= resume)
		ssa.ibm<- ibm.collapse(ibm)
		tpc.data<- ssa.run(ssa.ibm, ssa.ctime= 0, ssa.etime= 1,record.tpc= record.tpc, verbose= 1, save=paste(dir.name,"firstsampledtpc.R",sep='/'), resume= resume)		
		#load( paste(dir.name,"firstsampledtpc.R",sep='/') )
		cat(paste("\n total new infections",tpc.data[["attack.rate"]] * nrow(ibm[["curr.pop"]])))		
		tpc.data<- tpc.collapse(tpc.data)	
		#tmp<- tpc.init(ibm)
simu.time<- proc.time()		
		tpc.t<- tpc.tabulate(tpc.data)
print(	(proc.time()-simu.time)[3] )
		print(tpc.t)
simu.time<- proc.time()		
		
		sampled.individuals<- ibm.sample(ssa.ibm, 0.5)
#print(sampled.individuals)
print(	(proc.time()-simu.time)[3] )
simu.time<- proc.time()				
		tpc.t<- tpc.tabulate.after.sample(tpc.data, sampled.individuals)
print(	(proc.time()-simu.time)[3] )			

		print(tpc.t)
		stop()		
		
		tpc.agg.data<- tpc.tabulate(tpc.data)
		print(tpc.agg.data)

	}
	stop()	
}
###############################################################################
prj.acute.test.lkl.nosampling.bothUandT	<- function()
{
	require(RColorBrewer)
	fixed.susc			<- 0
	m.type				<- "Acute"
	loc.type			<- "Town II"
	m.popsize			<- NA
	m.known.states		<- 1
	theta0				<- c(8, 0.05, 0, 0)
	theta0				<- c(4, 0.058, 0, 0)
	theta0				<- c(2, 0.065, 0, 0)
	names(theta0)		<- c("acute","base","m.st1","m.st2")
	sample.prob			<- c(0.5,0.5)
	names(sample.prob)	<- c("Idx","E")
	cluster.tw			<- 3
	m.repeat			<- 1
	resume				<- 0
	verbose				<- 0	
	
	if(fixed.susc)
		dir.name		<- paste(DATA,"acutesimu_fxs",sep='/')
	else
		dir.name		<- paste(DATA,"acutesimu_nofxs",sep='/')
	f.name				<- paste(dir.name,paste("tpcdat_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta0[1],"_b",theta0[2],"_sIdx",sample.prob[1],"_sE",sample.prob[2],"_tw",cluster.tw,"_median",sep=''),sep='/')
	cat(paste("\nload ",paste(f.name,".R",sep='')))		
	load(paste(f.name,".R",sep=''))
	print(tpc.table.all.median)
	
	
	if(resume)
	{
		options(show.error.messages = FALSE)
		f.name			<- paste(dir.name,paste("tpcdat_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta0[1],"_b",theta0[2],"_sIdx",sample.prob[1],"_sE",sample.prob[2],"_tw",cluster.tw,"_lkl",sep=''),sep='/')		
		cat(paste("\ntry load ",paste(f.name,".R",sep='')))				
		readAttempt		<-try(suppressWarnings(load(paste(f.name,".R",sep=''))))
		options(show.error.messages = TRUE)
		if(!inherits(readAttempt, "try-error"))
			cat(paste("\nresumed file ",paste(f.name,".R",sep='')))
	}
	
	theta.acute		<- seq(1,20,0.25)
	theta.base		<- seq(0.03,0.085,0.002)
	theta			<- expand.grid(acute= theta.acute, base= theta.base)
	
	if(!resume || inherits(readAttempt, "try-error"))
	{
		
		#theta			<- matrix( c(2,4,6,8,    0.065, 0.058, 0.053, 0.05), 4,2, dimnames=list(c(),c("acute","base")) )		
		clu.n			<- clu.tipc.n(acute.MAX.TIPC.SIZE)
		
		ibm				<- ibm.collapse( ibm.init.model(m.type, loc.type, NA, theta0, save='', resume= 0, init.pop=0) )
	
		state.n			<- ibm$init.pop.distr$status * ibm$init.pop.distr$npop
		pop.n			<- ibm$init.pop.distr$npop
		print(pop.n)
		print(state.n)
		print(nrow(theta))
		#evaluate likelihood over parameter space
		tipc.lkl		<- sapply(seq_len(nrow(theta)),function(i)
				{				
					ibm[["beta"]][['i']][["status"]]['i']	<- theta[i,"acute"]
					ibm[["beta"]][["base"]]					<- theta[i,"base"]	
					rate.m									<- acute.get.rates(ibm[["beta"]], ibm.pop= NULL, pop.n=pop.n, state.n= as.matrix(state.n), per.capita.i= 1)
					sample.prob[]							<- c(1,1)	
					dT										<- cluster.tw
					tpc										<- tpc.table.all.median				
					if(1)
					{
						lkl				<- acute.loglkl(tpc.table.all.median, rate.m, dT)[["table.lkl"]]
					}
					else
					{
						tpc.n.mx		<- ncol(tpc)-1														#max number of transmissions in tip cluster  	
						tpc.closure		<- ifelse(all(sample.prob==1), tpc.n.mx, acute.MAX.TIPC.SIZE)		#without sampling, need to compute probabilities only up to the largest number of transmissions + 1 in the any tip cluster
						clu.n			<- clu.n[seq_len(tpc.closure+1),seq_len(tpc.closure+1)]
						clu.n.sum		<- apply(clu.n,2,sum)	#cheaper than 'seq_len(ncol(clu.n))^seq.int(-1,ncol(clu.n)-2)'
						clu.n			<- clu.n / matrix(clu.n.sum, nrow=nrow(clu.n), ncol=length(clu.n.sum), byrow=1)								
						part1			<- acute.subtree.lkl.PartNT(seq.int(0,ncol(clu.n)-1),rate.m['u','s'],rate.m['i','s'],dT, log=1)
						part1			<- matrix(part1,nrow(clu.n),length(part1),byrow=1,dimnames=list(paste('idx',seq.int(0,nrow(clu.n)-1),sep=''),paste('n',seq.int(0,length(part1)-1),sep='')))
						part2			<- acute.subtree.lkl.PartIdx(ncol(clu.n)-1,rate.m['u','s'],rate.m['i','s'], log=1)
						lkl.U			<- part1 + part2 + log(clu.n)
						dimnames(lkl.U)	<- dimnames(part2)
						lkl.U			<- apply(lkl.U,2,function(x) log(sum(exp(x), na.rm=1)))
						
						part1			<- acute.subtree.lkl.PartNT(seq.int(0,ncol(clu.n)-1),rate.m['t','s'],rate.m['i','s'],dT, log=1)
						part1			<- matrix(part1,nrow(clu.n),length(part1),byrow=1,dimnames=list(paste('idx',seq.int(0,nrow(clu.n)-1),sep=''),paste('n',seq.int(0,length(part1)-1),sep='')))
						part2			<- acute.subtree.lkl.PartIdx(ncol(clu.n)-1,rate.m['t','s'],rate.m['i','s'], log=1)
						lkl.T			<- part1 + part2 + log(clu.n)
						dimnames(lkl.T)	<- dimnames(part2)
						lkl.T			<- apply(lkl.T,2,function(x) log(sum(exp(x), na.rm=1)))
						
						init.freq.mle	<- apply(tpc.table.all.median,1,sum) / sum(tpc.table.all.median)
						lkl.U			<- lkl.U + log(init.freq.mle)['u']
						lkl.T			<- lkl.T + log(init.freq.mle)['t']
						
						lkl				<- rbind(rep(0,length(lkl.U)),lkl.T,lkl.U)
						lkl				<- sum( tpc.table.all.median * lkl )
					}
	#print(lkl); print(init.freq.mle); stop()
					lkl
				})	
		names(tipc.lkl)<- apply(theta,1,function(x)	paste(x,collapse='_',sep='') )
		f.name			<- paste(dir.name,paste("tpcdat_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta0[1],"_b",theta0[2],"_sIdx",sample.prob[1],"_sE",sample.prob[2],"_tw",cluster.tw,"_lkl",sep=''),sep='/')		
		cat(paste("\nsave tipc.lkl to file ",paste(f.name,".R",sep='')))
		save(tipc.lkl, file=paste(f.name,".R",sep=''))		
	}
	
	if(0)	#compute log likelihood of tip cluster table
	{
		#print(tipc.lkl)
		#print(tpc.table.all.median)
		
		x		<- theta.acute
		y		<- tipc.lkl
		plot(1,1,type='n',bty='n',xlim=range(x),ylim=range(y),xlab=expression(beta[EE]/beta[UE]),ylab="tipc log likelihood")
		lines(x,y)
		abline(v=theta0[1], col="red")
		theta.mle<- theta[which.max(y), ]
		abline(v=theta.mle[1])
		print(theta.mle)
		stop()
	}			
	if(1)
	{
		theta		<- sapply( strsplit(names(tipc.lkl),'_',fixed=1), as.numeric )
		theta.acute	<- unique(theta[1,])
		theta.base	<- unique(theta[2,])	
		print(theta.acute)
		print(theta.base)
		#stop()
		#note: theta.acute runs first in lkls
		
		#lkls.mean<- apply(lkl,2,mean)
		#lkls.mean<- apply(lkls,2,function(x) sd(x)/mean(x))
		tipc.lkl.m	<- matrix(tipc.lkl, nrow= length(theta.acute), ncol= length(theta.base), dimnames=list(theta.acute,theta.base))
		cat(paste("\nmax lkl is ",max(tipc.lkl.m)))
		theta.mle	<- as.numeric(strsplit(names(tipc.lkl)[which.max(tipc.lkl)],'_')[[1]])
		print(theta.mle)
		f.name		<- paste(dir.name,paste("tpc_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta0[1],"_b",theta0[2],"_s",m.known.states,sep=''),sep='/')
		cat(paste("\nplot",paste(f.name,"_2D.pdf",sep='')))
		pdf(paste(f.name,"_2D.pdf",sep=''),version="1.4",width=4,height=4)
		my.image(theta.acute,theta.base,tipc.lkl.m, xlab=expression(beta[EE]/beta[UE]),ylab=expression(beta[0]),nlevels=10)			
		abline(v=theta.mle[1],lty=2)
		abline(h=theta.mle[2],lty=2)
		points(theta0[1],theta0[2], pch=19, col="red")
		dev.off()		
	}
}
###############################################################################
prj.birth.test.lkl.wsampling<- function()
{
	require(RColorBrewer)
	m.type				<- "Acute"
	loc.type			<- "Town II"
	m.popsize			<- NA
	m.known.states		<- 1
	theta0				<- c(1, 0.05, 0, 0)
	names(theta0)		<- c("acute","base","m.st1","m.st2")
	sample.prob0		<- 0.2
	cluster.tw			<- 3
	m.repeat			<- 1
	resume				<- 0
	verbose				<- 0	
	#dir.name			<- paste(DATA,"acutesimu_fxs_onlyu",sep='/')
	dir.name			<- paste(DATA,"acutesimu_fxs_onlyu_tpcsample",sep='/')
	f.name				<- paste(dir.name,paste("tpcdat_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta0[1],"_b",theta0[2],"_sIdx",sample.prob0,"_sE",sample.prob0,"_tw",cluster.tw,"_median",sep=''),sep='/')
	cat(paste("\nload ",paste(f.name,".R",sep='')))		
	load(paste(f.name,".R",sep=''))	
	tpc.table.sample.median['u',]	<- tpc.table.sample.median['u',] + tpc.table.sample.median['i',]
	tpc.table.sample.median['i',]	<- 0
	
	print(tpc.table.all.median)
	#tmp					<- clu.sample(tpc.table.all.median, { tmp<- rep(sample.prob0,2); names(tmp)<- c("Idx","E"); tmp }, rtn.exp=1)
	#print(round(tmp))
	print(tpc.table.sample.median)	
#stop()	
	theta.base		<- seq(0.03,0.085,0.002)
	theta.sample	<- seq(0.1,1,0.1)		
	theta			<- expand.grid(sample= theta.sample, base= theta.base)
	#f.name.head		<- "tpcsampled_Uonly"
	f.name.head		<- "tpcsampled_Uonly_bdfun"
	if(resume)
	{
		options(show.error.messages = FALSE)
		f.name			<- paste(dir.name,paste(f.name.head,"_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta0[1],"_b",theta0[2],"_sIdx",sample.prob0,"_sE",sample.prob0,"_tw",cluster.tw,"_lkl",sep=''),sep='/')		
		cat(paste("\ntry load ",paste(f.name,".R",sep='')))				
		readAttempt		<-try(suppressWarnings(load(paste(f.name,".R",sep=''))))
		options(show.error.messages = TRUE)
		if(!inherits(readAttempt, "try-error"))
			cat(paste("\nresumed file ",paste(f.name,".R",sep='')))
	}
	if(!resume || inherits(readAttempt, "try-error"))
	{				
		ibm										<- ibm.collapse( ibm.init.model(m.type, loc.type, NA, theta0, save='', resume= 0, init.pop=0) )
		ibm[["beta"]][['i']][["status"]]['i']	<- 1
		ibm[["beta"]][['i']][["status"]]['t']	<- 1
		ibm[["init.pop.distr"]][["status"]]['u']<- ibm[["init.pop.distr"]][["status"]]['t'] + ibm[["init.pop.distr"]][["status"]]['u']
		ibm[["init.pop.distr"]][["status"]]['t']<- 0		
		state.n									<- ibm$init.pop.distr$status * ibm$init.pop.distr$npop
		pop.n									<- ibm$init.pop.distr$npop	
		print(pop.n)
		print(state.n)	

		acute.MAX.TIPC.SIZE	<<- 10 			
		clu.n			<- clu.tipc.n(acute.MAX.TIPC.SIZE)			
		clu.n.sum		<- apply(clu.n,2,sum)																	#faster than seq_len(ncol(clu.n))^seq.int(-1,ncol(clu.n)-2)
		lclu.n			<- log( clu.n / matrix(clu.n.sum, nrow=nrow(clu.n), ncol=length(clu.n.sum), byrow=1) ) 		
		#evaluate likelihood over parameter space
		tipc.lkl		<- sapply(seq_len(nrow(theta)),function(i)
				{				
#i<- 10
					print(theta[i,])				
					ibm[["beta"]][["base"]]					<- theta[i,"base"]	
					s										<- theta[i,"sample"]
					rate.m									<- acute.get.rates(ibm[["beta"]], ibm.pop= NULL, pop.n=pop.n, state.n= as.matrix(state.n), per.capita.i= 1)				
					dT										<- cluster.tw					
#					lkl										<- birthsampled.loglkl(tpc.table.sample.median['u',,drop=0], rate.m['u','s'], s, dT)
					lkl2									<- acutesampled.loglkl(tpc.table.sample.median, rate.m, s, dT, lclu.n=lclu.n)				
#print(lkl[[1]]); print(exp(lkl[[2]])); stop()
#print(lkl); print(lkl2); stop()
					lkl2[["table.lkl"]]
				})	
		names(tipc.lkl)<- apply(theta,1,function(x)	paste(x,collapse='_',sep='') )
		
		f.name			<- paste(dir.name,paste(f.name.head,"_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta0[1],"_b",theta0[2],"_sIdx",sample.prob0,"_sE",sample.prob0,"_tw",cluster.tw,"_lkl",sep=''),sep='/')		
		cat(paste("\nsave tipc.lkl to file ",paste(f.name,".R",sep='')))
		save(tipc.lkl, file=paste(f.name,".R",sep=''))				
	}
	#print(tipc.lkl)
	#stop()
	if(1)
	{
		theta		<- sapply( strsplit(names(tipc.lkl),'_',fixed=1), as.numeric )
		theta.sample<- unique(theta[1,])
		theta.sample	<- unique(theta[1,])
		theta.base	<- unique(theta[2,])	
		#print(theta.sample); print(theta.base)
		#stop()
		#note: theta.sample runs first in lkls
		
		#lkls.mean<- apply(lkl,2,mean)
		#lkls.mean<- apply(lkls,2,function(x) sd(x)/mean(x))
		tipc.lkl.m	<- matrix(tipc.lkl, nrow= length(theta.sample), ncol= length(theta.base), dimnames=list(theta.sample,theta.base))
		cat(paste("\nmax lkl is ",max(tipc.lkl.m)))
		theta.mle	<- as.numeric(strsplit(names(tipc.lkl)[which.max(tipc.lkl)],'_')[[1]])
		print(theta.mle)
		f.name		<- paste(dir.name,paste(f.name.head,"_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta0[1],"_b",theta0[2],"_s",sample.prob0,sep=''),sep='/')
		cat(paste("\nplot",paste(f.name,"_2D.pdf",sep='')))
		pdf(paste(f.name,"_2D.pdf",sep=''),version="1.4",width=4,height=4)
		my.image(theta.sample,theta.base,tipc.lkl.m, xlab="sample",ylab=expression(beta[0]),nlevels=20)			
		abline(v=theta.mle[1],lty=2)
		abline(h=theta.mle[2],lty=2)
		points(sample.prob0,theta0[2], pch=19, col="red")
		dev.off()		
	}
	if(0)	#compute log likelihood of tip cluster table
	{
		print(tipc.lkl)		
		x		<- theta.sample
		y		<- tipc.lkl
		plot(1,1,type='n',bty='n',xlim=range(x),ylim=range(y),xlab=expression(beta[EE]/beta[UE]),ylab="tipc log likelihood")
		lines(x,y)
		abline(v=theta0[1], col="red")
		theta.mle<- theta[which.max(y), ]
		abline(v=theta.mle[1])
		print(theta.mle)
		stop()
	}	
	stop()
}
###############################################################################
prj.acute.test.lkl.wsampling.onlyU<- function()
{
	require(RColorBrewer)
	m.type				<- "Acute"
	loc.type			<- "Town II"
	m.popsize			<- NA
	m.known.states		<- 1
	theta0				<- c(8, 0.05, 0, 0)
	theta0				<- c(2, 0.065, 0, 0)
	names(theta0)		<- c("acute","base","m.st1","m.st2")
	sample.prob			<- c(0.6,0.6)
	names(sample.prob)	<- c("Idx","E")
	cluster.tw			<- 3
	m.repeat			<- 1
	resume				<- 0
	verbose				<- 0	
	#dir.name			<- paste(DATA,"acutesimu_fxs_onlyu",sep='/')
	dir.name			<- paste(DATA,"acutesimu_fxs_onlyu_tpcsample",sep='/')
	f.name				<- paste(dir.name,paste("tpcdat_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta0[1],"_b",theta0[2],"_sIdx",sample.prob[1],"_sE",sample.prob[2],"_tw",cluster.tw,"_median",sep=''),sep='/')
	cat(paste("\nload ",paste(f.name,".R",sep='')))		
	load(paste(f.name,".R",sep=''))
	print(tpc.table.all.median)
	print(tpc.table.sample.median)
	
	if(any(diff(sample.prob)!=0))	stop("cannot handle unequal sampling probabilities")
	
	theta.acute		<- seq(1,20,1)
	theta.base		<- theta0[2]
	theta.acute		<- seq(1,20,0.25)
	theta.base		<- seq(0.03,0.085,0.002)	
	theta			<- expand.grid(acute= theta.acute, base= theta.base)
		
	if(resume)
	{
		options(show.error.messages = FALSE)
		f.name			<- paste(dir.name,paste("tpcsampled_Uonly_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta0[1],"_b",theta0[2],"_sIdx",sample.prob[1],"_sE",sample.prob[2],"_tw",cluster.tw,"_lkl",sep=''),sep='/')		
		cat(paste("\ntry load ",paste(f.name,".R",sep='')))				
		readAttempt		<-try(suppressWarnings(load(paste(f.name,".R",sep=''))))
		options(show.error.messages = TRUE)
		if(!inherits(readAttempt, "try-error"))
			cat(paste("\nresumed file ",paste(f.name,".R",sep='')))
	}
	if(!resume || inherits(readAttempt, "try-error"))
	{		
		
		ibm				<- ibm.collapse( ibm.init.model(m.type, loc.type, NA, theta0, save='', resume= 0, init.pop=0) )	
		state.n			<- ibm$init.pop.distr$status * ibm$init.pop.distr$npop
		pop.n			<- ibm$init.pop.distr$npop	
		print(pop.n)
		print(state.n)	
				
		acute.MAX.TIPC.SIZE	<<- 12 			
		clu.n			<- clu.tipc.n(acute.MAX.TIPC.SIZE)			
		clu.n.sum		<- apply(clu.n,2,sum)																	#faster than seq_len(ncol(clu.n))^seq.int(-1,ncol(clu.n)-2)
		lclu.n			<- log( clu.n / matrix(clu.n.sum, nrow=nrow(clu.n), ncol=length(clu.n.sum), byrow=1) ) 		
		#evaluate likelihood over parameter space
		tipc.lkl		<- sapply(seq_len(nrow(theta)),function(i)
				{				
#i<- 10
#					print(theta[i,])
					ibm[["beta"]][['i']][["status"]]['i']	<- theta[i,"acute"]
					ibm[["beta"]][["base"]]					<- theta[i,"base"]	
					rate.m									<- acute.get.rates(ibm[["beta"]], ibm.pop= NULL, pop.n=pop.n, state.n= as.matrix(state.n), per.capita.i= 1)	
					dT										<- cluster.tw
					lkl										<- acutesampled.loglkl(tpc.table.sample.median, rate.m, sample.prob[1], dT, lclu.n=lclu.n)				
#print(lkl); stop()
					lkl[["table.lkl"]]
				})	
		names(tipc.lkl)<- apply(theta,1,function(x)	paste(x,collapse='_',sep='') )
		
		f.name			<- paste(dir.name,paste("tpcsampled_Uonly_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta0[1],"_b",theta0[2],"_sIdx",sample.prob[1],"_sE",sample.prob[2],"_tw",cluster.tw,"_lkl",sep=''),sep='/')		
		cat(paste("\nsave tipc.lkl to file ",paste(f.name,".R",sep='')))
		save(tipc.lkl, file=paste(f.name,".R",sep=''))				
	}
	#print(tipc.lkl)	
	if(1)
	{
		theta		<- sapply( strsplit(names(tipc.lkl),'_',fixed=1), as.numeric )
		theta.acute	<- unique(theta[1,])
		theta.base	<- unique(theta[2,])	
		print(theta.acute)
		print(theta.base)
		#stop()
		#note: theta.acute runs first in lkls
		
		#lkls.mean<- apply(lkl,2,mean)
		#lkls.mean<- apply(lkls,2,function(x) sd(x)/mean(x))
		tipc.lkl.m	<- matrix(tipc.lkl, nrow= length(theta.acute), ncol= length(theta.base), dimnames=list(theta.acute,theta.base))
		cat(paste("\nmax lkl is ",max(tipc.lkl.m)))
		theta.mle	<- as.numeric(strsplit(names(tipc.lkl)[which.max(tipc.lkl)],'_')[[1]])
		print(theta.mle)
		f.name		<- paste(dir.name,paste("tpcsampled_Uonly_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta0[1],"_b",theta0[2],"_s",sample.prob[1],sep=''),sep='/')
		cat(paste("\nplot",paste(f.name,"_2D.pdf",sep='')))
		pdf(paste(f.name,"_2D.pdf",sep=''),version="1.4",width=4,height=4)
		my.image(theta.acute,theta.base,tipc.lkl.m, xlab=expression(beta[EE]/beta[UE]),ylab=expression(beta[0]),nlevels=20)			
		abline(v=theta.mle[1],lty=2)
		abline(h=theta.mle[2],lty=2)
		points(theta0[1],theta0[2], pch=19, col="red")
		dev.off()		
	}
	if(0)	#compute log likelihood of tip cluster table
	{
		print(tipc.lkl)		
		x		<- theta.acute
		y		<- tipc.lkl
		plot(1,1,type='n',bty='n',xlim=range(x),ylim=range(y),xlab=expression(beta[EE]/beta[UE]),ylab="tipc log likelihood")
		lines(x,y)
		abline(v=theta0[1], col="red")
		theta.mle<- theta[which.max(y), ]
		abline(v=theta.mle[1])
		print(theta.mle)
		stop()
	}	
	stop()
}
###############################################################################
prj.acute.test.lkl.wsampling.bothUandT	<- function()
{
	require(RColorBrewer)
	m.type				<- "Acute"
	loc.type			<- "Town II"
	m.popsize			<- NA
	m.known.states		<- 1
	theta0				<- c(8, 0.05, 0, 0)
	names(theta0)		<- c("acute","base","m.st1","m.st2")
	sample.prob			<- c(0.5,0.5)
	names(sample.prob)	<- c("Idx","E")
	cluster.tw			<- 3
	m.repeat			<- 1
	resume				<- 1
	verbose				<- 0	
	dir.name			<- paste(DATA,"acutesimu_fxs",sep='/')
	f.name				<- paste(dir.name,paste("tpcdat_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta0[1],"_b",theta0[2],"_sIdx",sample.prob[1],"_sE",sample.prob[2],"_tw",cluster.tw,"_median",sep=''),sep='/')
	cat(paste("\nload ",paste(f.name,".R",sep='')))		
	load(paste(f.name,".R",sep=''))
	print(tpc.table.all.median)
	print(tpc.table.sample.median)
			
	theta.acute		<- seq(1,20,0.25)
	theta.base		<- seq(0.03,0.085,0.002)
	theta			<- expand.grid(acute= theta.acute, base= theta.base)
				
	if(resume)
	{
		options(show.error.messages = FALSE)
		f.name			<- paste(dir.name,paste("tpcsampled_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta0[1],"_b",theta0[2],"_sIdx",sample.prob[1],"_sE",sample.prob[2],"_tw",cluster.tw,"_lkl",sep=''),sep='/')		
		cat(paste("\ntry load ",paste(f.name,".R",sep='')))				
		readAttempt		<-try(suppressWarnings(load(paste(f.name,".R",sep=''))))
		options(show.error.messages = TRUE)
		if(!inherits(readAttempt, "try-error"))
			cat(paste("\nresumed file ",paste(f.name,".R",sep='')))
	}
	if(!resume || inherits(readAttempt, "try-error"))
	{		
		acute.MAX.TIPC.SIZE	<<- 20 
		ibm				<- ibm.collapse( ibm.init.model(m.type, loc.type, NA, theta0, save='', resume= 0, init.pop=0) )
		state.n			<- ibm$init.pop.distr$status * ibm$init.pop.distr$npop
		pop.n			<- ibm$init.pop.distr$npop	
		print(pop.n)
		print(state.n)	
#stop()
		clu.n			<- clu.tipc.n(acute.MAX.TIPC.SIZE)			
		clu.n.sum		<- apply(clu.n,2,sum)																	#faster than seq_len(ncol(clu.n))^seq.int(-1,ncol(clu.n)-2)
		lclu.n			<- log( clu.n / matrix(clu.n.sum, nrow=nrow(clu.n), ncol=length(clu.n.sum), byrow=1) ) 		
		#evaluate likelihood over parameter space
		tipc.lkl		<- sapply(seq_len(nrow(theta)),function(i)
				{				
#i<- 2
print(theta[i,])
					ibm[["beta"]][['i']][["status"]]['i']	<- theta[i,"acute"]
					ibm[["beta"]][["base"]]					<- theta[i,"base"]	
					rate.m									<- acute.get.rates(ibm[["beta"]], ibm.pop= NULL, pop.n=pop.n, state.n= as.matrix(state.n), per.capita.i= 1)	
					dT										<- cluster.tw
					lkl										<- acutesampled.loglkl(tpc.table.sample.median, rate.m, sample.prob[1], dT, lclu.n=lclu.n)[["table.lkl"]]				
#print(lkl); stop()
					lkl
				})	
		names(tipc.lkl)<- apply(theta,1,function(x)	paste(x,collapse='_',sep='') )
		
		f.name			<- paste(dir.name,paste("tpcsampled_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta0[1],"_b",theta0[2],"_sIdx",sample.prob[1],"_sE",sample.prob[2],"_tw",cluster.tw,"_lkl",sep=''),sep='/')		
		cat(paste("\nsave tipc.lkl to file ",paste(f.name,".R",sep='')))
		save(tipc.lkl, file=paste(f.name,".R",sep=''))				
	}
	print(tipc.lkl)	
	if(1)
	{
		theta		<- sapply( strsplit(names(tipc.lkl),'_',fixed=1), as.numeric )
		theta.acute	<- unique(theta[1,])
		theta.base	<- unique(theta[2,])	
		print(theta.acute)
		print(theta.base)
		#stop()
		#note: theta.acute runs first in lkls
		
		#lkls.mean<- apply(lkl,2,mean)
		#lkls.mean<- apply(lkls,2,function(x) sd(x)/mean(x))
		tipc.lkl.m	<- matrix(tipc.lkl, nrow= length(theta.acute), ncol= length(theta.base), dimnames=list(theta.acute,theta.base))
		cat(paste("\nmax lkl is ",max(tipc.lkl.m)))
		theta.mle	<- as.numeric(strsplit(names(tipc.lkl)[which.max(tipc.lkl)],'_')[[1]])
		print(theta.mle)
		f.name		<- paste(dir.name,paste("tpc_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta0[1],"_b",theta0[2],"_s",m.known.states,sep=''),sep='/')
		cat(paste("\nplot",paste(f.name,"_2D.pdf",sep='')))
		#pdf(paste(f.name,"_2D.pdf",sep=''),version="1.4",width=4,height=4)
		my.image(theta.acute,theta.base,tipc.lkl.m, xlab=expression(beta[EE]/beta[UE]),ylab=expression(beta[0]),nlevels=10)			
		abline(v=theta.mle[1],lty=2)
		abline(h=theta.mle[2],lty=2)
		points(theta0[1],theta0[2], pch=19, col="red")
		#dev.off()		
	}
	stop()
}
###############################################################################
prj.acute.test.lkl.nosampling.onlyU	<- function()
{
		require(RColorBrewer)
		m.type				<- "Acute"
		loc.type			<- "Town II"
		m.popsize			<- NA
		m.known.states		<- 1
		theta0				<- c(8, 0.05, 0, 0)
		names(theta0)		<- c("acute","base","m.st1","m.st2")
		sample.prob			<- c(0.5,0.5)
		names(sample.prob)	<- c("Idx","E")
		cluster.tw			<- 3
		m.repeat			<- 1
		resume				<- 1
		verbose				<- 0	
		dir.name			<- paste(DATA,"acutesimu_fxs_onlyu",sep='/')
		f.name				<- paste(dir.name,paste("tpcdat_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta0[1],"_b",theta0[2],"_sIdx",sample.prob[1],"_sE",sample.prob[2],"_tw",cluster.tw,"_median",sep=''),sep='/')
		cat(paste("\nload ",paste(f.name,".R",sep='')))		
		load(paste(f.name,".R",sep=''))
		print(tpc.table.all.median)
		
		theta.acute		<- seq(1,20,1)
		theta.base		<- theta0[2]
		theta			<- expand.grid(acute= theta.acute, base= theta.base)
		
		#theta			<- matrix( c(2,4,6,8,    0.065, 0.058, 0.053, 0.05), 4,2, dimnames=list(c(),c("acute","base")) )		
		clu.n			<- clu.tipc.n(acute.MAX.TIPC.SIZE)
		
		ibm				<- ibm.collapse( ibm.init.model(m.type, loc.type, NA, theta0, save='', resume= 0, init.pop=0) )
		#simplify problem: INITIAL FREQUENCY OF U IS 1
		ibm$init.pop.distr$status['u']<- ibm$init.pop.distr$status['u']+ibm$init.pop.distr$status['t']
		ibm$init.pop.distr$status['t']<- 0
		state.n			<- ibm$init.pop.distr$status * ibm$init.pop.distr$npop
		pop.n			<- ibm$init.pop.distr$npop
		print(pop.n)
		print(state.n)
		log				<- 1
		#stop()
		#evaluate likelihood over parameter space
		chain.lkl		<- lapply(seq_len(nrow(theta)),function(i)
				{					
					ibm[["beta"]][['i']][["status"]]['i']	<- theta[i,"acute"]
					ibm[["beta"]][["base"]]					<- theta[i,"base"]	
					rate.m									<- acute.get.rates(ibm[["beta"]], ibm.pop= NULL, pop.n=pop.n, state.n= as.matrix(state.n), per.capita.i= 1)
					sample.prob[]							<- c(1,1)	
					dT										<- cluster.tw
					tpc										<- tpc.table.all.median															
					tpc.n.mx		<- ncol(tpc)-1														#max number of transmissions in tip cluster  	
					tpc.closure		<- ifelse(all(sample.prob==1), tpc.n.mx, acute.MAX.TIPC.SIZE)		#without sampling, need to compute probabilities only up to the largest number of transmissions + 1 in the any tip cluster
					clu.n			<- clu.n[seq_len(tpc.closure+1),seq_len(tpc.closure+1)]
					clu.n.sum		<- apply(clu.n,2,sum)	#cheaper than 'seq_len(ncol(clu.n))^seq.int(-1,ncol(clu.n)-2)'
					clu.n			<- clu.n / matrix(clu.n.sum, nrow=nrow(clu.n), ncol=length(clu.n.sum), byrow=1)					
					part1			<- acute.subtree.lkl.PartNT(seq.int(0,ncol(clu.n)-1),rate.m['u','s'],rate.m['i','s'],dT, log=log)
					part1			<- matrix(part1,nrow(clu.n),length(part1),byrow=1,dimnames=list(paste('idx',seq.int(0,nrow(clu.n)-1),sep=''),paste('n',seq.int(0,length(part1)-1),sep='')))
					part2			<- acute.subtree.lkl.PartIdx(ncol(clu.n)-1,rate.m['u','s'],rate.m['i','s'], log=log)
#print(part1); print(part2); print(log(clu.n)); stop()
					if(log)
						ans			<- part1 + part2 + log(clu.n)
					else
						ans			<- part1 * part2 * clu.n
					dimnames(ans)	<- dimnames(part2)
					#ans				<- ans / sum(ans)	-- DON T KNOW normalizing constant - this goes to n=Inf
					ans
				})
		#plot tipc likelihoods to see if they make sense - yes they do
		if(log)
		{
			tipc.lkl	<- sapply(seq_along(chain.lkl), function(i)   apply(chain.lkl[[i]],2,function(x) log(sum(exp(x), na.rm=1)))   )	
		}
		else
		{
			tipc.lkl	<- sapply(seq_along(chain.lkl), function(i)   apply(chain.lkl[[i]],2,function(x) sum(x[x!=-Inf], na.rm=1))   )
			tipc.lkl	<- log(tipc.lkl)
		}
		colnames(tipc.lkl)<- apply(theta,1,function(x)	paste(x,collapse='_',sep='') )	
		
		if(0)	#compute log likelihood of tip cluster table
		{
			print(tipc.lkl)
			print(tpc.table.all.median)
			
			cnts	<- matrix( tpc.table.all.median['u',], nrow=ncol(tpc.table.all.median), ncol=ncol(tipc.lkl), dimnames=dimnames(tipc.lkl) )
			lkl		<- apply( tipc.lkl[seq_len(nrow(cnts)),] * cnts,2,sum)
			
			x		<- theta.acute
			y		<- lkl
			plot(1,1,type='n',bty='n',xlim=range(x),ylim=range(y),xlab=expression(beta[EE]/beta[UE]),ylab="tipc log likelihood")
			lines(x,y)
			abline(v=theta0[1], col="red")
			theta.mle<- theta[which.max(lkl), ]
			abline(v=theta.mle[1])
			print(theta.mle)
			stop()
		}		
		if(0)	#x-axis is theta[,"acute"]
		{			
			m	<- tipc.lkl
			x	<- theta[,"acute"]
			xlim<- range(x)		
			ylim<- range(m)
			cols<- brewer.pal(nrow(m), "Set3")
			plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim,xlab=expression(beta[EE]/beta[UE]),ylab="log tipc probability")
			sapply(seq_len(nrow(m)),function(i)
					{
						lines(x,m[i,],col=cols[i])			
					})
			legend("bottomright", fill=c("transparent",cols), border=NA, legend= c("#transmissions",seq.int(0,nrow(m)-1)), bty='n')
			stop()
		}
		if(1)	#x-axis is total transmissions
		{			
			dT	<- cluster.tw
			ibm[["beta"]][['i']][["status"]]['i']	<- 1
			ibm[["beta"]][["base"]]					<- theta0["base"]	
			rate.m									<- acute.get.rates(ibm[["beta"]], ibm.pop= NULL, pop.n=pop.n, state.n= as.matrix(state.n), per.capita.i= 1)
			lambda									<- rate.m['u','s']
			m										<- tipc.lkl
			z										<- exp(-lambda*dT)*(1-exp(-lambda*dT))^seq.int(0,nrow(m)-1)
						
			x	<- seq.int(0,nrow(m)-1)
			xlim<- range(x)		
			ylim<- range(m)
			cols<- colorRampPalette(c("red", "orange", "blue"), space = "Lab")( ncol(m) ) 
			plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim,xlab="total transmissions in tip cluster (one U, others E)",ylab="log tipc probability")
			sapply(seq_len(ncol(m)),function(j)
					{
						lines(x,m[,j],col=cols[j])			
					})
			#lines(x, log(z))
			lines(seq.int(0,ncol(tpc.table.all.median)-1),log( tpc.table.all.median['u',]/sum(tpc.table.all.median['u',] )))
			legend("bottomleft", fill=c("transparent",cols), border=NA, legend= c(expression(beta[EE]/beta[UE]),theta[,"acute"]), bty='n')
			stop()
		}	
}
###############################################################################
prj.acute.test	<- function()
{
	require(RColorBrewer)
	if(0)
	{
		print("here 1,3")		
		print( clu.subtrees.find(1,3) )		
		print("here 2,1")
		print( clu.subtrees.find(2, 1) )
		print("here 2,2")
		print( clu.subtrees.find(2, 2) )
		print("here 2,4")
		print( clu.subtrees.find(2, 4) )
		print("here 4,2")
		print( clu.subtrees.find(4, 2) )								
	}
	if(0)
	{				
		print("here 2,3")
		print( clu.subtrees.n(clu.subtrees.find(2, 3)) )
				
		print("here 2,0 -- should be 1")
		print( clu.subtrees.n(clu.subtrees.find(2, 0)) )				
		
		print("here 3,3 -- should be 48 54  6")
		print( clu.subtrees.n(clu.subtrees.find(3, 3)) )
		
		print("here 4,2 -- should be 12 12")
		print( clu.subtrees.n(clu.subtrees.find(4,2)) )
		
		print("here 2,4 -- should be 250 128 54")
		print( clu.subtrees.n(clu.subtrees.find(2,4)) )		
		stop()
	}
	if(0)	#test clu.n
	{			
		clu.n<- clu.tipc.n(10)
		tmp<- cbind(	as.matrix(apply(clu.n,2,sum)),
						as.matrix(seq.int(1,ncol(clu.n))^( seq.int(1,ncol(clu.n))-2 ))	)
		colnames(tmp)	<- c("incode","correct")
		print(tmp)	
		
		clu.n<- clu.tipc.n(15)
		print(clu.n)
		tmp<- cbind(	as.matrix(apply(clu.n,2,sum)),
				as.matrix(seq.int(1,ncol(clu.n))^( seq.int(1,ncol(clu.n))-2 ))	)
		colnames(tmp)	<- c("incode","correct")
		print(tmp)	
		stop()
	}
	if(0)
	{
		prj.acute.test.lkl.nosampling.onlyU()		
	}
	if(0)
	{
		prj.acute.test.lkl.nosampling.bothUandT()
	}
	if(0)
	{
		prj.birth.test.lkl.wsampling()
	}
	if(1)
	{
		prj.acute.test.lkl.wsampling.onlyU()
	}
	if(0)
	{
		prj.acute.test.lkl.wsampling.bothUandT()
	}
	
}
###############################################################################
prj.plotfisherhettransm<- function()
{
	cex.axis		<- 0.75
	
	r.overall		<- c(0.66)
	names(r.overall)<- c("Overall")
	r.age			<- c(NA,1.86,0.39,0.14)
	names(r.age)	<- c("Age group","<35","35-44",">45")
	r.cd4			<- c(NA,0.34,0.71,0.77,0.71,0.38)
	names(r.cd4)	<- c("CD4 count (cell/mcl)","<200","201-350","351-500",">500","Not known")
	r.vl			<- c(NA,0.06,0.42,1.17,1.59,2.29,0.56)
	names(r.vl)		<- c("Viral load (copies/ml)","<50","50-1000","1001-10000","10001-100001",">100000","Not known")
	r.cat			<- c(NA,5.67,1.28,0.06,0.96)
	names(r.cat)	<- c("Infection category","Recently infected","Chronic, untreated","Chronic, treated","Chronic, interrupted")
	r.aids			<- c(NA,0.81,0.08)
	names(r.aids)	<- c("AIDS","No","Yes")
	r.std			<- c(NA,0.52,6.25)
	names(r.std)	<- c("STD diagnosis","No","Yes")
	r.yr			<- c(NA,0.41,0.54,0.85)
	names(r.yr)		<- c("Year","2000-2001","2002-2003",">2004")
	
	r				<- rev(c(r.overall,r.age,r.cd4,r.vl,r.cat,r.aids,r.std,r.yr))
	names(r)		<- rev(names(c(r.overall,r.age,r.cd4,r.vl,r.cat,r.aids,r.std,r.yr)))
	
	col.ticks			<- rep("black",length(r))
	col.ticks[is.na(r)]	<- "white" 
	print(col.ticks)
	dir.name<- DATA
	f.name<- paste(dir.name,paste("fisher_transmhet.pdf",sep=''),sep='/')
	cat(paste("print to",f.name))
	pdf(paste(f.name),version="1.4",width=5,height=8)
	par(mar=c(5,9,1,1))
	plot(1,1,type='n',bty='n',xlim=range(r, na.rm=1),ylim=range(seq_along(r)),cex.lab=cex.axis,cex.axis=cex.axis,xlab="transmission rate/100PYFU",ylab='',log='x',yaxt='n')
	
	sapply(list(which(is.na(r)),which(!is.na(r))),function(x)
			{
				tmp<- ifelse(is.na(r[x[1]]),1.2,1)*cex.axis
				axis(2,at=x,col.ticks=col.ticks[x],labels=names(r)[x],cex.axis=tmp,las=1)				
			})	
	points(r,seq_along(r),pch=15)
	abline(v=r.overall,lty=3,col="gray50")
	dev.off()
	#print(r)
}
###############################################################################
prj.plotgatesfigs<- function()
{
	loc.type<- "ZA-C"
	m.popsize<- 3e4
	theta<- c(2, 0.112)
	#lambda_u 0.0969920 lambda_I is multiple
	
	#vary acute from 1 to 8 for different size of tip clusters
	tipc.p1<- c(0.9075633, 0.08389218, 0.00387736, 0.0001194702, 0.00454772)
	tipc.p2<- c(0.9075633, 0.08001482, 0.00587872, 0.0003174552, 0.006225735)
	tipc.p4<- c(0.9075633, 0.07296036, 0.008798088, 0.000830085, 0.0098482)
	tipc.p8<- c(0.9075633, 0.06122958, 0.01170424, 0.001814424, 0.01768849)
	
	tipc.p<- matrix(c(tipc.p1,tipc.p2,tipc.p4,tipc.p8),nrow= 5, ncol= 4, dimnames= list(1:5,paste("l",c(1,2,4,8),sep='')))
	print(tipc.p)
	
	dir.name<- DATA
	f.name<- paste(dir.name,paste("tipcprobs_1.pdf",sep=''),sep='/')
	pdf(paste(f.name),version="1.4",width=4,height=4)
	par(mar=c(4.5,4,.5,.5))
	plot(1,1,type='n',xlim=c(1,8),ylim=range(tipc.p[-1,]),ylab="prob",xlab=expression(r[I]))
	x<- c(1,2,4,8)
	lines(x,tipc.p[1,])
	lines(x,tipc.p[2,],col="darkblue")
	points(x,tipc.p[2,],col="darkblue",pch=19)
	#lines(x,tipc.p[3,],col="darkblue")
	#lines(x,tipc.p[4,],col="darkmagenta")
	lines(x,tipc.p[5,],col="red")
	points(x,tipc.p[5,],col="red",pch=19)
	legend(x=0.8,y=0.06,bty='n',fill=c("darkblue","red"),legend=c("tip cluster size 2","tip cluster size 3"))
	dev.off()
}
###############################################################################
prj.acutesampling.rejabc<- function()
{
	m.type<- "Acute"
	loc.type<- "ZA-C"
	m.popsize<- 3e4
	abc.sit<- NA
	theta0<- c(2, 0.112, 0.75, 0.5)#c(8, 0.09, 0.75, 0.5)#c(2, 0.112, 0.75, 0.5)
	theta.n<- 20
	names(theta0)<- c("acute","base","m.st1","m.st2")
	my.mkdir(DATA,"acutesamplingabc")
	
	if(exists("args"))
	{		
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,2),
									n= return(as.numeric(substr(arg,3,nchar(arg)))),NA)	}))
		if(length(tmp)>0) m.popsize<- tmp[1]
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,2),
									l= return(substr(arg,3,nchar(arg))),NA)	}))
		if(length(tmp)>0) loc.type<- tmp[1]	
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,2),
									a= return(as.numeric(substr(arg,3,nchar(arg)))),NA)	}))
		if(length(tmp)>0) theta0[1]<- tmp[1]
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,2),
									b= return(as.numeric(substr(arg,3,nchar(arg)))),NA)	}))
		if(length(tmp)>0) theta0[2]<- tmp[1]
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,2),
									c= return(as.numeric(substr(arg,3,nchar(arg)))),NA)	}))
		if(length(tmp)>0) theta0[3]<- tmp[1]
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,2),
									s= return(as.numeric(substr(arg,3,nchar(arg)))),NA)	}))
		if(length(tmp)>0) abc.sit<- tmp[1]
	}
	
	cat(paste("\nm.popsize ",m.popsize))
	cat(paste("\nloc.type ",loc.type))
	cat(paste("\ntheta0.acute ",theta0[1]))
	cat(paste("\ntheta0.base ",theta0[2]))
	cat(paste("\ntheta0.missing.st1 ",theta0[3]))
	cat(paste("\ntheta0.missing.st2 ",theta0[4]))
	cat(paste("\nabc.sit ",abc.sit))
						
	theta.acute.prior<- 	c(1,10)
	theta.base.prior<-		c(0.05,0.16)
	abc.nsample<- 			1.5e3
	abc.nit<-				abc.sit + abc.nsample
	
	
	dir.name<- paste(DATA,"acutesamplingabc",sep='/')
	if(is.na(abc.sit))	#load precomputed abc.simus
	{
		#collect all abc simus
		f.name<- paste(dir.name,paste("abcsimu_",m.type,"_",loc.type,"_n",m.popsize,"_mst1",theta0[3],"_all.R",sep=''),sep='/')
		cat(paste("\n try load",f.name))		
		options(show.error.messages = FALSE)	
		readAttempt<-try(suppressWarnings(load(f.name)))
		options(show.error.messages = TRUE)
		if(!inherits(readAttempt, "try-error"))
			cat(paste("\nresumed file ",f.name))
		else
		{
			f.name<- list.files(dir.name,pattern=paste("^abcsimu_",m.type,"_",loc.type,"_n",m.popsize,"_mst1",theta0[3],sep=''),full.names=1)		
			abc.simus<- unlist(lapply(seq_along(f.name),function(i)
					{
						cat(paste("\nload",f.name[i]))
						load(f.name[i])						
						abc.simus					
					}), recursive=0)
			f.name<- paste(dir.name,paste("abcsimu_",m.type,"_",loc.type,"_n",m.popsize,"_mst1",theta0[3],"_all.R",sep=''),sep='/')
			cat(paste("\n try save abc.simus to",f.name))					
			save(abc.simus,file= f.name)
		}
		cat(paste("\nnumber of loaded abc.simus",length(abc.simus)))		
	}
	else				#precompute abc.simus
	{
		f.name<- paste(dir.name,paste("abcsimu_",m.type,"_",loc.type,"_n",m.popsize,"_mst1",theta0[3],"_sit",abc.sit,".R",sep=''),sep='/')
		cat(paste("\ntry load",f.name))
		options(show.error.messages = FALSE)	
		readAttempt<-try(suppressWarnings(load(f.name)))
		options(show.error.messages = TRUE)	
		if(!inherits(readAttempt, "try-error"))
			cat(paste("\nloaded",f.name))
		else
		{
			#for abc.sit to abc.nit, draw theta and 1) compute lkl from full table and 2) compute mean tpc table under sampling
			abc.simus<- lapply(abc.sit:abc.nit,function(abc.i)
				{
					cat(paste("\nabc iteration",abc.i))
					ans<- vector("list",6)
					names(ans)<- c("lkl.complete","tpc.complete","tpc.meansampled","simu.time","simu.attackrate","simu.theta")
					
					simu.time<- proc.time()
					ans[["simu.theta"]]<- theta0
					ans[["simu.theta"]]["acute"]<- runif(1,theta.acute.prior[1],theta.acute.prior[2])
					ans[["simu.theta"]]["base"]<-  runif(1,theta.base.prior[1],theta.base.prior[2])
					cat(paste("\ntheta.acute",ans[["simu.theta"]]["acute"],"theta.base",ans[["simu.theta"]]["base"]))
					
					ibm<- ibm.init.model(m.type, loc.type, m.popsize, ans[["simu.theta"]], resume= 0)
					ibm<- ibm.collapse(ibm)
					lkl.rates<- acute.get.rates(ibm, per.capita.i=1)
					trees.complete<- ssa.run(ibm, ssa.ctime= 0, ssa.etime= 1,record.tpc= 1, verbose= 0, resume= 0)
					ans[["simu.attackrate"]]<- trees.complete[["attack.rate"]]
					trees.complete<- tpc.collapse(trees.complete)
					ans[["tpc.complete"]]<- tpc.tabulate(trees.complete)				
					ans[["lkl.complete"]]<- acute.loglkl(ans[["tpc.complete"]], lkl.rates, ibm[["beta"]][["dT"]])
					cat(paste("\ntrees simulated, tip cluster likelihood is",sum(ans[["lkl.complete"]])))
					tpc.sampled<- lapply(seq_len(abc.nsample), function(i)
							{				
								sampled.individuals<- ibm.sample(ibm, 1-ans[["simu.theta"]]["m.st1"])							
								tpc.tabulate.after.sample(trees.complete, sampled.individuals)
							})				
					ans[["tpc.meansampled"]]<- tpc.mean(tpc.sampled)
					ans[["simu.time"]]<- (proc.time()-simu.time)[3]	
					cat(paste("\nsimu time is",ans[["simu.time"]]))
					ans
				})
			dir.name<- paste(DATA,"acutesamplingabc",sep='/')
			f.name<- paste(dir.name,paste("abcsimu_",m.type,"_",loc.type,"_n",m.popsize,"_mst1",theta0[3],"_sit",abc.sit,".R",sep=''),sep='/')
			cat(paste("\nsave abc.simus to",f.name))
			save(abc.simus,file=f.name)
		}
		if(1)		#add mle s
		{
			#load observed data
			dir.name<- paste(DATA,"acutesamplingabc",sep='/')
			f.name<- paste(dir.name,paste("tpc_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta0[1],"_b",theta0[2],"_mst1",theta0[3],sep=''),sep='/')		
			cat(paste("\ntry load",paste(f.name,"_tpctables.R",sep='')))
			options(show.error.messages = FALSE)	
			readAttempt<-try(suppressWarnings(load(paste(f.name,"_tpctables.R",sep=''))))
			if(!inherits(readAttempt, "try-error")) 
				cat(paste("\nloaded",paste(f.name,"_tpctables.R",sep='')))
			else
				stop("prj.acutesampling.rejabc: failed to load observed data")
			obs.ibm<- ibm.as.data.table(obs.ibm)
		
			cat(paste("\n add simu.mle to abc.simus",f.name))
			abc.simus<- lapply(seq_along(abc.simus), function(i)
					{
						abc.simus[[i]][["simu.mle"]]<- tipc.mle(round(abc.simus[[i]][["tpc.meansampled"]]),obs.ibm,exp(seq(log(.01),log(10),length.out=theta.n)),seq(0.01,0.16,length.out=theta.n))
						abc.simus[[i]]
					})
			dir.name<- paste(DATA,"acutesamplingabc",sep='/')
			f.name<- paste(dir.name,paste("abcsimu_",m.type,"_",loc.type,"_n",m.popsize,"_mst1",theta0[3],"_sit",abc.sit,".R",sep=''),sep='/')
			cat(paste("\nsave abc.simus to",f.name))
			save(abc.simus,file=f.name)
			stop()
		}
		stop()
	}

	#get observed tpc
	dir.name<- paste(DATA,"acutesamplingabc",sep='/')
	f.name<- paste(dir.name,paste("tpc_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta0[1],"_b",theta0[2],"_mst1",theta0[3],sep=''),sep='/')		
	cat(paste("\ntry load",paste(f.name,"_tpctables.R",sep='')))
	options(show.error.messages = FALSE)	
	readAttempt<-try(suppressWarnings(load(paste(f.name,"_tpctables.R",sep=''))))
	options(show.error.messages = TRUE)
	
	if(!inherits(readAttempt, "try-error"))
		cat(paste("\nloaded",paste(f.name,"_tpctables.R",sep='')))
	else
	{
		dir.name<- paste(DATA,"simudata",sep='/')
		f.name<- paste(dir.name,paste("tpc_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta0[1],"_b",theta0[2],sep=''),sep='/')
		cat(paste("\ntry load",paste(f.name,"_tpc.R",sep='')))
		load(file=paste(f.name,"_tpc.R",sep=''))
		
		tpc.t<- lapply(seq_along(tpc), function(i)
				{				
					tpc[[i]][["ibm"]]<- ibm.as.data.table(tpc[[i]][["ibm"]])
					sampled.individuals<- ibm.sample(tpc[[i]][["ibm"]], 1-theta0["m.st1"])
					tpc.tabulate.after.sample(tpc[[i]][["tpc"]], sampled.individuals)
				})
		obs.mean.tpc<- tpc.mean(tpc.t)
		obs.ibm<- ibm.as.data.table(tpc[[1]][["ibm"]])
		obs.mle<- tipc.mle(round(obs.mean.tpc),obs.ibm,exp(seq(log(.01),log(10),length.out=theta.n)),seq(0.01,0.16,length.out=theta.n))
		dir.name<- paste(DATA,"acutesamplingabc",sep='/')
		f.name<- paste(dir.name,paste("tpc_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta0[1],"_b",theta0[2],"_mst1",theta0[3],sep=''),sep='/')		
		cat(paste("\nsave tpc.t to",paste(f.name,"_tpctables.R",sep='')))	
		save(obs.mean.tpc,obs.mle,obs.ibm, file=paste(f.name,"_tpctables.R",sep=''))
	}
	cat(paste("\nobserved mean tpc is\n"))
	print(obs.mean.tpc)
	#cat(paste("\nobserved mle tpc is\n"))
	#print(obs.mle)
	#print(obs.ibm)
	
	#compute mle for each abc simu
	if(0)
	{
		best.n<- seq_len(6)
		obs.mle<- c(apply( obs.mle[1:2,best.n], 1, mean),0.012)
		cat(paste("\nobserved mle tpc is\n"))
		print(obs.mle)		
		sim.mle<-	sapply(seq_along(abc.simus),function(i)
			{
				c( apply( abc.simus[[i]][["simu.mle"]][1:2,best.n], 1, mean), abc.simus[[i]][["simu.attackrate"]] )
			})
		rownames(sim.mle)<- c("acute","base","attack")
		print(sim.mle[,1:10])
		diff.mle<- sim.mle - obs.mle
		#diff.mle<- diff.mle[1:2,]
		#print(diff.mle[,1:10])	
		#diff.mle.sd<- apply(diff.mle,1,sd)
		#print(diff.mle.sd)
		#diff.mle<- diff.mle/diff.mle.sd
		#print(diff.mle[,1:10])
		diff.tolerance<- apply(abs(diff.mle), 1, function(x) quantile( x, probs=c(0.05,0.15,0.25)  ) )
		print(diff.tolerance)	
		diff.accept<- lapply(seq_len(nrow(diff.tolerance)),function(i)
				{
					which( apply( abs(diff.mle), 2, function(x) all(x<=diff.tolerance[i,])   ) )				
				})
		names(diff.accept)<- rownames( diff.tolerance )
		print(sapply(diff.accept,length))
		#stop()
	}
	if(1)
	{
		#compute distances between observed and simulated mean tpcs, and accepted abc.simus
		max.col<-	max(sapply(seq_along(abc.simus),function(i)  ncol( abc.simus[[i]][["tpc.meansampled"]] )  ))
		#print(max.col)	
		sim.mean.tpc<-	sapply(seq_along(abc.simus),function(i)
				{
					c( 	abc.simus[[i]][["tpc.meansampled"]]['t',]/sum(abc.simus[[i]][["tpc.meansampled"]]['t',]),rep(0,max.col-ncol(abc.simus[[i]][["tpc.meansampled"]])),
						abc.simus[[i]][["tpc.meansampled"]]['u',]/sum(abc.simus[[i]][["tpc.meansampled"]]['u',]),rep(0,max.col-ncol(abc.simus[[i]][["tpc.meansampled"]]))	)
				})
		sim.mean.tpc.sd<- apply(sim.mean.tpc, 1, sd)
		sim.mean.tpc.sd[sim.mean.tpc.sd==0]<- 1
		#sim.mean.tpc<- sim.mean.tpc/sim.mean.tpc.sd
		#print(sim.mean.tpc[,1:10])
		obs.mean.tpc<- c(obs.mean.tpc['t',]/sum(obs.mean.tpc['t',]),rep(0,max.col-ncol(obs.mean.tpc)),obs.mean.tpc['u',]/sum(obs.mean.tpc['u',]),rep(0,max.col-ncol(obs.mean.tpc)))
		#obs.mean.tpc<- obs.mean.tpc/sim.mean.tpc.sd
		#print(obs.mean.tpc)	
		diff.mean.tpc<- sim.mean.tpc - obs.mean.tpc
		print(diff.mean.tpc[,1:10])
		
		diff.mean.tpc<- apply(abs(diff.mean.tpc),2,sum)
		#print( c( c(1,rep(2,nrow(diff.mean.tpc)/2-1)),c(1,rep(2,nrow(diff.mean.tpc)/2-1)) ) )
		
		#diff.mean.tpc<- apply(abs(diff.mean.tpc),2,function(x)  x^c(1,rep(2,ncol(diff.mean.tpc)/2-1)) )
		diff.tolerance<- quantile( diff.mean.tpc, probs=c(0.005,0.01,0.02)  )
		diff.accept<- lapply(seq_along(diff.tolerance),function(i)
				{
					which(diff.mean.tpc<diff.tolerance[i])
				})
		names(diff.accept)<- diff.tolerance
	}
	#print(diff.accept)
	abc.accepted<- lapply(seq_along(diff.accept),function(i)
			{
				lapply(diff.accept[[i]],function(j)		abc.simus[[j]]	)					
			})
	names(abc.accepted)<- names(diff.accept)
	#print(abc.accepted[[1]])
	#stop()
	#print abc posterior based on difference in sampled tipc tables
	abc.stdposterior<- lapply(seq_along(abc.accepted),function(i)
			{
				tmp<- sapply( abc.accepted[[i]],function(x)	x[["simu.theta"]][c("acute","base")] )
				rownames(tmp)<- c("acute","base")
				tmp
			})
	names(abc.stdposterior)<- names(diff.accept)
	#tmp<- rev(sort(abc.stdposterior[[2]]["lkl.complete",],index.return=1)$ix)
	#print(tmp[1:20])
	k<- 2
	#print( abc.stdposterior[[k]][,tmp[1:20]] )
	#print( exp(abc.stdposterior[[k]][3,tmp[1:20]] ))
	
	
	f.name<- paste(dir.name,paste("abcpost_",m.type,"_",loc.type,"_n",m.popsize,"_mst1",theta0[3],".pdf",sep=''),sep='/')
	pdf(f.name,version="1.4",width=5,height=5)
	par(mar=c(4.5,4.5,.5,.5))
	my.plot.2D.dens(abc.stdposterior[[k]]["acute",],abc.stdposterior[[k]]["base",],expression(r[I]),expression(beta[0]), palette= "gray",method="ash",xlim=c(1,10),ylim=c(0.05,0.16))
	points(theta0[1],theta0[2],col="red",pch=19)
	dev.off()
	stop()
	my.image.smooth(abc.stdposterior[[k]]["acute",],abc.stdposterior[[k]]["base",],exp(abc.stdposterior[[k]]["lkl.complete",]), xlab=expression(r[a]), ylab=expression(beta[0]) )
	#abline(v=theta0[1],col="red",lty=2)
	#abline(h=theta0[2],col="red",lty=2)	
}
###############################################################################
prj.acute.clusterdistribution<- function()
{
	m.type				<- "Acute"
	loc.type			<- "ZA-C"
	m.popsize			<- 6e4
	m.known.states		<- 0
	theta0				<- c(8, 0.09, 0, 0)
	names(theta0)		<- c("acute","base","m.st1","m.st2")
	sample.prob			<- c(0.7,0.7)
	names(sample.prob)	<- c("Idx","E")	
	m.repeat			<- 1
	resume				<- 1
	verbose				<- 1	
	
	cat(paste("\nm.popsize ",m.popsize))
	cat(paste("\nloc.type ",loc.type))	
	cat(paste("\nm.known.states ",m.known.states))
	cat(paste("\nm.repeat ",m.repeat))
	cat(paste("\nacute.MAX.TIPC.SIZE ",acute.MAX.TIPC.SIZE))	
	cat(paste("\ntheta0.acute ",theta0[1]))
	cat(paste("\ntheta0.base ",theta0[2]))
	cat(paste("\nsample.prob.Idx ",sample.prob[1]))
	cat(paste("\nsample.prob.E ",sample.prob[2]))
	cat(paste("\nverbose ",verbose))
	cat(paste("\nresume ",resume))
	my.mkdir(DATA,"acutecld")
	
	dir.name<- paste(DATA,"simudata",sep='/')
	f.name<- paste(dir.name,paste("tpcdat_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta0[1],"_b",theta0[2],"_sIdx",sample.prob[1],"_sE",sample.prob[2],sep=''),sep='/')					
	#tmp<- paste(f.name,"_example_contingencytable_upto3.R",sep='')	
	tmp<- paste(f.name,".R",sep='')
	cat(paste("\nload ",tmp))		
	load(tmp)
	
	theta.acute		<- seq(1,10,1)
	theta			<- expand.grid(acute= theta.acute, base= theta0[2])
	print(theta)
	
	sample.prob2	<- sample.prob
	sample.prob2[]	<- c(1,1)
	#compare tip cluster distribution as 'acute' varies
	cld		<- lapply(seq_len(nrow(theta)),function(i)
			{																		
				ans	<- lapply(c(1), function(r)
						{							
							ibm		<- tpc[[r]][["ibm"]]
							if(!m.known.states)
							{
								ibm	<- ibm.init.model(m.type, loc.type, m.popsize, theta0, resume= 0)
								ibm	<- ibm.collapse(ibm)								
							}	
							ibm[["beta"]][['i']][["status"]]['i']	<- theta[i,"acute"]
							ibm[["beta"]][["base"]]					<- theta[i,"base"]	
							beta.stratified							<- acute.get.rates(ibm, per.capita.i=1)
							tipc.table								<- tpc[[r]][["tpc.table.all"]]
							#print(beta.stratified)
							#print(tipc.table)
							tipc.lkl								<- acute.loglkl(tipc.table, beta.stratified, tpc[[r]][["ibm"]][["beta"]][["dT"]])[["tipc.lkl"]]
							tipc.lkl #/ sum(tipc.lkl)							
						})
				names(ans)<- paste('r',seq_along(ans),sep='')
				ans[[1]]
			})
	
	tmp<- sapply(cld, function(x)   apply(x,2,sum)		)
	plot(1,1,type='n',xlim=range(theta.acute),ylim=range(tmp))
	cols<- rainbow(ncol(tmp))
	sapply(seq_len(ncol(tmp)),function(j)   lines(tmp[,j],col=cols[j]) )
	print(tmp)
	stop()
}
###############################################################################
prj.acute.loglklsurface<- function()
{
	#call syntax pkg/misc/phdes.startme.R -exeACUTE.LKL -acute=8 -baseline=0.09 -r=1 -v=1 -l=ZA-C -sIdx=0.7 -sE=0.7 -knownstates=1
	
	#how well can we estimate the true parameter of the acute model
	#- for different theta
	#- for different popsizes ?
	#- unknown population states ?
	
	m.type				<- "Acute"
	loc.type			<- "Town II"
	m.popsize			<- NA
	m.known.states		<- 1
	theta0				<- c(8, 0.05, 0, 0)
	names(theta0)		<- c("acute","base","m.st1","m.st2")
	sample.prob			<- c(0.5,0.5)
	names(sample.prob)	<- c("Idx","E")
	cluster.tw			<- 3
	m.repeat			<- 1
	resume				<- 1
	verbose				<- 0	
	if(exists("args"))
	{		
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,2),
									r= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) m.repeat<- tmp[1]
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,2),
									l= return(substr(arg,4,nchar(arg))),NA)	}))
		if(length(tmp)>0) loc.type<- tmp[1]
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,2),
									n= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) m.popsize<- tmp[1]
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,11),
									cluster.tw= return(as.numeric(substr(arg,13,nchar(arg)))),NA)	}))
		if(length(tmp)>0) cluster.tw<- tmp[1]
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,6),
									acute= return(as.numeric(substr(arg,8,nchar(arg)))),NA)	}))
		if(length(tmp)>0) theta0[1]<- tmp[1]
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,9),
									baseline= return(as.numeric(substr(arg,11,nchar(arg)))),NA)	}))
		if(length(tmp)>0) theta0[2]<- tmp[1]
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,5),
									sIdx= return(as.numeric(substr(arg,7,nchar(arg)))),NA)	}))
		if(length(tmp)>0) sample.prob[1]<- tmp[1]
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,3),
									sE= return(as.numeric(substr(arg,5,nchar(arg)))),NA)	}))
		if(length(tmp)>0) sample.prob[2]<- tmp[1]				
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]								
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,12),
									knownstates= return(as.numeric(substr(arg,14,nchar(arg)))),NA)	}))
		if(length(tmp)>0) m.known.states<- tmp[1]
	}
	cat(paste("\nm.popsize ",m.popsize))
	cat(paste("\nloc.type ",loc.type))	
	cat(paste("\nm.known.states ",m.known.states))
	cat(paste("\nm.repeat ",m.repeat))
	cat(paste("\nacute.MAX.TIPC.SIZE ",acute.MAX.TIPC.SIZE))	
	cat(paste("\ntheta0.acute ",theta0[1]))
	cat(paste("\ntheta0.base ",theta0[2]))
	cat(paste("\nsample.prob.Idx ",sample.prob[1]))
	cat(paste("\nsample.prob.E ",sample.prob[2]))
	cat(paste("\nverbose ",verbose))
	cat(paste("\nresume ",resume))
	my.mkdir(DATA,"acutelkl")		
	if(resume)
	{
		options(show.error.messages = FALSE)
		dir.name	<- paste(DATA,"acutelkl",sep='/')
		f.name		<- paste(dir.name,paste("tpclkl_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta0[1],"_b",theta0[2],"_sIdx",sample.prob[1],"_sE",sample.prob[2],"_st",m.known.states,"_si",acute.MAX.TIPC.SIZE,sep=''),sep='/')
				
		cat(paste("\ntry load ",paste(f.name,".R",sep='')))				
		readAttempt	<-try(suppressWarnings(load(paste(f.name,".R",sep=''))))
		options(show.error.messages = TRUE)
		if(!inherits(readAttempt, "try-error"))
			cat(paste("\nresumed file ",paste(f.name,".R",sep='')))
	}
	if(!resume || inherits(readAttempt, "try-error"))
	{
		dir.name	<- paste(DATA,"acutesimu_fxs_onlyu",sep='/')
		f.name		<- paste(dir.name,paste("tpcdat_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta0[1],"_b",theta0[2],"_sIdx",sample.prob[1],"_sE",sample.prob[2],"_tw",cluster.tw,"_median",sep=''),sep='/')
		cat(paste("\nload ",paste(f.name,".R",sep='')))		
		load(paste(f.name,".R",sep=''))

		theta.acute		<- seq(1,16,1)
		theta.base		<- theta0[2]
		theta			<- expand.grid(acute= theta.acute, base= theta.base)
		tipc.n			<- clu.tipc.n(acute.MAX.TIPC.SIZE)
		print(theta)
		
		ibm				<- ibm.init.model(m.type, loc.type, NA, theta0, save='', resume= 0, init.pop=0)
		ibm				<- ibm.collapse(ibm)
		state.n			<- ibm$init.pop.distr$status * ibm$init.pop.distr$npop
		pop.n			<- ibm$init.pop.distr$npop
		print(pop.n)
		print(as.matrix(state.n))

		#evaluate likelihood over parameter space
		loglkl			<- sapply(seq_len(nrow(theta)),function(i)
				{
					if(verbose) cat(paste("\nprocess theta",theta[i,"acute"], theta[i,"base"]))														
					if(!m.known.states)
					{
						ibm	<- ibm.init.model(m.type, loc.type, m.popsize, theta0, resume= 0)
						ibm	<- ibm.collapse(ibm)
						stop("not fully implemented")
					}	
					ibm[["beta"]][['i']][["status"]]['i']	<- theta[i,"acute"]
					ibm[["beta"]][["base"]]					<- theta[i,"base"]								
					beta.stratified							<- acute.get.rates(ibm[["beta"]], ibm.pop= NULL, pop.n=pop.n, state.n= as.matrix(state.n), per.capita.i= 1)					
					ans										<- acute.loglkl(tpc.table.all.median, beta.stratified, cluster.tw, clu.n=tipc.n)
					ans[["table.lkl"]]
				})
		names(loglkl)	<- apply(theta,1,function(x) paste(x,collapse='_'))
		print(loglkl)
		stop()
		
		dir.name<- paste(DATA,"acutelkl",sep='/')
		f.name<- paste(dir.name,paste("tpc_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta0[1],"_b",theta0[2],"_s",m.known.states,sep=''),sep='/')		
		cat(paste("\nsave ",paste(f.name,"_tpc_lkl.R",sep='')))
		save(theta,lkl, file= paste(f.name,"_tpc_lkl.R",sep=''))
	}
	theta<- sapply( strsplit(colnames(lkl),'_',fixed=1), as.numeric )
	theta.acute<- unique(theta[1,])
	theta.base<- unique(theta[2,])		
	#note: theta.acute runs first in lkls
	lkls.mean<- apply(lkl,2,function(x) x[40] )#4 for 5000	 #39 for 1e4
	#lkls.mean<- apply(lkl,2,mean)
	#lkls.mean<- apply(lkls,2,function(x) sd(x)/mean(x))
	lkls.mean<- matrix(lkls.mean, nrow= length(theta.acute), ncol= length(theta.base), dimnames=list(theta.acute,theta.base))
	cat(paste("\nmax lkl is ",max(lkls.mean)))
	
	f.name<- paste(dir.name,paste("tpc_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta0[1],"_b",theta0[2],"_s",m.known.states,sep=''),sep='/')
	cat(paste("\nplot",paste(f.name,"_2D.pdf",sep='')))
	pdf(paste(f.name,"_2D.pdf",sep=''),version="1.4",width=4,height=4)
	my.image(theta.acute,theta.base,lkls.mean, xlab=expression(r[a]),ylab=expression(beta[0]))	
	points(theta0[1],theta0[2], pch=19, col="red")
	dev.off()
	
	cat(paste("\nplot",paste(f.name,"_3D.pdf",sep='')))
	pdf(paste(f.name,"_3D.pdf",sep=''),version="1.4",width=5,height=5)	
	out<- my.plot.persplocfit(theta.acute,theta.base,exp(lkls.mean), theta= 30, phi= 10, palette= "gray", xlab="acute multiplier",ylab="base",zlab="prob")
	dev.off()
	
	#compute log likelihood ratio
	cat(paste("\nplot",paste(f.name,"_lklr.pdf",sep='')))
	pdf(paste(f.name,"_lklr.pdf",sep=''),version="1.4",width=5,height=5)
	loglkl.ratio<- -2*apply(lkls.mean,1,max)
	loglkl.ratio<- loglkl.ratio[as.character(theta0[1])] - loglkl.ratio
	plot(as.numeric(names(loglkl.ratio)),loglkl.ratio,xlab=expression(r[a]),ylab="log lkl ratio", type='l')
	abline(v=theta0[1],lty=2,col="blue")
	dev.off()
	cat(print.v(loglkl.ratio, print.char=0, as.R=1))
	
	#print(lkls.mean)
	
}
###############################################################################
prj.wh.sleeper<- function()
{
	require(ape)
	d.name<- "/Users/Oliver/duke/2012_HIVLythgoe"
	f.name<- paste(d.name,"g26312_gp41_njtree.nwk",sep='/')
	
	ph<- read.tree(f.name)
	plot(ph)
}
###############################################################################
prj.pipeline<- function()
{
	if(0)	#simulate debug tip cluster data sets
	{
		dir.name			<- CODE.HOME
		#acute				<- c(2,4,6,8)
		#base				<- c(0.065, 0.058, 0.053, 0.05)
		acute				<- c(1.7,	11,		2.7,	15)					#first two: attack 0.014 last two: attack 0.01
		base				<- c(0.06,	0.03, 	0.042, 	0.03)				#two in row: E2E 0.1 and 0.4
		sIdx				<- 0.2
		sE					<- 0.2
		cluster.tw			<- 3
		debug.susc.const	<- 1
		debug.only.u		<- 1
		cmd			<-	sapply(seq_along(acute),function(i)
							{
								cmd			<- prj.simudata.cmd(dir.name, "Town II", acute[i], base[i], 50, sIdx, sE, cluster.tw, save=1, debug.susc.const=debug.susc.const, debug.only.u=debug.only.u)
								cmd			<- prj.hpcwrapper(cmd, hpc.walltime=8, hpc.mem="1600mb", hpc.load="module load R/2.15",hpc.nproc=1, hpc.q=NA)
								cat(cmd)								
								signat		<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
								outdir		<- paste(CODE.HOME,"misc",sep='/')
								outfile		<- paste("phd",signat,"qsub",sep='.')
								prj.hpccaller(outdir, outfile, cmd)
							})
		#cmd			<- paste(cmd,sep='',collapse='')		
	}
	if(0)	#simulate tip cluster data sets
	{
		dir.name			<- CODE.HOME
		#acute				<- c(2,4,6,8)
		#base				<- c(0.065, 0.058, 0.053, 0.05)
		acute				<- c(2.2,	15 )
		base				<- c(0.052,	0.032 )
		sIdx				<- 0.2
		sE					<- 0.2
		cluster.tw			<- 3
		debug.susc.const	<- 0
		debug.only.u		<- 0
		cmd			<-	sapply(seq_along(acute),function(i)
				{
					cmd			<- prj.simudata.cmd(dir.name, "Town II", acute[i], base[i], 50, sIdx, sE, cluster.tw, save=1, debug.susc.const=debug.susc.const, debug.only.u=debug.only.u)
					cmd			<- prj.hpcwrapper(cmd, hpc.walltime=8, hpc.mem="1600mb", hpc.load="module load R/2.15",hpc.nproc=1, hpc.q=NA)
					cat(cmd)
					signat		<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
					outdir		<- paste(CODE.HOME,"misc",sep='/')
					outfile		<- paste("phd",signat,"qsub",sep='.')
					prj.hpccaller(outdir, outfile, cmd)
				})
		#cmd			<- paste(cmd,sep='',collapse='')		
	}
	if(1)	#start auxiliary job
	{
		dir.name	<- CODE.HOME
		cmd			<- paste("\n",dir.name,"/misc/phdes.startme.R -exeSIMU.MATCH",sep='')
		cmd			<- prj.hpcwrapper(cmd, hpc.walltime=71, hpc.mem="1600mb", hpc.load="module load R/2.15",hpc.nproc=1, hpc.q="pqeph")
		cat(cmd)
		signat		<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
		outdir		<- paste(CODE.HOME,"misc",sep='/')
		outfile		<- paste("phd",signat,"qsub",sep='.')
		prj.hpccaller(outdir, outfile, cmd)
	}
	if(0)
	{
		dir.name	<- CODE.HOME
		acute		<- 1
		base		<- c(rep(0.07,5), rep(0.05,5))
		sIdx		<- c(seq(0.2,1,0.2),seq(0.2,1,0.2))
		sE			<- sIdx
		cluster.tw	<- 3
		cmd			<-	sapply(seq_along(base),function(i)
				{
					cmd			<- prj.simudata.cmd(dir.name, "Town II", acute, base[i], 50, sIdx[i], sE[i], cluster.tw, save=1, debug.susc.const=debug.susc.const, debug.only.u=debug.only.u)
					cmd			<- prj.hpcwrapper(cmd, hpc.walltime=8, hpc.mem="1600mb", hpc.load="module load R/2.15",hpc.nproc=1, hpc.q=NA)
					cat(cmd)
					
					signat		<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
					outdir		<- paste(CODE.HOME,"misc",sep='/')
					outfile		<- paste("phd",signat,"qsub",sep='.')
					prj.hpccaller(outdir, outfile, cmd)
				})
		#cmd			<- paste(cmd,sep='',collapse='')		
	}
}
###############################################################################
prj.hpcwrapper<- function(cmd, hpc.walltime=3, hpc.mem="400mb", hpc.load="",hpc.nproc=1, hpc.q=NA)
{
	wrap<- "#!/bin/sh"	
	if(1)
	{				
		tmp	<- paste("#PBS -l walltime=",hpc.walltime,":59:59,pcput=",hpc.walltime,":45:00",sep='')
		wrap<- paste(wrap, tmp, sep='\n')		
		tmp	<- paste("#PBS -l select=1:ncpus=",hpc.nproc,":mem=",hpc.mem,sep='')
		wrap<- paste(wrap, tmp, sep='\n')
		wrap<- paste(wrap, "#PBS -j oe", sep='\n')
		if(!is.na(hpc.q))
			wrap<- paste(wrap, paste("#PBS -q",hpc.q), sep='\n\n')
		wrap<- paste(wrap, hpc.load, sep='\n')
		#tmp	<- paste("export PATH=$PATH:",HPC.BIN,sep='')
		#wrap<- paste(wrap, tmp, sep='\n')
		
	}
	else if(tmp=='')
		cat(paste("\ndetected no HPC system and no hpcwrapper generated, domain name is",tmp))
	else
		stop(paste("unknown hpc system with domain name",tmp))
	
	cmd<- lapply(seq_along(cmd),function(i){	paste(wrap,cmd[[i]],sep='\n')	})
	if(length(cmd)==1)
		cmd<- unlist(cmd)
	cmd	
}
###############################################################################
prj.hpccaller<- function(outdir, outfile, cmd)
{
	file<- paste(outdir,outfile,sep='/')
	cat(paste("\nwrite cmd to",file,"\n"))
	cat(cmd,file=file)
	cmd<- paste("qsub",file)
	cat( cmd )
	cat( system(cmd, intern=TRUE) )	
	Sys.sleep(1)
}
###############################################################################
prj.simudata.cmd<- function(dir.name, loc, acute, base, rep, sIdx, sE, cluster.tw, debug.susc.const=0, debug.only.u=0, save=1, resume=1, verbose=1)
{		
	cmd<- paste("\n",dir.name,"/misc/phdes.startme.R -exeSIMU.DATA ",sep='')
	#cmd<- paste(cmd, " -l=",loc,sep='')
	cmd<- paste(cmd, " -acute=",acute,sep='')
	cmd<- paste(cmd, " -baseline=",base,sep='')
	cmd<- paste(cmd, " -rep=",rep,sep='')
	cmd<- paste(cmd, " -sIdx=",sIdx,sep='')
	cmd<- paste(cmd, " -sE=",sE,sep='')
	cmd<- paste(cmd, " -cluster.tw=",cluster.tw,sep='')
	cmd<- paste(cmd, " -save=",save,sep='')
	cmd<- paste(cmd, " -resume=",resume,sep='')
	cmd<- paste(cmd, " -debug.only.u=",debug.only.u,sep='')
	cmd<- paste(cmd, " -debug.susc.const=",debug.susc.const,sep='')	
	cmd<- paste(cmd, " -v=",verbose,sep='')
	cmd
}
###############################################################################
prj.simudata<- function()
{		
	require(phylodesign)
	#call with eg pkg/misc/phdes.startme.R -exeSIMU.DATA -v1 -acute=8 -baseline=0.05 -r=1 -sIdx=0.5 -sE=0.5 -cluster.tw=3
	m.type				<- "Acute"
	loc.type			<- "Town II"
	m.popsize			<- NA
	resume				<- 1
	verbose				<- 0
	record.tpc			<- 1
	tpc.repeat			<- 2
	cluster.tw			<- 3
	save				<- 1
	#theta				<- c(8, 0.05, 0, 0)
	theta				<- c(2, 0.065, 0, 0)
	names(theta)		<- c("acute","base","m.st1","m.st2")
	sample.prob			<- c(0.6,0.6)
	names(sample.prob)	<- c("Idx","E")
	debug.susc.const	<- 1
	debug.only.u		<- 1
	dir.name			<- paste("acutesimu_fxs",debug.susc.const,"_onlyu",debug.only.u,sep='')

	if(exists("args"))
	{
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,4),
									rep= return(as.numeric(substr(arg,6,nchar(arg)))),NA)	}))
		if(length(tmp)>0) tpc.repeat<- tmp[1]
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,2),
									l= return(substr(arg,4,nchar(arg))),NA)	}))
		if(length(tmp)>0) loc.type<- tmp[1]
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,2),
									n= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) m.popsize<- tmp[1]
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,11),
									cluster.tw= return(as.numeric(substr(arg,13,nchar(arg)))),NA)	}))
		if(length(tmp)>0) cluster.tw<- tmp[1]				
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,6),
									acute= return(as.numeric(substr(arg,8,nchar(arg)))),NA)	}))
		if(length(tmp)>0) theta[1]<- tmp[1]
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,9),
									baseline= return(as.numeric(substr(arg,11,nchar(arg)))),NA)	}))
		if(length(tmp)>0) theta[2]<- tmp[1]
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,5),
									sIdx= return(as.numeric(substr(arg,7,nchar(arg)))),NA)	}))
		if(length(tmp)>0) sample.prob[1]<- tmp[1]
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,3),
									sE= return(as.numeric(substr(arg,5,nchar(arg)))),NA)	}))
		if(length(tmp)>0) sample.prob[2]<- tmp[1]				
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,5),
									save= return(as.numeric(substr(arg,7,nchar(arg)))),NA)	}))
		if(length(tmp)>0) save<- tmp[1]
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,8),
									resume= return(as.numeric(substr(arg,10,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,17),
									debug.susc.const= return(as.numeric(substr(arg,19,nchar(arg)))),NA)	}))
		if(length(tmp)>0) debug.susc.const<- tmp[1]
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,13),
									debug.only.u= return(as.numeric(substr(arg,15,nchar(arg)))),NA)	}))
		if(length(tmp)>0) debug.only.u<- tmp[1]
		
		dir.name			<- paste("acutesimu_fxs",debug.susc.const,"_onlyu",debug.only.u,sep='')
	}
	if(verbose)
	{
		cat(paste("\nm.popsize ",m.popsize))
		cat(paste("\nloc.type ",loc.type))
		cat(paste("\ncluster.tw ",cluster.tw))	
		cat(paste("\ntpc.repeat ",tpc.repeat))
		cat(paste("\ntheta0.acute ",theta[1]))
		cat(paste("\ntheta0.base ",theta[2]))
		cat(paste("\nsample.prob.Idx ",sample.prob[1]))
		cat(paste("\nsample.prob.E ",sample.prob[2]))
		cat(paste("\ndebug.susc.const ",debug.susc.const))
		cat(paste("\ndebug.only.u ",debug.only.u))
		cat(paste("\nverbose ",verbose))
		cat(paste("\nresume ",resume))
		cat(paste("\nsave ",save))
	}
	my.mkdir(DATA,dir.name)
	dir.name	<- paste(DATA,dir.name,sep='/')
	f.name		<- paste(dir.name,paste("tpcdat_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta[1],"_b",theta[2],"_sIdx",sample.prob[1],"_sE",sample.prob[2],"_tw",cluster.tw,sep=''),sep='/')
	if(resume)
	{
		options(show.error.messages = FALSE)		
		if(verbose)
			cat(paste("\nprj.simudata: try to resume file ",paste(f.name,".R",sep='')))
		readAttempt<-	try(suppressWarnings(load(paste(f.name,".R",sep=''))))
		options(show.error.messages = TRUE)
		if(!inherits(readAttempt, "try-error") && verbose)
			cat(paste("\nprj.simudata: resumed file ",paste(f.name,".R",sep='')))
	}
	if(!resume || inherits(readAttempt, "try-error"))
	{			
		tpc<- lapply(seq_len(tpc.repeat),function(i)
			{
				if(verbose)	
					cat(paste("\nprj.simudata: replicate",i))
				ans			<- vector("list",4)
				names(ans)	<- c("ibm","tpc.internal","tpc.table.all","tpc.table.sample")				
				ibm			<- ibm.init.model( m.type, loc.type, m.popsize, theta, resume= 0, debug=debug.only.u )
				ans[["ibm"]]<- ibm.collapse( ibm )
				tpc.data	<- ssa.run(ans[["ibm"]], ssa.ctime= 0, ssa.etime= cluster.tw,record.tpc= record.tpc, verbose= 0, resume= 0, debug=debug.susc.const)
				if(record.tpc)
					ans[["tpc.internal"]]	<- tpc.collapse(tpc.data)
				else
					ans[["tpc.internal"]]	<- tpc.data				
				ans[["tpc.table.all"]]		<- tpc.tabulate( ans[["tpc.internal"]] )								
				ans[["tpc.table.sample"]]	<- tpc.tabulate( tpc.sample( ans[["tpc.internal"]], sample.prob ) )				
				ans
			})		
		if(save)
		{
			cat(paste("\nprj.simudata: write tpc data to file",paste(f.name,".R",sep='')))
			save(tpc,file=paste(f.name,".R",sep=''))	
		}		
	}
	if(verbose)
	{
		sum.attack	<- summary( sapply(seq_along(tpc), function(i) tpc[[i]][["tpc.internal"]][["attack.rate"]]) )	
		sum.E2E		<- summary( sapply(seq_along(tpc), function(i)
						{
							tmp	<- tpc.proportion.E2E(tpc[[i]][["tpc.internal"]])
							tmp["E2E"] / tmp["X2E"]
						}) )
		print(sum.attack)
		print(sum.E2E)
	}
	if(save)
	{
		table.name				<- "tpc.table.all"
		max.ntr					<- max(sapply(seq_along(tpc), function(i) ncol(tpc[[i]][[table.name]]) ))
		tpc.table				<- sapply(seq_along(tpc), function(i)
				{
					c( as.vector(tpc[[i]][[table.name]]), rep(0, nrow(tpc[[i]][[table.name]]) * (max.ntr - ncol(tpc[[i]][[table.name]]))) )			
				})
		tpc.table				<- matrix( round(apply(tpc.table, 1, median )), ncol=max.ntr )
		tpc.table				<- tpc.table[ , apply(tpc.table,2,function(x)  any(x!=0) ) ]		
		dimnames(tpc.table)		<- list(rownames(tpc[[1]][[table.name]]), paste("n",seq.int(0,ncol(tpc.table)-1),sep=''))
		tpc.table.all.median	<- tpc.table
		print(tpc.table.all.median)
		
		table.name				<- "tpc.table.sample"
		max.ntr					<- max(sapply(seq_along(tpc), function(i) ncol(tpc[[i]][[table.name]]) ))
		tpc.table				<- sapply(seq_along(tpc), function(i)
				{
					c( as.vector(tpc[[i]][[table.name]]), rep(0, nrow(tpc[[i]][[table.name]]) * (max.ntr - ncol(tpc[[i]][[table.name]]))) )			
				})
		tpc.table				<- matrix( round(apply(tpc.table, 1, median )), ncol=max.ntr )
		tpc.table				<- tpc.table[ , apply(tpc.table,2,function(x)  any(x!=0) ) ]		
		dimnames(tpc.table)		<- list(rownames(tpc[[1]][[table.name]]), paste("ns",seq.int(0,ncol(tpc.table)-1),sep=''))
		tpc.table.sample.median	<- tpc.table
		print(tpc.table.sample.median)
		
		f.name<- paste(dir.name,paste("tpcdat_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta[1],"_b",theta[2],"_sIdx",sample.prob[1],"_sE",sample.prob[2],"_tw",cluster.tw,"_median",sep=''),sep='/')
		if(verbose)
			cat(paste("\nprj.simudata: write tpc data medians to file",paste(f.name,".R",sep='')))
		save(tpc.table.all.median,tpc.table.sample.median,sum.attack,sum.E2E,file=paste(f.name,".R",sep=''))				
	}
	tpc	
}
###############################################################################
prj.simudata.match.theta.to.Inc.E2E<- function()
{
	require(data.table)
	#require(multicore)
		
	m.type				<- "Acute"
	loc.type			<- "Town II"
	m.popsize			<- NA
	resume				<- 1
	verbose				<- 1	
	cluster.tw			<- 3
	save				<- 1
	prior.acute			<- c(1,20)
	prior.base			<- c(0.015, 0.09)		
	sample.prob			<- 1	
	debug.susc.const	<- 0
	debug.only.u		<- 0
	target.su.INC		<- 0.01
	target.su.E2E		<- 0.1
	abc.nit				<- 5e4
	abc.cores			<- 8
	dir.name			<- paste("acutesimu_fxs",debug.susc.const,"_onlyu",debug.only.u,sep='')
	
	f.name		<- paste(DATA,'/',dir.name,'/',"match_INC_E2E",'_',m.type,'_',loc.type,'_',cluster.tw,"_acute_",prior.acute[1],'_',prior.acute[2],"_base_",prior.base[1],'_',prior.base[2],".R",sep='')
	if(resume)
	{
		options(show.error.messages = FALSE)		
		if(verbose)
			cat(paste("\nprj.simudata.match.theta.to.Inc.E2E: try to resume file ",f.name))
		readAttempt<-	try(suppressWarnings(load(f.name)))
		options(show.error.messages = TRUE)
		if(!inherits(readAttempt, "try-error") && verbose)
			cat(paste("\nprj.simudata.match.theta.to.Inc.E2E: resumed file ",f.name))
	}
	if(!resume || inherits(readAttempt, "try-error"))
	{	
		abc.struct	<- c( runif(abc.nit,prior.acute[1],prior.acute[2]), runif(abc.nit,prior.base[1],prior.base[2]) )
		abc.struct	<- as.data.table( matrix(abc.struct, ncol=2, nrow= abc.nit, dimnames=list(c(),c("acute","base"))) )		
		#print(abc.struct)
		#abc.struct	<- mclapply(seq_len(nrow(abc.struct)),function(i)
		abc.struct	<- lapply(seq_len(nrow(abc.struct)),function(i)
				{		
					print(unlist(c(i,abc.struct[i,])))
					args<<- prj.simudata.cmd(CODE.HOME, loc.type, abc.struct[i,acute], abc.struct[i,base], 1, sample.prob, sample.prob, cluster.tw, save=0, resume=0, verbose=0, debug.susc.const=debug.susc.const, debug.only.u=debug.only.u)
					args<<- unlist(strsplit(args,' '))
					tpc	<- prj.simudata()
					sum.attack	<- median( sapply(seq_along(tpc), function(i) tpc[[i]][["tpc.internal"]][["attack.rate"]]) )					
					sum.E2E		<- median( sapply(seq_along(tpc), function(i)
									{
										tmp	<- tpc.proportion.E2E(tpc[[i]][["tpc.internal"]])
										tmp["E2E"] / tmp["X2E"]
									}) )
					ans			<- c(abc.struct[i,acute], abc.struct[i,base], sum.attack, sum.E2E)				
					ans
				})
		#}, mc.cores= abc.cores)
		#print(unlist(abc.struct))
		abc.struct	<- as.data.table( matrix(unlist(abc.struct),byrow=1,ncol=4,nrow=abc.nit, dimnames=list(c(),c("acute","base","INC","E2E"))) )
		cat(paste("\nprj.simudata.match.theta.to.Inc.E2E: write data to file",f.name))
		save(abc.struct, file=f.name)		
	}	
	#print(tmp)	
}
###############################################################################

