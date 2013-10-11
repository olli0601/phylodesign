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
#simulate high acute and low acute tip clusters for each community	
prj.popart.powercalc.by.acutelklratio.tpcobs<- function(theta.EE.H0, theta.EE.H1, cohort.dur, p.consent.coh, p.consent.clu, p.lab, p.vhcc.prev.AB, p.vhcc.inc.AB, p.vhcc.prev.C, p.vhcc.inc.C, opt.sampling, pooled.n, dir.name=DATA, verbose=1, resume=1, standalone=0)
{
	m.type			<- "Acute"	
	theta.model.Hx	<- NULL
	f.name			<- paste(dir.name,'/',"tpcobs_",m.type,'_',theta.EE.H0,'_',theta.EE.H1,'_',opt.sampling,'_',"central",'_',p.lab,'_',p.consent.coh,sep='')	
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
		#get sites and add target effects to 'sites'
		sites						<- popart.getdata.randomized.arm( pooled.n, rtn.fixed=debug, rtn.phylostudy=1 )
		#TODO this would change for each arm and depending on pess, central and opt target
		sites[,"mu.inc.rate.H0"]	<- sites[,"inc.rate"]
		sites[,"mu.inc.rate.H1"]	<- 0.013
		sites[,"mu.pE2E.H0"]		<- theta.EE.H0
		sites[,"mu.pE2E.H1"]		<- theta.EE.H1
		#print(sites)
		samples.CD4					<- popart.predicted.firstCD4()
		samples.seq					<- popart.predicted.sequences(sites,  samples.CD4, p.consent.coh, p.consent.clu, p.lab, p.vhcc.prev.AB, p.vhcc.inc.AB, p.vhcc.prev.C, p.vhcc.inc.C, method= opt.sampling)
		#print(samples.seq)
		#print(apply(samples.seq,2,sum))
		samples.seq					<- cbind( samples.seq, samples.seq[,"%avg",drop=0] * 0.1 / 2 )		#add simple prior on seq.cov
		colnames(samples.seq)[ncol(samples.seq)]	<- c("sigma")
		sites						<- cbind(sites, samples.CD4, samples.seq)
		
		#each site was calibrated to match 		"mu.inc.rate.H0","mu.pE2E.H0"		"mu.inc.rate.H1","mu.pE2E.H1"
		if(verbose)	cat(paste("\nset up default theta corresponding to %E2E 10 vs 40 -- should only use these for arm C as these theta don t reflect the induced %E2E in A,B"))
		theta.model.Hx	<- matrix(	c(		1.25, 0.062, 4.4, 0.09,
											1.1, 0.088, 5.8, 0.073,
											1, 0.1, 7, 0.065,
											1.8, 0.045, 5.8, 0.069,
											1.6, 0.063, 6, 0.062,
											0.7, 0.125, 5.6, 0.08,
											3.8, 0.027, 10.8, 0.045,
											2.2, 0.053, 11, 0.045,
											1.7, 0.067, 12, 0.043,
											2, 0.048, 5.9, 0.071,
											0.7, 0.125, 3.8, 0.099,
											1, 0.09, 7.4, 0.059		),byrow=T, ncol=4,nrow=nrow(sites), dimnames=list(c(),c("acute.H0","base.H0","acute.H1","base.H1")))
		theta.model.Hx	<- as.data.table(theta.model.Hx)
		theta.model.Hx[,comid_old:=sites[,"comid_old"]]		
			
		if(1)		#take matching parameters from 95% cloud around target parameters
		{						
			if(verbose)	cat(paste("\nfind theta corresponding to incidence and %E2E under H0 and H1 "))
			theta.model.Hx	<- lapply(seq_len(nrow(sites)),function(i)
							{					
								site		<- sites[i,"comid_old"]
								if(verbose)	cat(paste("\nprocess ",site))
								f.name		<- paste(DATA,'/',paste("acutesimu_fxs0_onlyu0",sep=''),'/',"match_INC_E2E",'_',m.type,'_',site,'_',cohort.dur,"_acute_0.5_13_base_0.01_0.13.R",sep='')								
								readAttempt	<-	try(suppressWarnings(load(f.name)))
								if(!inherits(readAttempt, "try-error"))	
								{	
									if(verbose)	cat(paste("\nloaded precomputed simulations from file",f.name))
									setnames(abc.struct, "INC", "Inc")
									abc.struct[,comid_old:=site]
									tmp			<- prj.popart.powercalc.by.acutelklratio.matchpa( abc.struct, sites[,c("comid_old","mu.inc.rate.H0","mu.inc.rate.H1","mu.pE2E.H0","mu.pE2E.H1")], hetclu.scale=1, return.top=25 )						
#									print(tmp)
									if(verbose)	cat(paste("\nfound params matching H0, n with nonzero weights=",nrow(subset(tmp[["H0"]],weight>0))))
									if(verbose)	cat(paste("\nfound params matching H1, n with nonzero weights=",nrow(subset(tmp[["H1"]],weight>0))))
									#take weighted mean of matching params					
									ans			<- c(	acute.H0=weighted.mean(tmp[["H0"]][, acute], tmp[["H0"]][, weight]), base.H0=weighted.mean(tmp[["H0"]][, base], tmp[["H0"]][, weight]), 
														acute.H1=weighted.mean(tmp[["H1"]][, acute], tmp[["H1"]][, weight]), base.H1=weighted.mean(tmp[["H1"]][, base], tmp[["H1"]][, weight])	)
									#use default param if not sensible			
									if(nrow(subset(tmp[["H0"]],weight>0))<1)
									{
										if(verbose)	cat(paste("\ncould not find params matching H0, use default for arm C"))
										ans[1:2]<- as.numeric(subset(theta.model.Hx, comid_old==site,select=c(acute.H0,base.H0)))
									}
									if(nrow(subset(tmp[["H1"]],weight>0))<1)
									{
										if(verbose)	cat(paste("\ncould not find params matching H1, use default for arm C"))
										ans[3:4]<- as.numeric(subset(theta.model.Hx, comid_old==site,select=c(acute.H1,base.H1)))
									}
								}
								else
								{
									options(warn=1)
									warning("could not load file",f.name,"\nusing inappropriate default parameters")
									options(warn=2)
									ans			<- c(acute.H0=2.2, base.H0=0.052, acute.H1=15, base.H1=0.032)
								}									
								data.table(comid_old=site,acute.H0=ans[1], base.H0=ans[2], acute.H1=ans[3], base.H1=ans[4])								
							})
			theta.model.Hx	<- rbindlist(theta.model.Hx)
			if(verbose)	cat(paste("\nfind theta corresponding to incidence and %E2E under H0 and H1"))
			print(theta.model.Hx)			
		}
		#have data table with model parameters for each site to simulate from
		#	-- take subset from data table and simulate 50				
		tpc.obs			<- lapply(seq_len(nrow(sites)),function(i)
				{					
					if(verbose)
						cat(paste("\nprocess ",sites[i,"comid_old"]))
					tmp	<- subset(theta.model.Hx, comid_old==sites[i,"comid_old"])
					if(!standalone)
					{
						args		<<- prj.simudata.cmd(CODE.HOME, sites[i,"comid_old"], tmp[1,acute.H0], tmp[1,base.H0], 50, round(sites[i,"%avg"],d=2), round(sites[i,"%avg"],d=2), cohort.dur, save=1, resume=1, verbose=1, debug.susc.const=0, debug.only.u=0)
						args		<<- unlist(strsplit(args,' '))		#print(args)					
						tpc			<- prj.simudata()					#print(tpc[["sum.attack"]]["Median"]); print(tpc[["sum.E2E"]]["Median"])										
						tpc.H0		<- tpc[["tpc.table.sample.median"]]
						
						args		<<- prj.simudata.cmd(CODE.HOME, sites[i,"comid_old"], tmp[1,acute.H1], tmp[1,base.H1], 50, round(sites[i,"%avg"],d=2), round(sites[i,"%avg"],d=2), cohort.dur, save=1, resume=1, verbose=1, debug.susc.const=0, debug.only.u=0)
						args		<<- unlist(strsplit(args,' '))					
						tpc			<- prj.simudata()
						tpc.H1		<- tpc[["tpc.table.sample.median"]]
						ans			<- list(H0=tpc.H0, H1=tpc.H1)
					}
					else	#precompute low and high acute scenarios
					{
						#low acute
						cmd			<- prj.simudata.cmd(CODE.HOME, sites[i,"comid_old"], tmp[1,acute.H0], tmp[1,base.H0], 50, round(sites[i,"%avg"],d=2), round(sites[i,"%avg"],d=2), cohort.dur, save=1, resume=1, verbose=1, debug.susc.const=0, debug.only.u=0)
						cmd			<- prj.hpcwrapper(cmd, hpc.walltime=8, hpc.mem="1600mb", hpc.load="module load R/2.15",hpc.nproc=1, hpc.q="pqeph")
						cat(cmd)								
						signat		<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
						outdir		<- paste(CODE.HOME,"misc",sep='/')
						outfile		<- paste("phd_ct0",signat,"qsub",sep='.')
						prj.hpccaller(outdir, outfile, cmd)
						#high acute
						cmd			<- prj.simudata.cmd(CODE.HOME, sites[i,"comid_old"], tmp[1,acute.H1], tmp[1,base.H1], 50, round(sites[i,"%avg"],d=2), round(sites[i,"%avg"],d=2), cohort.dur, save=1, resume=1, verbose=1, debug.susc.const=0, debug.only.u=0)
						cmd			<- prj.hpcwrapper(cmd, hpc.walltime=8, hpc.mem="1600mb", hpc.load="module load R/2.15",hpc.nproc=1, hpc.q="pqeph")
						cat(cmd)								
						signat		<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
						outdir		<- paste(CODE.HOME,"misc",sep='/')
						outfile		<- paste("phd_ct1",signat,"qsub",sep='.')
						prj.hpccaller(outdir, outfile, cmd)
						ans			<- NULL
					}	
					ans
				})
		if(!standalone)
		{
			names(tpc.obs)	<- sites[,"comid_old"]		
			f.name			<- paste(dir.name,'/',"tpcobs_",m.type,'_',theta.EE.H0,'_',theta.EE.H1,'_',opt.sampling,'_',"central",'_',p.lab,'_',p.consent.coh,sep='')
			if(verbose)	cat(paste("\nwrite tpc.obs to",paste(f.name,".R",sep='')))	
			save(tpc.obs, theta.model.Hx, sites, file=paste(f.name,".R",sep=''))
		}
	}
	
	if(!standalone)
		ans	<- list(tpc.obs=tpc.obs, theta.model.Hx=theta.model.Hx, sites=sites)
	else
		ans	<- NULL
	ans
}	
###############################################################################
#	mlkl.n= 1e3; cohort.dur=3; replace<- 1; remote<- 0;
prj.popart.powercalc.by.acutelklratio.lkl4Precomputed<- function(sites=NULL, tpc.obs=NULL, cohort.dur=3, f.name=NA, replace=1, resume=1, verbose=1, remote=0, remote.signat=NA)
{	
	m.type			<- "Acute"
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
		#
		#	load precomputed data tables for each site and add 'sample' from prior determined by '%avg' and 'sigma'
		#
		tmp			<- seq_len(nrow(sites))[-which(sites$comid_old=="Ikhwezi")]
		lkl.theta	<- rbindlist( lapply( tmp, function(i)
						{					
							site		<- sites[i,"comid_old"]
							if(verbose)	cat(paste("\nprocess ",site))
							file		<- paste(DATA,'/',paste("acutesimu_fxs0_onlyu0",sep=''),'/',"match_INC_E2E",'_',m.type,'_',site,'_',cohort.dur,"_acute_0.5_13_base_0.01_0.13.R",sep='')								
							readAttempt	<-	try(suppressWarnings(load(file)))
							if(!inherits(readAttempt, "try-error"))	
							{	
								if(verbose)	cat(paste("\nloaded precomputed simulations from file",file))
								setnames(abc.struct, "INC", "Inc")
								abc.struct[,comid_old:=site]						
							}
							else
							{						
								stop("could not load file",file,"\nusing inappropriate default parameters")					
							}
							abc.struct[,sample:= rnorm(nrow(abc.struct),mean=sites[i,"%avg"], sd=sites[i,"sigma"])]
							abc.struct								
						}) )
		if(remote)
		{
			#save job specifics to temporary file		
			f.name.remote	<- paste(f.name,'_tmp_',remote.signat,sep='')
			if(verbose)	cat(paste("\nremote job. saving lkl.theta to file", f.name.remote))			
			save(lkl.theta, sites,  tpc.obs, cohort.dur, file=paste(f.name.remote,".R",sep=''))
			
			#call acute.loglkl.batch remotely
			dummy			<- lapply(seq_len(nrow(sites)),function(i)
					{		
						#print( mlkl.theta[[i]] )
						#H0.H0
						cmd			<- prog.acute.loglkl.batch.cmd(CODE.HOME, paste(f.name.remote,".R",sep=''),	paste(f.name.remote,'_',sites[i,"comid_old"],"_H0",".R",sep=''),	sites[i,"comid_old"],	"H0")
						cat(cmd)
						cmd			<- prj.hpcwrapper(cmd, hpc.walltime=8, hpc.mem="1600mb", hpc.load="module load R/2.15",hpc.nproc=1, hpc.q="pqeph")
						prj.hpccaller(paste(CODE.HOME,"misc",sep='/'), paste("phd_l0",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),"qsub",sep='.'), cmd)
						#H0.H1
						cmd			<- prog.acute.loglkl.batch.cmd(CODE.HOME, paste(f.name.remote,".R",sep=''),	paste(f.name.remote,'_',sites[i,"comid_old"],"_H1",".R",sep=''),	sites[i,"comid_old"],	"H1")
						cat(cmd)
						cmd			<- prj.hpcwrapper(cmd, hpc.walltime=8, hpc.mem="1600mb", hpc.load="module load R/2.15",hpc.nproc=1, hpc.q="pqeph")
						prj.hpccaller(paste(CODE.HOME,"misc",sep='/'), paste("phd_l1",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),"qsub",sep='.'), cmd)						
					})
		}
		else
		{	
			lkl.theta.comid		<- seq_len(nrow(sites))[-which(sites$comid_old=="Ikhwezi")]			
			lkl.theta			<- rbindlist( lapply(lkl.theta.comid[1:2],function(i)
					{		
						#print( mlkl.theta[[i]] )
						site		<- sites[i,"comid_old"]
						tmp			<- subset(lkl.theta, comid_old==site)
						if(verbose)	cat(paste("\nacute.loglkl.batch, process ",site,"\nnumber lkl evaluations: ",nrow(tmp)))
						#tmp		<- tmp[1:10,]
						ans			<- acute.loglkl.batch(site, tpc.obs[[i]][["H0"]], cohort.dur, tmp, clu.closure= 12, verbose=verbose)
						ans[, h:="H0"]
						tmp		<- acute.loglkl.batch(site, tpc.obs[[i]][["H1"]], cohort.dur, tmp, clu.closure= 12, verbose=verbose)
						tmp[, h:="H1"]
						ans			<- rbind(ans, tmp)
						#				print(ans); stop()
						ans
					}) )
			if(verbose)
				cat(paste("\nwrite tpc.lkl to",paste(f.name,".R",sep='')))	
			save(tpc.obs, sites, lkl.theta, file=paste(f.name,".R",sep=''))				
		}					
	}
	lkl.theta
}
###############################################################################
#	mlkl.n= 1e3; cohort.dur=3; replace<- 1; remote<- 0;
prj.popart.powercalc.by.acutelklratio.lklH0H1<- function(sites=NULL, tpc.obs=NULL, mlkl.theta.model=NULL, mlkl.n= 1e3, cohort.dur=3, f.name=NA, replace=1, resume=1, verbose=1, remote=0, remote.signat=NA)
{
	if(resume && remote && !is.null(sites) )
	{
		mlkl.theta	<- NULL
		mlkl.theta	<- for(site in sites[,"comid_old"])
						{					
							f.name.remote			<- paste(paste(f.name,'_tmp_',remote.signat,sep=''),'_',site,"_H0H0",".R",sep='')
							readAttempt				<- try(suppressWarnings(load(paste(f.name.remote,".R",sep=''))))
							if(!inherits(readAttempt, "try-error"))
								lkl.tpcH0.thetaH0	<- ans
							else	
								break								
							f.name.remote			<- paste(paste(f.name,'_tmp_',remote.signat,sep=''),'_',site,"_H0H1",".R",sep='')
							readAttempt				<- try(suppressWarnings(load(paste(f.name.remote,".R",sep=''))))
							if(!inherits(readAttempt, "try-error"))
								lkl.tpcH0.thetaH1	<- ans
							else	
								break
							f.name.remote			<- paste(paste(f.name,'_tmp_',remote.signat,sep=''),'_',site,"_H1H0",".R",sep='')
							readAttempt				<- try(suppressWarnings(load(paste(f.name.remote,".R",sep=''))))
							if(!inherits(readAttempt, "try-error"))
								lkl.tpcH1.thetaH0	<- ans
							else	
								break
							f.name.remote			<- paste(paste(f.name,'_tmp_',remote.signat,sep=''),'_',site,"_H1H1",".R",sep='')
							readAttempt				<- try(suppressWarnings(load(paste(f.name.remote,".R",sep=''))))
							if(!inherits(readAttempt, "try-error"))
								lkl.tpcH1.thetaH1	<- ans
							else	
								break
							mlkl.theta				<- c(mlkl.theta, list(tpcH0.thetaH0=lkl.tpcH0.thetaH0, tpcH0.thetaH1=lkl.tpcH0.thetaH1, tpcH1.thetaH0=lkl.tpcH1.thetaH0, tpcH1.thetaH1=lkl.tpcH1.thetaH1))							
						}
		#print(mlkl.theta)
		if(length(mlkl.theta)==nrow(sites))
		{
			if(verbose)
				cat(paste("\nwrite tpc.lkl to",paste(f.name,".R",sep='')))	
			save(tpc.obs, sites, mlkl.theta, file=paste(f.name,".R",sep=''))
		}
	}
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
		#get a table with param combinations for each site to sample from
		mlkl.theta			<- lapply(seq_len(nrow(sites)),function(i)
								{									
									tmp.H0	<- subset( mlkl.theta.model, comid_old==sites[i, "comid_old"] & h=="H0" )
									if(all(tmp.H0[,sample.sigma]==0))
										tmp.H0	<- tmp.H0[ sample.int(nrow(tmp.H0), min(nrow(tmp.H0), mlkl.n), replace=0), ]										 									
									else
										tmp.H0	<- tmp.H0[ sample.int(nrow(tmp.H0), mlkl.n, replace=1), ]
									tmp.H0[,sample:= rnorm( nrow(tmp.H0), tmp.H0[,sample.mu], tmp.H0[,sample.sigma] )]
									tmp.H1	<- subset( mlkl.theta.model, comid_old==sites[i, "comid_old"] & h=="H1" )
									if(all(tmp.H1[,sample.sigma]==0))
										tmp.H1	<- tmp.H1[ sample.int(nrow(tmp.H1), min(nrow(tmp.H0), mlkl.n), replace=0), ]
									else
										tmp.H1	<- tmp.H1[ sample.int(nrow(tmp.H1), mlkl.n, replace=1), ]
									tmp.H1[,sample:= rnorm( nrow(tmp.H1), tmp.H1[,sample.mu], tmp.H1[,sample.sigma] )]									
									list(	H0= subset(tmp.H0, select=c(acute,base,sample)), H1= subset(tmp.H1, select=c(acute,base,sample)) 	) 	
								})
		names(mlkl.theta)	<- sites[,"comid_old"]
		
		if(remote)
		{
			#save job specifics to temporary file		
			f.name.remote	<- paste(f.name,'_tmp_',remote.signat,sep='')			
			print(paste(f.name.remote,".R",sep=''))
			save(mlkl.theta, sites,  tpc.obs, cohort.dur, file=paste(f.name.remote,".R",sep=''))
			
			#call acute.loglkl.batch remotely
			mlkl.theta		<- lapply(seq_len(nrow(sites)),function(i)
								{		
									#print( mlkl.theta[[i]] )
									#H0.H0
									cmd			<- prog.acute.loglkl.batch.cmd(CODE.HOME, paste(f.name.remote,".R",sep=''),	paste(f.name.remote,'_',sites[i,"comid_old"],"_H0H0",".R",sep=''),	sites[i,"comid_old"],	"H0",	"H0")
									cat(cmd)
									cmd			<- prj.hpcwrapper(cmd, hpc.walltime=8, hpc.mem="1600mb", hpc.load="module load R/2.15",hpc.nproc=1, hpc.q="pqeph")
									prj.hpccaller(paste(CODE.HOME,"misc",sep='/'), paste("phd",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),"qsub",sep='.'), cmd)
									#H0.H1
									cmd			<- prog.acute.loglkl.batch.cmd(CODE.HOME, paste(f.name.remote,".R",sep=''),	paste(f.name.remote,'_',sites[i,"comid_old"],"_H0H1",".R",sep=''),	sites[i,"comid_old"],	"H0",	"H1")
									cat(cmd)
									cmd			<- prj.hpcwrapper(cmd, hpc.walltime=8, hpc.mem="1600mb", hpc.load="module load R/2.15",hpc.nproc=1, hpc.q="pqeph")
									prj.hpccaller(paste(CODE.HOME,"misc",sep='/'), paste("phd",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),"qsub",sep='.'), cmd)
									#H1.H0
									cmd			<- prog.acute.loglkl.batch.cmd(CODE.HOME, paste(f.name.remote,".R",sep=''),	paste(f.name.remote,'_',sites[i,"comid_old"],"_H1H0",".R",sep=''),	sites[i,"comid_old"],	"H1",	"H0")
									cat(cmd)
									cmd			<- prj.hpcwrapper(cmd, hpc.walltime=8, hpc.mem="1600mb", hpc.load="module load R/2.15",hpc.nproc=1, hpc.q="pqeph")
									prj.hpccaller(paste(CODE.HOME,"misc",sep='/'), paste("phd",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),"qsub",sep='.'), cmd)
									#H1.H1
									cmd			<- prog.acute.loglkl.batch.cmd(CODE.HOME, paste(f.name.remote,".R",sep=''),	paste(f.name.remote,'_',sites[i,"comid_old"],"_H1H1",".R",sep=''),	sites[i,"comid_old"],	"H1",	"H1")
									cat(cmd)
									cmd			<- prj.hpcwrapper(cmd, hpc.walltime=8, hpc.mem="1600mb", hpc.load="module load R/2.15",hpc.nproc=1, hpc.q="pqeph")
									prj.hpccaller(paste(CODE.HOME,"misc",sep='/'), paste("phd",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),"qsub",sep='.'), cmd)
								})
		}
		else
		{			 
			mlkl.theta			<- lapply(seq_len(nrow(sites)),function(i)
									{		
										#print( mlkl.theta[[i]] )
										if(verbose)	cat(paste("\nacute.loglkl.batch, process ",sites[i,"comid_old"],"\nnumber lkl evaluations: 4x ",nrow(mlkl.theta[[i]][["H0"]])))
										lkl.tpcH0.thetaH0	<- acute.loglkl.batch(sites[i,"comid_old"], tpc.obs[[i]][["H0"]], cohort.dur, mlkl.theta[[i]][["H0"]], clu.closure= 12, verbose=verbose)
										lkl.tpcH0.thetaH1	<- acute.loglkl.batch(sites[i,"comid_old"], tpc.obs[[i]][["H0"]], cohort.dur, mlkl.theta[[i]][["H1"]], clu.closure= 12, verbose=verbose)
										lkl.tpcH1.thetaH0	<- acute.loglkl.batch(sites[i,"comid_old"], tpc.obs[[i]][["H1"]], cohort.dur, mlkl.theta[[i]][["H0"]], clu.closure= 12, verbose=verbose)
										lkl.tpcH1.thetaH1	<- acute.loglkl.batch(sites[i,"comid_old"], tpc.obs[[i]][["H1"]], cohort.dur, mlkl.theta[[i]][["H1"]], clu.closure= 12, verbose=verbose)
										ans					<- list(tpcH0.thetaH0=lkl.tpcH0.thetaH0, tpcH0.thetaH1=lkl.tpcH0.thetaH1, tpcH1.thetaH0=lkl.tpcH1.thetaH0, tpcH1.thetaH1=lkl.tpcH1.thetaH1)
										#				print(ans); stop()
										ans
									})
			if(verbose)
				cat(paste("\nwrite tpc.lkl to",paste(f.name,".R",sep='')))	
			save(tpc.obs, sites, mlkl.theta, file=paste(f.name,".R",sep=''))				
		}					
	}
	mlkl.theta
}
###############################################################################
prog.acute.loglkl.batch.cmd<- function(dir.name,infile,outfile,site,tpcHx)
{
	cmd<- paste("\n",dir.name,"/misc/phdes.startme.R -exeACUTE.LKL.BATCH ",sep='')
	cmd<- paste(cmd, " -infile=",infile,sep='')
	cmd<- paste(cmd, " -outfile=",outfile,sep='')
	cmd<- paste(cmd, " -site=",site,sep='')
	cmd<- paste(cmd, " -tpcHx=",tpcHx,sep='')
	cmd
}
###############################################################################
prog.acute.loglkl.batch<- function()
{
	require(data.table)
	require(phylodesign)
	
	infile	<- NA
	outfile	<- NA
	site	<- NA
	verbose	<- 1
	tpcHx	<- "H0"
	
	if(exists("args"))
	{
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,7),
									infile= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile<- tmp[1]
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,8),
									outfile= return(substr(arg,10,nchar(arg))),NA)	}))
		if(length(tmp)>0) outfile<- tmp[1]
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,5),
									site= return(substr(arg,7,nchar(arg))),NA)	}))
		if(length(tmp)>0) site<- tmp[1]			
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,7),
									tpcHx= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) tpcHx<- tmp[1]	
	}
	if(verbose)
	{
		print(infile)
		print(outfile)
		print(site)
		print(tpcHx)
	}
	if(verbose)	cat(paste("\nloading file",infile))
	load(infile)
	
	i	<- which( sites[,"comid_old"]==site )
	tmp	<- subset(lkl.theta, comid_old==site)
	if(verbose)	cat(paste("\nprocessing"))
	if(verbose)	print(tmp)
	ans	<- acute.loglkl.batch(sites[i,"comid_old"], tpc.obs[[i]][[tpcHx]], cohort.dur, tmp, clu.closure= 12, verbose=verbose)
	if(verbose)	cat(paste("\nsaving output to file",outfile))
	save(ans, file=outfile)
	
}
###############################################################################
prj.popart.powercalc.by.acutelklratio.mlkl.plot<- function(mlkl.criterion, f.name)
{
	require(RColorBrewer)	
	mlkl.criterion[,pch:=15]		
	set(mlkl.criterion, which(mlkl.criterion[,arm=="B"]), "pch", 16)
	set(mlkl.criterion, which(mlkl.criterion[,arm=="C"]), "pch", 17)
	setkey(mlkl.criterion, "H1/H0|H1")		
	xlim					<- c(0,1)
	ylim					<- range( c(12,mlkl.criterion[,"H1/H0|H1", with=F]) )
	ylim[1]					<- ylim[1]*0.9
	ylim[2]					<- ylim[2]*1.5
	cols					<- brewer.pal(9, "PuBuGn")[1:4]
	cols.bf					<- matrix( c( min(0,ylim[1]),2,2,6,6,10,10,ylim[2] ), nrow=2, ncol=4 )
	cat(paste("plot centraltarget_BFpoint to file", f.name))
	pdf(file=f.name, width=4,height=8)
	par(mar=c(0,4,2,0))
	plot(1,1,type='n',bty='n',xlim=xlim, ylim=ylim, ylab="2 ln BF", xlab="", xaxt='n', log='y')	
	sapply(seq_len(ncol(cols.bf)), function(i){		polygon(	c(xlim, rev(xlim)), c(cols.bf[1,i],cols.bf[1,i],cols.bf[2,i],cols.bf[2,i]), border=NA, col=cols[i] )			})
	points(rep(0.5,nrow(mlkl.criterion)),unlist(mlkl.criterion[,"H1/H0|H1", with=F]), pch=mlkl.criterion[,pch])		
	text(rep(0.5,nrow(mlkl.criterion)) + rep(c(0.25,-0.25),nrow(mlkl.criterion)/2),unlist(mlkl.criterion[,"H1/H0|H1", with=F]), labels=gsub(' ','',mlkl.criterion[,comid_old]))
	legend("topright",pch=unique(mlkl.criterion[,pch]), legend=paste("arm",unique(mlkl.criterion[,arm])), bg="white", box.col="white")
	dev.off()
}		
###############################################################################
prj.popart.powercalc.by.acutelklratio.mlkl.increasingcoverage<- function(mlkl.theta)
{
	#get marginal likelihoods
	#	-- take mean of the above likelihood values
	#print(names(mlkl.theta))
	mlkl.criterion			<- lapply(seq_len(length(mlkl.theta)),function(i)
			{
				tmp				<- merge( mlkl.theta[[i]][["tpcH0.thetaH0"]], mlkl.theta[[i]][["tpcH0.thetaH1"]], by="sample")
				setnames(tmp, c("acute.x","base.x","lkl.x","acute.y","base.y","lkl.y"), c("acute.H0","base.H0","lkl.H0","acute.H1","base.H1","lkl.H1"))
				mlkl.ratio.H0	<- tmp[,	{
												scale.max<- max(c(lkl.H0, lkl.H1));
												list("H0/H1|H0"=2*log( mean( exp(lkl.H0-scale.max) ) / mean( exp(lkl.H1-scale.max) ) ) ) 
											},by=sample]
				#print(tmp)
				#print(mlkl.ratio.H0)				
				tmp				<- merge( mlkl.theta[[i]][["tpcH1.thetaH0"]], mlkl.theta[[i]][["tpcH1.thetaH1"]], by="sample")
				setnames(tmp, c("acute.x","base.x","lkl.x","acute.y","base.y","lkl.y"), c("acute.H0","base.H0","lkl.H0","acute.H1","base.H1","lkl.H1"))
				mlkl.ratio.H1	<- tmp[,	{
												scale.max	<- max(c(lkl.H0, lkl.H1));
												tmp			<- 2*log( mean( exp(lkl.H1-scale.max) ) / mean( exp(lkl.H0-scale.max) ) );												
												list("H1/H0|H1"=tmp, "comment"= ifelse( tmp<=2, "nil", ifelse( tmp<=6, "pos", ifelse( tmp<=10, "strong", "vstrong" ) ) )  ) 
											},by=sample]
				
				ans				<- merge(mlkl.ratio.H0,mlkl.ratio.H1, by="sample")
				ans[,comid_old:=names(mlkl.theta)[i]]
				ans				
			})	
	mlkl.criterion			<- do.call("rbind", mlkl.criterion)
	mlkl.criterion
}
###############################################################################
prj.popart.powercalc.by.acutelklratio.mlkl.increasingcoverage.plot<- function(mlkl.theta.cov, f.name)
{
	require(RColorBrewer)
	sapply( unique(mlkl.theta.cov[,arm]), function(a)
			{
				#print(a)
				mlkl.theta.cov.arm	<- subset(mlkl.theta.cov,arm==a)
				#print(mlkl.theta.cov.arm)
				file				<- paste(f.name,'_',a,".pdf",sep='')
				cols				<- brewer.pal(9, "PuBuGn")[1:4]
				xlim				<- c(0.1,1)	#range(unique(mlkl.theta.cov.arm[,sample]))		
				ylim				<- range(c(12,unlist( mlkl.theta.cov.arm[,"H1/H0|H1",with=F] )))		
				cols.bf				<- matrix( c( min(0,ylim[1]),2,2,6,6,10,10,ylim[2] ), nrow=2, ncol=4 )
				pdf(file=file, width=5,height=5)
				plot(1,1,type='n',bty='n',xlim=xlim, ylim=ylim, ylab="2 ln BF", xlab="sampling coverage")	
				sapply(seq_len(ncol(cols.bf)), function(i){		polygon(	c(xlim, rev(xlim)), c(cols.bf[1,i],cols.bf[1,i],cols.bf[2,i],cols.bf[2,i]), border=NA, col=cols[i] )			})
				mlkl.theta.cov.arm.s<- unique(mlkl.theta.cov.arm[,comid_old])
				sapply(seq_along(mlkl.theta.cov.arm.s), function(i)
						{				
							tmp		<- subset(mlkl.theta.cov.arm, comid_old==mlkl.theta.cov.arm.s[i])
							#print(tmp)
							x		<- tmp[,sample]
							y		<- unlist( tmp[,"H1/H0|H1",with=F] )
							lines(x,y, lty=i)
							points(x[1],y[1], pch=18)							
						})
				legend("topleft", bty='n', lty=seq_along(mlkl.theta.cov.arm.s), legend=mlkl.theta.cov.arm.s)
				dev.off()
			})
}
###############################################################################
#sim<- abc.struct; target<- sites[,c("comid_old","mu.inc.rate.H0","mu.inc.rate.H1","mu.pE2E.H0","mu.pE2E.H1")]; hetclu.scale<- 1; hetclu.cov<- ( hetclu.scale*matrix(c(0.0015/2,0.006,0.006,0.1/2),2,2) )^2
prj.popart.powercalc.by.acutelklratio.matchpa<- function(sim, target, hetclu.scale= 1, hetclu.cov= ( hetclu.scale*matrix(c(0.015/2,0.006,0.006,0.1/2),2,2) )^2, return.top= 10)
{
	hetclu.icov	<- solve(hetclu.cov)
	tmp			<- merge(sim, target, by="comid_old", all.x=1)
	x			<- as.matrix(subset(tmp,select=c(Inc,E2E)) - subset(tmp,select=c(mu.inc.rate.H0,mu.pE2E.H0)))
	d.H0		<- rowSums((x %*% hetclu.icov) * x)	
	x			<- as.matrix(subset(tmp,select=c(Inc,E2E)) - subset(tmp,select=c(mu.inc.rate.H1,mu.pE2E.H1)))
	d.H1		<- rowSums((x %*% hetclu.icov) * x)
	ans.H0		<- copy(sim)
	ans.H0[,dist:= d.H0]	
	ans.H0[,weight:= dnorm(ans.H0[,dist])]	
	ans.H1		<- sim
	ans.H1[,dist:= d.H1]
	ans.H1[,weight:= dnorm(ans.H1[,dist])]
	
	setkey(ans.H1, dist)
	setkey(ans.H0, dist)	
	ans.H0		<- ans.H0[seq_len(return.top),]
	ans.H1		<- ans.H1[seq_len(return.top),]
	
	ans			<- list(H0= ans.H0, H1= ans.H1)
	ans
}
###############################################################################
prj.popart.powercalc.by.acutelklratio	<- function()
{
	require(data.table)
	
	dir.name		<- "popartpower_acute"
	my.mkdir(DATA,dir.name)
	dir.name		<- paste(DATA,dir.name,sep='/')	
	resume			<- 1
	verbose			<- 1
	plot.increment	<- 0.05
	
	m.type			<- "Acute"	
	cohort.size		<- 2500	
	cohort.dur		<- 3	
	theta.EE.H0		<- 0.1
	theta.EE.H1		<- 0.4
	theta.UE		<- 0.3
	theta.TE		<- theta.UE / 5
	test.alpha		<- 0.05
	p.nocontam		<- 0.95
	debug			<- 1
	pooled.n		<- 1
	opt.pooled		<- "no pooling"#"pooled across ZA"#"pooled across trial"#"no pooling"
	opt.pooled		<- "pooled across SA"
	opt.clu.closure	<- 14
	opt.sampling	<- "PC12+HCC"	#"PC and HCC"#"only HCC"	#"PC and HCC"	#
	opt.power		<- "All"
	
	p.lab			<- 0.75*0.9			#set lower as discussed	70% from CD4 90% from sequencing
	p.consent.coh	<- 0.9*0.9			#90% consent to main study and of those 90% consent to phylo study				
	p.consent.clu	<- 1				#waiver
	p.vhcc.prev.AB	<- 1				#already in PopART model estimate
	p.vhcc.inc.AB	<- 1				#already in PopART model estimate
	p.vhcc.prev.C	<- 1				#already in PopART model estimate
	p.vhcc.inc.C	<- 1				#already in PopART model estimate
	#p.contam		<- seq(0.05,0.2,0.025)
	
	if(verbose)
	{
		cat(paste("\ncohort.size",cohort.size))
		cat(paste("\ncohort.dur",cohort.dur))
		cat(paste("\ntheta.EE.H0",theta.EE.H0))
		cat(paste("\ntheta.EE.H1",theta.EE.H1))
		cat(paste("\ntest.alpha",test.alpha))	
		cat(paste("\np.nocontam",p.nocontam))
		cat(paste("\npooled.n",pooled.n))
		cat(paste("\nopt.pooled",opt.pooled))
		cat(paste("\nopt.power",opt.power))
		cat(paste("\nopt.sampling",opt.sampling))
	}

	#simulate high acute and low acute tip clusters for each community
	tmp				<- prj.popart.powercalc.by.acutelklratio.tpcobs(theta.EE.H0, theta.EE.H1, cohort.dur, p.consent.coh, p.consent.clu, p.lab, p.vhcc.prev.AB, p.vhcc.inc.AB, p.vhcc.prev.C, p.vhcc.inc.C, opt.sampling, pooled.n, dir.name=dir.name, verbose=1, resume=1, standalone=0)
	tpc.obs			<- tmp$tpc.obs
	theta.model.Hx	<- tmp$theta.model.Hx 
	sites			<- tmp$sites	
	print(theta.model.Hx)
	print(sites)
	#print(tpc.obs)
	if(0)
	{
		require(xtable)		
		tpc.obs.latex	<- lapply(seq_along(tpc.obs), function(i)
								{
									tmp<- lapply(seq_along(tpc.obs[[i]]), function(j)
											{
												txt.start	<- "\\begin{table}[htbp]\n\\centering\n{\\footnotesize\n"
												txt.table	<- print( xtable(tpc.obs[[i]][[j]],digits=0), floating=FALSE, print.results=FALSE )
												txt.caption	<- paste( names(tpc.obs)[i], '_', names(tpc.obs[[i]])[j], sep='')
												txt.end		<- paste( "}\n\\caption{", txt.caption ,"}\n\\end{table}\n", sep='')
												paste(txt.start, txt.table, txt.end, sep='')
											})	
									tmp	<- paste(unlist(tmp), collapse="\\n\n", sep='')
									tmp
								})
		tpc.obs.latex	<- paste(tmp, collapse="\\n\n", sep='')
		#cat(tmp)
	}

	#TODO get acute parameters that match target effects under H0, H1. 
	#	this accounts for between cluster heterogeneity by adding a site specific random effect whose magnitude is specified by 'hetclu.scale'
	if(0)
	{
		require(MASS)	
		tmp	<- mvrnorm(1e1, c(sites[1,"mu.inc.rate.H0"], sites[1,"mu.pE2E.H0"]), hetclu.cov)
		tmp	<- rbind(tmp, mvrnorm(1e1, c(sites[1,"mu.inc.rate.H1"], sites[1,"mu.pE2E.H1"]), hetclu.cov))
		tmp	<- data.table(comid_old="Ndeke", Inc= tmp[,1], E2E=tmp[,2])	
		#TODO load matches.R and create data.table
		tmp<- prj.popart.powercalc.by.acutelklratio.matchpa( tmp, sites[,c("comid_old","mu.inc.rate.H0","mu.inc.rate.H1","mu.pE2E.H0","mu.pE2E.H1")], hetclu.scale= 0.25 )
		#returns list of H0 data.table with comid_old and H1 data.table with comid_old
	}
#print(sites)

	#get likelihood values
	#	-- use the predicted seq coverage as stored in 'sites' to build a prior on the sampling prob
	#	-- get likelihood values of the simulated tip cluster for this set of model parameters
	
	mlkl.theta.model	<- rbindlist( lapply(seq_len(nrow(sites)),function(i)
							{
								tmp	<- subset(theta.model.Hx, comid_old==sites[i,"comid_old"])								
								data.table( comid_old= sites[i,"comid_old"], h="H0", acute=tmp[1, acute.H0], base= tmp[1, base.H0], sample.mu=sites[i,"%avg"], sample.sigma=sites[i,"sigma"]  )
							}))
	mlkl.theta.model	<- rbind(mlkl.theta.model, rbindlist( lapply(seq_len(nrow(sites)),function(i)
							{
								tmp	<- subset(theta.model.Hx, comid_old==sites[i,"comid_old"])
								data.table( comid_old= sites[i,"comid_old"], h="H1", acute=tmp[1, acute.H1], base= tmp[1, base.H1], sample.mu=sites[i,"%avg"], sample.sigma=sites[i,"sigma"]  )
							})))	
	remote.signat		<- "Fri_Oct_11_10:30:19_2013"
	#print( mlkl.theta.model.H0 )
	#stop()			
	f.name				<- paste(dir.name,'/',"tpclkl_",m.type,'_',theta.EE.H0,'_',theta.EE.H1,'_',opt.sampling,'_',"central",'_',p.lab,'_',p.consent.coh,sep='')
	mlkl.theta			<- prj.popart.powercalc.by.acutelklratio.lklH0H1(sites, tpc.obs, mlkl.theta.model, mlkl.n= 1e3, cohort.dur=cohort.dur, f.name=f.name, resume=1, verbose=1, remote=0, remote.signat=remote.signat)

	if(0)	#coverage required calculations
	{
		#get likelihood values for increasing sampling coverage
		mlkl.theta.model.H0	<- do.call("rbind", lapply(seq_len(nrow(sites)),function(i)
								{
									tmp	<- subset(theta.model.Hx, comid_old==sites[i,"comid_old"])
									data.table( comid_old= sites[i,3], acute=tmp[1, acute.H0], base= tmp[1, base.H0], sample.mu=seq(round(sites[i,"%avg"],d=2),1,by=0.025), sample.sigma=0  )
								}))
		mlkl.theta.model.H1	<- do.call("rbind", lapply(seq_len(nrow(sites)),function(i)
								{
									tmp	<- subset(theta.model.Hx, comid_old==sites[i,"comid_old"])
									data.table( comid_old= sites[i,3], acute=tmp[1, acute.H1], base= tmp[1, base.H1], sample.mu=seq(round(sites[i,"%avg"],d=2),1,by=0.025), sample.sigma=0  )
								}))
		f.name				<- paste(dir.name,'/',"tpclkl_",m.type,'_',theta.EE.H0,'_',theta.EE.H1,'_',opt.sampling,'_',"increasingcoverage",'_',p.lab,'_',p.consent.coh,sep='')	
		mlkl.theta.cov		<- prj.popart.powercalc.by.acutelklratio.lklH0H1(sites, tpc.obs, mlkl.theta.model.H0, mlkl.theta.model.H1, mlkl.n= 1e3, cohort.dur=3, f.name=f.name, resume=1, verbose=1, remote=0, remote.signat=remote.signat)
		names(mlkl.theta.cov)	<- sites[,"comid_old"]	
		mlkl.theta.cov		<- prj.popart.powercalc.by.acutelklratio.mlkl.increasingcoverage(mlkl.theta.cov)	
		mlkl.theta.cov		<- merge(mlkl.theta.cov, as.data.table(sites[,c("comid_old","arm","%avg","popsize")]), by="comid_old")
		prj.popart.powercalc.by.acutelklratio.mlkl.increasingcoverage.plot(mlkl.theta.cov, f.name)
		#print( mlkl.theta.cov )
	}	
	if(1)	#plot integrated Bayes factor for central targets
	{
		qu<-	0.9 
		#compute marginal likelihood values on grid and select best 10% as those corresponding to H0 and H1.
		mlkl.criterion	<- sapply(seq_len(nrow(sites)),function(i)
				{
					args			<<- prj.acute.test.lkl.wsampling.bothUandT.cmd(CODE.HOME, sites[i,"comid_old"], theta.model.Hx[i,acute.H0], theta.model.Hx[i,base.H0], round(sites[i,"%avg"],d=2), round(sites[i,"%avg"],d=2), cohort.dur, by.sample=2, save=1, resume=1, verbose=0, debug.susc.const=0, debug.only.u=0)		 
					args			<<- unlist(strsplit(args,' '))							
					tipc.lkl		<- prj.acute.test.lkl.wsampling.bothUandT()				
					ans				<- t( sapply( strsplit(names(tipc.lkl),'_',fixed=1), as.numeric ) )
					colnames(ans)	<- c("acute","base","sample")
					ans				<- as.data.table(ans)
					ans[,comid_old:= sites[i,"comid_old"]]
					ans[,lkl.H0:=tipc.lkl]								
					args			<<- prj.acute.test.lkl.wsampling.bothUandT.cmd(CODE.HOME, sites[i,"comid_old"], theta.model.Hx[i,acute.H1], theta.model.Hx[i,base.H1], round(sites[i,"%avg"],d=2), round(sites[i,"%avg"],d=2), cohort.dur, by.sample=2, save=1, resume=1, verbose=0, debug.susc.const=0, debug.only.u=0)		 
					args			<<- unlist(strsplit(args,' '))							
					tipc.lkl		<- prj.acute.test.lkl.wsampling.bothUandT()				
					ans[,lkl.H1:=tipc.lkl]
					
					#subset of all those theta that explain H0 data best		
					thetaH0			<- subset(ans, lkl.H0>=quantile(lkl.H0, probs=qu))
					#subset of all those theta that explain H1 data best
					thetaH1			<- subset(ans, lkl.H1>=quantile(lkl.H1, probs=qu))					
					#lkl for H0 data for the thetaH0 para
					tpcH0.thetaH0	<- subset(thetaH0, select=c(comid_old, acute, base, sample, lkl.H0))					
					#lkl for H0 data for the thetaH1 para
					tpcH0.thetaH1	<- subset(thetaH1, select=c(comid_old, acute, base, sample, lkl.H0))
					#lkl for H1 data for the thetaH0 para
					tpcH1.thetaH0	<- subset(thetaH0, select=c(comid_old, acute, base, sample, lkl.H1))
					#lkl for H1 data for the thetaH1 para
					tpcH1.thetaH1	<- subset(thetaH1, select=c(comid_old, acute, base, sample, lkl.H1))										
					#print(ans); print(tpcH0.thetaH0); print(tpcH0.thetaH1); print(tpcH1.thetaH0); print(tpcH1.thetaH1)
					
					scale.max			<- max( tpcH0.thetaH0[,lkl.H0], tpcH0.thetaH1[,lkl.H0] )	
					mlkl.scaled			<- do.call("cbind", list( H0=tpcH0.thetaH0[,lkl.H0] - scale.max, H1=tpcH0.thetaH1[,lkl.H0] - scale.max ))	#norm constant for both is now unknown C_x times exp(scale.max)
					#print( exp(mlkl.scaled ) ) 
					mlkl.ratio			<- apply(exp(mlkl.scaled),2,mean)
					mlkl.ratio.H0		<- mlkl.ratio["H0"]/mlkl.ratio["H1"]
					
					scale.max			<- max( tpcH1.thetaH0[,lkl.H1], tpcH1.thetaH1[,lkl.H1] )	
					mlkl.scaled			<- do.call("cbind", list( H0=tpcH1.thetaH0[,lkl.H1] - scale.max, H1=tpcH1.thetaH1[,lkl.H1] - scale.max ))	#norm constant for both is now unknown C_x times exp(scale.max)
					#print( exp(mlkl.scaled ) ) 
					mlkl.ratio			<- apply(exp(mlkl.scaled),2,mean)
					mlkl.ratio.H1		<- mlkl.ratio["H1"]/mlkl.ratio["H0"]
					
					ans					<- 2 * log( c( mlkl.ratio.H0, mlkl.ratio.H1 ) )
					ans					<- c(ans, ifelse( ans[2]<=3, 0, ifelse( ans[2]<=20, 1, ifelse( ans[2]<=150, 2, 3 ) ) ) )  
					names(ans)			<- c("H0/H1|H0","H1/H0|H1","comment")
					ans					
				})		
		mlkl.criterion			<- as.data.table( t(mlkl.criterion) )		
print(mlkl.criterion)		
		set(mlkl.criterion, NULL, "comment", factor(mlkl.criterion[,comment], levels= 0:3, labels= c("nil","pos","strong","vstrong") ) )
		setnames(sites, "%avg", "perc.avg.seq.cov")
		mlkl.criterion			<- cbind(subset(sites, select= c(comid_old, popsize, hivcomb, artadjust, arm, perc.avg.seq.cov)), mlkl.criterion)
		print( mlkl.criterion )
		mlkl.criterion			<- as.data.table(mlkl.criterion)
		f.name					<- paste(dir.name,'/',"tpclkl_",m.type,'_',theta.EE.H0,'_',theta.EE.H1,'_',opt.sampling,'_',"centraltarget_BFavg",'_',p.lab,'_',p.consent.coh,".pdf",sep='')
		prj.popart.powercalc.by.acutelklratio.mlkl.plot(mlkl.criterion, f.name)
	
	}
	if(0)	#plot point Bayes factor for central targets
	{
		#get marginal likelihoods
		#	-- take mean of the abovelikelihood values
		mlkl.criterion			<- sapply(seq_len(nrow(sites)),function(i)
									{
										scale.max			<- max( mlkl.theta[[i]][["tpcH0.thetaH0"]][,lkl], mlkl.theta[[i]][["tpcH0.thetaH1"]][,lkl] )	
										mlkl.scaled			<- do.call("cbind", list( H0=mlkl.theta[[i]][["tpcH0.thetaH0"]][,lkl] - scale.max, H1=mlkl.theta[[i]][["tpcH0.thetaH1"]][,lkl] - scale.max ))	#norm constant for both is now unknown C_x times exp(scale.max)
										#print( exp(mlkl.scaled ) ) 
										mlkl.ratio			<- apply(exp(mlkl.scaled),2,mean)
										mlkl.ratio.H0		<- mlkl.ratio["H0"]/mlkl.ratio["H1"]
										
										scale.max			<- max( mlkl.theta[[i]][["tpcH1.thetaH0"]][,lkl], mlkl.theta[[i]][["tpcH1.thetaH1"]][,lkl] )	
										mlkl.scaled			<- do.call("cbind", list( H0=mlkl.theta[[i]][["tpcH1.thetaH0"]][,lkl] - scale.max, H1=mlkl.theta[[i]][["tpcH1.thetaH1"]][,lkl] - scale.max ))	#norm constant for both is now unknown C_x times exp(scale.max)
										#print( exp(mlkl.scaled ) ) 
										mlkl.ratio			<- apply(exp(mlkl.scaled),2,mean)
										mlkl.ratio.H1		<- mlkl.ratio["H1"]/mlkl.ratio["H0"]
										
										ans					<- 2 * log( c( mlkl.ratio.H0, mlkl.ratio.H1 ) )
										ans					<- c(ans, ifelse( ans[2]<=3, 0, ifelse( ans[2]<=20, 1, ifelse( ans[2]<=150, 2, 3 ) ) ) )  
										names(ans)			<- c("H0/H1|H0","H1/H0|H1","comment")
										ans
									})		
		#colnames(mlkl.criterion)<- sites[,"comid_old"]
		#print(sites)	
		mlkl.criterion			<- as.data.table( t(mlkl.criterion) )		
		set(mlkl.criterion, NULL, "comment", factor(mlkl.criterion[,comment], levels= 0:3, labels= c("nil","pos","strong","vstrong") ) )
		setnames(sites, "%avg", "perc.avg.seq.cov")
		mlkl.criterion			<- cbind(subset(sites, select= c(comid_old, popsize, hivcomb, artadjust, arm, perc.avg.seq.cov)), mlkl.criterion)
		print( mlkl.criterion )
		
		mlkl.criterion			<- as.data.table(mlkl.criterion)
		f.name					<- paste(dir.name,'/',"tpclkl_",m.type,'_',theta.EE.H0,'_',theta.EE.H1,'_',opt.sampling,'_',"centraltarget_BFpoint",'_',p.lab,'_',p.consent.coh,".pdf",sep='')
		prj.popart.powercalc.by.acutelklratio.mlkl.plot(mlkl.criterion, f.name)				
	}
	stop()
	
}
###############################################################################
prj.acute.test.lkl.wsampling.bothUandT.cmd<- function(dir.name, loc, acute, base, sIdx, sE, cluster.tw, by.sample=1, debug.susc.const=0, debug.only.u=0, save=1, resume=1, verbose=1, plot.2D=0)
{		
	cmd<- paste("\n",dir.name,"/misc/phdes.startme.R -exeWSAMPLINGBOTHUANDT ",sep='')
	cmd<- paste(cmd, " -loc=",loc,sep='')
	cmd<- paste(cmd, " -acute=",acute,sep='')
	cmd<- paste(cmd, " -baseline=",base,sep='')
	cmd<- paste(cmd, " -sIdx=",sIdx,sep='')
	cmd<- paste(cmd, " -sE=",sE,sep='')
	cmd<- paste(cmd, " -cluster.tw=",cluster.tw,sep='')
	cmd<- paste(cmd, " -by.sample=",by.sample,sep='')
	cmd<- paste(cmd, " -save=",save,sep='')
	cmd<- paste(cmd, " -resume=",resume,sep='')
	cmd<- paste(cmd, " -debug.only.u=",debug.only.u,sep='')
	cmd<- paste(cmd, " -debug.susc.const=",debug.susc.const,sep='')	
	cmd<- paste(cmd, " -v=",verbose,sep='')
	cmd<- paste(cmd, " -plot=",plot.2D,sep='')
	cmd
}
###############################################################################
prj.acute.test.lkl.wsampling.bothUandT	<- function()
{
	require(RColorBrewer)
	m.type				<- "Acute"
	loc.type			<- "TownII"
	loc.type			<- "Ndeke"
	m.popsize			<- NA
	m.known.states		<- 1	
	theta0				<- c(4.4, 0.09, 0, 0)		#1.25   0.062      4.4   0.090      Ndeke
	#theta0				<- c(1.25, 0.062, 0, 0)
	sample.prob			<- 0.62
	#theta0				<- c(15, 0.032, 0, 0)		#TownII
	#theta0				<- c(2.2, 0.052, 0, 0)
	#sample.prob		<- 0.4
	
	names(theta0)		<- c("acute","base","m.st1","m.st2")
	sample.prob			<- c(sample.prob,sample.prob)
	names(sample.prob)	<- c("Idx","E")
	cluster.tw			<- 3
	m.repeat			<- 1
	resume				<- 1
	verbose				<- 0	
	debug.onlyu			<- 0
	debug.suscconst		<- 0
	by.sample			<- 1
	plot.2D				<- 1
	if(exists("args"))
	{
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,4),
									loc= return(substr(arg,6,nchar(arg))),NA)	}))
		if(length(tmp)>0) loc.type<- tmp[1]
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
						{	switch(substr(arg,2,10),
									by.sample= return(as.numeric(substr(arg,12,nchar(arg)))),NA)	}))
		if(length(tmp)>0) by.sample<- tmp[1]							
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,2),
									v= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) verbose<- tmp[1]
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
		if(length(tmp)>0) resume<- tmp[1]
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,5),
									plot= return(as.numeric(substr(arg,7,nchar(arg)))),NA)	}))
		if(length(tmp)>0) plot.2D<- tmp[1]		
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,17),
									debug.susc.const= return(as.numeric(substr(arg,19,nchar(arg)))),NA)	}))
		if(length(tmp)>0) debug.susc.const<- tmp[1]
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,13),
									debug.only.u= return(as.numeric(substr(arg,15,nchar(arg)))),NA)	}))
		if(length(tmp)>0) debug.only.u<- tmp[1]				
	}
	if(verbose)
	{
		cat(paste("\nloc.type ",loc.type))
		cat(paste("\ncluster.tw ",cluster.tw))	
		cat(paste("\ntheta0.acute ",theta0[1]))
		cat(paste("\ntheta0.base ",theta0[2]))
		cat(paste("\nsample.prob.Idx ",sample.prob[1]))
		cat(paste("\nsample.prob.E ",sample.prob[2]))
		cat(paste("\ndebug.susc.const ",debug.susc.const))
		cat(paste("\ndebug.only.u ",debug.only.u))
		cat(paste("\nby.sample ",by.sample))		
		cat(paste("\nverbose ",verbose))
		cat(paste("\nresume ",resume))
		cat(paste("\nplot.2D ",plot.2D))		
	}
	dir.name			<- paste(DATA,"/acutesimu_fxs",debug.suscconst,"_onlyu",debug.onlyu,sep='')
	f.name				<- paste(dir.name,paste("tpcdat_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta0[1],"_b",theta0[2],"_sIdx",sample.prob[1],"_sE",sample.prob[2],"_tw",cluster.tw,"_median",sep=''),sep='/')
	cat(paste("\nload ",paste(f.name,".R",sep='')))		
	load(paste(f.name,".R",sep=''))
	if(verbose)
	{
		print(tpc.table.all.median)
		print(tpc.table.sample.median)
	}
	if(by.sample==0)
	{
		theta			<- expand.grid(acute= seq(0.5,13,0.25), base= theta0["base"], sample= seq(0.1,0.9,0.05))
	}
	else if(by.sample==1)
	{
		theta			<- expand.grid(acute= seq(0.5,13,0.25), base= seq(0.01,0.13,0.0025), sample=sample.prob[1] )		
	}
	else if(by.sample==2)
	{
		theta			<- expand.grid(acute= seq(0.5,13,0.25), base= seq(0.01,0.13,0.0025), sample= seq(0.1,0.9,0.05) )		
	}
	else
		stop("unknown by.sample")
	if(resume)
	{
		options(show.error.messages = FALSE)
		f.name			<- paste(dir.name,paste("tpcsampled_bysample",by.sample,"_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta0[1],"_b",theta0[2],"_sIdx",sample.prob[1],"_sE",sample.prob[2],"_tw",cluster.tw,"_lkl",sep=''),sep='/')		
		cat(paste("\ntry load ",paste(f.name,".R",sep='')))				
		readAttempt		<-try(suppressWarnings(load(paste(f.name,".R",sep=''))))
		options(show.error.messages = TRUE)
		if(!inherits(readAttempt, "try-error"))
			cat(paste("\nresumed file ",paste(f.name,".R",sep='')))
	}
	if(!resume || inherits(readAttempt, "try-error"))
	{		
		acute.MAX.TIPC.SIZE	<<- 12 
		ibm				<- ibm.collapse( ibm.init.model(m.type, loc.type, NA, theta0, save='', resume= 0, init.pop=0) )
		state.n			<- ibm$init.pop.distr$status * ibm$init.pop.distr$npop
		pop.n			<- ibm$init.pop.distr$npop
		if(verbose)
		{
			print(pop.n)
			print(state.n)	
		}
		clu.n			<- clu.tipc.n(acute.MAX.TIPC.SIZE)			
		clu.n.sum		<- apply(clu.n,2,sum)																	#faster than seq_len(ncol(clu.n))^seq.int(-1,ncol(clu.n)-2)
		lclu.n			<- log( clu.n / matrix(clu.n.sum, nrow=nrow(clu.n), ncol=length(clu.n.sum), byrow=1) ) 		
		#evaluate likelihood over parameter space
		tipc.lkl		<- sapply(seq_len(nrow(theta)),function(i)
				{				
#i<- 10
#print(theta[i,])
					ibm[["beta"]][['i']][["status"]]['i']	<- theta[i,"acute"]
					ibm[["beta"]][["base"]]					<- theta[i,"base"]	
					rate.m									<- acute.get.rates(ibm[["beta"]], ibm.pop= NULL, pop.n=pop.n, state.n= as.matrix(state.n), per.capita.i= 1)	
					dT										<- cluster.tw
					lkl										<- acutesampled.loglkl(tpc.table.sample.median, rate.m, theta[i,"sample"], dT, lclu.n=lclu.n)				
					lkl[["table.lkl"]]
				})	
		names(tipc.lkl)<- apply(theta,1,function(x)	paste(x,collapse='_',sep='') )					
		cat(paste("\nsave tipc.lkl to file ",paste(f.name,".R",sep='')))
		save(tipc.lkl, file=paste(f.name,".R",sep=''))				
	}
	#print(tipc.lkl)	
	if(plot.2D)
	{
		theta.mle		<- as.numeric(strsplit(names(tipc.lkl)[which.max(tipc.lkl)],'_')[[1]])
		print(theta.mle)
		
		theta			<- t( sapply( strsplit(names(tipc.lkl),'_',fixed=1), as.numeric ) )
		colnames(theta)	<- c("acute","base","sample")
		if(by.sample)
		{
			theta		<- lapply(c("acute","base"), function(x) unique( theta[, x] ) )
			theta.mle	<- theta.mle[c(1,2)]
			theta.true	<- c(theta0[1],theta0[2])
			plot.labels	<- list(acute=expression(beta[EE]/beta[UE]), base=expression(beta[0]))
		}
		else
		{
			theta		<- lapply(c("acute","sample"), function(x) unique( theta[, x] ) )
			theta.mle	<- theta.mle[c(1,3)]
			theta.true	<- c(theta0[1],sample.prob[1])
			plot.labels	<- list(acute=expression(beta[EE]/beta[UE]), sample=expression(s))
		}		
		print(theta)
		tipc.lkl.m	<- matrix(tipc.lkl, nrow= length(theta[[1]]), ncol= length(theta[[2]]), dimnames=list(theta[[1]],theta[[2]]))
		cat(paste("\nmax lkl is ",max(tipc.lkl.m)))
		f.name		<- paste(dir.name,paste("tpc_bysample",by.sample,"_fxs",debug.suscconst,"_onlyu",debug.onlyu,"_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta0[1],"_b",theta0[2],"_s",sample.prob[1],sep=''),sep='/')
		cat(paste("\nplot",paste(f.name,"_2D.pdf",sep='')))
		pdf(paste(f.name,"_2D.pdf",sep=''),version="1.4",width=4,height=4)
		my.image(theta[[1]],theta[[2]],tipc.lkl.m, xlab=plot.labels[[1]],ylab=plot.labels[[2]],nlevels=10)			
		abline(v=theta.mle[1],lty=2)
		abline(h=theta.mle[2],lty=2)
		points(theta.true[1],theta.true[2], pch=19, col="red")
		dev.off()		
	}
	tipc.lkl
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
	if(0)
	{
		prj.acute.test.lkl.wsampling.onlyU()
	}
	if(1)
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
	#
	#	start 'prj.simudata.match.theta.to.Inc.E2E'	
	#	this simulates tip clusters for prior acute, beta and records INC and %E2E
	#	
	if(0)	
	{
		sites		<- popart.getdata.randomized.arm( 1, rtn.fixed=debug, rtn.phylostudy=1 )
		nit			<- 25e3
		sites		<- subset(sites, comid_old=="Ikhwezi")
		dir.name	<- CODE.HOME
		sapply(sites$comid_old, function(loc)
				{
					cmd			<- paste("\n",dir.name,"/misc/phdes.startme.R -exeSIMU.MATCH",sep='')
					cmd			<- paste(cmd, " -loc=",loc,sep='')
					cmd			<- paste(cmd, " -nit=",nit,sep='')
					
					cmd			<- prj.hpcwrapper(cmd, hpc.walltime=500, hpc.mem="1600mb", hpc.load="module load R/2.15",hpc.nproc=1, hpc.q="pqeph")
					cat(cmd)
					signat		<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
					outdir		<- paste(CODE.HOME,"misc",sep='/')
					outfile		<- paste("phd",signat,"qsub",sep='.')
					prj.hpccaller(outdir, outfile, cmd)			
				})		
	}
	#
	#	compute representative theta corresponding to H0 and H1 and simulate tip cluster table for this theta
	#
	if(0)	
	{
		require(data.table)		
		dir.name		<- "popartpower_acute"
		my.mkdir(DATA,dir.name)
		dir.name		<- paste(DATA,dir.name,sep='/')	
		resume			<- 1
		verbose			<- 1
		plot.increment	<- 0.05
		
		m.type			<- "Acute"	
		cohort.size		<- 2500	
		cohort.dur		<- 3	
		theta.EE.H0		<- 0.1
		theta.EE.H1		<- 0.4
		theta.UE		<- 0.3
		theta.TE		<- theta.UE / 5
		test.alpha		<- 0.05
		p.nocontam		<- 0.95
		debug			<- 1
		
		p.lab			<- 0.75*0.9			#set lower as discussed	70% from CD4 90% from sequencing
		p.consent.coh	<- 0.9*0.9			#90% consent to main study and of those 90% consent to phylo study				
		p.consent.clu	<- 1				#waiver
		p.vhcc.prev.AB	<- 1				#already in PopART model estimate
		p.vhcc.inc.AB	<- 1				#already in PopART model estimate
		p.vhcc.prev.C	<- 1				#already in PopART model estimate
		p.vhcc.inc.C	<- 1				#already in PopART model estimate
		
		opt.sampling	<- "PC12+HCC"
		pooled.n		<- 1
		
		if(verbose)
		{
			cat(paste("\ncohort.size",cohort.size))
			cat(paste("\ncohort.dur",cohort.dur))
			cat(paste("\ntheta.EE.H0",theta.EE.H0))
			cat(paste("\ntheta.EE.H1",theta.EE.H1))
			cat(paste("\ntest.alpha",test.alpha))	
			cat(paste("\np.nocontam",p.nocontam))			
		}	
		dummy			<- prj.popart.powercalc.by.acutelklratio.tpcobs(theta.EE.H0, theta.EE.H1, cohort.dur, p.consent.coh, p.consent.clu, p.lab, p.vhcc.prev.AB, p.vhcc.inc.AB, p.vhcc.prev.C, p.vhcc.inc.C, opt.sampling, pooled.n, dir.name=dir.name, verbose=verbose, resume=resume, standalone=1)
	}
	if(1)
	{
		require(data.table)		
		dir.name		<- "popartpower_acute"
		my.mkdir(DATA,dir.name)
		dir.name		<- paste(DATA,dir.name,sep='/')	
		resume			<- 1
		verbose			<- 1
		plot.increment	<- 0.05
		
		m.type			<- "Acute"	
		cohort.size		<- 2500	
		cohort.dur		<- 3	
		theta.EE.H0		<- 0.1
		theta.EE.H1		<- 0.4
		theta.UE		<- 0.3
		theta.TE		<- theta.UE / 5
		test.alpha		<- 0.05
		p.nocontam		<- 0.95
		debug			<- 1
		pooled.n		<- 1		
		opt.sampling	<- "PC12+HCC"		#"PC and HCC"#"only HCC"	#"PC and HCC"	#
		
		p.lab			<- 0.75*0.9			#set lower as discussed	70% from CD4 90% from sequencing
		p.consent.coh	<- 0.9*0.9			#90% consent to main study and of those 90% consent to phylo study				
		p.consent.clu	<- 1				#waiver
		p.vhcc.prev.AB	<- 1				#already in PopART model estimate
		p.vhcc.inc.AB	<- 1				#already in PopART model estimate
		p.vhcc.prev.C	<- 1				#already in PopART model estimate
		p.vhcc.inc.C	<- 1				#already in PopART model estimate				
			
		#	load high acute and low acute tip clusters for each community
		tmp				<- prj.popart.powercalc.by.acutelklratio.tpcobs(theta.EE.H0, theta.EE.H1, cohort.dur, p.consent.coh, p.consent.clu, p.lab, p.vhcc.prev.AB, p.vhcc.inc.AB, p.vhcc.prev.C, p.vhcc.inc.C, opt.sampling, pooled.n, dir.name=dir.name, verbose=1, resume=1, standalone=0)
		tpc.obs			<- tmp$tpc.obs
		theta.model.Hx	<- tmp$theta.model.Hx 
		sites			<- tmp$sites
		#	if not sensitivity analysis, reset 'sigma' to zero
		if(1)
			sites$sigma	<- 0	
		#	get likelihood values for precomputed theta-E2E/Inc values
		remote.signat		<- "Fri_Oct_11_10:30:19_2013"
		f.name				<- paste(dir.name,'/',"tpclkl_",m.type,'_',theta.EE.H0,'_',theta.EE.H1,'_',opt.sampling,'_',"central",'_',p.lab,'_',p.consent.coh,sep='')
		dummy				<- prj.popart.powercalc.by.acutelklratio.lkl4Precomputed(sites, tpc.obs, cohort.dur=cohort.dur, f.name=f.name, resume=1, verbose=1, remote=1, remote.signat=remote.signat)
		#	if remote, need to collect results
		
	}	
	if(0)	#start 'prj.popart.powercalc.by.acutelklratio' for all locations
	{		
		
		dir.name	<- CODE.HOME
		cmd			<- paste("\n",dir.name,"/misc/phdes.startme.R -exePOPART.POWER.ACUTELKLRATIO",sep='')
		cmd			<- prj.hpcwrapper(cmd, hpc.walltime=71, hpc.mem="3800mb", hpc.load="module load R/2.15",hpc.nproc=1, hpc.q="pqeph")
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
	cmd<- paste(cmd, " -loc=",loc,sep='')
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
	resume				<- 0
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
						{	switch(substr(arg,2,4),
									loc= return(substr(arg,6,nchar(arg))),NA)	}))
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
						{	switch(substr(arg,2,7),
									resume= return(as.numeric(substr(arg,9,nchar(arg)))),NA)	}))
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
	
	f.name<- paste(dir.name,paste("tpcdat_",m.type,"_",loc.type,"_n",m.popsize,"_rI",theta[1],"_b",theta[2],"_sIdx",sample.prob[1],"_sE",sample.prob[2],"_tw",cluster.tw,"_median",sep=''),sep='/')
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
		sum.attack	<- summary( sapply(seq_along(tpc), function(i) tpc[[i]][["tpc.internal"]][["attack.rate"]]) )	
		sum.E2E		<- summary( sapply(seq_along(tpc), function(i)
						{
							tmp	<- tpc.proportion.E2E(tpc[[i]][["tpc.internal"]])
							tmp["E2E"] / tmp["X2E"]
						}) )		
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
		#print(tpc.table.all.median)
		
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
		#print(tpc.table.sample.median)
				
		if(verbose)
			cat(paste("\nprj.simudata: write tpc data medians to file",paste(f.name,".R",sep='')))
		if(save)
			save(tpc.table.all.median,tpc.table.sample.median,sum.attack,sum.E2E,file=paste(f.name,".R",sep=''))				
	}
	if(verbose)
	{
		print(sum.attack)
		print(sum.E2E)
	}
	ans<- list(tpc=tpc, tpc.table.all.median=tpc.table.all.median, tpc.table.sample.median=tpc.table.sample.median, sum.attack=sum.attack, sum.E2E=sum.E2E )
	ans
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
	prior.acute			<- c(0.5,13)
	prior.base			<- c(0.01, 0.13)		
	sample.prob			<- 1	
	debug.susc.const	<- 0
	debug.only.u		<- 0
	target.su.INC		<- 0.01
	target.su.E2E		<- 0.1
	abc.nit				<- 5e4
	abc.cores			<- 8
	dir.name			<- paste("acutesimu_fxs",debug.susc.const,"_onlyu",debug.only.u,sep='')
	
	if(exists("args"))
	{
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,4),
									loc= return(substr(arg,6,nchar(arg))),NA)	}))
		if(length(tmp)>0) loc.type<- tmp[1]
		tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,4),
									nit= return(as.numeric(substr(arg,6,nchar(arg)))),NA)	}))
		if(length(tmp)>0) abc.nit<- tmp[1]
	}
	if(verbose)
	{
		cat(paste("\nabc.nit ",abc.nit))
		cat(paste("\nloc.type ",loc.type))		
	}
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
		print(abc.struct)
		#abc.struct	<- mclapply(seq_len(nrow(abc.struct)),function(i)
		abc.struct	<- lapply(seq_len(nrow(abc.struct)),function(i)
				{		
					print(unlist(c(i,abc.struct[i,])))
					args<<- prj.simudata.cmd(CODE.HOME, loc.type, abc.struct[i,acute], abc.struct[i,base], 1, sample.prob, sample.prob, cluster.tw, save=0, resume=0, verbose=0, debug.susc.const=debug.susc.const, debug.only.u=debug.only.u)
					args<<- unlist(strsplit(args,' '))
					tpc	<- prj.simudata()[["tpc"]]
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

