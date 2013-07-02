###############################################################################
#' compute the power in differentiating H0: acute to acute is <10\% vs H0: acute to acute is >40\% with the linkage approach
#' 
#' This script varies the probability that HIV+ individuals consent to participation in the phylogenetic study at HCC
#' @param p.phylosignal				Probability that a true (potentially indirect) transmission is identified with a phylogenetic method. 
#' @param p.nocontam				Frequency with which transmissions occur within the community
#' @param p.prev.instudy.clu.armC	Probability that HIV+ individuals visit an HCC in arm C (sensitivity analysis).
#' @param opt.pooled				Pooling option for power analysis
#' @param opt.sampling				Sampling option for trial
prj.popart.powercalc_link_consenting<- function(p.phylosignal=0.7,p.nocontam=0.85, p.prev.instudy.clu.armC= 0.2, opt.pooled= "no pooling", opt.sampling= "PC and HCC")
{
	require(binom)
	require(phylodesign)
	my.mkdir(DATA,"popartpowercalc")
	dir.name<- paste(DATA,"popartpowercalc",sep='/')	
	resume<- 0
	verbose<- 1
	plot.increment<- 0.05
	
	m.type		<- "Acute"	
	cohort.size	<- 2500
	pc24.size	<- 6000
	cohort.dur	<- 3	
	theta.EE.H0	<- 0.10
	theta.EE.H1	<- 0.4
	test.alpha	<- 0.05		
	pooled.n	<- 1
	opt.pooled	<- "no pooling"#"pooled across trial"#"no pooling"#"pooled across ZA"#"pooled across trial"#"no pooling"
	#opt.sampling<- "PC and HCC"#"only HCC"#"PC and HCC"
	opt.power<-	"All"	
	if(!opt.pooled%in%c("pooled across country","pooled across ZA","pooled across SA","pooled across trial","no pooling"))
		stop("prj.popart.powercalc_medsize: unknown method to pool")
	if(!opt.power%in%c("All","PonlyPC","IonlyPC","PonlyPCandIonlyPC"))
		stop("prj.popart.powercalc_medsize: unknown method to pool")
	cat(paste("\ncohort.size",cohort.size))
	cat(paste("\ncohort.dur",cohort.dur))
	cat(paste("\ntheta.EE.H0",theta.EE.H0))
	cat(paste("\ntheta.EE.H1",theta.EE.H1))
	cat(paste("\ntest.alpha",test.alpha))	
	cat(paste("\np.nocontam",p.nocontam))
	cat(paste("\np.prev.instudy.clu.armC",p.prev.instudy.clu.armC))
	cat(paste("\npooled.n",pooled.n))
	cat(paste("\nopt.pooled",opt.pooled))
	cat(paste("\nopt.power",opt.power))
	cat(paste("\nopt.sampling",opt.sampling))
	
	sites<-	popart.getdata.randomized.arm( pooled.n )
	sites<-	popart.getdata.randomized.n(sites, cohort.size, cohort.dur)
	print(sites)
	###############################################################################
	#vary %consenting in HCC
	###############################################################################
	cat("\ncompute sampled transmissions for %consenting")
	s.consent<- seq(0.5,0.75,plot.increment)		
	x2i<- sapply(s.consent,function(x)
			{
				p.lab			<- 0.9			
				p.consent.coh	<- 0.9
				p.consent.clu	<- x
				p.vhcc.prev.AB<- 0.95
				p.vhcc.inc.AB<- 0.8
				p.vhcc.prev.C<- p.prev.instudy.clu.armC
				p.vhcc.inc.C<- p.vhcc.prev.C/2
				
				popart.get.sampled.transmissions(	sites, 
						opt.sampling,
						p.vhcc.prev.AB=p.vhcc.prev.AB, p.vhcc.prev.C=p.vhcc.prev.C, p.vhcc.inc.AB=p.vhcc.inc.AB, p.vhcc.inc.C=p.vhcc.inc.C, 
						consent.PC=p.consent.coh, consent.HCC=p.consent.clu, p.lab=p.lab, p.community=p.nocontam
				)																				
			})
	colnames(x2i)<- s.consent
	
	cat("\ncomputed linked transmissions")
	x2i<- phdes.get.linked.transmissions(x2i,p.phylosignal)
#print(x2i)
#stop()
	cat("\npool linked transmissions")
	print(sites)
	tmp		<- popart.pool(sites, x2i, method=opt.pooled)
	x2i		<- tmp[["transm"]]
	idx.A	<- tmp[["idx.A"]]
	idx.B	<- tmp[["idx.B"]]
	idx.C	<- tmp[["idx.C"]]
#print(x2i); print(idx.A); print(idx.B); print(idx.C)	
#stop()	
	cat("\ncompute power for dectecting differences in sampled acute to acute transmissions")
	confint.lw			<- lapply(	list(idx.A,idx.B,idx.C),	
									function(arm)	phdes.binom.power(	x2i[arm,,drop=0], round(x2i[arm,,drop=0]*theta.EE.H0), theta.EE.H0, theta.EE.H1, test.alpha, verbose=0, method.pw="cloglog", method.ci="asymptotic")[["conf"]]		
									)
	names(confint.lw)	<- c("A","B","C")
	tmp					<- lapply(	list(idx.A,idx.B,idx.C),
									function(arm)	phdes.binom.power(	x2i[arm,,drop=0], round(x2i[arm,,drop=0]*theta.EE.H1), theta.EE.H0, theta.EE.H1, test.alpha, verbose=0, method.pw="cloglog", method.ci="asymptotic")
									)
	confint.hg			<- lapply(tmp, function(x) x$conf	)							
	is.conf.hg			<- lapply(tmp, function(x) x$is.conf	)
	power.hg			<- lapply(tmp, function(x) x$power	)
	names(confint.hg)	<- c("A","B","C")
	names(is.conf.hg)	<- c("A","B","C")
	names(power.hg)		<- c("A","B","C")
	
	
	if(p.nocontam>=0.85)	
		cols<- c("deepskyblue","dodgerblue4")
	else					
		cols<- c("firebrick1","firebrick4")
	
	#plot power
	f.name<- paste(dir.name,paste("CFLINK_consent",p.nocontam,"phsig",p.phylosignal,"power",opt.power,"pool",opt.pooled,"sample",opt.sampling,"pwcalc",theta.EE.H0,theta.EE.H1,test.alpha,".pdf",sep='_'),sep='/')	
	cat(paste("\nplot power to\n",f.name))
	pdf(paste(f.name),version="1.4",width=6,height=6)
	phdes.plot.power(	power.hg[["A"]], power.hg[["C"]], is.conf.hg[["A"]], is.conf.hg[["C"]],
						xlab="% consenting to ph study at HCC", 
						ylab=paste("power to distinguish acute < ",theta.EE.H0*100,"% vs > ",theta.EE.H1*100,"%\n",opt.pooled,sep=''),
						legend.txt=c("arm A","arm C",paste("contamination",(1-p.nocontam)*100,"%, linked to HCC/C",p.prev.instudy.clu.armC*100,"%")),
						cols=cols,
						legend.loc="bottomright", 
						verbose= 0	)
	dev.off()
	
	#plot conf intervals
	lapply(names(confint.hg),function(arm)
			{
				f.name<- paste(dir.name,paste("CFLINK_consent",p.nocontam,"phsig",p.phylosignal,"power",opt.power,"pool",opt.pooled,"sample",opt.sampling,arm,"confint",theta.EE.H0,theta.EE.H1,test.alpha,".pdf",sep='_'),sep='/')
				cat(paste("\nplot binomial confidence intervals to\n",f.name))
				pdf(paste(f.name),version="1.4",width=6,height=6)
				phdes.plot.confint(	rep(theta.EE.H0,nrow(confint.lw[[arm]])), rep(theta.EE.H1,nrow(confint.lw[[arm]])),
						confint.lw[[arm]], confint.hg[[arm]],			
						xlab="% consenting to ph study at HCC",
						ylab=paste("estimated proportion acute to acute transmission\n arm",arm,",",opt.pooled),
						legend.loc="topright",
						legend.txt=c(paste("true prop",theta.EE.H0*100,"%",sep=' '), paste("true prop",theta.EE.H1*100,"%",sep=' ')),
						cols=cols		)													
				dev.off()
			})
		
	f.name<- paste(dir.name,paste("CFLINK_consent",p.nocontam,"phsig",p.phylosignal,"power",opt.power,"pool",opt.pooled,"sample",opt.sampling,"pwcalc",theta.EE.H0,theta.EE.H1,test.alpha,".R",sep='_'),sep='/')
	cat(paste("\nsave R objects to\n",f.name))
	save(sites, x2i, idx.A, idx.B, idx.C, power.hg, is.conf.hg, confint.lw, confint.hg, file=f.name)
	stop()													
}
###############################################################################
#' compuare the power in differentiating H0: acute to acute is <10\% vs H0: acute to acute is >40\% between the tipc cluster approach and the linkage approach
#' 
#' This script varies the probability that a true (potentially indirect) transmission is identified with a phylogenetic method.
#' @param p.nocontam				Frequency with which transmissions occur within the community
#' @param p.prev.instudy.clu.armC	Probability that HIV+ individuals visit an HCC in arm C (sensitivity analysis).
#' @param opt.pooled				Pooling option for power analysis
#' @param opt.sampling				Sampling option for trial
prj.popart.powercalc_cmp_link_tipc<- function(p.nocontam=0.85, p.prev.instudy.clu.armC= 0.2, opt.pooled= "no pooling", opt.sampling= "PC and HCC")
{
	require(binom)
	require(phylodesign)
	my.mkdir(DATA,"popartpowercalc")
	dir.name<- paste(DATA,"popartpowercalc",sep='/')	
	resume<- 0
	verbose<- 1
	plot.increment<- 0.05
	
	m.type<- "Acute"	
	cohort.size	<- 2500
	pc24.size	<- 6000
	cohort.dur	<- 3	
	theta.EE.H0	<- 0.1
	theta.EE.H1	<- 0.4
	test.alpha	<- 0.05		
	
	pooled.n<- 500
	#opt.pooled<- "pooled across ZA"#"pooled across trial"#"no pooling"
	#opt.sampling<- "PC and HCC"#"only HCC"#"PC and HCC"
	opt.power<-	"All"
	###############################################################################
	#set arguments
	###############################################################################	
	if(!opt.pooled%in%c("pooled across country","pooled across ZA","pooled across SA","pooled across trial","no pooling"))
		stop("prj.popart.powercalc_medsize: unknown method to pool")
	if(!opt.power%in%c("All","PonlyPC","IonlyPC","PonlyPCandIonlyPC"))
		stop("prj.popart.powercalc_medsize: unknown method to pool")
	debug<- 0
	
	sites<-	popart.getdata.randomized.arm( pooled.n, rtn.fixed= debug )
	sites<-	popart.getdata.randomized.n(sites, cohort.size, cohort.dur, rtn.exp= debug)
	#print(sites)
	data("popart.tipcp.highacute")
	data("popart.tipcp.lowacute")	
	#fixed parameters here:
	p.lab			<- 0.9			
	p.consent.coh	<- 0.9
	p.consent.clu	<- 0.5
	p.vhcc.prev.AB	<- 0.95
	p.vhcc.inc.AB	<- 0.8
	p.vhcc.prev.C	<- p.prev.instudy.clu.armC
	p.vhcc.inc.C	<- p.vhcc.prev.C/2	
	if(0)
	{
		p.prev.instudy.clu.armC <- 1
		p.lab			<- 1			
		p.consent.coh	<- 1
		p.consent.clu	<- 1
		p.vhcc.prev.AB	<- 1
		p.vhcc.inc.AB	<- 1
		p.vhcc.prev.C	<- 1
		p.vhcc.inc.C	<- 1
		p.nocontam		<- 1
	}
	
	cat(paste("\ncohort.size",cohort.size))
	cat(paste("\ncohort.dur",cohort.dur))
	cat(paste("\ntheta.EE.H0",theta.EE.H0))
	cat(paste("\ntheta.EE.H1",theta.EE.H1))
	cat(paste("\ntest.alpha",test.alpha))	
	cat(paste("\np.nocontam",p.nocontam))
	cat(paste("\np.prev.instudy.clu.armC",p.prev.instudy.clu.armC))
	cat(paste("\npooled.n",pooled.n))
	cat(paste("\nopt.pooled",opt.pooled))
	cat(paste("\nopt.power",opt.power))
	cat(paste("\nopt.sampling",opt.sampling))	
	###############################################################################
	#vary p.phylosignal with link approach
	###############################################################################
	cat("\ncompute sampled transmissions for p.phylosignal with link approach")	
	p.phylosignal<- seq(0.3,1,plot.increment)		
	x2i<- sapply(p.phylosignal,function(x)
			{				
				x2i<- popart.get.sampled.transmissions(	sites, 
														opt.sampling,
														p.vhcc.prev.AB=p.vhcc.prev.AB, p.vhcc.prev.C=p.vhcc.prev.C, p.vhcc.inc.AB=p.vhcc.inc.AB, p.vhcc.inc.C=p.vhcc.inc.C, 
														consent.PC=p.consent.coh, consent.HCC=p.consent.clu, p.lab=p.lab, p.community=p.nocontam,
														rtn.int=!debug
														)
				x2i<- phdes.get.linked.transmissions(x2i,x, rtn.exp=debug)
				x2i
			})
	colnames(x2i)<- p.phylosignal
	###############################################################################
	#compute acute to acute with tipc approach
	###############################################################################
	cat("\ncompute sampled acute to acute transmissions with tipc approach")
	tipc.p			<- phdes.get.hyp.tipc.probs(popart.tipcp.highacute,1-p.nocontam)
	a2a.hg			<- popart.get.sampled.acute2acute( 	sites, tipc.p, opt.sampling, 
														p.vhcc.prev.AB=p.vhcc.prev.AB, p.vhcc.prev.C=p.vhcc.prev.C, p.vhcc.inc.AB=p.vhcc.inc.AB, p.vhcc.inc.C=p.vhcc.inc.C, 
														consent.PC=p.consent.coh, consent.HCC=p.consent.clu, p.lab=p.lab,
														rtn.int= !debug)
	tipc.p			<- phdes.get.hyp.tipc.probs(popart.tipcp.lowacute,1-p.nocontam)
	a2a.lw			<- popart.get.sampled.acute2acute( 	sites, tipc.p, opt.sampling, 
														p.vhcc.prev.AB=p.vhcc.prev.AB, p.vhcc.prev.C=p.vhcc.prev.C, p.vhcc.inc.AB=p.vhcc.inc.AB, p.vhcc.inc.C=p.vhcc.inc.C, 
														consent.PC=p.consent.coh, consent.HCC=p.consent.clu, p.lab=p.lab,
														rtn.int= !debug)																
	x2i.lw	<- as.matrix(a2a.lw["x2i.s",])
	i2i.lw	<- as.matrix(a2a.lw["i2i.s",])	
	x2i.hg	<- as.matrix(a2a.hg["x2i.s",])
	i2i.hg	<- as.matrix(a2a.hg["i2i.s",])
	###############################################################################
	#pool transmissions
	###############################################################################	
	cat(paste("\npool transmissions",opt.pooled))
	x2i.lw						<- 		popart.pool(sites, x2i.lw, method=opt.pooled)[["transm"]]
	i2i.lw						<- 		popart.pool(sites, i2i.lw, method=opt.pooled)[["transm"]]
	x2i.hg						<- 		popart.pool(sites, x2i.hg, method=opt.pooled)[["transm"]]
	i2i.hg						<- 		popart.pool(sites, i2i.hg, method=opt.pooled)[["transm"]]		
	g(x2i, idx.A, idx.B, idx.C)	%<-%	popart.pool(sites, x2i, method=opt.pooled)
	
	#compute test.biased.H0 as fraction over means to see a pattern	
	test.biased.H0		<- apply(i2i.lw,2,mean)/apply(x2i.lw,2,mean)	
	test.biased.H1		<- apply(i2i.hg,2,mean)/apply(x2i.hg,2,mean)	
	test.biased.H1.arm	<- lapply(	list(idx.A, idx.B, idx.C),
									function(arm)	apply(i2i.hg[arm,,drop=0],2,mean)/apply(x2i.hg[arm,,drop=0],2,mean)
									)				
	test.biased.H0.arm	<- lapply(	list(idx.A, idx.B, idx.C),
									function(arm)	apply(i2i.lw[arm,,drop=0],2,mean)/apply(x2i.lw[arm,,drop=0],2,mean)
									)
	names(test.biased.H0.arm)<- c("A","B","C")
	names(test.biased.H1.arm)<- c("A","B","C")										
	###############################################################################
	#plot number of acute to acute transmissions under high scenario
	#adjust for sampling bias
	###############################################################################
	cat("\nplot number of acute to acute transmissions under high scenario\nadjust for sampling bias")
	f.name	<- paste(dir.name,paste("VARYLINKAGE_LINK_a2a_nocontam",p.nocontam,"power",opt.power,"pool",opt.pooled,"sample",opt.sampling,"pwcalc",theta.EE.H0,theta.EE.H1,test.alpha,".pdf",sep='_'),sep='/')
	cat(paste("\nplot to\n",f.name))
	pdf(paste(f.name),version="1.4",width=6,height=6)
	par(mar=c(5,5.5,0.5,2))
	ylim	<- range(sapply(test.biased.H1.arm, function(x)	range(round(x2i*x)) ))
	plot(1,1,type='n',xlim=range(p.phylosignal),ylim=ylim,xlab="% linked w phylogenetic method",ylab=paste("acute 2 acute, high scenario\n",opt.pooled))
	lines(p.phylosignal, apply(round(x2i*test.biased.H1),2,mean))
	abline(h=apply(i2i.hg,2,mean), lty=2)		
	lines(p.phylosignal, apply(round(x2i[idx.A,,drop=0]*test.biased.H1.arm[["A"]]),2,mean),col="red")
	abline(h=apply(i2i.hg[idx.A,,drop=0],2,mean), lty=2,col="red")
	lines(p.phylosignal, apply(round(x2i[idx.B,,drop=0]*test.biased.H1.arm[["B"]]),2,mean),col="green")
	abline(h=apply(i2i.hg[idx.B,,drop=0],2,mean), lty=2,col="green")	
	lines(p.phylosignal, apply(round(x2i[idx.C,,drop=0]*test.biased.H1.arm[["C"]]),2,mean),col="blue")
	abline(h=apply(i2i.hg[idx.C,,drop=0],2,mean), lty=2,col="blue")
	legend("topleft",fill=c("black","red","green","blue"),legend=c("overall","arm A","arm B","arm C"),bty='n',border=NA)
	legend("topright",lty=c(1,2),legend=c("linked","tipc"),bty='n')
	dev.off()
	
	stop()													
}
###########################################################################
prj.popart.powercalc_tipc_test_ukhivrdb<- function(p.phylosignal=0.7,p.nocontam=1, opt.pooled= "no pooling", opt.sampling= "PC and HCC")
{
	require(binom)
	require(phylodesign)
	my.mkdir(DATA,"popartpowercalc_test")
	dir.name<- paste(DATA,"popartpowercalc_test",sep='/')	
	resume<- 0
	verbose<- 0
	plot.increment<- 0.05
	
	m.type		<- "Acute"	
	cohort.size	<- 0
	cohort.dur	<- 5	
	theta.EE.H0	<- 0.10
	theta.EE.H1	<- 0.4
	theta.UE	<- 0.3
	theta.TE	<- theta.UE / 5
	test.alpha	<- 0.05		 
	debug		<- 1
	pooled.n	<- 1
	opt.pooled	<- "no pooling"#"pooled across ZA"#"pooled across trial"#"no pooling"
	opt.pooled	<- "pooled across SA"
	opt.clu.closure	<- 14
	opt.sampling<- "only HCC"
	#opt.sampling<- "PC after yr 1 and HCC"
	#opt.sampling<- "PC only incident and HCC"
	opt.power	<-	"All"
	
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
	
	sites			<- popart.getdata.randomized.arm( pooled.n, rtn.fixed=debug )
	sites$arm		<- 'C'
	sites$inc.rate	<- seq(1.6,2.3,len=nrow(sites))/100
	sites$p.adults	<- 1
	sites$hivcomb	<- 20
	sites$hivsero	<- 20
	sites$artadjust	<- 85
	sites$artcrude	<- 85
	sites$popsize	<- 6e3 * 0.8 / (1-seq(0.2,0.3,len=nrow(sites))) / 0.2  	#number in UK HIV RDB is 6k and 80% is MSM	assume %undiagnosed 20%-30%	-> total prevalence. prevalence 20%,-> total population		
			
	sites<-	popart.getdata.randomized.n(sites, cohort.size, cohort.dur, rtn.exp=debug)	
	print(sites)

	#total annual incidence between 315-600
	print(range(sites$n.inc/cohort.dur))
	
	#compute complete tip cluster probs under H0 and H1
	clu.n		<- clu.n.of.tchain(opt.clu.closure)	
	theta.H0	<- clu.p.init(theta.EE.H0, theta.UE, theta.TE, 1-p.nocontam)
	theta.H1	<- clu.p.init(theta.EE.H1, theta.UE, theta.TE, 1-p.nocontam)
	tipc.p.H0	<- clu.probabilities(clu.n, theta.H0, with.ntr.weight= 1)
	tipc.p.H1	<- clu.probabilities(clu.n, theta.H1, with.ntr.weight= 1)
	
	#for H1 (high E->E), get n(E->E) and n(x->E) for each arm under baseline scenario
	if(1)
	{
		p.lab			<- seq(0.1,0.7,0.05)		
		p.consent.coh	<- 0						#there is no cohort
		p.consent.clu	<- 1
		p.vhcc.prev.AB	<- 0
		p.vhcc.inc.AB	<- 0
		p.vhcc.prev.C	<- 0.75						#this setting gives 300-550 clusters for 30% sampling
		p.vhcc.inc.C	<- 0.75 / 1.5					#	
		p.baseline		<- 0.5						#sequences that cluster are 50%
		xlab			<- "coverage"
		
		NCLU<<- matrix(NA,0,2)
		ans			<- lapply(p.lab,function(x)
				{		
					cat(paste("\nprocess",x))
					#theta.H0	<- clu.p.init(theta.EE.H0, theta.UE, theta.TE, p.contam)
					#theta.H1	<- clu.p.init(theta.EE.H1, theta.UE, theta.TE, p.contam)
					#tipc.p.H0	<- clu.probabilities(clu.n, theta.H0, with.ntr.weight=1)
					#tipc.p.H1	<- clu.probabilities(clu.n, theta.H1, with.ntr.weight=1)
					sampling<- popart.sampling.init(sites,  p.consent.coh, p.consent.clu, x, p.vhcc.prev.AB, p.vhcc.inc.AB, p.vhcc.prev.C, p.vhcc.inc.C, method= opt.sampling)
					sampling$baseline<- p.baseline
					print(sampling)
					
					ntr.hg.H0<- popart.get.sampled.transmissions.from.tipc(sites, tipc.p.H0, clu.n, theta.H0, sampling, opt.sampling, mx.sampled.ntr=6, exclude.O= 1, rtn.int=!debug, verbose=0)
					ntr.hg.H1<- popart.get.sampled.transmissions.from.tipc(sites, tipc.p.H1, clu.n, theta.H1, sampling, opt.sampling, mx.sampled.ntr=6, exclude.O= 1, rtn.int=!debug, verbose=1)
					#print(ntr.hg)						
					list(i2i.s.H0= ntr.hg.H0["i2i.s",], x2i.s.H0=ntr.hg.H0["x2i.s",],i2i.s.H1= ntr.hg.H1["i2i.s",], x2i.s.H1=ntr.hg.H1["x2i.s",]) 
				})
		names(ans)<- p.lab	
		#print(ans)
		
		ans.x2E.H0<- sapply(seq_along(ans),function(i) ans[[i]][["x2i.s.H0"]])
		ans.E2E.H0<- sapply(seq_along(ans),function(i) ans[[i]][["i2i.s.H0"]])
		ans.x2E.H1<- sapply(seq_along(ans),function(i) ans[[i]][["x2i.s.H1"]])
		ans.E2E.H1<- sapply(seq_along(ans),function(i) ans[[i]][["i2i.s.H1"]])
		
		#print(ans.x2E.H1)		
		#print(ans.E2E.H1 / ans.x2E.H1 )
		arms<- c("C")
		ans	<- lapply(seq_along(arms),function(i)
				{
					cat(paste("\n arm ",arms[i]))
					#power with t-test per arm
					idx			<- sites$arm==arms[i]
					p.H0		<- ans.E2E.H0[idx,,drop=0] / ans.x2E.H0[idx,,drop=0]
					p.H1		<- ans.E2E.H1[idx,,drop=0] / ans.x2E.H1[idx,,drop=0]															
					
					#power with simple binom test, averaging over different assumptions on %diagnosed and cum true incidence/year
					p.H0		<- apply(ans.E2E.H0[idx,,drop=0] / ans.x2E.H0[idx,,drop=0],2,mean)
					p.H1		<- apply(ans.E2E.H1[idx,,drop=0] / ans.x2E.H1[idx,,drop=0],2,mean)
					e.H0		<- apply(floor(ans.E2E.H0[idx,,drop=0]),2,mean)
					e.H1		<- apply(floor(ans.E2E.H1[idx,,drop=0]),2,mean)
					n.H0		<- apply(floor(ans.x2E.H0[idx,,drop=0]),2,mean)
					n.H1		<- apply(floor(ans.x2E.H1[idx,,drop=0]),2,mean) 
					#print(c(n.H0,n.H1,e.H0,e.H1))
					
					pw.b		<- .phdes.binom.power(n.H1, p.H1, p.H0, test.alpha, method="asymp")
					names(pw.b)	<- names(ans)
					#print(pw.b)					
					#stop()
					list(pw.b=pw.b, p.H0.m=p.H0, p.H1.m=p.H1)
				})
		#plot computed power for each arm	
		pw.b<- sapply(seq_along(ans),function(i)
				{
					ans[[i]][["pw.b"]]
				})		
		colnames(pw.b)<- arms
		print(pw.b)		
		print(NCLU)

		
		
		clr<- c("red","blue","green")
		f.name	<- paste(dir.name,"/",xlab,"_power.pdf",sep='')
		cat(paste("\n plot to ",f.name))
		pdf(f.name,version="1.4",width=5,height=5)		
		plot(1,1,type='n',bty='n',xlim=range(p.lab),ylim=range(c(pw.b)),xlab=xlab,ylab="power")
		sapply(seq_len(ncol(pw.b)),function(j)
				{
					lines(p.lab,pw.b[,j])
				})		
		legend("bottomright",bty='n',legend=c("prob(acute/early)>40% vs <10%"))
		dev.off()
		
		f.name	<- paste(dir.name,"/",xlab,"_nclu.pdf",sep='')
		cat(paste("\n plot to ",f.name))
		pdf(f.name,version="1.4",width=5,height=5)		
		plot(1,1,type='n',bty='n',xlim=range(p.lab),ylim=range(c(NCLU[,2])),xlab=xlab,ylab="number of clusters")
		sapply(seq_len(ncol(pw.b)),function(j)
				{
					lines(p.lab,NCLU[,2])
				})				
		dev.off()				
	}
	stop()
}
###########################################################################
prj.popart.powercalc_tipc_test_residual<- function(p.phylosignal=0.7,p.nocontam=0.85, opt.pooled= "no pooling", opt.sampling= "PC and HCC")
{
	require(binom)
	require(phylodesign)
	my.mkdir(DATA,"popartpowercalc_test")
	dir.name<- paste(DATA,"popartpowercalc_test",sep='/')	
	resume<- 0
	verbose<- 0
	plot.increment<- 0.05
	
	m.type		<- "Acute"	
	cohort.size	<- 2500
	pc24.size	<- 6000
	cohort.dur	<- 3	
	theta.EE.H0	<- 0.10
	theta.EE.H1	<- 0.4
	theta.UE	<- 0.3
	theta.TE	<- theta.UE / 5
	test.alpha	<- 0.05		 
	debug		<- 1
	pooled.n	<- 1
	opt.pooled	<- "no pooling"#"pooled across ZA"#"pooled across trial"#"no pooling"
	opt.pooled	<- "pooled across SA"
	opt.clu.closure	<- 14
	opt.sampling<- "PC and HCC"#"only HCC"	#"PC and HCC"	#
	#opt.sampling<- "PC after yr 1 and HCC"
	#opt.sampling<- "PC only incident and HCC"
	opt.power	<-	"All"
	
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
	
	sites<-	popart.getdata.randomized.arm( pooled.n, rtn.fixed=debug, rtn.phylostudy=1 )
	sites<-	popart.getdata.randomized.n(sites, cohort.size, cohort.dur, rtn.exp=debug)	
	print(sites)
	 
	#compute complete tip cluster probs under H0 and H1
	#clu.n		<- clu.n.of.tchain(opt.clu.closure)	
	#theta.H0	<- clu.p.init(theta.EE.H0, theta.UE, theta.TE, 1-p.nocontam)
	#theta.H1	<- clu.p.init(theta.EE.H1, theta.UE, theta.TE, 1-p.nocontam)
	#tipc.p.H0	<- clu.probabilities(clu.n, theta.H0, with.ntr.weight= 1)
	#tipc.p.H1	<- clu.probabilities(clu.n, theta.H1, with.ntr.weight= 1)
	
	#for H1 (high E->E), get n(E->E) and n(x->E) for each arm under residual sampling scenario
	if(1)
	{
		p.lab			<- 0.7*0.9			#set lower as discussed		
		p.consent.coh	<- 0.9				
		p.consent.clu	<- 1
		p.vhcc.prev.AB	<- 1				#already in PopART model estimate
		p.vhcc.inc.AB	<- 1				#already in PopART model estimate
		p.vhcc.prev.C	<- 1				#already in PopART model estimate
		p.vhcc.inc.C	<- 1				#already in PopART model estimate		
		p.contam		<- seq(0.05,0.2,0.025)
		opt.sampling	<- "PC and HCC"	
		opt.sampling	<- "PC after yr 1 and HCC"		
		ans			<- lapply(p.contam,function(x)
				{
					sampling<- popart.sampling.init.PopARTmodel(sites,  p.consent.coh, p.consent.clu, p.lab, p.vhcc.prev.AB, p.vhcc.inc.AB, p.vhcc.prev.C, p.vhcc.inc.C, method= opt.sampling)
					print(sampling[c("Ndeke","Chimwemwe","Ngungu","Maramba","Dambwa","Shampande","Kuyasa","Luvuyo","Town II","Ikhwezi","Bloekombos","Delft South"),])
					print(sampling[,"total prev"]*sampling[,"total inc"]*(1-x))
					stop()
										
					theta.H0	<- clu.p.init(theta.EE.H0, theta.UE, theta.TE, x)
					theta.H1	<- clu.p.init(theta.EE.H1, theta.UE, theta.TE, x)
					tipc.p.H0	<- clu.probabilities(clu.n, theta.H0, with.ntr.weight=1)
					tipc.p.H1	<- clu.probabilities(clu.n, theta.H1, with.ntr.weight=1)
					
					#sampling<- popart.sampling.init(sites,  p.consent.coh, p.consent.clu, p.lab, p.vhcc.prev.AB, p.vhcc.inc.AB, p.vhcc.prev.C, p.vhcc.inc.C, method= opt.sampling)
					#print(xtable(sampling, digits=2), floating=FALSE)
					ntr.hg<- popart.get.sampled.transmissions.from.tipc(sites, tipc.p.H1, clu.n, theta.H1, sampling, opt.sampling, mx.sampled.ntr=6, exclude.O= 1, rtn.int=!debug)
					#print(ntr.hg)						
					list(i2i.s= ntr.hg["i2i.s",], x2i.s=ntr.hg["x2i.s",]) 
				})
		names(ans)<- p.contam	
		#print(ans)
		ans.x2E<- sapply(seq_along(ans),function(i) ans[[i]][["x2i.s"]])
		ans.E2E<- sapply(seq_along(ans),function(i) ans[[i]][["i2i.s"]])
		
		print(ans.E2E)
		clr<- c("red","blue","green")
		arms<- c("A","B","C")
		p.vary<- as.numeric(names(ans))
		xlab<- "proportion transmission from outside cluster"
		sapply(seq_along(arms),function(j)
				{
					plotmat	<- ans.E2E[sites$arm==arms[j],]
					ylim	<- c(0,max( sapply(seq_along(arms),function(k) apply( ans.E2E[sites$arm==arms[k],],2,sum ) ) ))
					
					f.name	<- paste(dir.name,"/residual_",opt.sampling,p.vhcc.prev.C,p.vhcc.inc.C,p.lab,xlab,"_vs_E2E_arm",arms[j],".pdf",sep='')
					cat(paste("\n plot to ",f.name))
					pdf(f.name,version="1.4",width=5,height=5)
					plot(1,1,type='n',bty='n',xlim=range(p.vary),ylim=ylim,xlab=xlab,ylab="n(E->E)")
					z		<- rep(0,ncol(plotmat))
					cols	<- colorRampPalette(c(clr[j],"yellow"))(nrow(plotmat))
					for(i in seq_len(nrow(plotmat)))
					{
						polygon( c(p.vary,rev(p.vary)), c(z+plotmat[i,], rev(z)), col=cols[i] )
						z	<- z+plotmat[i,]
					}										
					legend("topleft",bty='n',legend=paste("arm",arms[j]))
					dev.off()
				})
		clr<- c("black","black","black")
		arms<- c("A","B","C")
		sapply(seq_along(arms),function(j)
				{
					plotmat	<- ans.x2E[sites$arm==arms[j],]
					ylim	<- c(0,max( sapply(seq_along(arms),function(k) apply( ans.x2E[sites$arm==arms[k],],2,sum ) ) ))
					
					f.name	<- paste(dir.name,"/residual_",opt.sampling,p.vhcc.prev.C,p.vhcc.inc.C,p.lab,xlab,"_vs_x2E_arm",arms[j],".pdf",sep='')
					cat(paste("\n plot to ",f.name))
					pdf(f.name,version="1.4",width=5,height=5)
					plot(1,1,type='n',bty='n',xlim=range(p.vary),ylim=ylim,xlab=xlab,ylab="n(E->E)")
					z		<- rep(0,ncol(plotmat))
					cols	<- colorRampPalette(c(clr[j],"gray50"))(nrow(plotmat))
					for(i in seq_len(nrow(plotmat)))
					{
						polygon( c(p.vary,rev(p.vary)), c(z+plotmat[i,], rev(z)), col=cols[i] )
						z	<- z+plotmat[i,]
					}										
					legend("topleft",bty='n',legend=paste("arm",arms[j]))
					dev.off()
				})
	}
	#if samples available from multiple visits to HCC, how would this change p.lab ?
	if(0)
	{
		p.lab<- 0.5
		n<- seq.int(1,3)
		print( 1-pbinom(0,n,p.lab) )
		p.lab<- 0.6
		print( 1-pbinom(0,n,p.lab) )
	}	
	#for H1 (high E->E), get n(E->E) and n(x->E) for each arm under residual baseline scenario
	if(0)
	{
		p.lab			<- 0.7		
		p.consent.coh	<- 0.9
		p.consent.clu	<- 1
		p.vhcc.prev.AB	<- 0.95
		p.vhcc.inc.AB	<- 0.8
		p.vhcc.prev.C	<- 0.4
		p.vhcc.inc.C	<- p.vhcc.prev.C/2					
		p.contam		<- seq(0.05,0.2,0.025)
		opt.sampling	<- "PC and HCC"	
		#opt.sampling	<- "only HCC"
		xlab			<- "proportion of transmissions from outside cluster"
		
		ans			<- lapply(p.contam,function(x)
				{
					theta.H0	<- clu.p.init(theta.EE.H0, theta.UE, theta.TE, x)
					theta.H1	<- clu.p.init(theta.EE.H1, theta.UE, theta.TE, x)
					tipc.p.H0	<- clu.probabilities(clu.n, theta.H0, with.ntr.weight=1)
					tipc.p.H1	<- clu.probabilities(clu.n, theta.H1, with.ntr.weight=1)
					sampling<- popart.sampling.init(sites,  p.consent.coh, p.consent.clu, p.lab, p.vhcc.prev.AB, p.vhcc.inc.AB, p.vhcc.prev.C, p.vhcc.inc.C, method= opt.sampling)
					#print(sampling)
					ntr.hg.H0<- popart.get.sampled.transmissions.from.tipc(sites, tipc.p.H0, clu.n, theta.H0, sampling, opt.sampling, mx.sampled.ntr=6, exclude.O= 1, rtn.int=!debug)
					ntr.hg.H1<- popart.get.sampled.transmissions.from.tipc(sites, tipc.p.H1, clu.n, theta.H1, sampling, opt.sampling, mx.sampled.ntr=6, exclude.O= 1, rtn.int=!debug)
					#print(ntr.hg)						
					list(i2i.s.H0= ntr.hg.H0["i2i.s",], x2i.s.H0=ntr.hg.H0["x2i.s",],i2i.s.H1= ntr.hg.H1["i2i.s",], x2i.s.H1=ntr.hg.H1["x2i.s",]) 
				})
		names(ans)<- p.contam	
		#print(ans)
		ans.x2E.H0<- sapply(seq_along(ans),function(i) ans[[i]][["x2i.s.H0"]])
		ans.E2E.H0<- sapply(seq_along(ans),function(i) ans[[i]][["i2i.s.H0"]])
		ans.x2E.H1<- sapply(seq_along(ans),function(i) ans[[i]][["x2i.s.H1"]])
		ans.E2E.H1<- sapply(seq_along(ans),function(i) ans[[i]][["i2i.s.H1"]])
		
		#print(ans.x2E.H1)		
		#print(ans.E2E.H1 / ans.x2E.H1 )
		arms<- c("A","B","C")
		ans	<- lapply(seq_along(arms),function(i)
				{
					cat(paste("\n arm ",arms[i]))
					#power with t-test per arm
					idx			<- sites$arm==arms[i]
					p.H0		<- ans.E2E.H0[idx,] / ans.x2E.H0[idx,]
					print(p.H0)
					p.H1		<- ans.E2E.H1[idx,] / ans.x2E.H1[idx,]					
					var.pH0		<- (apply(p.H0,2,mean)*0.5)^2 	#10*10*apply(ans.E2E.H0[idx,] / ans.x2E.H0[idx,],2,var)
					var.pH1		<- (apply(p.H1,2,mean)*0.5)^2	#10*10*apply(ans.E2E.H0[idx,] / ans.x2E.H0[idx,],2,var)
					
					pw			<- sapply(seq_along(ans),function(j)
							{
								n.H1	<- floor(ans.x2E.H1[idx,j])										
								phdes.power.ttest.cl(n.H1, p.H0[,j], p.H1[,j], var.pH0[j], var.pH1[j], alpha= test.alpha)
							})
					names(pw)	<- names(ans)
					#print(pw)
					
					#design effect: 1+(harmonic mean(n.H1)-1) * within clu corr
					des.effect	<- 1 + (1/apply(1/ans.x2E.H1[idx,],2,mean)-1) * var.pH1 / (apply(p.H1,2,mean)*(1-apply(p.H1,2,mean)))
					des.effect	<- 1 + (apply(ans.x2E.H1[idx,],2,mean)-1) * var.pH1 / (apply(p.H1,2,mean)*(1-apply(p.H1,2,mean)))
					#print(des.effect)
					
					#power with simple binom test, adjusting for design effect
					p.H0		<- apply(ans.E2E.H0[idx,] / ans.x2E.H0[idx,],2,mean)
					p.H1		<- apply(ans.E2E.H1[idx,] / ans.x2E.H1[idx,],2,mean)
					e.H0		<- apply(floor(ans.E2E.H0[idx,]),2,sum) / des.effect
					e.H1		<- apply(floor(ans.E2E.H1[idx,]),2,sum) / des.effect
					n.H0		<- apply(floor(ans.x2E.H0[idx,]),2,sum) / des.effect
					n.H1		<- apply(floor(ans.x2E.H1[idx,]),2,sum) / des.effect
					pw.b		<- .phdes.binom.power(n.H1, p.H1, p.H0, test.alpha, method="asymp")
					names(pw.b)	<- names(ans)
					#print(pw.b)
					
					#95% conf interval with simple binom test, adjusting for design effect
					conf.H1		<- binom.confint(e.H1, n.H1, conf.level = 0.95, methods="cloglog")[,c("lower","upper")]
					conf.H0		<- binom.confint(e.H0, n.H0, conf.level = 0.95, methods="cloglog")[,c("lower","upper")]
					rownames(conf.H0)<- rownames(conf.H1)<- names(ans)
					#print( conf.H0 ); print( conf.H1 )
					
					list(pw.n=pw, pw.b=pw.b,conf.b.H0=conf.H0, conf.b.H1=conf.H1, p.H0.m=p.H0, p.H1.m=p.H1)
				})
		print(ans)		
		#plot computed power for each arm	
		pw.n<- sapply(seq_along(ans),function(i)
				{
					ans[[i]][["pw.n"]]
				})		
		colnames(pw.n)<- arms
		pw.b<- sapply(seq_along(ans),function(i)
				{
					ans[[i]][["pw.b"]]
				})		
		colnames(pw.b)<- arms
		
		clr<- c("red","blue","green")
		f.name	<- paste(dir.name,"/residual_",opt.sampling,p.vhcc.prev.C,p.vhcc.inc.C,p.lab,xlab,"_power.pdf",sep='')
		cat(paste("\n plot to ",f.name))
		pdf(f.name,version="1.4",width=5,height=5)		
		plot(1,1,type='n',bty='n',xlim=range(p.contam),ylim=range(c(0,1,pw.b,pw.n)),xlab=xlab,ylab="power")
		sapply(seq_len(ncol(pw.n)),function(j)
				{
					#lines(p.contam,pw.n[,j],col=clr[j],lty=1)
				})
		sapply(seq_len(ncol(pw.b)),function(j)
				{
					lines(p.contam,pw.b[,j],col=clr[j],lty=2)
				})
		#legend("bottomleft",bty='n',lty=c(1,2),legend=c("normal approx (Hayes1999)","Binomial using effective size"))
		legend("bottomright",bty='n',fill=clr,legend=c("arm A","arm B","arm C"))
		dev.off()
		
		#plot confidence intervals for each arm
		sapply(seq_along(ans),function(i)
				{
					ylim<- range(c(ans[[i]][["conf.b.H0"]],ans[[i]][["conf.b.H1"]]))
					ylim<- ylim*c(1,1.3)
					
					f.name	<- paste(dir.name,"/residual_",xlab,"_confint_",arms[i],".pdf",sep='')
					cat(paste("\n plot to ",f.name))
					pdf(f.name,version="1.4",width=5,height=5)		
					
					plot(1,1,type='n',bty='n',xlim=range(p.contam),ylim=ylim,xlab=xlab,ylab="95% Binomial confidence interval")
					x<- ans[[i]][["conf.b.H0"]]
					#polygon( c(p.contam,rev(p.contam)), c(x[,1],rev(x[,2])), col=my.fade.col(clr[i],0.4), border=NA )
					polygon( c(p.contam,rev(p.contam)), c(x[,1],rev(x[,2])), col="gray50", border=NA )
					x<- ans[[i]][["conf.b.H1"]]
					polygon( c(p.contam,rev(p.contam)), c(x[,1],rev(x[,2])), col=my.fade.col(clr[i],0.7), border=NA )					
					lines(p.contam,ans[[i]][["p.H0.m"]], lty=2)
					lines(p.contam,ans[[i]][["p.H1.m"]], lty=3)
					legend("topleft",bty='n',legend=paste("arm",arms[i]))
					dev.off()
				})  
		
	}
}
###########################################################################
prj.popart.powercalc_tipc_test<- function(p.phylosignal=0.7,p.nocontam=0.85, opt.pooled= "no pooling", opt.sampling= "PC and HCC")
{
	require(binom)
	require(phylodesign)
	my.mkdir(DATA,"popartpowercalc_test")
	dir.name<- paste(DATA,"popartpowercalc_test",sep='/')	
	resume<- 0
	verbose<- 0
	plot.increment<- 0.05
	
	m.type		<- "Acute"	
	cohort.size	<- 2500
	pc24.size	<- 6000
	cohort.dur	<- 3	
	theta.EE.H0	<- 0.10
	theta.EE.H1	<- 0.4
	theta.UE	<- 0.3
	theta.TE	<- theta.UE / 5
	test.alpha	<- 0.05		 
	debug		<- 1
	pooled.n	<- 1
	opt.pooled	<- "no pooling"#"pooled across ZA"#"pooled across trial"#"no pooling"
	opt.pooled	<- "pooled across SA"
	opt.clu.closure	<- 14
	opt.sampling<- "PC and HCC"#"only HCC"	#"PC and HCC"	#
	#opt.sampling<- "PC after yr 1 and HCC"
	#opt.sampling<- "PC only incident and HCC"
	opt.power	<-	"All"
	
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
	
	sites<-	popart.getdata.randomized.arm( pooled.n, rtn.fixed=debug )
	sites<-	popart.getdata.randomized.n(sites, cohort.size, cohort.dur, rtn.exp=debug)	
	print(sites)

	#compute complete tip cluster probs under H0 and H1
	clu.n		<- clu.n.of.tchain(opt.clu.closure)	
	theta.H0	<- clu.p.init(theta.EE.H0, theta.UE, theta.TE, 1-p.nocontam)
	theta.H1	<- clu.p.init(theta.EE.H1, theta.UE, theta.TE, 1-p.nocontam)
	tipc.p.H0	<- clu.probabilities(clu.n, theta.H0, with.ntr.weight= 1)
	tipc.p.H1	<- clu.probabilities(clu.n, theta.H1, with.ntr.weight= 1)
	
	
	#for H1 (high E->E), get sampling biases
	if(0)
	{
		s.intensity	<- seq(0.3,1,0.1)
		ans			<- sapply(s.intensity,function(x)
						{
							p.lab			<- x		
							p.consent.coh	<- p.consent.clu<- p.vhcc.prev.AB	<- p.vhcc.inc.AB	<- p.vhcc.prev.C	<- p.vhcc.inc.C	<- 1	
							sampling<- popart.sampling.init(sites,  p.consent.coh, p.consent.clu, p.lab, p.vhcc.prev.AB, p.vhcc.inc.AB, p.vhcc.prev.C, p.vhcc.inc.C, method= opt.sampling)
							sampling$baseline<- x
							#print(sampling)
							ntr.hg<- popart.get.sampled.transmissions.from.tipc(sites, tipc.p.H1, clu.n, theta.H1, sampling, opt.sampling, mx.sampled.ntr=6, exclude.O= 1, rtn.int=!debug)
							#print(ntr.hg)						
							ntr.hg["i2i.s",] / ntr.hg["x2i.s",] 
						})
		colnames(ans)<- s.intensity		
		
		f.name	<- paste(dir.name,"/propH1_tipcluster.pdf",sep='')
		cat(paste("\n plot to ",f.name))
		pdf(f.name,version="1.4",width=5,height=5)		
		plot(1,1,type='n',bty='n',xlim=range(s.intensity),ylim=range(c(tipc.p.H1,ans)),xlab="sampling intensity",ylab="p(E->E)")
		lines(s.intensity, ans[1,])
		abline(h=theta.H1["E2E"], lty=2)
		legend("bottomright",bty='n',lty=c(2,1),legend=c("true proportion","estimate from tip clusters"))
		dev.off()
	}
	
	#compare delta for power analyses under sampling bias
	if(0)
	{
		s.intensity	<- seq(0.3,1,0.1)
		ans			<- sapply(s.intensity,function(x)
				{
					p.lab			<- x		
					p.consent.coh	<- p.consent.clu<- p.vhcc.prev.AB	<- p.vhcc.inc.AB	<- p.vhcc.prev.C	<- p.vhcc.inc.C	<- 1	
					sampling<- popart.sampling.init(sites,  p.consent.coh, p.consent.clu, p.lab, p.vhcc.prev.AB, p.vhcc.inc.AB, p.vhcc.prev.C, p.vhcc.inc.C, method= opt.sampling)
					sampling$baseline<- x
					#print(sampling)
					ntr.hg.H0<- popart.get.sampled.transmissions.from.tipc(sites, tipc.p.H0, clu.n, theta.H0, sampling, opt.sampling, mx.sampled.ntr=6, exclude.O= 1, rtn.int=!debug)
					ntr.hg.H1<- popart.get.sampled.transmissions.from.tipc(sites, tipc.p.H1, clu.n, theta.H1, sampling, opt.sampling, mx.sampled.ntr=6, exclude.O= 1, rtn.int=!debug)
					#print(ntr.hg)						
					ntr.hg.H1["i2i.s",] / ntr.hg.H1["x2i.s",] - ntr.hg.H0["i2i.s",] / ntr.hg.H0["x2i.s",] 
				})
		colnames(ans)<- s.intensity		

		f.name	<- paste(dir.name,"/delta_tipcluster.pdf",sep='')
		cat(paste("\n plot to ",f.name))
		pdf(f.name,version="1.4",width=5,height=5)				
		plot(1,1,type='n',bty='n',xlim=range(s.intensity),ylim=range(c(theta.H1["E2E"]-theta.H0["E2E"],ans)),xlab="sampling intensity",ylab="p(E->E | H1) - p(E->E | H0)")
		lines(s.intensity, ans[1,])
		abline(h=theta.H1["E2E"]-theta.H0["E2E"], lty=2)
		legend("bottomright",bty='n',lty=c(2,1),legend=c("betw true proportions","betw estimates from tip clusters"))
		dev.off()
	}
	
	#for H1 (high E->E), get n(E->E) and n(x->E) for each arm under increasing sampling intensity
	if(0)
	{
		s.intensity	<- seq(0.15,0.4,0.05)
		ans			<- lapply(s.intensity,function(x)
				{
					p.lab			<- x		
					p.consent.coh	<- p.consent.clu<- p.vhcc.prev.AB	<- p.vhcc.inc.AB	<- p.vhcc.prev.C	<- p.vhcc.inc.C	<- 1	
					sampling<- popart.sampling.init(sites,  p.consent.coh, p.consent.clu, p.lab, p.vhcc.prev.AB, p.vhcc.inc.AB, p.vhcc.prev.C, p.vhcc.inc.C, method= opt.sampling)
					sampling$baseline<- x
					ntr.hg<- popart.get.sampled.transmissions.from.tipc(sites, tipc.p.H1, clu.n, theta.H1, sampling, opt.sampling, mx.sampled.ntr=6, exclude.O= 1, rtn.int=!debug)
					#print(ntr.hg)						
					list(i2i.s= ntr.hg["i2i.s",], x2i.s=ntr.hg["x2i.s",]) 
				})
		names(ans)<- s.intensity	
		#print(ans)
		ans.x2E<- sapply(seq_along(ans),function(i) ans[[i]][["x2i.s"]])
		ans.E2E<- sapply(seq_along(ans),function(i) ans[[i]][["i2i.s"]])
		
		print(ans.E2E)
				
		clr<- c("red","blue","green")
		arms<- c("A","B","C")
		sapply(seq_along(arms),function(j)
				{
					plotmat	<- ans.E2E[sites$arm==arms[j],]
					ylim	<- c(0,max( sapply(seq_along(arms),function(k) apply( ans.E2E[sites$arm==arms[k],],2,sum ) ) ))
					
					f.name	<- paste(dir.name,"/sampling.intensity_vs_E2E_arm",arms[j],".pdf",sep='')
					cat(paste("\n plot to ",f.name))
					pdf(f.name,version="1.4",width=5,height=5)
					plot(1,1,type='n',bty='n',xlim=range(s.intensity),ylim=ylim,xlab="sampling intensity",ylab="n(E->E)")
					z		<- rep(0,ncol(plotmat))
					cols	<- colorRampPalette(c(clr[j],"yellow"))(nrow(plotmat))
					for(i in seq_len(nrow(plotmat)))
					{
						polygon( c(s.intensity,rev(s.intensity)), c(z+plotmat[i,], rev(z)), col=cols[i] )
						z	<- z+plotmat[i,]
					}										
					legend("topleft",bty='n',legend=paste("arm",arms[j]))
					dev.off()
				})
		clr<- c("black","black","black")
		arms<- c("A","B","C")
		sapply(seq_along(arms),function(j)
				{
					plotmat	<- ans.x2E[sites$arm==arms[j],]
					ylim	<- c(0,max( sapply(seq_along(arms),function(k) apply( ans.x2E[sites$arm==arms[k],],2,sum ) ) ))
					
					f.name	<- paste(dir.name,"/sampling.intensity_vs_x2E_arm",arms[j],".pdf",sep='')
					cat(paste("\n plot to ",f.name))
					pdf(f.name,version="1.4",width=5,height=5)
					plot(1,1,type='n',bty='n',xlim=range(s.intensity),ylim=ylim,xlab="sampling intensity",ylab="n(E->E)")
					z		<- rep(0,ncol(plotmat))
					cols	<- colorRampPalette(c(clr[j],"gray50"))(nrow(plotmat))
					for(i in seq_len(nrow(plotmat)))
					{
						polygon( c(s.intensity,rev(s.intensity)), c(z+plotmat[i,], rev(z)), col=cols[i] )
						z	<- z+plotmat[i,]
					}										
					legend("topleft",bty='n',legend=paste("arm",arms[j]))
					dev.off()
				})
	}
	
	#for H1 (high E->E), get n(E->E) and n(x->E) for each arm under baseline scenario
	if(0)
	{
		p.lab			<- 0.9		
		p.consent.coh	<- 0.9
		p.consent.clu	<- 0.5
		p.vhcc.prev.AB	<- 0.95
		p.vhcc.inc.AB	<- 0.8
		p.vhcc.prev.C	<- 0.4
		p.vhcc.inc.C	<- 0.4/2					
		p.contam		<- seq(0.1,0.4,0.05)
		
		ans			<- lapply(p.contam,function(x)
				{
					theta.H0	<- clu.p.init(theta.EE.H0, theta.UE, theta.TE, x)
					theta.H1	<- clu.p.init(theta.EE.H1, theta.UE, theta.TE, x)
					tipc.p.H0	<- clu.probabilities(clu.n, theta.H0)
					tipc.p.H1	<- clu.probabilities(clu.n, theta.H1)
					
					sampling<- popart.sampling.init(sites,  p.consent.coh, p.consent.clu, p.lab, p.vhcc.prev.AB, p.vhcc.inc.AB, p.vhcc.prev.C, p.vhcc.inc.C, method= opt.sampling)
					print(sampling)
					print(xtable(sampling, digits=2), floating=FALSE)
					stop()
					ntr.hg<- popart.get.sampled.transmissions.from.tipc(sites, tipc.p.H1, clu.n, theta.H1, sampling, opt.sampling, mx.sampled.ntr=6, exclude.O= 1, rtn.int=!debug)
					#print(ntr.hg)						
					list(i2i.s= ntr.hg["i2i.s",], x2i.s=ntr.hg["x2i.s",]) 
				})
		names(ans)<- p.contam	
		#print(ans)
		ans.x2E<- sapply(seq_along(ans),function(i) ans[[i]][["x2i.s"]])
		ans.E2E<- sapply(seq_along(ans),function(i) ans[[i]][["i2i.s"]])
		
		print(ans.E2E)
		
		clr<- c("red","blue","green")
		arms<- c("A","B","C")
		p.vary<- as.numeric(names(ans))
		xlab<- "proportion transmission from outside cluster"
		sapply(seq_along(arms),function(j)
				{
					plotmat	<- ans.E2E[sites$arm==arms[j],]
					ylim	<- c(0,max( sapply(seq_along(arms),function(k) apply( ans.E2E[sites$arm==arms[k],],2,sum ) ) ))
					
					f.name	<- paste(dir.name,"/",xlab,"_vs_E2E_arm",arms[j],".pdf",sep='')
					cat(paste("\n plot to ",f.name))
					pdf(f.name,version="1.4",width=5,height=5)
					plot(1,1,type='n',bty='n',xlim=range(p.vary),ylim=ylim,xlab=xlab,ylab="n(E->E)")
					z		<- rep(0,ncol(plotmat))
					cols	<- colorRampPalette(c(clr[j],"yellow"))(nrow(plotmat))
					for(i in seq_len(nrow(plotmat)))
					{
						polygon( c(p.vary,rev(p.vary)), c(z+plotmat[i,], rev(z)), col=cols[i] )
						z	<- z+plotmat[i,]
					}										
					legend("topleft",bty='n',legend=paste("arm",arms[j]))
					dev.off()
				})
		clr<- c("black","black","black")
		arms<- c("A","B","C")
		sapply(seq_along(arms),function(j)
				{
					plotmat	<- ans.x2E[sites$arm==arms[j],]
					ylim	<- c(0,max( sapply(seq_along(arms),function(k) apply( ans.x2E[sites$arm==arms[k],],2,sum ) ) ))
					
					f.name	<- paste(dir.name,"/",xlab,"_vs_x2E_arm",arms[j],".pdf",sep='')
					cat(paste("\n plot to ",f.name))
					pdf(f.name,version="1.4",width=5,height=5)
					plot(1,1,type='n',bty='n',xlim=range(p.vary),ylim=ylim,xlab=xlab,ylab="n(E->E)")
					z		<- rep(0,ncol(plotmat))
					cols	<- colorRampPalette(c(clr[j],"gray50"))(nrow(plotmat))
					for(i in seq_len(nrow(plotmat)))
					{
						polygon( c(p.vary,rev(p.vary)), c(z+plotmat[i,], rev(z)), col=cols[i] )
						z	<- z+plotmat[i,]
					}										
					legend("topleft",bty='n',legend=paste("arm",arms[j]))
					dev.off()
				})
	}
	
	#for H1 (high E->E), get n(E->E) and n(x->E) for each arm under baseline scenario
	if(0)
	{
		p.lab			<- 0.9		
		p.consent.coh	<- 0.9
		p.consent.clu	<- 0.5
		p.vhcc.prev.AB	<- 0.95
		p.vhcc.inc.AB	<- 0.8
		p.vhcc.prev.C	<- 0.4
		p.vhcc.inc.C	<- 0.4/2					
		p.contam		<- seq(0.1,0.4,0.05)
		
		ans			<- lapply(p.contam,function(x)
				{
					theta.H0	<- clu.p.init(theta.EE.H0, theta.UE, theta.TE, x)
					theta.H1	<- clu.p.init(theta.EE.H1, theta.UE, theta.TE, x)
					tipc.p.H0	<- clu.probabilities(clu.n, theta.H0)
					tipc.p.H1	<- clu.probabilities(clu.n, theta.H1)
					
					sampling<- popart.sampling.init(sites,  p.consent.coh, p.consent.clu, p.lab, p.vhcc.prev.AB, p.vhcc.inc.AB, p.vhcc.prev.C, p.vhcc.inc.C, method= opt.sampling)
					#print(sampling)
					#print(xtable(sampling, digits=2), floating=FALSE)
					#stop()
					ntr.hg<- popart.get.sampled.transmissions.from.tipc(sites, tipc.p.H1, clu.n, theta.H1, sampling, opt.sampling, mx.sampled.ntr=6, exclude.O= 1, rtn.int=!debug)
					#print(ntr.hg)						
					list(i2i.s= ntr.hg["i2i.s",], x2i.s=ntr.hg["x2i.s",]) 
				})
		names(ans)<- p.contam	
		#print(ans)
		ans.x2E<- sapply(seq_along(ans),function(i) ans[[i]][["x2i.s"]])
		ans.E2E<- sapply(seq_along(ans),function(i) ans[[i]][["i2i.s"]])
		
		print(ans.E2E)
		clr<- c("red","blue","green")
		arms<- c("A","B","C")
		p.vary<- as.numeric(names(ans))
		xlab<- "proportion transmission from outside cluster"
		sapply(seq_along(arms),function(j)
				{
					plotmat	<- ans.E2E[sites$arm==arms[j],]
					ylim	<- c(0,max( sapply(seq_along(arms),function(k) apply( ans.E2E[sites$arm==arms[k],],2,sum ) ) ))
					
					f.name	<- paste(dir.name,"/",xlab,"_vs_E2E_arm",arms[j],".pdf",sep='')
					cat(paste("\n plot to ",f.name))
					pdf(f.name,version="1.4",width=5,height=5)
					plot(1,1,type='n',bty='n',xlim=range(p.vary),ylim=ylim,xlab=xlab,ylab="n(E->E)")
					z		<- rep(0,ncol(plotmat))
					cols	<- colorRampPalette(c(clr[j],"yellow"))(nrow(plotmat))
					for(i in seq_len(nrow(plotmat)))
					{
						polygon( c(p.vary,rev(p.vary)), c(z+plotmat[i,], rev(z)), col=cols[i] )
						z	<- z+plotmat[i,]
					}										
					legend("topleft",bty='n',legend=paste("arm",arms[j]))
					dev.off()
				})
		clr<- c("black","black","black")
		arms<- c("A","B","C")
		sapply(seq_along(arms),function(j)
				{
					plotmat	<- ans.x2E[sites$arm==arms[j],]
					ylim	<- c(0,max( sapply(seq_along(arms),function(k) apply( ans.x2E[sites$arm==arms[k],],2,sum ) ) ))
					
					f.name	<- paste(dir.name,"/",xlab,"_vs_x2E_arm",arms[j],".pdf",sep='')
					cat(paste("\n plot to ",f.name))
					pdf(f.name,version="1.4",width=5,height=5)
					plot(1,1,type='n',bty='n',xlim=range(p.vary),ylim=ylim,xlab=xlab,ylab="n(E->E)")
					z		<- rep(0,ncol(plotmat))
					cols	<- colorRampPalette(c(clr[j],"gray50"))(nrow(plotmat))
					for(i in seq_len(nrow(plotmat)))
					{
						polygon( c(p.vary,rev(p.vary)), c(z+plotmat[i,], rev(z)), col=cols[i] )
						z	<- z+plotmat[i,]
					}										
					legend("topleft",bty='n',legend=paste("arm",arms[j]))
					dev.off()
				})
	}
	
	#for H1 (high E->E), get n(E->E) and n(x->E) for each arm under baseline scenario
	if(0)
	{
		p.lab			<- 0.9		
		p.consent.coh	<- 0.9
		p.consent.clu	<- 0.5
		p.vhcc.prev.AB	<- 0.95
		p.vhcc.inc.AB	<- 0.8
		p.vhcc.prev.C	<- 0.4
		p.vhcc.inc.C	<- 0.4/2					
		p.contam		<- seq(0.1,0.4,0.1)
		xlab			<- "proportion of transmissions from outside cluster"

		ans			<- lapply(p.contam,function(x)
				{
					theta.H0	<- clu.p.init(theta.EE.H0, theta.UE, theta.TE, x)
					theta.H1	<- clu.p.init(theta.EE.H1, theta.UE, theta.TE, x)
					tipc.p.H0	<- clu.probabilities(clu.n, theta.H0, with.ntr.weight=1)
					tipc.p.H1	<- clu.probabilities(clu.n, theta.H1, with.ntr.weight=1)
					sampling<- popart.sampling.init(sites,  p.consent.coh, p.consent.clu, p.lab, p.vhcc.prev.AB, p.vhcc.inc.AB, p.vhcc.prev.C, p.vhcc.inc.C, method= opt.sampling)
					#print(sampling)
					ntr.hg.H0<- popart.get.sampled.transmissions.from.tipc(sites, tipc.p.H0, clu.n, theta.H0, sampling, opt.sampling, mx.sampled.ntr=6, exclude.O= 1, rtn.int=!debug)
					ntr.hg.H1<- popart.get.sampled.transmissions.from.tipc(sites, tipc.p.H1, clu.n, theta.H1, sampling, opt.sampling, mx.sampled.ntr=6, exclude.O= 1, rtn.int=!debug)
					#print(ntr.hg)						
					list(i2i.s.H0= ntr.hg.H0["i2i.s",], x2i.s.H0=ntr.hg.H0["x2i.s",],i2i.s.H1= ntr.hg.H1["i2i.s",], x2i.s.H1=ntr.hg.H1["x2i.s",]) 
				})
		names(ans)<- p.contam	
		#print(ans)
		ans.x2E.H0<- sapply(seq_along(ans),function(i) ans[[i]][["x2i.s.H0"]])
		ans.E2E.H0<- sapply(seq_along(ans),function(i) ans[[i]][["i2i.s.H0"]])
		ans.x2E.H1<- sapply(seq_along(ans),function(i) ans[[i]][["x2i.s.H1"]])
		ans.E2E.H1<- sapply(seq_along(ans),function(i) ans[[i]][["i2i.s.H1"]])
		
		#print(ans.x2E.H1)		
		#print(ans.E2E.H1 / ans.x2E.H1 )
		arms<- c("A","B","C")
		ans	<- lapply(seq_along(arms),function(i)
				{
					cat(paste("\n arm ",arms[i]))
					#power with t-test per arm
					idx			<- sites$arm==arms[i]
					p.H0		<- ans.E2E.H0[idx,] / ans.x2E.H0[idx,]
					print(p.H0)
					p.H1		<- ans.E2E.H1[idx,] / ans.x2E.H1[idx,]					
					var.pH0		<- (apply(p.H0,2,mean)*0.5)^2 	#10*10*apply(ans.E2E.H0[idx,] / ans.x2E.H0[idx,],2,var)
					var.pH1		<- (apply(p.H1,2,mean)*0.5)^2	#10*10*apply(ans.E2E.H0[idx,] / ans.x2E.H0[idx,],2,var)
					
					pw			<- sapply(seq_along(ans),function(j)
									{
										n.H1	<- floor(ans.x2E.H1[idx,j])										
										phdes.power.ttest.cl(n.H1, p.H0[,j], p.H1[,j], var.pH0[j], var.pH1[j], alpha= test.alpha)
									})
					names(pw)	<- names(ans)
					#print(pw)
					
					#design effect: 1+(harmonic mean(n.H1)-1) * within clu corr
					des.effect	<- 1 + (1/apply(1/ans.x2E.H1[idx,],2,mean)-1) * var.pH1 / (apply(p.H1,2,mean)*(1-apply(p.H1,2,mean)))
					des.effect	<- 1 + (apply(ans.x2E.H1[idx,],2,mean)-1) * var.pH1 / (apply(p.H1,2,mean)*(1-apply(p.H1,2,mean)))
					#print(des.effect)
					
					#power with simple binom test, adjusting for design effect
					p.H0		<- apply(ans.E2E.H0[idx,] / ans.x2E.H0[idx,],2,mean)
					p.H1		<- apply(ans.E2E.H1[idx,] / ans.x2E.H1[idx,],2,mean)
					e.H0		<- apply(floor(ans.E2E.H0[idx,]),2,sum) / des.effect
					e.H1		<- apply(floor(ans.E2E.H1[idx,]),2,sum) / des.effect
					n.H0		<- apply(floor(ans.x2E.H0[idx,]),2,sum) / des.effect
					n.H1		<- apply(floor(ans.x2E.H1[idx,]),2,sum) / des.effect
					pw.b		<- .phdes.binom.power(n.H1, p.H1, p.H0, test.alpha, method="asymp")
					names(pw.b)	<- names(ans)
					#print(pw.b)
					
					#95% conf interval with simple binom test, adjusting for design effect
					conf.H1		<- binom.confint(e.H1, n.H1, conf.level = 0.95, methods="cloglog")[,c("lower","upper")]
					conf.H0		<- binom.confint(e.H0, n.H0, conf.level = 0.95, methods="cloglog")[,c("lower","upper")]
					rownames(conf.H0)<- rownames(conf.H1)<- names(ans)
					#print( conf.H0 ); print( conf.H1 )
					
					list(pw.n=pw, pw.b=pw.b,conf.b.H0=conf.H0, conf.b.H1=conf.H1, p.H0.m=p.H0, p.H1.m=p.H1)
				})
		print(ans)		
		#plot computed power for each arm	
		pw.n<- sapply(seq_along(ans),function(i)
				{
					ans[[i]][["pw.n"]]
				})		
		colnames(pw.n)<- arms
		pw.b<- sapply(seq_along(ans),function(i)
				{
					ans[[i]][["pw.b"]]
				})		
		colnames(pw.b)<- arms
		
		clr<- c("red","blue","green")
		f.name	<- paste(dir.name,"/",xlab,"_power.pdf",sep='')
		cat(paste("\n plot to ",f.name))
		pdf(f.name,version="1.4",width=5,height=5)		
		plot(1,1,type='n',bty='n',xlim=range(p.contam),ylim=range(c(pw.b,pw.n)),xlab=xlab,ylab="power")
		sapply(seq_len(ncol(pw.n)),function(j)
				{
					lines(p.contam,pw.n[,j],col=clr[j],lty=1)
				})
		sapply(seq_len(ncol(pw.b)),function(j)
				{
					lines(p.contam,pw.b[,j],col=clr[j],lty=2)
				})
		legend("bottomleft",bty='n',lty=c(1,2),legend=c("normal approx (Hayes1999)","Binomial using effective size"))
		legend("bottomright",bty='n',fill=clr,legend=c("arm A","arm B","arm C"))
		dev.off()
		
		#plot confidence intervals for each arm
		sapply(seq_along(ans),function(i)
				{
					ylim<- range(c(ans[[i]][["conf.b.H0"]],ans[[i]][["conf.b.H1"]]))
					ylim<- ylim*c(1,1.3)
					
					f.name	<- paste(dir.name,"/",xlab,"_confint_",arms[i],".pdf",sep='')
					cat(paste("\n plot to ",f.name))
					pdf(f.name,version="1.4",width=5,height=5)		
					
					plot(1,1,type='n',bty='n',xlim=range(p.contam),ylim=ylim,xlab=xlab,ylab="95% Binomial confidence interval")
					x<- ans[[i]][["conf.b.H0"]]
					#polygon( c(p.contam,rev(p.contam)), c(x[,1],rev(x[,2])), col=my.fade.col(clr[i],0.4), border=NA )
					polygon( c(p.contam,rev(p.contam)), c(x[,1],rev(x[,2])), col="gray50", border=NA )
					x<- ans[[i]][["conf.b.H1"]]
					polygon( c(p.contam,rev(p.contam)), c(x[,1],rev(x[,2])), col=my.fade.col(clr[i],0.7), border=NA )					
					lines(p.contam,ans[[i]][["p.H0.m"]], lty=2)
					lines(p.contam,ans[[i]][["p.H1.m"]], lty=3)
					legend("topleft",bty='n',legend=paste("arm",arms[i]))
					dev.off()
				})  
		
	}
	
	
	stop()
	
	a2a.hg<- popart.get.sampled.transmissions.from.tipc(sites, tipc.p.H0, clu.n, theta.H0, sampling, opt.sampling, mx.sampled.ntr=6, exclude.O= 1, rtn.int=!debug)
	print(a2a.hg)
	
	p.lab			<- 0.9		
	p.consent.coh	<- 0.9
	p.consent.clu	<- 0.5
	p.vhcc.prev.AB	<- 0.95
	p.vhcc.inc.AB	<- 0.8
	p.vhcc.prev.C	<- 0.4
	p.vhcc.inc.C	<- 0.4/2
	#compute sampling probabilities at baseline and during trial
	sampling<- popart.sampling.init(sites,  p.consent.coh, p.consent.clu, p.lab, p.vhcc.prev.AB, p.vhcc.inc.AB, p.vhcc.prev.C, p.vhcc.inc.C, method= opt.sampling)
	print(sampling)
	a2a.hg<- popart.get.sampled.transmissions.from.tipc(sites, tipc.p.H0, clu.n, theta.H0, sampling, opt.sampling, mx.sampled.ntr=6, exclude.O= 1, rtn.int=!debug)
	print(a2a.hg)
#update 'popart.get.sampled.acute2acute'
}
	
###############################################################################
#' compute the power in differentiating H0: acute to acute is <10\% vs H0: acute to acute is >40\% with the tipc cluster approach
#' 
#' This script varies the probability that HIV+ individuals consent to participation in the phylogenetic study at HCC
#' @param p.phylosignal				Probability that a true (potentially indirect) transmission is identified with a phylogenetic method. (Here not used)
#' @param p.nocontam				Frequency with which transmissions occur within the community
#' @param opt.pooled				Pooling option for power analysis
#' @param opt.sampling				Sampling option for trial
prj.popart.powercalc_tipc_consenting<- function(p.phylosignal=0.7,p.nocontam=0.85, opt.pooled= "no pooling", opt.sampling= "PC and HCC")
{
	require(binom)
	require(phylodesign)
	my.mkdir(DATA,"popartpowercalc")
	dir.name<- paste(DATA,"popartpowercalc",sep='/')	
	resume<- 0
	verbose<- 0
	plot.increment<- 0.05
	
	m.type		<- "Acute"	
	cohort.size	<- 2500
	pc24.size	<- 6000
	cohort.dur	<- 3	
	theta.EE.H0	<- 0.10
	theta.EE.H1	<- 0.4
	theta.UE	<- 0.3
	theta.TE	<- theta.UE / 5
	test.alpha	<- 0.05		 
	debug		<- 1
	pooled.n	<- 1
	opt.pooled	<- "no pooling"#"pooled across ZA"#"pooled across trial"#"no pooling"
	opt.pooled	<- "pooled across SA"
	opt.clu.closure	<- 14
	opt.sampling<- "PC and HCC"#"only HCC"	#"PC and HCC"	#
	#opt.sampling<- "PC after yr 1 and HCC"
	#opt.sampling<- "PC only incident and HCC"
	opt.power	<-	"All"
	if(!opt.pooled%in%c("pooled across country","pooled across ZA","pooled across SA","pooled across trial","no pooling"))
		stop("prj.popart.powercalc_medsize: unknown method to pool")
	if(!opt.power%in%c("All","PonlyPC","IonlyPC","PonlyPCandIonlyPC"))
		stop("prj.popart.powercalc_medsize: unknown method to pool")
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
	
	sites<-	popart.getdata.randomized.arm( pooled.n, rtn.fixed=debug )
	sites<-	popart.getdata.randomized.n(sites, cohort.size, cohort.dur, rtn.exp=debug)
	#data("popart.tipcp.highacute")
	#data("popart.tipcp.lowacute")
		
	print(sites)
	#compute complete tip cluster probs under H0 and H1
	clu.n		<- clu.n.of.tchain(opt.clu.closure)	
	theta.H0	<- clu.p.init(theta.EE.H0, theta.UE, theta.TE, 1-p.nocontam)
	theta.H1	<- clu.p.init(theta.EE.H1, theta.UE, theta.TE, 1-p.nocontam)
	tipc.p.H0	<- clu.probabilities(clu.n, theta.H0)
	tipc.p.H1	<- clu.probabilities(clu.n, theta.H1)
	
	p.lab			<- 0.9		
	p.consent.coh	<- 0.9
	p.consent.clu	<- 0.5
	p.vhcc.prev.AB	<- 0.95
	p.vhcc.inc.AB	<- 0.8
	p.vhcc.prev.C	<- 0.4
	p.vhcc.inc.C	<- 0.4/2
	#compute sampling probabilities at baseline and during trial
	sampling<- popart.sampling.init(sites,  p.consent.coh, p.consent.clu, p.lab, p.vhcc.prev.AB, p.vhcc.inc.AB, p.vhcc.prev.C, p.vhcc.inc.C, method= opt.sampling)
	print(sampling)
	a2a.hg<- popart.get.sampled.transmissions.from.tipc(sites, tipc.p.H0, clu.n, theta.H0, sampling, opt.sampling, mx.sampled.ntr=6, exclude.O= 1, rtn.int=!debug)
	print(a2a.hg)
	#update 'popart.get.sampled.acute2acute'

	stop()
	###############################################################################
	#vary %consenting in HCC
	###############################################################################
	cat("\ncompute sampled acute to acute transmissions for %consenting")
	p.vhcc.prev.Cs<- seq(0.15,0.4,plot.increment)
	p.consent.clus<- seq(0.5,0.8,plot.increment)
	out<- lapply(p.vhcc.prev.Cs,function(x)
			{
				inc<- lapply(p.consent.clus,function(y)
						{
							p.lab			<- 0.9		
							p.consent.coh	<- 0.9
							p.consent.clu	<- y
							p.vhcc.prev.AB	<- 0.95
							p.vhcc.inc.AB	<- 0.8
							p.vhcc.prev.C	<- x
							p.vhcc.inc.C	<- x/2
							tipc.p			<- phdes.get.hyp.tipc.probs(popart.tipcp.highacute,1-p.nocontam)
							a2a.hg			<- popart.get.sampled.acute2acute( 	sites, tipc.p, opt.sampling, 
																				p.vhcc.prev.AB=p.vhcc.prev.AB, p.vhcc.prev.C=p.vhcc.prev.C, p.vhcc.inc.AB=p.vhcc.inc.AB, p.vhcc.inc.C=p.vhcc.inc.C, 
																				consent.PC=p.consent.coh, consent.HCC=p.consent.clu, p.lab=p.lab,
																				rtn.int=!debug)
							tipc.p			<- phdes.get.hyp.tipc.probs(popart.tipcp.lowacute,1-p.nocontam)
							a2a.lw			<- popart.get.sampled.acute2acute( 	sites, tipc.p, opt.sampling, 
																				p.vhcc.prev.AB=p.vhcc.prev.AB, p.vhcc.prev.C=p.vhcc.prev.C, p.vhcc.inc.AB=p.vhcc.inc.AB, p.vhcc.inc.C=p.vhcc.inc.C, 
																				consent.PC=p.consent.coh, consent.HCC=p.consent.clu, p.lab=p.lab,
																				rtn.int=!debug)																
							list(a2a.lw=a2a.lw, a2a.hg=a2a.hg)				
						})
				names(inc)<- p.consent.clus										
				x2i.lw	<- sapply(inc,function(x)		x[["a2a.lw"]]["x2i.s",])
				i2i.lw	<- sapply(inc,function(x)		x[["a2a.lw"]]["i2i.s",])	
				x2i.hg	<- sapply(inc,function(x)		x[["a2a.hg"]]["x2i.s",])
				i2i.hg	<- sapply(inc,function(x)		x[["a2a.hg"]]["i2i.s",])					
				# pool transmissions
				x2i.lw							<- 		popart.pool(sites, x2i.lw, method=opt.pooled)[["transm"]]
				i2i.lw							<- 		popart.pool(sites, i2i.lw, method=opt.pooled)[["transm"]]
				x2i.hg							<- 		popart.pool(sites, x2i.hg, method=opt.pooled)[["transm"]]
				g(i2i.hg, idx.A, idx.B, idx.C)	%<-% 	popart.pool(sites, i2i.hg, method=opt.pooled)
				
				test.biased.H0	<- lapply(	list(idx.A, idx.B, idx.C),
												function(arm)	apply(i2i.lw[arm,,drop=0],2,mean)/apply(x2i.lw[arm,,drop=0],2,mean)
												)
				test.biased.H1	<- lapply(	list(idx.A, idx.B, idx.C),
												function(arm)	apply(i2i.hg[arm,,drop=0],2,mean)/apply(x2i.hg[arm,,drop=0],2,mean)
												)																
				names(test.biased.H0)<- c("A","B","C")
				names(test.biased.H1)<- c("A","B","C")
				
				# compute power and confidence intervals
				conf.lw.A		<- phdes.binom.power(x2i.lw[idx.A,,drop=0], i2i.lw[idx.A,,drop=0], test.biased.H0[["A"]], test.biased.H1[["A"]], test.alpha, verbose=0, method.ci="asymptotic",method.pw="cloglog")[["conf"]]
				conf.lw.B		<- phdes.binom.power(x2i.lw[idx.B,,drop=0], i2i.lw[idx.B,,drop=0], test.biased.H0[["B"]], test.biased.H1[["B"]], test.alpha, verbose=0, method.ci="asymptotic",method.pw="cloglog")[["conf"]]
				conf.lw.C		<- phdes.binom.power(x2i.lw[idx.C,,drop=0], i2i.lw[idx.C,,drop=0], test.biased.H0[["C"]], test.biased.H1[["C"]], test.alpha, verbose=0, method.ci="asymptotic",method.pw="cloglog")[["conf"]]
				  				 	 
				g(conf.hg.A,is.conf.hg.A,power.hg.A)	%<-% phdes.binom.power(x2i.hg[idx.A,,drop=0], i2i.hg[idx.A,,drop=0], test.biased.H0[["A"]], test.biased.H1[["A"]], test.alpha, verbose=0, method.ci="asymptotic",method.pw="cloglog")
				g(conf.hg.B,is.conf.hg.B,power.hg.B)	%<-% phdes.binom.power(x2i.hg[idx.B,,drop=0], i2i.hg[idx.B,,drop=0], test.biased.H0[["B"]], test.biased.H1[["B"]], test.alpha, verbose=0, method.ci="asymptotic",method.pw="cloglog")
				g(conf.hg.C,is.conf.hg.C,power.hg.C)	%<-% phdes.binom.power(x2i.hg[idx.C,,drop=0], i2i.hg[idx.C,,drop=0], test.biased.H0[["C"]], test.biased.H1[["C"]], test.alpha, verbose=0, method.ci="asymptotic",method.pw="cloglog")
												
				tmp<- list(	i2i.hg.A= apply(i2i.hg[idx.A,,drop=0],2,mean),
							i2i.hg.B= apply(i2i.hg[idx.B,,drop=0],2,mean),
							i2i.hg.C= apply(i2i.hg[idx.C,,drop=0],2,mean),
							x2i.hg.A= apply(x2i.hg[idx.A,,drop=0],2,mean),
							x2i.hg.B= apply(x2i.hg[idx.B,,drop=0],2,mean),
							x2i.hg.C= apply(x2i.hg[idx.C,,drop=0],2,mean),
							is.conf.A= is.conf.hg.A,
							is.conf.B= is.conf.hg.B,
							is.conf.C= is.conf.hg.C,
							power.A= power.hg.A,
							power.B= power.hg.B,
							power.C= power.hg.C,
							conf.hg.A= conf.hg.A,
							conf.hg.B= conf.hg.B,
							conf.hg.C= conf.hg.C,
							conf.lw.A= conf.lw.A,
							conf.lw.B= conf.lw.B,
							conf.lw.C= conf.lw.C
							)	
				tmp
			})
	names(out)<- p.vhcc.prev.Cs
	i2i.hg			<- lapply(c("i2i.hg.A","i2i.hg.B","i2i.hg.C"), function(arm)	sapply(out,function(x)		x[[arm]]	)	)
	names(i2i.hg)	<- c("A","B","C")
	conf.hg.u		<- lapply(c("conf.hg.A","conf.hg.B","conf.hg.C"), function(arm)	sapply(out,function(x)		x[[arm]][,"upper"]	)	)
	names(conf.hg.u)<- c("A","B","C")
	conf.lw.u		<- lapply(c("conf.lw.A","conf.lw.B","conf.lw.C"), function(arm)	sapply(out,function(x)		x[[arm]][,"upper"]	)	)
	names(conf.lw.u)<- c("A","B","C")
	conf.hg.l		<- lapply(c("conf.hg.A","conf.hg.B","conf.hg.C"), function(arm)	sapply(out,function(x)		x[[arm]][,"lower"]	)	)
	names(conf.hg.l)<- c("A","B","C")
	conf.lw.l		<- lapply(c("conf.lw.A","conf.lw.B","conf.lw.C"), function(arm)	sapply(out,function(x)		x[[arm]][,"lower"]	)	)
	names(conf.lw.l)<- c("A","B","C")
	
	###############################################################################
	#plot sampled a2a
	###############################################################################
	if(0)
	{
		require(fields)
		f.name<- paste(dir.name,paste("VARYCONSENT_TIPC_a2a_",p.nocontam,"power",opt.power,"pool",opt.pooled,"sample",opt.sampling,"pwcalc",theta.EE.H0,theta.EE.H1,test.alpha,".pdf",sep='_'),sep='/')
		cat(paste("\nplot a2a to\n",f.name))
		pdf(paste(f.name),version="1.4",width=12,height=6)		
		breaks		<- diff(range(c(i2i.hg[["A"]],i2i.hg[["B"]],i2i.hg[["C"]])))/50
		breaks		<- seq(min(c(i2i.hg[["A"]],i2i.hg[["B"]],i2i.hg[["C"]])),by=breaks, len=51)	
		#def.par <- par(no.readonly = TRUE)
		layout.m<- matrix(data= seq_len(3),ncol=3,nrow=1,byrow=1)
		layout(layout.m)				
		sapply(c("A","B","C"),function(arm)
			{
				if(arm!="C")
					image(main=paste("arm",arm),p.consent.clus,p.vhcc.prev.Cs,i2i.hg[[arm]], breaks=breaks, col=head( rev(gray(seq(0,.95,len=trunc(50*1.4)))), 50))
				else
					image.plot(main=paste("arm",arm),p.consent.clus,p.vhcc.prev.Cs,i2i.hg[[arm]], breaks=breaks, col=head( rev(gray(seq(0,.95,len=trunc(50*1.4)))), 50),zlim= range(i2i.hg))
			})		
		dev.off()
	}
	###############################################################################
	#plot panel of confidence intervals
	###############################################################################
	if(1)
	{
		cols<- c("deepskyblue","dodgerblue4")
		cat(paste("\nplot confidence panels\n"))
		sapply(names(conf.lw.u),function(arm)
				{
					f.name<- paste(dir.name,paste("VARYCONSENT_TIPC_confint",arm,p.nocontam,"power",opt.power,"pool",opt.pooled,"sample",opt.sampling,"pwcalc",theta.EE.H0,theta.EE.H1,test.alpha,".pdf",sep='_'),sep='/')
					cat(paste("\nplot confint",arm,"to\n",f.name))
					pdf(paste(f.name),version="1.4",width=6,height=12)
					phdes.plot.confint.panel(t(conf.lw.l[[arm]]),t(conf.lw.u[[arm]]),t(conf.hg.l[[arm]]),t(conf.hg.u[[arm]]),p.consent.clus,p.vhcc.prev.Cs,"p.vhcc.prev.Cs","p.consent.clus", cols=cols)
					dev.off()					
				})		
	}
	stop()	
}
###############################################################################
#' compute the power in differentiating H0: acute to acute is <10\% vs H0: acute to acute is >40\% with the tipc cluster approach
#' 
#' This script varies the probability that transmission occurs outside the community.
#' @param p.prev.instudy.clu.armC	Probability that HIV+ individuals visit an HCC in arm C (sensitivity analysis).
#' @param p.phylosignal				Probability that a true (potentially indirect) transmission is identified with a phylogenetic method. (Here not used)
#' @param opt.pooled				Pooling option for power analysis
#' @param opt.sampling				Sampling option for trial
prj.popart.powercalc_tipc_contam<- function(p.prev.instudy.clu.armC=0.4, p.phylosignal=0.7, opt.pooled= "no pooling", opt.sampling= "PC and HCC")
{
	require(binom)
	require(phylodesign)
	my.mkdir(DATA,"popartpowercalc")
	dir.name<- paste(DATA,"popartpowercalc",sep='/')	
	resume<- 0
	verbose<- 0
	plot.increment<- 0.05
	
	m.type<- "Acute"	
	cohort.size<- 2500
	pc24.size<- 6000
	cohort.dur<- 3	
	theta.EE.H0<- 0.10
	theta.EE.H1<- 0.4
	test.alpha<- 0.05		 
	debug<- 0
	pooled.n<- 1
	opt.pooled<- "pooled across T5"#"pooled across ZA"#"pooled across trial"#"no pooling"
	opt.sampling<- "PC and HCC"#"only HCC"#"PC and HCC"
	opt.power<-	"All"
	if(!opt.pooled%in%c("pooled across country","pooled across ZA","pooled across SA","pooled across trial","pooled across T5","no pooling"))
		stop("prj.popart.powercalc_medsize: unknown method to pool")
	if(!opt.power%in%c("All","PonlyPC","IonlyPC","PonlyPCandIonlyPC"))
		stop("prj.popart.powercalc_medsize: unknown method to pool")
	cat(paste("\ncohort.size",cohort.size))
	cat(paste("\ncohort.dur",cohort.dur))
	cat(paste("\ntheta.EE.H0",theta.EE.H0))
	cat(paste("\ntheta.EE.H1",theta.EE.H1))
	cat(paste("\ntest.alpha",test.alpha))	
	cat(paste("\np.prev.instudy.clu.armC",p.prev.instudy.clu.armC))
	cat(paste("\npooled.n",pooled.n))
	cat(paste("\nopt.pooled",opt.pooled))
	cat(paste("\nopt.power",opt.power))
	cat(paste("\nopt.sampling",opt.sampling))
	
	sites<-	popart.getdata.randomized.arm( pooled.n, rtn.fixed=debug )
	sites<-	popart.getdata.randomized.n(sites, cohort.size, cohort.dur, rtn.exp=debug)
	data("popart.tipcp.highacute")
	data("popart.tipcp.lowacute")
	#print(sites)
	###############################################################################
	#vary %consenting in HCC
	###############################################################################
	cat("\ncompute sampled acute to acute transmissions for %contamination")
	p.nocontams		<- seq(0.5,0.95,0.1) 	
	p.consent.clus	<- seq(0.4,0.8,plot.increment)
	out<- lapply(p.nocontams,function(x)
			{
				inc	<- lapply(p.consent.clus,function(y)
						{
							p.lab			<- 0.95		
							p.consent.coh	<- 0.9
							p.consent.clu	<- y
							p.vhcc.prev.AB	<- 0.95
							p.vhcc.inc.AB	<- 0.8
							p.vhcc.prev.C	<- p.prev.instudy.clu.armC
							p.vhcc.inc.C	<- p.prev.instudy.clu.armC/2
							p.nocontam		<- x
							tipc.p			<- phdes.get.hyp.tipc.probs(popart.tipcp.highacute,1-p.nocontam)
							a2a.hg			<- popart.get.sampled.acute2acute( 	sites, tipc.p, opt.sampling, 
																				p.vhcc.prev.AB=p.vhcc.prev.AB, p.vhcc.prev.C=p.vhcc.prev.C, p.vhcc.inc.AB=p.vhcc.inc.AB, p.vhcc.inc.C=p.vhcc.inc.C, 
																				consent.PC=p.consent.coh, consent.HCC=p.consent.clu, p.lab=p.lab,
																				rtn.int=!debug)
							tipc.p			<- phdes.get.hyp.tipc.probs(popart.tipcp.lowacute,1-p.nocontam)
							a2a.lw			<- popart.get.sampled.acute2acute( 	sites, tipc.p, opt.sampling, 
																				p.vhcc.prev.AB=p.vhcc.prev.AB, p.vhcc.prev.C=p.vhcc.prev.C, p.vhcc.inc.AB=p.vhcc.inc.AB, p.vhcc.inc.C=p.vhcc.inc.C, 
																				consent.PC=p.consent.coh, consent.HCC=p.consent.clu, p.lab=p.lab,
																				rtn.int=!debug)																
							list(a2a.lw=a2a.lw, a2a.hg=a2a.hg)				
						})
				names(inc)<- p.consent.clus										
				x2i.lw	<- sapply(inc,function(x)		x[["a2a.lw"]]["x2i.s",])
				i2i.lw	<- sapply(inc,function(x)		x[["a2a.lw"]]["i2i.s",])	
				x2i.hg	<- sapply(inc,function(x)		x[["a2a.hg"]]["x2i.s",])
				i2i.hg	<- sapply(inc,function(x)		x[["a2a.hg"]]["i2i.s",])					
				# pool transmissions
				x2i.lw							<- 		popart.pool(sites, x2i.lw, method=opt.pooled)[["transm"]]
				i2i.lw							<- 		popart.pool(sites, i2i.lw, method=opt.pooled)[["transm"]]
				x2i.hg							<- 		popart.pool(sites, x2i.hg, method=opt.pooled)[["transm"]]
				g(i2i.hg, idx.A, idx.B, idx.C)	%<-% 	popart.pool(sites, i2i.hg, method=opt.pooled)
				
				test.biased.H0	<- lapply(	list(idx.A, idx.B, idx.C),
						function(arm)	apply(i2i.lw[arm,,drop=0],2,mean)/apply(x2i.lw[arm,,drop=0],2,mean)
				)
				test.biased.H1	<- lapply(	list(idx.A, idx.B, idx.C),
						function(arm)	apply(i2i.hg[arm,,drop=0],2,mean)/apply(x2i.hg[arm,,drop=0],2,mean)
				)																
				names(test.biased.H0)<- c("A","B","C")
				names(test.biased.H1)<- c("A","B","C")
				
				# compute power and confidence intervals
				conf.lw.A		<- phdes.binom.power(x2i.lw[idx.A,,drop=0], i2i.lw[idx.A,,drop=0], test.biased.H0[["A"]], test.biased.H1[["A"]], test.alpha, verbose=0, method.ci="asymptotic",method.pw="cloglog")[["conf"]]
				conf.lw.B		<- phdes.binom.power(x2i.lw[idx.B,,drop=0], i2i.lw[idx.B,,drop=0], test.biased.H0[["B"]], test.biased.H1[["B"]], test.alpha, verbose=0, method.ci="asymptotic",method.pw="cloglog")[["conf"]]
				conf.lw.C		<- phdes.binom.power(x2i.lw[idx.C,,drop=0], i2i.lw[idx.C,,drop=0], test.biased.H0[["C"]], test.biased.H1[["C"]], test.alpha, verbose=0, method.ci="asymptotic",method.pw="cloglog")[["conf"]]
				
				g(conf.hg.A,is.conf.hg.A,power.hg.A)	%<-% phdes.binom.power(x2i.hg[idx.A,,drop=0], i2i.hg[idx.A,,drop=0], test.biased.H0[["A"]], test.biased.H1[["A"]], test.alpha, verbose=0, method.ci="asymptotic",method.pw="cloglog")
				g(conf.hg.B,is.conf.hg.B,power.hg.B)	%<-% phdes.binom.power(x2i.hg[idx.B,,drop=0], i2i.hg[idx.B,,drop=0], test.biased.H0[["B"]], test.biased.H1[["B"]], test.alpha, verbose=0, method.ci="asymptotic",method.pw="cloglog")
				g(conf.hg.C,is.conf.hg.C,power.hg.C)	%<-% phdes.binom.power(x2i.hg[idx.C,,drop=0], i2i.hg[idx.C,,drop=0], test.biased.H0[["C"]], test.biased.H1[["C"]], test.alpha, verbose=0, method.ci="asymptotic",method.pw="cloglog")
				
				tmp<- list(	i2i.hg.A= apply(i2i.hg[idx.A,,drop=0],2,mean),
							i2i.hg.B= apply(i2i.hg[idx.B,,drop=0],2,mean),
							i2i.hg.C= apply(i2i.hg[idx.C,,drop=0],2,mean),
							x2i.hg.A= apply(x2i.hg[idx.A,,drop=0],2,mean),
							x2i.hg.B= apply(x2i.hg[idx.B,,drop=0],2,mean),
							x2i.hg.C= apply(x2i.hg[idx.C,,drop=0],2,mean),
							is.conf.A= is.conf.hg.A,
							is.conf.B= is.conf.hg.B,
							is.conf.C= is.conf.hg.C,
							power.A= power.hg.A,
							power.B= power.hg.B,
							power.C= power.hg.C,
							conf.hg.A= conf.hg.A,
							conf.hg.B= conf.hg.B,
							conf.hg.C= conf.hg.C,
							conf.lw.A= conf.lw.A,
							conf.lw.B= conf.lw.B,
							conf.lw.C= conf.lw.C
							)	
				tmp						
			})
	names(out)<- p.nocontams
	i2i.hg			<- lapply(c("i2i.hg.A","i2i.hg.B","i2i.hg.C"), function(arm)	sapply(out,function(x)		x[[arm]]	)	)
	names(i2i.hg)	<- c("A","B","C")
	conf.hg.u		<- lapply(c("conf.hg.A","conf.hg.B","conf.hg.C"), function(arm)	sapply(out,function(x)		x[[arm]][,"upper"]	)	)
	names(conf.hg.u)<- c("A","B","C")
	conf.lw.u		<- lapply(c("conf.lw.A","conf.lw.B","conf.lw.C"), function(arm)	sapply(out,function(x)		x[[arm]][,"upper"]	)	)
	names(conf.lw.u)<- c("A","B","C")
	conf.hg.l		<- lapply(c("conf.hg.A","conf.hg.B","conf.hg.C"), function(arm)	sapply(out,function(x)		x[[arm]][,"lower"]	)	)
	names(conf.hg.l)<- c("A","B","C")
	conf.lw.l		<- lapply(c("conf.lw.A","conf.lw.B","conf.lw.C"), function(arm)	sapply(out,function(x)		x[[arm]][,"lower"]	)	)
	names(conf.lw.l)<- c("A","B","C")		
	###############################################################################
	#plot sampled a2a
	###############################################################################
	if(1)
	{
		require(fields)
		f.name<- paste(dir.name,paste("VARYCONTAM_TIPC_a2a_visitHCC.C",p.prev.instudy.clu.armC,"power",opt.power,"pool",opt.pooled,"sample",opt.sampling,"pwcalc",theta.EE.H0,theta.EE.H1,test.alpha,".pdf",sep='_'),sep='/')
		cat(paste("\nplot a2a to\n",f.name))
		pdf(paste(f.name),version="1.4",width=12,height=6)		
		breaks		<- diff(range(unlist(i2i.hg)))/50
		breaks		<- seq(min(unlist(i2i.hg)),by=breaks, len=51)		
		#def.par <- par(no.readonly = TRUE)
		layout.m<- matrix(data= seq_len(3),ncol=3,nrow=1,byrow=1)
		layout(layout.m)				
		sapply(c("A","B","C"),function(arm)
				{
					if(arm!="C")
						image(main=paste("arm",arm),p.consent.clus,p.nocontams,i2i.hg[[arm]], breaks=breaks, col=head( rev(gray(seq(0,.95,len=trunc(50*1.4)))), 50))
					else
						image.plot(main=paste("arm",arm),p.consent.clus,p.nocontams,i2i.hg[[arm]], breaks=breaks, col=head( rev(gray(seq(0,.95,len=trunc(50*1.4)))), 50),zlim= range(i2i.hg))
				})
		dev.off()	
	}
	###############################################################################
	#plot panel of confidence intervals
	###############################################################################
	if(1)
	{
		cols<- c("deepskyblue","dodgerblue4")
		cat(paste("\nplot confidence panels\n"))
		sapply(names(conf.lw.u),function(arm)
				{
					f.name<- paste(dir.name,paste("VARYCONTAM_TIPC_confint",arm,"visitHCC.C",p.prev.instudy.clu.armC,"power",opt.power,"pool",opt.pooled,"sample",opt.sampling,"pwcalc",theta.EE.H0,theta.EE.H1,test.alpha,".pdf",sep='_'),sep='/')
					cat(paste("\nplot confint C to\n",f.name))
					pdf(paste(f.name),version="1.4",width=6,height=12)
					phdes.plot.confint.panel(t(conf.lw.l[[arm]]),t(conf.lw.u[[arm]]),t(conf.hg.l[[arm]]),t(conf.hg.u[[arm]]),p.consent.clus,p.nocontams,"p.vhcc.prev.Cs","p.nocontams", cols=cols)
					dev.off()					
				})				
	}
	stop()	
}
###############################################################################
prj.popart.tchain_test<- function()
{
	clu.clo	<- 10
	clu.n	<- clu.n.of.tchain(clu.clo)
	my.mkdir(DATA,"popartpowercalc_test")
	dir.name<- paste(DATA,"popartpowercalc_test",sep='/')	
	
	if(0)
	{
		require(xtable)
		print(xtable(clu.n, digits=0), floating=FALSE)		
	}
	
	#densities for given tip cluster size
	if(0)
	{
		p.U2E	<- 0.2
		r.E2E	<- c(0.5,1,2,4)
		clu.p	<- lapply(r.E2E,function(x)
				{
					theta<- p.U2E*c(x,1)
					names(theta)<- c("E2E","U2E")				
					clu.p.of.tchain(clu.n, theta )								
				})
		
		clu.p	<- lapply(seq_along(clu.p),function(i)		
				{
					norm<- apply(clu.p[[i]],2,function(x)  sum(x,na.rm=1))
					clu.p[[i]] / matrix( rep(norm,each=nrow(clu.p[[i]])), nrow(clu.p[[i]]), ncol(clu.p[[i]]) )
				})
		print(clu.p)
		#plot
		col.idx<- 4
		clu.dens<- sapply(clu.p,function(x) x[,col.idx])
		#print(clu.dens)
		xlim<- range(c(1,col.idx))
		ylim<- range(clu.dens,na.rm=1)
		par(mar=c(4,5,1,1))
		plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim,xlab="i",ylab=expression("p("*U^i*" "*E^(n-i)*")"))
		sapply(seq_len(ncol(clu.dens)),function(j)
				{
					lines(1:col.idx, clu.dens[1:col.idx,j],lty=j)
				})
		legend("topright",lty=seq_len(ncol(clu.dens)), legend=r.E2E, bty='n')
	}
	
	#log prob of tip clusters, normalized after a closure at 'clu.clo'
	if(0)
	{
		ltys	<- c(4,1,2,3)
		p.U2E	<- 0.2
		r.E2E	<- c(0.5,1,2,4)
		clu.p	<- lapply(r.E2E,function(x)
				{
					theta<- p.U2E*c(x,1)
					theta<- theta / sum(theta)
					names(theta)<- c("E2E","U2E")				
					clu.p.of.tchain(clu.n, theta )								
				})
		
		clu.nchain	<- apply( clu.n, 2, sum )				
		clu.logp	<- sapply(seq_along(clu.p),function(i)
				{
					tmp<- apply(clu.p[[i]],2,function(x)	sum(x,na.rm=1)	)
					tmp<- tmp / clu.nchain
					log( tmp / sum(tmp) )
				})
		print(clu.logp)
		
		#plot
		f.name	<- paste(dir.name,"/tipcluster_loglkl_wsize.pdf",sep='')
		cat(paste("\n plot to ",f.name))
		pdf(f.name,version="1.4",width=5,height=5)		
		
		xlim<- c(1,nrow(clu.logp))
		ylim<- range(clu.logp)
		plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim,xlab="number of transmissions in tip cluster",ylab=expression("log likelihood"))
		sapply(seq_len(ncol(clu.logp)),function(j)
				{
					lines(1:nrow(clu.logp), clu.logp[,j], lty=ltys[j])
				})
		legend("bottomleft",lty=ltys, legend=r.E2E, bty='n')
		dev.off()
		
		#same but no division by number of possible spanning trees
		clu.logp	<- sapply(seq_along(clu.p),function(i)
				{
					tmp<- apply(clu.p[[i]],2,function(x)	sum(x,na.rm=1)	)				
					log( tmp / sum(tmp) )
				})
		print(clu.logp)
		
		#plot
		f.name	<- paste(dir.name,"/tipcluster_loglkl_nosize.pdf",sep='')
		cat(paste("\n plot to ",f.name))
		pdf(f.name,version="1.4",width=5,height=5)		
		
		xlim<- c(1,nrow(clu.logp))
		ylim<- range(clu.logp)
		plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim,xlab="number of transmissions in tip cluster",ylab=expression("log likelihood"))
		sapply(seq_len(ncol(clu.logp)),function(j)
				{
					lines(1:nrow(clu.logp), clu.logp[,j], lty=ltys[j])
				})
		legend("bottomright",lty=ltys, legend=r.E2E, bty='n')
		dev.off()
	}
	
	#cluster probabilities under sampling and explore loss under smaller closures
	if(0)
	{
		sampling.pa	<- c(0.05,0.1,0.2,0.3, 0.4, 0.6)
		mx.s.ntr	<- 4		
		clu.clo		<- c(6,8,10,12,14,16)
		theta		<- { tmp<- c(0.6,0.3); names(tmp)<- c("E2E","U2E"); tmp }
		
		mse<- sapply(sampling.pa, function(y)
				{
					sampling<- { tmp<- c(y,y); names(tmp)<- c("Idx","E"); tmp }									
					clu.ps	<- lapply(clu.clo,function(x)
							{
								clu.n	<- clu.n.of.tchain(x)
								clu.p	<- clu.p.of.tchain(clu.n, theta )
								clu.ps	<- clu.p.of.tchain.rnd.sampling(clu.p, sampling, mx.s.ntr, rtn.only.closure.sum=1 )
								clu.ps
							})
					clu.lps	<- clu.ps[[length(clu.ps)]]	
					mse		<- sapply(seq_along(clu.ps),function(i)
							{
								tmp<- abs(clu.ps[[i]]-clu.lps )#*( clu.ps[[i]]-clu.lps )
								sum( tmp, na.rm=1 )
							})
					names(mse)<- clu.clo
					mse
				})
		colnames(mse)<- sampling.pa
		cat("\nmean square error is\n")
		print(mse)
		
		#plot
		xlim<- range(clu.clo)
		ylim<- range(mse, na.rm=1)
		plot(1,1,type='n',bty='n',xlab="closure",ylab="total absolute error",xlim=xlim,ylim=ylim)
		sapply(seq_len(ncol(mse)),function(j)
				{
					lines(clu.clo, mse[,j],lty=j)
				})
		legend("topright",bty='n',lty=seq_len(ncol(mse)), legend=sampling.pa)
	}
	
	#cluster probabilities under sampling and explore how large tip cluster table should be
	if(0)
	{
		sampling.pa	<- c(0.05,0.1,0.2,0.3,0.4)
		mx.s.ntr	<- c(4,5,6,7,8,9)		
		clu.clo		<- 18
		theta		<- { tmp<- c(0.6,0.3); names(tmp)<- c("E2E","U2E"); tmp }
		clu.n		<- clu.n.of.tchain(clu.clo)
		clu.nchain	<- apply( clu.n, 2, sum )
		clu.p		<- clu.p.of.tchain(clu.n, theta )
		clu.p		<- apply(clu.p,2,function(x)	sum(x,na.rm=1)	) / clu.nchain
		clu.p		<- clu.p / sum( clu.p ) 
		print(clu.p)
		
		mse<- sapply(sampling.pa, function(y)
				{
					sampling<- { tmp<- c(y,y); names(tmp)<- c("Idx","E"); tmp }									
					clu.ps	<- lapply(seq_along(mx.s.ntr),function(i)
							{								
								tmp			<- clu.p.of.tipc.rnd.sampling(clu.p, sampling, mx.s.ntr[i], rtn.only.closure.sum=0 )
								c(tmp, rep(0,length(mx.s.ntr)-i))								
							})
					clu.lps	<- clu.ps[[length(clu.ps)]]	
					mse		<- sapply(seq_along(clu.ps),function(i)
							{
								tmp<- abs(clu.ps[[i]]-clu.lps )#*( clu.ps[[i]]-clu.lps )
								sum( tmp, na.rm=1 )
							})
					names(mse)<- mx.s.ntr
					mse
				})		
		colnames(mse)<- sampling.pa
		cat("\nmean square error is\n")
		print(mse)
		
		#plot
		xlim<- range(mx.s.ntr)
		ylim<- range(mse, na.rm=1)
		plot(1,1,type='n',bty='n',xlab="number of transmissions in tip cluster",ylab="total absolute error",xlim=xlim,ylim=ylim)
		sapply(seq_len(ncol(mse)),function(j)
				{
					lines(mx.s.ntr, mse[,j],lty=j)
				})
		legend("topright",bty='n',lty=seq_len(ncol(mse)), legend=sampling.pa)
	}
	
	#exp number x->E
	if(0)
	{
		ltys	<- c(4,1,2,3)
		p.U2E	<- 0.2
		r.E2E	<- c(0.5,1,2,4)
		clu.p	<- lapply(r.E2E,function(x)
				{
					theta<- p.U2E*c(x,1)
					theta<- theta / sum(theta)
					names(theta)<- c("E2E","U2E")				
					clu.p.of.tchain(clu.n, theta )								
				})
		
		print( clu.p[[1]] )
		clu.E2E<-	sapply(seq_along(clu.p),function(i)	clu.exp.X2E(clu.p[[i]])[["n.E2E"]]	)
		clu.x2E<-	sapply(seq_along(clu.p),function(i)	clu.exp.X2E(clu.p[[i]])[["n.Idx2E"]]	)
		print( clu.E2E )
		print( clu.x2E )
		#plot
		xlim<- c(1,nrow(clu.x2E))
		ylim<- range(c(0,3,clu.x2E))
		f.name	<- paste(dir.name,"/tipcluster_entr.pdf",sep='')
		cat(paste("\n plot to ",f.name))
		pdf(f.name,version="1.4",width=5,height=5)		
		
		plot(1,1,type='n',bty='n',xlim=xlim,ylim=ylim,xlab="number of transmissions in tip cluster",ylab="expected U->E ")
		sapply(seq_len(ncol(clu.x2E)),function(j)
				{
					lines(1:nrow(clu.x2E), clu.x2E[,j], lty=ltys[j])
				})
		abline(h=1,col="blue",lty=2)
		legend("bottomleft",lty=ltys, legend=r.E2E, bty='n')
		dev.off()
	}
	
	#simulate tip cluster
	if(0)
	{
		p.U2E		<- 0.4
		p.T2E		<- 0.1
		p.O2E		<- 0.1
		p.E2E		<- 0.4		
		theta		<- { tmp<- c(p.E2E, p.U2E, p.T2E, p.O2E); names(tmp)<- c("E2E","U2E","T2E","O2E"); tmp }
		incidence	<- 1000
		
		clu.n		<- clu.n.of.tchain(15)
		clu.p		<- clu.probabilities(clu.n, theta)
		clu.simulate( clu.p, incidence )
	}
	
	#simulate tip cluster under sampling
	if(1)
	{
		p.U2E		<- 0.4
		p.T2E		<- 0.1
		p.O2E		<- 0.1
		p.E2E		<- 0.4		
		theta		<- { tmp<- c(p.E2E, p.U2E, p.T2E, p.O2E); names(tmp)<- c("E2E","U2E","T2E","O2E"); tmp }
		sampling	<- { tmp<- c(0.3,0.3); names(tmp)<- c("Idx","E"); tmp }
		incidence	<- 1000		
		clu.n		<- clu.n.of.tchain(15)
		mx.s.ntr	<- ncol(clu.n)#6
		print(clu.n)
		clu.p		<- clu.probabilities(clu.n, theta)
		print(clu.p)
		clu.sim		<- clu.simulate(clu.p, incidence)
		print(clu.sim)
		stop()
		clu.sim		<- clu.sample(clu.sim, sampling, mx.s.ntr, rtn.exp=1)		
		print(clu.sim)
		stop()
		clu.p	<- lapply(r.E2E,function(x)
				{
					theta<- p.U2E*c(x,1)
					names(theta)<- c("E2E","U2E")				
					clu.p.of.tchain(clu.n, theta )								
				})
		stop()
	}
	
	#compute expected number of x->E and E->E transmissions for sampled tip cluster, and compare pi^s as sampling varies
	if(0)
	{
		p.U2E		<- 0.4
		p.T2E		<- 0.3
		p.O2E		<- 0.3
		p.E2E		<- 0		
		theta		<- { tmp<- c(p.E2E, p.U2E, p.T2E, p.O2E); names(tmp)<- c("E2E","U2E","T2E","O2E"); tmp }
		sampling.pa	<- seq(1,1,0.05)		
		incidence	<- 1000		
		clu.n		<- clu.n.of.tchain(14)
		mx.s.ntr	<- 14
		clu.p		<- clu.probabilities(clu.n, theta)		
		clu.sim		<- clu.simulate( clu.p, incidence )
		print(apply(clu.sim,1,sum))

		clu.trans	<- sapply(sampling.pa,function(x)
						{
							sampling	<- { tmp<- c(x,x); names(tmp)<- c("Idx","E"); tmp }
							clu.sim		<- clu.sample(clu.sim, sampling, rtn.exp=1)	
							clu.trans	<- clu.exp.transmissions(clu.sim, clu.n, theta, sampling, mx.s.ntr, exclude.O= 0)
							clu.trans
						})
		colnames(clu.trans)	<- sampling.pa	
		clu.E2E		<- apply(clu.trans,2,function(x)  x["E2E"]/sum(x) )		
		print(clu.trans)
		print(sum(clu.trans))
		print(clu.E2E)
	#print( clu.trans["E2E"] / sum(clu.trans) )
	}
	
	#compute expected number of x->E and E->E transmissions for sampled tip cluster,
	#and compute power based on cluster randomized t - test
	if(0)
	{
		theta.EE.H0	<- 0.15
		theta.EE.H1	<- 0.35
		theta.UE	<- 0.3
		theta.TE	<- theta.UE / 5
		theta.OE	<- 0.1
		theta.H0	<- clu.p.init(theta.EE.H0, theta.UE, theta.TE, theta.OE)
		theta.H1	<- clu.p.init(theta.EE.H1, theta.UE, theta.TE, theta.OE)				
		sampling	<- { tmp<- c(0.4,0.4); names(tmp)<- c("Idx","E"); tmp }		
		cl.inc		<- c(50,50,50)		
		clu.n		<- clu.n.of.tchain(12)
		mx.s.ntr	<- 6
		
		clu.trans.H0<- sapply(cl.inc, function(incidence)
						{
							clu.p		<- clu.probabilities(clu.n, theta.H0)
							clu.sim		<- clu.simulate(clu.p, incidence)							
							clu.sim		<- clu.sample(clu.sim, sampling, rtn.exp=1)
							tmp			<- clu.exp.transmissions(clu.sim, clu.n, theta.H0, sampling, mx.s.ntr, exclude.O= 0)
							c( tmp["E2E"] / sum(tmp), tmp["E2E"], sum(tmp) ) 
						})
		colnames(clu.trans.H0)<- cl.inc		
		rownames(clu.trans.H0)<- c("p.E2E","n.E2E","n.x2E")
		
		clu.trans.H1<- sapply(cl.inc, function(incidence)
						{
							clu.p		<- clu.probabilities(clu.n, theta.H1)
							clu.sim		<- clu.simulate(clu.p, incidence)								
							clu.sim		<- clu.sample(clu.sim, sampling, rtn.exp=1)			
							tmp			<- clu.exp.transmissions(clu.sim, clu.n, theta.H1, sampling, mx.s.ntr, exclude.O= 0)
							c( tmp["E2E"] / sum(tmp), tmp["E2E"], sum(tmp) )
						})
		colnames(clu.trans.H1)<- cl.inc
		rownames(clu.trans.H1)<- c("p.E2E","n.E2E","n.x2E")
		
		print(clu.trans.H0)
		print(clu.trans.H1)
		phdes.power.ttest.cl(clu.trans.H1["n.x2E",], clu.trans.H0["p.E2E",], clu.trans.H1["p.E2E",], 0.6*0.6 )
	}
	stop()
}
###############################################################################
prj.popart.power_test<- function()
{
	require(binom)
	require(phylodesign)
	my.mkdir(DATA,"popartpowercalc")
	dir.name<- paste(DATA,"popartpowercalc",sep='/')	
	resume<- 0
	verbose<- 1
	plot.increment<- 0.05
	
	m.type<- "Acute"	
	cohort.size<- 2500
	pc24.size<- 6000
	cohort.dur<- 3	
	theta.EE.H0<- 0.10
	theta.EE.H1<- 0.4
	test.alpha<- 0.05		 
	p.nocontam<- 	1	
	pooled.n<- 1
	p.prev.instudy.clu.armC<- 0.2
	opt.pooled<- "no pooling"#"pooled across ZA"#"pooled across trial"#"no pooling"
	opt.sampling<- "PC and HCC"#"only HCC"#"PC and HCC"
	opt.power<-	"All"
	
	debug<- 1	#for checking / debugging
	sites<-	popart.getdata.randomized.arm( pooled.n, rtn.fixed=debug)
	sites<-	popart.getdata.randomized.n(sites, cohort.size, cohort.dur, rtn.exp=debug)
	data("popart.tipcp.highacute")
	data("popart.tipcp.lowacute")	
	#print(sites)	
	###############################################################################
	#plot the biased proportions for baseline parameters
	###############################################################################
	cat("\nplot the biased proportions for baseline parameters")
	s.consent<- seq(0.4,0.8,0.05)		
	inc<- lapply(s.consent,function(x)
			{
				p.lab			<- 0.9		
				p.consent.coh	<- 0.9
				p.consent.clu	<- x
				p.vhcc.prev.AB	<- 0.95
				p.vhcc.inc.AB	<- 0.8
				p.vhcc.prev.C	<- p.prev.instudy.clu.armC
				p.vhcc.inc.C	<- p.prev.instudy.clu.armC#*0.9
				p.nocontam	 	<- 0.85
				tipc.p			<- phdes.get.hyp.tipc.probs(popart.tipcp.highacute,1-p.nocontam)
				a2a.hg			<- popart.get.sampled.acute2acute( 	sites, tipc.p, opt.sampling, rtn.int=!debug,
																	p.vhcc.prev.AB=p.vhcc.prev.AB, p.vhcc.prev.C=p.vhcc.prev.C, p.vhcc.inc.AB=p.vhcc.inc.AB, p.vhcc.inc.C=p.vhcc.inc.C, 
																	consent.PC=p.consent.coh, consent.HCC=p.consent.clu, p.lab=p.lab	)
				tipc.p			<- phdes.get.hyp.tipc.probs(popart.tipcp.lowacute,1-p.nocontam)
				a2a.lw			<- popart.get.sampled.acute2acute( 	sites, tipc.p, opt.sampling, rtn.int=!debug,
																	p.vhcc.prev.AB=p.vhcc.prev.AB, p.vhcc.prev.C=p.vhcc.prev.C, p.vhcc.inc.AB=p.vhcc.inc.AB, p.vhcc.inc.C=p.vhcc.inc.C, 
																	consent.PC=p.consent.coh, consent.HCC=p.consent.clu, p.lab=p.lab	)																
				list(a2a.lw=a2a.lw, a2a.hg=a2a.hg)				
			})
	names(inc)<- s.consent
	x2i.lw	<- sapply(inc,function(x)		x[["a2a.lw"]]["x2i.s",])
	i2i.lw	<- sapply(inc,function(x)		x[["a2a.lw"]]["i2i.s",])	
	x2i.hg	<- sapply(inc,function(x)		x[["a2a.hg"]]["x2i.s",])
	i2i.hg	<- sapply(inc,function(x)		x[["a2a.hg"]]["i2i.s",])
	#print( i2i.lw )
	#print( x2i.lw )	
	# pool transmissions
	x2i.lw	<- popart.pool(sites, x2i.lw, method=opt.pooled)[["transm"]]
	i2i.lw	<- popart.pool(sites, i2i.lw, method=opt.pooled)[["transm"]]
	x2i.hg	<- popart.pool(sites, x2i.hg, method=opt.pooled)[["transm"]]
	tmp		<- popart.pool(sites, i2i.hg, method=opt.pooled)
	i2i.hg	<- tmp[["transm"]]
	idx.A	<- tmp[["idx.A"]]
	idx.C	<- tmp[["idx.C"]]
	#print( i2i.lw[idx.C,,drop=0] )
	#print( x2i.lw[idx.C,,drop=0] )
	#compute biased proportions
	test.biased.H0.A<-	apply(i2i.lw[idx.A,,drop=0],2,mean)/apply(x2i.lw[idx.A,,drop=0],2,mean)	
	test.biased.H1.A<-	apply(i2i.hg[idx.A,,drop=0],2,mean)/apply(x2i.hg[idx.A,,drop=0],2,mean)
	test.biased.H0.C<-	apply(i2i.lw[idx.C,,drop=0],2,mean)/apply(x2i.lw[idx.C,,drop=0],2,mean)	
	test.biased.H1.C<-	apply(i2i.hg[idx.C,,drop=0],2,mean)/apply(x2i.hg[idx.C,,drop=0],2,mean)			
	#print( test.biased.H0.C )
	#stop()
	#plot proportion I->I sampled against sampling fraction for arm A
	f.name<- paste(dir.name,paste("SAMPLINGBIAS_TIPC_A_visitHCC.C",p.prev.instudy.clu.armC,"power",opt.power,"pool",opt.pooled,"sample",opt.sampling,"pwcalc",theta.EE.H0,theta.EE.H1,test.alpha,".pdf",sep='_'),sep='/')
	cat(paste("\nplot sampling bias A to\n",f.name))
	pdf(paste(f.name),version="1.4",width=6,height=6)
	ylim<- range(c(test.biased.H0.A,test.biased.H1.A,theta.EE.H0,theta.EE.H1))
	plot(1,1,type='n',xlim=range(s.consent),ylim=ylim, xlab="%consenting to phyl study at HCC",ylab="proportion acute to acute" )
	abline(h=theta.EE.H0, col="red")
	abline(h=theta.EE.H1, col="blue")
	lines(s.consent, test.biased.H0.A, col="red", lty=2)
	lines(s.consent, test.biased.H1.A, col="blue", lty=2)
	legend(x=0.38,y=0.4,bty='n',legend=c("a2a 40%","a2a 10%"), fill=c("blue","red"), border=NA)
	dev.off()
	#plot proportion I->I sampled against sampling fraction for arm C
	f.name<- paste(dir.name,paste("SAMPLINGBIAS_TIPC_C_visitHCC.C",p.prev.instudy.clu.armC,"power",opt.power,"pool",opt.pooled,"sample",opt.sampling,"pwcalc",theta.EE.H0,theta.EE.H1,test.alpha,".pdf",sep='_'),sep='/')
	cat(paste("\nplot sampling bias C to\n",f.name))
	pdf(paste(f.name),version="1.4",width=6,height=6)
	ylim<- range(c(test.biased.H0.C,test.biased.H1.C,theta.EE.H0,theta.EE.H1))
	plot(1,1,type='n',xlim=range(s.consent),ylim=ylim, xlab="%consenting to phyl study at HCC",ylab="proportion acute to acute" )
	abline(h=theta.EE.H0, col="red")
	abline(h=theta.EE.H1, col="blue")
	lines(s.consent, test.biased.H0.C, col="red", lty=2)
	lines(s.consent, test.biased.H1.C, col="blue", lty=2)
	legend(x=0.38,y=0.4,bty='n',legend=c("a2a 40%","a2a 10%"), fill=c("blue","red"), border=NA)
	dev.off()
	
	###############################################################################
	#check that the biased proportions approach the correct ones as the sampling approaches one.
	###############################################################################
	cat("\ncheck that the biased proportions approach the correct ones as the sampling approaches one")
	s.consent<- seq(0.95,1,0.05)		
	inc<- lapply(s.consent,function(x)
			{
				p.lab			<- 1		
				p.consent.coh	<- 1
				p.consent.clu	<- x
				p.vhcc.prev.AB	<- 1
				p.vhcc.inc.AB	<- 1
				p.vhcc.prev.C	<- 1
				p.vhcc.inc.C	<- 1
				p.nocontam	 	<- 1
				tipc.p			<- phdes.get.hyp.tipc.probs(popart.tipcp.highacute,1-p.nocontam)
				a2a.hg			<- popart.get.sampled.acute2acute( 	sites, tipc.p, opt.sampling, rtn.int=!debug,
						p.vhcc.prev.AB=p.vhcc.prev.AB, p.vhcc.prev.C=p.vhcc.prev.C, p.vhcc.inc.AB=p.vhcc.inc.AB, p.vhcc.inc.C=p.vhcc.inc.C, 
						consent.PC=p.consent.coh, consent.HCC=p.consent.clu, p.lab=p.lab	)
				tipc.p			<- phdes.get.hyp.tipc.probs(popart.tipcp.lowacute,1-p.nocontam)
				a2a.lw			<- popart.get.sampled.acute2acute( 	sites, tipc.p, opt.sampling, rtn.int=!debug,
						p.vhcc.prev.AB=p.vhcc.prev.AB, p.vhcc.prev.C=p.vhcc.prev.C, p.vhcc.inc.AB=p.vhcc.inc.AB, p.vhcc.inc.C=p.vhcc.inc.C, 
						consent.PC=p.consent.coh, consent.HCC=p.consent.clu, p.lab=p.lab	)																
				list(a2a.lw=a2a.lw, a2a.hg=a2a.hg)				
			})
	names(inc)<- s.consent
	x2i.lw	<- sapply(inc,function(x)		x[["a2a.lw"]]["x2i.s",])
	i2i.lw	<- sapply(inc,function(x)		x[["a2a.lw"]]["i2i.s",])	
	x2i.hg	<- sapply(inc,function(x)		x[["a2a.hg"]]["x2i.s",])
	i2i.hg	<- sapply(inc,function(x)		x[["a2a.hg"]]["i2i.s",])		
	test.biased.H0	<- apply(i2i.lw,2,mean)/apply(x2i.lw,2,mean)
	test.biased.H1	<- apply(i2i.hg,2,mean)/apply(x2i.hg,2,mean)
	#plot proportion I->I sampled against sampling fraction
	f.name<- paste(dir.name,paste("SAMPLINGBIAS_TIPC_approacges_one.pdf",sep='_'),sep='/')
	cat(paste("\nplot sampling bias approaches one to\n",f.name))
	pdf(paste(f.name),version="1.4",width=6,height=6)
	ylim<- range(c(test.biased.H0,test.biased.H1,theta.EE.H0,theta.EE.H1))
	plot(1,1,type='n',xlim=range(s.consent),ylim=ylim, xlab="%consenting to phyl study at HCC",ylab="proportion acute to acute" )
	abline(h=theta.EE.H0, col="red")
	abline(h=theta.EE.H1, col="blue")
	lines(s.consent, test.biased.H0, col="red", lty=2)
	lines(s.consent, test.biased.H1, col="blue", lty=2)
	legend(x=0.38,y=0.4,bty='n',legend=c("a2a 40%","a2a 10%"), fill=c("blue","red"), border=NA)
	dev.off()
	###############################################################################
	#check that the x2i's are all consistent -- we have 3 calculations
	###############################################################################
	cat("\nccheck that the x2i's are all consistently calculated")
	x2i<- sapply(s.consent,function(x)
			{
				p.lab			<- 1			
				p.consent.coh	<- 1
				p.consent.clu	<- x
				p.vhcc.prev.AB	<- 1
				p.vhcc.inc.AB	<- 1
				p.vhcc.prev.C	<- 1
				p.vhcc.inc.C	<- 1				
				popart.get.sampled.transmissions(	sites, opt.sampling, rtn.int=!debug,
						p.vhcc.prev.AB=p.vhcc.prev.AB, p.vhcc.prev.C=p.vhcc.prev.C, p.vhcc.inc.AB=p.vhcc.inc.AB, p.vhcc.inc.C=p.vhcc.inc.C, 
						consent.PC=p.consent.coh, consent.HCC=p.consent.clu, p.lab=p.lab, p.community=p.nocontam
				)																				
			})
	colnames(x2i)<- s.consent
	
	sapply(seq_len(nrow(x2i)),function(i)
			{
				if(any(round(x2i.lw[i,],d=3)!=round(x2i.hg[i,],d=3)))
					stop("x2i.lw and x2i.hg different ?")
				if(any(round(x2i[i,],d=3)!=round(x2i.hg[i,],d=3)))
					stop("x2i and x2i.hg different ?")				
			})
	stop()	
}
###############################################################################