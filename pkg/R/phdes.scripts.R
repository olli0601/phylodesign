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
	
	m.type<- "Acute"	
	cohort.size<- 2500
	pc24.size<- 6000
	cohort.dur<- 3	
	test.prop0<- 0.10
	test.prop1<- 0.4
	test.alpha<- 0.05		
	pooled.n<- 200
	#opt.pooled<- "no pooling"#"pooled across ZA"#"pooled across trial"#"no pooling"
	#opt.sampling<- "PC and HCC"#"only HCC"#"PC and HCC"
	opt.power<-	"All"	
	if(!opt.sampling%in%c("PC and HCC","only HCC"))
		stop("prj.popart.powercalc_medsize: unknown method to sample")
	if(!opt.pooled%in%c("pooled across country","pooled across ZA","pooled across SA","pooled across trial","no pooling"))
		stop("prj.popart.powercalc_medsize: unknown method to pool")
	if(!opt.power%in%c("All","PonlyPC","IonlyPC","PonlyPCandIonlyPC"))
		stop("prj.popart.powercalc_medsize: unknown method to pool")
	cat(paste("\ncohort.size",cohort.size))
	cat(paste("\ncohort.dur",cohort.dur))
	cat(paste("\ntest.prop0",test.prop0))
	cat(paste("\ntest.prop1",test.prop1))
	cat(paste("\ntest.alpha",test.alpha))	
	cat(paste("\np.nocontam",p.nocontam))
	cat(paste("\np.prev.instudy.clu.armC",p.prev.instudy.clu.armC))
	cat(paste("\npooled.n",pooled.n))
	cat(paste("\nopt.pooled",opt.pooled))
	cat(paste("\nopt.power",opt.power))
	cat(paste("\nopt.sampling",opt.sampling))
	
	sites<-	popart.getdata.randomized.arm( pooled.n )
	sites<-	popart.getdata.randomized.n(sites, cohort.size, cohort.dur)
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
	cat("\npool linked transmissions")
	tmp		<- popart.pool(sites, x2i, method=opt.pooled)
	x2i		<- tmp[["transm"]]
	idx.A	<- tmp[["idx.A"]]
	idx.C	<- tmp[["idx.C"]]
#print(x2i); print(idx.A); print(idx.C)	
	
	cat("\ncompute power for dectecting differences in sampled acute to acute transmissions")
	tmp<- phdes.binom.power(	x2i, round(x2i*test.prop0), idx.A, idx.C, test.prop0, test.prop1, test.alpha, verbose=0)
	conf.lw.med.armA	<- tmp[["conf.A"]] 
	conf.lw.med.armC	<- tmp[["conf.C"]] 
	
	tmp<- phdes.binom.power(	x2i, round(x2i*test.prop1), idx.A, idx.C, test.prop0, test.prop1, test.alpha, verbose=0)
	conf.hg.med.armA	<- tmp[["conf.A"]]
	conf.hg.med.armC	<- tmp[["conf.C"]]
	is.conf.hg.med.armA	<- tmp[["is.conf.A"]] 
	is.conf.hg.med.armC	<- tmp[["is.conf.C"]] 
	power.hg.med.armA	<- tmp[["power.A"]]
	power.hg.med.armC	<- tmp[["power.C"]] 
		
	if(p.nocontam>=0.85)	
		cols<- c("deepskyblue","dodgerblue4")
	else					
		cols<- c("firebrick1","firebrick4")
	
	
	f.name<- paste(dir.name,paste("CFLINK_consent",p.nocontam,"phsig",p.phylosignal,"power",opt.power,"pool",opt.pooled,"sample",opt.sampling,"pwcalc",test.prop0,test.prop1,test.alpha,".pdf",sep='_'),sep='/')	
	cat(paste("\nplot power to\n",f.name))
	print("HERE")
	pdf(paste(f.name),version="1.4",width=6,height=6)
	phdes.plot.power(	power.hg.med.armA, power.hg.med.armC, is.conf.hg.med.armA, is.conf.hg.med.armC,
			f.name, "% consenting to ph study at HCC", 
			paste("power to distinguish acute < ",test.prop0*100,"% vs > ",test.prop1*100,"%\n",opt.pooled,sep=''),
			c("arm A","arm C",paste("contamination",(1-p.nocontam)*100,"%, linked to HCC/C",p.prev.instudy.clu.armC*100,"%")),
			cols=cols,"bottomright", verbose= 0	)
	dev.off()
	
	f.name<- paste(dir.name,paste("CFLINK_consent",p.nocontam,"phsig",p.phylosignal,"power",opt.power,"pool",opt.pooled,"sample",opt.sampling,"A_confint",test.prop0,test.prop1,test.alpha,".pdf",sep='_'),sep='/')
	cat(paste("\nplot binomial confidence intervals to\n",f.name))
	pdf(paste(f.name),version="1.4",width=6,height=6)
	phdes.plot.confint(	rep(test.prop0,nrow(conf.lw.med.armA)), rep(test.prop1,nrow(conf.lw.med.armA)),
			conf.lw.med.armA, conf.hg.med.armA,			
			f.name=f.name, xlab="% consenting to ph study at HCC",ylab=paste("estimated proportion acute to acute transmission\n arm A,",opt.pooled),
			legend.loc="topright",legend.txt=c(paste("true prop",test.prop0*100,"%",sep=' '), paste("true prop",test.prop1*100,"%",sep=' ')),cols=cols		)													
	dev.off()
	
	f.name<- paste(dir.name,paste("CFLINK_consent",p.nocontam,"phsig",p.phylosignal,"power",opt.power,"pool",opt.pooled,"sample",opt.sampling,"C_confint",test.prop0,test.prop1,test.alpha,".pdf",sep='_'),sep='/')
	pdf(paste(f.name),version="1.4",width=6,height=6)
	phdes.plot.confint(	rep(test.prop0,nrow(conf.lw.med.armC)), rep(test.prop1,nrow(conf.lw.med.armC)),
			conf.lw.med.armC, conf.hg.med.armC,			
			f.name=f.name, xlab="% consenting to ph study at HCC",ylab=paste("estimated proportion acute to acute transmission\n arm C,",opt.pooled),
			legend.loc="topright",legend.txt=c(paste("true prop",test.prop0*100,"%",sep=' '), paste("true prop",test.prop1*100,"%",sep=' ')),cols=cols		)
	dev.off()
	
	f.name<- paste(dir.name,paste("CFLINK_consent",p.nocontam,"phsig",p.phylosignal,"power",opt.power,"pool",opt.pooled,"sample",opt.sampling,"pwcalc",test.prop0,test.prop1,test.alpha,".R",sep='_'),sep='/')
	cat(paste("\nsave R objects to\n",f.name))
	save(sites, x2i, idx.A, idx.C, power.hg.med.armA, power.hg.med.armC, is.conf.hg.med.armA, is.conf.hg.med.armC, conf.lw.med.armA, conf.lw.med.armC, conf.hg.med.armA, conf.hg.med.armC, file=f.name)
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
	test.prop0	<- 0.1
	test.prop1	<- 0.4
	test.alpha	<- 0.05		
	
	pooled.n<- 500
	#opt.pooled<- "pooled across ZA"#"pooled across trial"#"no pooling"
	#opt.sampling<- "PC and HCC"#"only HCC"#"PC and HCC"
	opt.power<-	"All"
	###############################################################################
	#set arguments
	###############################################################################
	
	if(!opt.sampling%in%c("PC and HCC","only HCC"))
		stop("prj.popart.powercalc_medsize: unknown method to sample")
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
	cat(paste("\ntest.prop0",test.prop0))
	cat(paste("\ntest.prop1",test.prop1))
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
	x2i.lw	<- popart.pool(sites, x2i.lw, method=opt.pooled)[["transm"]]
	i2i.lw	<- popart.pool(sites, i2i.lw, method=opt.pooled)[["transm"]]
	x2i.hg	<- popart.pool(sites, x2i.hg, method=opt.pooled)[["transm"]]
	i2i.hg	<- popart.pool(sites, i2i.hg, method=opt.pooled)[["transm"]]	
	tmp		<- popart.pool(sites, x2i, method=opt.pooled)
	x2i		<- tmp[["transm"]]
	idx.A	<- tmp[["idx.A"]]
	idx.C	<- tmp[["idx.C"]]
	#compute test.biased.H0 as fraction over means to see a pattern	
	test.biased.H0		<- apply(i2i.lw[c(idx.A,idx.C),,drop=0],2,mean)/apply(x2i.lw[c(idx.A,idx.C),,drop=0],2,mean)	
	test.biased.H1		<- apply(i2i.hg[c(idx.A,idx.C),,drop=0],2,mean)/apply(x2i.hg[c(idx.A,idx.C),,drop=0],2,mean)
	test.biased.H0.A	<- apply(i2i.lw[idx.A,,drop=0],2,mean)/apply(x2i.lw[idx.A,,drop=0],2,mean)	
	test.biased.H1.A	<- apply(i2i.hg[idx.A,,drop=0],2,mean)/apply(x2i.hg[idx.A,,drop=0],2,mean)
	test.biased.H0.C	<- apply(i2i.lw[idx.C,,drop=0],2,mean)/apply(x2i.lw[idx.C,,drop=0],2,mean)	
	test.biased.H1.C	<- apply(i2i.hg[idx.C,,drop=0],2,mean)/apply(x2i.hg[idx.C,,drop=0],2,mean)	
	###############################################################################
	#plot number of acute to acute transmissions under high scenario
	#adjust for sampling bias
	###############################################################################
	cat("\nplot number of acute to acute transmissions under high scenario\nadjust for sampling bias")
	f.name<- paste(dir.name,paste("VARYLINKAGE_LINK_a2a_nocontam",p.nocontam,"power",opt.power,"pool",opt.pooled,"sample",opt.sampling,"pwcalc",test.prop0,test.prop1,test.alpha,".pdf",sep='_'),sep='/')
	cat(paste("\nplot to\n",f.name))
	pdf(paste(f.name),version="1.4",width=6,height=6)
	par(mar=c(5,5.5,0.5,2))
	plot(1,1,type='n',xlim=range(p.phylosignal),ylim=range(round(x2i*test.biased.H1.A)),xlab="% linked w phylogenetic method",ylab=paste("acute 2 acute, high scenario\n",opt.pooled))
	lines(p.phylosignal, apply(round(x2i[c(idx.A,idx.C),,drop=0]*test.biased.H1),2,mean))
	abline(h=apply(i2i.hg[c(idx.A,idx.C),,drop=0],2,mean), lty=2)
	lines(p.phylosignal, apply(round(x2i[idx.A,,drop=0]*test.biased.H1.A),2,mean),col="red")
	abline(h=apply(i2i.hg[idx.A,,drop=0],2,mean), lty=2,col="red")
	lines(p.phylosignal, apply(round(x2i[idx.C,,drop=0]*test.biased.H1.C),2,mean),col="blue")
	abline(h=apply(i2i.hg[idx.C,,drop=0],2,mean), lty=2,col="blue")
	legend("topleft",fill=c("black","red","blue"),legend=c("overall","arm A","arm C"),bty='n',border=NA)
	legend("topright",lty=c(1,2),legend=c("linked","tipc"),bty='n')
	dev.off()
	stop()
	###############################################################################
	#compare power - adjust for sampling bias
	###############################################################################
	if(p.nocontam>=0.85)	
		cols<- c("deepskyblue","dodgerblue4")
	else					
		cols<- c("firebrick1","firebrick4")
	
	
	cat("\ncompute power for dectecting differences in sampled acute to acute transmissions")	
	tmp<- phdes.binom.power(x2i, round(x2i*test.biased.H1.A), idx.A, idx.C, test.biased.H0.A, test.biased.H1.A, test.biased.H0.C, test.biased.H1.C, test.alpha, verbose=0)
	conf.hg.med.armA	<- tmp[["conf.A"]]
	is.conf.hg.med.armA	<- tmp[["is.conf.A"]]
	power.hg.med.armA	<- tmp[["power.A"]]
	tmp<- phdes.binom.power(x2i, round(x2i*test.biased.H1.C), idx.A, idx.C, test.biased.H0.A, test.biased.H1.A, test.biased.H0.C, test.biased.H1.C, test.alpha, verbose=0)
	conf.hg.med.armC	<- tmp[["conf.C"]]	 
	is.conf.hg.med.armC	<- tmp[["is.conf.C"]] 	
	power.hg.med.armC	<- tmp[["power.C"]] 		
	tmp<- phdes.binom.power(x2i.hg, i2i.hg, idx.A, idx.C, test.biased.H0.A, test.biased.H1.A, test.biased.H0.C, test.biased.H1.C, test.alpha, verbose=0)	 
		
	f.name<- paste(dir.name,paste("VARYLINKAGE_LINK_power_nocontam",p.nocontam,"power",opt.power,"pool",opt.pooled,"sample",opt.sampling,"pwcalc",test.prop0,test.prop1,test.alpha,".pdf",sep='_'),sep='/')	
	cat(paste("\nplot power to\n",f.name))
	pdf(paste(f.name),version="1.4",width=6,height=6)
	phdes.plot.power(	power.hg.med.armA, power.hg.med.armC, is.conf.hg.med.armA, is.conf.hg.med.armC,
			"% linked w phylogenetic method", 
			paste("power to distinguish acute < ",test.prop0*100,"% vs > ",test.prop1*100,"%\n",opt.pooled,sep=''),
			c("arm A","arm C",paste("contamination",(1-p.nocontam)*100,"%, linked to HCC/C",p.prev.instudy.clu.armC*100,"%")),
			cols=cols,"bottomright", verbose= 0	)
	abline(h=tmp[["power.A"]], col=cols[1])
	abline(h=tmp[["power.C"]], col=cols[2])
	dev.off()
	
	stop()													
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
	
	m.type<- "Acute"	
	cohort.size<- 2500
	pc24.size<- 6000
	cohort.dur<- 3	
	test.prop0<- 0.10
	test.prop1<- 0.4
	test.alpha<- 0.05		 
	debug<- 0
	pooled.n<- 1
	opt.pooled<- "pooled across trial"#"pooled across ZA"#"pooled across trial"#"no pooling"
	opt.sampling<- "only HCC"#"PC and HCC"
	opt.power<-	"All"
	if(!opt.sampling%in%c("PC and HCC","only HCC"))
		stop("prj.popart.powercalc_medsize: unknown method to sample")
	if(!opt.pooled%in%c("pooled across country","pooled across ZA","pooled across SA","pooled across trial","no pooling"))
		stop("prj.popart.powercalc_medsize: unknown method to pool")
	if(!opt.power%in%c("All","PonlyPC","IonlyPC","PonlyPCandIonlyPC"))
		stop("prj.popart.powercalc_medsize: unknown method to pool")
	cat(paste("\ncohort.size",cohort.size))
	cat(paste("\ncohort.dur",cohort.dur))
	cat(paste("\ntest.prop0",test.prop0))
	cat(paste("\ntest.prop1",test.prop1))
	cat(paste("\ntest.alpha",test.alpha))	
	cat(paste("\np.nocontam",p.nocontam))
	cat(paste("\npooled.n",pooled.n))
	cat(paste("\nopt.pooled",opt.pooled))
	cat(paste("\nopt.power",opt.power))
	cat(paste("\nopt.sampling",opt.sampling))
	
	sites<-	popart.getdata.randomized.arm( pooled.n, rtn.fixed=debug )
	sites<-	popart.getdata.randomized.n(sites, cohort.size, cohort.dur, rtn.exp=debug)
	data("popart.tipcp.highacute")
	data("popart.tipcp.lowacute")
		
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
				x2i.lw	<- popart.pool(sites, x2i.lw, method=opt.pooled)[["transm"]]
				i2i.lw	<- popart.pool(sites, i2i.lw, method=opt.pooled)[["transm"]]
				x2i.hg	<- popart.pool(sites, x2i.hg, method=opt.pooled)[["transm"]]
				tmp		<- popart.pool(sites, i2i.hg, method=opt.pooled)
				i2i.hg	<- tmp[["transm"]]
				idx.A	<- tmp[["idx.A"]]
				idx.C	<- tmp[["idx.C"]]
				
				test.biased.H0.A<-	apply(i2i.lw[idx.A,,drop=0],2,mean)/apply(x2i.lw[idx.A,,drop=0],2,mean)	
				test.biased.H1.A<-	apply(i2i.hg[idx.A,,drop=0],2,mean)/apply(x2i.hg[idx.A,,drop=0],2,mean)
				test.biased.H0.C<-	apply(i2i.lw[idx.C,,drop=0],2,mean)/apply(x2i.lw[idx.C,,drop=0],2,mean)	
				test.biased.H1.C<-	apply(i2i.hg[idx.C,,drop=0],2,mean)/apply(x2i.hg[idx.C,,drop=0],2,mean)	
				
				# compute power and confidence intervals
				tmp<- phdes.binom.power(x2i.lw, i2i.lw, idx.A, idx.C, test.biased.H0.A, test.biased.H1.A, test.biased.H0.C, test.biased.H1.C, test.alpha, verbose=0, method.ci="asymptotic",method.pw="cloglog")
				conf.lw.A				<- tmp[["conf.A"]] 
				conf.lw.C				<- tmp[["conf.C"]] 				
				tmp<- phdes.binom.power(x2i.hg, i2i.hg, idx.A, idx.C, test.biased.H0.A, test.biased.H1.A, test.biased.H0.C, test.biased.H1.C, test.alpha, verbose=0, method.ci="asymptotic",method.pw="cloglog")	 
				is.conf.hg.med.armA	<- tmp[["is.conf.A"]] 
				is.conf.hg.med.armC	<- tmp[["is.conf.C"]] 
				power.hg.med.armA	<- tmp[["power.A"]]
				power.hg.med.armC	<- tmp[["power.C"]]
				conf.hg.A				<- tmp[["conf.A"]] 
				conf.hg.C				<- tmp[["conf.C"]] 
								
				tmp<- list(	i2i.hg.A= apply(i2i.hg[idx.A,,drop=0],2,mean), 
							i2i.hg.C= apply(i2i.hg[idx.C,,drop=0],2,mean),
							x2i.hg.A= apply(x2i.hg[idx.A,,drop=0],2,mean), 
							x2i.hg.C= apply(x2i.hg[idx.C,,drop=0],2,mean),
							test.biased.H0.A= test.biased.H0.A,
							test.biased.H1.A= test.biased.H1.A,
							test.biased.H0.C= test.biased.H0.C,
							test.biased.H1.C= test.biased.H1.C,							
							is.conf.A= is.conf.hg.med.armA, 
							is.conf.C= is.conf.hg.med.armC,
							power.A= power.hg.med.armA, 
							power.C= power.hg.med.armC,
							conf.hg.A= conf.hg.A,
							conf.hg.C= conf.hg.C,
							conf.lw.A= conf.lw.A,
							conf.lw.C= conf.lw.C
							)	
				tmp
			})
	names(out)<- p.vhcc.prev.Cs
	i2i.hg.A	<- sapply(out,function(x)		x[["i2i.hg.A"]]	)
	i2i.hg.C	<- sapply(out,function(x)		x[["i2i.hg.C"]]	)
	power.A		<- sapply(out,function(x)		x[["power.A"]]	)
	power.C		<- sapply(out,function(x)		x[["power.C"]]	)
	conf.hg.A.u	<- sapply(out,function(x)		x[["conf.hg.A"]][,"upper"]	)
	conf.hg.C.u	<- sapply(out,function(x)		x[["conf.hg.C"]][,"upper"]	)
	conf.hg.A.l	<- sapply(out,function(x)		x[["conf.hg.A"]][,"lower"]	)
	conf.hg.C.l	<- sapply(out,function(x)		x[["conf.hg.C"]][,"lower"]	)
	conf.lw.A.u	<- sapply(out,function(x)		x[["conf.lw.A"]][,"upper"]	)
	conf.lw.C.u	<- sapply(out,function(x)		x[["conf.lw.C"]][,"upper"]	)
	conf.lw.A.l	<- sapply(out,function(x)		x[["conf.lw.A"]][,"lower"]	)
	conf.lw.C.l	<- sapply(out,function(x)		x[["conf.lw.C"]][,"lower"]	)	
	###############################################################################
	#plot sampled a2a
	###############################################################################
	if(1)
	{
		require(fields)
		f.name<- paste(dir.name,paste("VARYCONSENT_TIPC_a2a_",p.nocontam,"power",opt.power,"pool",opt.pooled,"sample",opt.sampling,"pwcalc",test.prop0,test.prop1,test.alpha,".pdf",sep='_'),sep='/')
		cat(paste("\nplot a2a to\n",f.name))
		pdf(paste(f.name),version="1.4",width=12,height=6)		
		breaks		<- diff(range(c(i2i.hg.A,i2i.hg.C)))/50
		breaks		<- seq(min(c(i2i.hg.A,i2i.hg.C)),by=breaks, len=51)	
		def.par <- par(no.readonly = TRUE)
		layout.m<- matrix(data= c(1,2),ncol=2,nrow=1,byrow=1)
		layout(layout.m)		
		image(main="arm C",p.consent.clus,p.vhcc.prev.Cs,i2i.hg.C, breaks=breaks, col=head( rev(gray(seq(0,.95,len=trunc(50*1.4)))), 50))
		image.plot(p.consent.clus,p.vhcc.prev.Cs,i2i.hg.A, breaks=breaks, col=head( rev(gray(seq(0,.95,len=trunc(50*1.4)))), 50),main="arm A")	
		par(def.par)
		dev.off()
	}
	###############################################################################
	#plot panel of confidence intervals
	###############################################################################
	if(1)
	{
		cols<- c("deepskyblue","dodgerblue4")
		f.name<- paste(dir.name,paste("VARYCONSENT_TIPC_confint_C_",p.nocontam,"power",opt.power,"pool",opt.pooled,"sample",opt.sampling,"pwcalc",test.prop0,test.prop1,test.alpha,".pdf",sep='_'),sep='/')
		cat(paste("\nplot confint C to\n",f.name))
		pdf(paste(f.name),version="1.4",width=6,height=12)
		phdes.plot.confint.panel(conf.lw.C.l,conf.lw.C.u,conf.hg.C.l,conf.hg.C.u,p.vhcc.prev.Cs,p.consent.clus,"p.vhcc.prev.Cs","p.consent.clus", cols=cols)
		dev.off()
		
		f.name<- paste(dir.name,paste("VARYCONSENT_TIPC_confint_A_",p.nocontam,"power",opt.power,"pool",opt.pooled,"sample",opt.sampling,"pwcalc",test.prop0,test.prop1,test.alpha,".pdf",sep='_'),sep='/')
		cat(paste("\nplot confint A to\n",f.name))
		pdf(paste(f.name),version="1.4",width=6,height=12)
		phdes.plot.confint.panel(conf.lw.A.l,conf.lw.A.u,conf.hg.A.l,conf.hg.A.u,p.vhcc.prev.Cs,p.consent.clus,"p.vhcc.prev.Cs","p.consent.clus", cols=cols)
		dev.off()
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
	test.prop0<- 0.10
	test.prop1<- 0.4
	test.alpha<- 0.05		 
	debug<- 0
	pooled.n<- 200
	opt.pooled<- "pooled across T5"#"pooled across ZA"#"pooled across trial"#"no pooling"
	opt.sampling<- "PC and HCC"#"only HCC"#"PC and HCC"
	opt.power<-	"All"
	if(!opt.sampling%in%c("PC and HCC","only HCC"))
		stop("prj.popart.powercalc_medsize: unknown method to sample")
	if(!opt.pooled%in%c("pooled across country","pooled across ZA","pooled across SA","pooled across trial","pooled across T5","no pooling"))
		stop("prj.popart.powercalc_medsize: unknown method to pool")
	if(!opt.power%in%c("All","PonlyPC","IonlyPC","PonlyPCandIonlyPC"))
		stop("prj.popart.powercalc_medsize: unknown method to pool")
	cat(paste("\ncohort.size",cohort.size))
	cat(paste("\ncohort.dur",cohort.dur))
	cat(paste("\ntest.prop0",test.prop0))
	cat(paste("\ntest.prop1",test.prop1))
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
				x2i.lw	<- popart.pool(sites, x2i.lw, method=opt.pooled)[["transm"]]
				i2i.lw	<- popart.pool(sites, i2i.lw, method=opt.pooled)[["transm"]]
				x2i.hg	<- popart.pool(sites, x2i.hg, method=opt.pooled)[["transm"]]
				tmp		<- popart.pool(sites, i2i.hg, method=opt.pooled)
				i2i.hg	<- tmp[["transm"]]
				idx.A	<- tmp[["idx.A"]]
				idx.C	<- tmp[["idx.C"]]
				test.biased.H0.A<-	apply(i2i.lw[idx.A,,drop=0],2,mean)/apply(x2i.lw[idx.A,,drop=0],2,mean)	
				test.biased.H1.A<-	apply(i2i.hg[idx.A,,drop=0],2,mean)/apply(x2i.hg[idx.A,,drop=0],2,mean)
				test.biased.H0.C<-	apply(i2i.lw[idx.C,,drop=0],2,mean)/apply(x2i.lw[idx.C,,drop=0],2,mean)	
				test.biased.H1.C<-	apply(i2i.hg[idx.C,,drop=0],2,mean)/apply(x2i.hg[idx.C,,drop=0],2,mean)	
				
				# compute power and confidence intervals
				tmp<- phdes.binom.power(x2i.lw, i2i.lw, idx.A, idx.C, test.biased.H0.A, test.biased.H1.A, test.biased.H0.C, test.biased.H1.C, test.alpha, verbose=0, method.ci="asymptotic",method.pw="cloglog")
				conf.lw.A				<- tmp[["conf.A"]] 
				conf.lw.C				<- tmp[["conf.C"]] 				
				tmp<- phdes.binom.power(x2i.hg, i2i.hg, idx.A, idx.C, test.biased.H0.A, test.biased.H1.A, test.biased.H0.C, test.biased.H1.C, test.alpha, verbose=0, method.ci="asymptotic",method.pw="cloglog")	 
				is.conf.hg.med.armA	<- tmp[["is.conf.A"]] 
				is.conf.hg.med.armC	<- tmp[["is.conf.C"]] 
				power.hg.med.armA	<- tmp[["power.A"]]
				power.hg.med.armC	<- tmp[["power.C"]]
				conf.hg.A				<- tmp[["conf.A"]] 
				conf.hg.C				<- tmp[["conf.C"]] 
				
				tmp<- list(	i2i.hg.A= apply(i2i.hg[idx.A,,drop=0],2,mean), 
						i2i.hg.C= apply(i2i.hg[idx.C,,drop=0],2,mean),
						x2i.hg.A= apply(x2i.hg[idx.A,,drop=0],2,mean), 
						x2i.hg.C= apply(x2i.hg[idx.C,,drop=0],2,mean),
						test.biased.H0.A= test.biased.H0.A,
						test.biased.H1.A= test.biased.H1.A,
						test.biased.H0.C= test.biased.H0.C,
						test.biased.H1.C= test.biased.H1.C,							
						is.conf.A= is.conf.hg.med.armA, 
						is.conf.C= is.conf.hg.med.armC,
						power.A= power.hg.med.armA, 
						power.C= power.hg.med.armC,
						conf.hg.A= conf.hg.A,
						conf.hg.C= conf.hg.C,
						conf.lw.A= conf.lw.A,
						conf.lw.C= conf.lw.C
				)	
				tmp
			})
	names(out)<- p.nocontams
	i2i.hg.A	<- sapply(out,function(x)		x[["i2i.hg.A"]]	)
	i2i.hg.C	<- sapply(out,function(x)		x[["i2i.hg.C"]]	)
	power.A		<- sapply(out,function(x)		x[["power.A"]]	)
	power.C		<- sapply(out,function(x)		x[["power.C"]]	)
	conf.hg.A.u	<- sapply(out,function(x)		x[["conf.hg.A"]][,"upper"]	)
	conf.hg.C.u	<- sapply(out,function(x)		x[["conf.hg.C"]][,"upper"]	)
	conf.hg.A.l	<- sapply(out,function(x)		x[["conf.hg.A"]][,"lower"]	)
	conf.hg.C.l	<- sapply(out,function(x)		x[["conf.hg.C"]][,"lower"]	)
	conf.lw.A.u	<- sapply(out,function(x)		x[["conf.lw.A"]][,"upper"]	)
	conf.lw.C.u	<- sapply(out,function(x)		x[["conf.lw.C"]][,"upper"]	)
	conf.lw.A.l	<- sapply(out,function(x)		x[["conf.lw.A"]][,"lower"]	)
	conf.lw.C.l	<- sapply(out,function(x)		x[["conf.lw.C"]][,"lower"]	)	
	###############################################################################
	#plot sampled a2a
	###############################################################################
	if(1)
	{
		require(fields)
		f.name<- paste(dir.name,paste("VARYCONTAM_TIPC_a2a_visitHCC.C",p.prev.instudy.clu.armC,"power",opt.power,"pool",opt.pooled,"sample",opt.sampling,"pwcalc",test.prop0,test.prop1,test.alpha,".pdf",sep='_'),sep='/')
		cat(paste("\nplot a2a to\n",f.name))
		pdf(paste(f.name),version="1.4",width=12,height=6)		
		breaks		<- diff(range(c(i2i.hg.A,i2i.hg.C)))/50
		breaks		<- seq(min(c(i2i.hg.A,i2i.hg.C)),by=breaks, len=51)	
		def.par <- par(no.readonly = TRUE)
		layout.m<- matrix(data= c(1,2),ncol=2,nrow=1,byrow=1)
		layout(layout.m)		
		image(main="arm C",p.consent.clus,p.nocontams,i2i.hg.C, breaks=breaks, col=head( rev(gray(seq(0,.95,len=trunc(50*1.4)))), 50))
		image.plot(p.consent.clus,p.nocontams,i2i.hg.A, breaks=breaks, col=head( rev(gray(seq(0,.95,len=trunc(50*1.4)))), 50),main="arm A")	
		par(def.par)
		dev.off()
	}
	###############################################################################
	#plot panel of confidence intervals
	###############################################################################
	if(1)
	{
		cols<- c("deepskyblue","dodgerblue4")
		f.name<- paste(dir.name,paste("VARYCONTAM_TIPC_confint_C_visitHCC.C",p.prev.instudy.clu.armC,"power",opt.power,"pool",opt.pooled,"sample",opt.sampling,"pwcalc",test.prop0,test.prop1,test.alpha,".pdf",sep='_'),sep='/')
		cat(paste("\nplot confint C to\n",f.name))
		pdf(paste(f.name),version="1.4",width=6,height=12)
		phdes.plot.confint.panel(t(conf.lw.C.l),t(conf.lw.C.u),t(conf.hg.C.l),t(conf.hg.C.u),p.consent.clus,p.nocontams,"p.vhcc.prev.Cs","p.nocontams", cols=cols)
		dev.off()
		
		f.name<- paste(dir.name,paste("VARYCONTAM_TIPC_confint_A_visitHCC.C",p.prev.instudy.clu.armC,"power",opt.power,"pool",opt.pooled,"sample",opt.sampling,"pwcalc",test.prop0,test.prop1,test.alpha,".pdf",sep='_'),sep='/')
		cat(paste("\nplot confint A to\n",f.name))
		pdf(paste(f.name),version="1.4",width=6,height=12)
		phdes.plot.confint.panel(t(conf.lw.A.l),t(conf.lw.A.u),t(conf.hg.A.l),t(conf.hg.A.u),p.consent.clus,p.nocontams,"p.vhcc.prev.Cs","p.nocontams", cols=cols)
		dev.off()
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
	test.prop0<- 0.10
	test.prop1<- 0.4
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
	f.name<- paste(dir.name,paste("SAMPLINGBIAS_TIPC_A_visitHCC.C",p.prev.instudy.clu.armC,"power",opt.power,"pool",opt.pooled,"sample",opt.sampling,"pwcalc",test.prop0,test.prop1,test.alpha,".pdf",sep='_'),sep='/')
	cat(paste("\nplot sampling bias A to\n",f.name))
	pdf(paste(f.name),version="1.4",width=6,height=6)
	ylim<- range(c(test.biased.H0.A,test.biased.H1.A,test.prop0,test.prop1))
	plot(1,1,type='n',xlim=range(s.consent),ylim=ylim, xlab="%consenting to phyl study at HCC",ylab="proportion acute to acute" )
	abline(h=test.prop0, col="red")
	abline(h=test.prop1, col="blue")
	lines(s.consent, test.biased.H0.A, col="red", lty=2)
	lines(s.consent, test.biased.H1.A, col="blue", lty=2)
	legend(x=0.38,y=0.4,bty='n',legend=c("a2a 40%","a2a 10%"), fill=c("blue","red"), border=NA)
	dev.off()
	#plot proportion I->I sampled against sampling fraction for arm C
	f.name<- paste(dir.name,paste("SAMPLINGBIAS_TIPC_C_visitHCC.C",p.prev.instudy.clu.armC,"power",opt.power,"pool",opt.pooled,"sample",opt.sampling,"pwcalc",test.prop0,test.prop1,test.alpha,".pdf",sep='_'),sep='/')
	cat(paste("\nplot sampling bias C to\n",f.name))
	pdf(paste(f.name),version="1.4",width=6,height=6)
	ylim<- range(c(test.biased.H0.C,test.biased.H1.C,test.prop0,test.prop1))
	plot(1,1,type='n',xlim=range(s.consent),ylim=ylim, xlab="%consenting to phyl study at HCC",ylab="proportion acute to acute" )
	abline(h=test.prop0, col="red")
	abline(h=test.prop1, col="blue")
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
	ylim<- range(c(test.biased.H0,test.biased.H1,test.prop0,test.prop1))
	plot(1,1,type='n',xlim=range(s.consent),ylim=ylim, xlab="%consenting to phyl study at HCC",ylab="proportion acute to acute" )
	abline(h=test.prop0, col="red")
	abline(h=test.prop1, col="blue")
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

