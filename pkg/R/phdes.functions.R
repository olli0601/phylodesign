


###############################################################################
#' randomize allocation of arms to triplets 
#' @param randomize.n 	number of times the randomization is repeated
#' @param rtn.fixed		indicator if a fixed allocation of arms is returned (for debugging)
#' @return data frame of trial sites with allocated arms
#' @examples	popart.getdata.randomized.arm(1)
popart.getdata.randomized.arm<- function( randomize.n, r.inc.A= 0.0056, r.inc.B=0.01, r.inc.C=0.013, p.adults.ZA=0.53, p.adults.SA=0.7, rtn.fixed=0, rtn.phylostudy=1 )
{	
	data(popart.triplets.130207)					
	# select studies for phylogenetics	
	if(rtn.phylostudy)
		simul	<- popart.triplets.130207[rep(which(popart.triplets.130207$triplet %in% c(1,4,5,6)),randomize.n),]
	else
		simul	<- popart.triplets.130207[rep(seq_len(nrow(popart.triplets.130207)),randomize.n),]
	simul$comid_old	<- gsub(' ','',simul$comid_old)
	# randomise communities arms
	#if(!rtn.fixed)
	#	simul$random	<- runif(nrow(simul))
	#else
	#	simul$random	<- seq_len(nrow(simul)) / nrow(simul)
	#simul$arm.code	<- unlist(tapply(simul$random,rep(1:(nrow(simul)/3),each=3),rank))
	#simul$arm		<- ifelse(simul$arm.code==1,"A",simul$arm.code)
	#simul$arm		<- ifelse(simul$arm.code==2,"B",simul$arm)
	#simul$arm		<- ifelse(simul$arm.code==3,"C",simul$arm)
	#drops 			<- c("random","arm.code")
	#simul			<- simul[,!names(simul) %in% drops]
	simul$inc.rate	<- sapply(simul$arm,function(x){(x %in% "A")*r.inc.A+(x %in% "B")*r.inc.B+(x %in% "C")*r.inc.C})
	simul$p.adults	<- ifelse(simul$country==1,p.adults.ZA,p.adults.SA)
	#simul$rid		<- rep(1:randomize.n, each=3*4)
	simul[ order(simul$triplet.id, simul$arm), ]	
}
###############################################################################
#' randomize the number of individuals with HIV, on ART at each trial site and in the population cohort 
#' @param x	data frame of trial sites with allocated arms
#' @param pc.size size of the population cohort
#' @param d.study duration of the trial in years
#' @param rtn.exp indicator if expectations are returned instead of draws from binomial (for debugging)
#' @return data frame of trial sites with allocated arms
#' @examples	sites<- popart.getdata.randomized.arm(1)
#' 				popart.getdata.randomized.n(sites, 2500, 3)
popart.getdata.randomized.n<- function(x, pc.size, d.study, rtn.exp=0)
{
	n			<- nrow(x)
	if(!rtn.exp)
	{
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
	}
	else
	{
		x$n.adults	<- x$popsize*x$p.adults
		x$n.prev	<- x$n.adults*x$hivcomb/100
		x$n.on.art	<- x$n.prev*x$artadjust/100
		x$n.not.art	<- x$n.prev-x$n.on.art	
		x$n.inc		<- (x$n.adults-x$n.prev)*d.study*x$inc.rate
		x$n.neg		<- with(x,n.adults-n.prev-n.inc)	
		x$p.in.PC	<- pc.size/x$n.adults	
		x$PC.on.art	<- x$n.on.art*x$p.in.PC
		x$PC.not.art<- x$n.not.art*x$p.in.PC
		x$PC.inc	<- x$n.inc*x$p.in.PC
		x$PC.neg	<- pc.size-x$PC.on.art-x$PC.not.art-x$PC.inc
	}
	x
}
###############################################################################
#' plot binomial confidence intervals
phdes.plot.confint<- function(p0,p1,c0,c1,xlab='s.consent',ylab='conf.int',cols=c("red","blue"),legend.loc="topright",legend.txt=c('',''))
{
	#print(p0); print(p1); print(c0); print(c1)
	par(mar=c(5,5.5,0.5,2))
	x<- as.numeric(rownames(c0))
	plot(1,1,type='n',xlim=range(x),ylim=range(rbind(c0,c1)),xlab=xlab,ylab=ylab)
	
	lines(x,p0,col=cols[1], lty=4)		
	lines(x,p1,col=cols[2], lty=4)		
	polygon( c(x,rev(x)),c(c0[,1],rev(c0[,2])), border=NA, col=my.fade.col(cols[1],0.3) )
	polygon( c(x,rev(x)),c(c1[,1],rev(c1[,2])), border=NA, col=my.fade.col(cols[2],0.3) )
	legend(legend.loc,fill=c(cols,"transparent","transparent"),legend=c(legend.txt,"","mean and 95% CI"),bty='n',border=NA)		
}
###############################################################################
#' plot binomial confidence intervals
phdes.plot.confint.panel<- function(lw.l,lw.u,hg.l,hg.u,x,y,xlab,ylab,cols=c("red","blue"))
{
	xlim<- range(x)
	ylim<- range(c(hg.l,hg.u,lw.l,lw.u))
	
	def.par <- par(no.readonly = TRUE)
	layout.m<- matrix(data= seq_len(length(y)),ncol=1,nrow=length(y))
	layout(layout.m)
	sapply(seq_len(length(y)),function(i)
			{												
				z<- cbind(lw.l[i,],lw.u[i,])
				z2<- cbind(hg.l[i,],hg.u[i,])
				par(mar=c(2,2,.5,.5))
				plot(1,1,type='n',xlim=xlim,ylim=ylim, xlab='', ylab="")
				polygon( c(x,rev(x)),c(z2[,1],rev(z2[,2])), border=NA, col=my.fade.col(cols[1],0.3) )
				polygon( c(x,rev(x)),c(z[,1],rev(z[,2])), border=NA, col=my.fade.col(cols[2],0.3) )
				legend("topleft",legend=paste(ylab,"=",round(y[i],2)),bty='n')				
			})
	par(def.par)
	
}
###############################################################################
#' plot power under a binomial test
phdes.plot.power<- function(m, m2, c, c2, xlab="", ylab="", legend.txt="", legend.loc="bottomright", cols=c("red","blue"), verbose= 0)
{	
	x<- as.numeric(names(m))
	x2<- as.numeric(names(m2))
	legend.txt<- c(legend.txt,"lines w/o dots cannot be trusted")	
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
}
###############################################################################
phdes.get.hyp.tipc.probs<- function(tipc.p,f.contam)
{
	tipc.p<- tipc.p / sum(tipc.p)								#tipc.p may not be standardized
	if(rownames(tipc.p)[nrow(tipc.p)]!='O')	stop("phdes.get.hyp.tipc.probs: error at 1a")
	
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
#' compute the number of acute to acute transmissions from tip clusters up to order 3
#' this is replaced by "popart.get.sampled.acute2acute"
popart.get.sampled.acute2acute.p<- function(tipc.p,pop.size,f.coh,f.susc,f.untreated,cohort.dur,perc.inc,p.lab,p.nocontam,p.consent.coh,p.consent.clu,p.prev.instudy.clu,p.inc.instudy.clu, opt.power="All")
{	
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
	s.prev.pooled<-	(p.consent.coh*f.coh	+  p.consent.clu*(1-f.coh)*p.prev.instudy.clu)	*	p.lab*f.untreated/(1-f.susc)	
	s.inc.nonPC<-	p.lab*p.consent.clu*p.inc.instudy.clu
	s.inc.PC<-	p.lab*p.consent.coh		
	if(length(s.inc.nonPC)!=1)	stop("popart.getsampled.tipcseq: cannot vectorize s.inc.nonPC")
	if(length(s.inc.PC)!=1)	stop("popart.getsampled.tipcseq: cannot vectorize s.inc.PC")
	
	s.clu<- t(matrix(	c(	s.inc.nonPC,		2*s.inc.nonPC*(1-s.inc.nonPC),		3*s.inc.nonPC*(1-s.inc.nonPC)^2,	
							0,			s.inc.nonPC*s.inc.nonPC,				3*s.inc.nonPC*s.inc.nonPC*(1-s.inc.nonPC),
							0,			0,							s.inc.nonPC*s.inc.nonPC*s.inc.nonPC 	),
						nrow=ncol(tipc.p),ncol=ncol(tipc.p),byrow=1))
	s.coh<- t(matrix(	c(	s.inc.PC,		2*s.inc.PC*(1-s.inc.PC),		3*s.inc.PC*(1-s.inc.PC)^2,	
							0,			s.inc.PC*s.inc.PC,				3*s.inc.PC*s.inc.PC*(1-s.inc.PC),
							0,			0,							s.inc.PC*s.inc.PC*s.inc.PC 	),
						nrow=ncol(tipc.p),ncol=ncol(tipc.p),byrow=1))
	if(verbose){ print(s.prev.pooled*s.clu); print(s.prev.pooled*s.coh) }
	
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
				
				tipc<- tipc.e[i]*s.prev.pooled[i]*( f.coh[i]*tipc.p%*%s.coh	+ (1-f.coh[i])*tipc.p%*%s.clu)	#tip clusters sampled in PC + tip clusters sampled in HCC
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
#' compute the number of incident sequences that can be linked to a prevalent case with an imperfect phylogenetic method
#' @param x2i 		number of incident sequences that could be linked to a prevalent case
#' @param p.linkage proportion of incident sequences that are linked with a phylogenetic method
#' @param rtn.exp indicator if expectations are returned instead of draws from binomial (for debugging)
#' @return linked.x2i	number of incident sequences that are linked to a prevalent case under a Binomial link model
#' @examples 	sites<- popart.getdata.randomized.arm(1)
#' 	sites<- popart.getdata.randomized.n(sites, 2500, 3)
#' 	x2i<- popart.get.sampled.transmissions(sites)
#'	phdes.get.linked.transmissions(x2i,0.8)			
phdes.get.linked.transmissions<- function(x2i, p.linkage, rtn.exp=0)
{
	if(rtn.exp)
		linked.x2i<- x2i*p.linkage
	else
	{
		linked.x2i<- rbinom(length(x2i),as.vector(x2i),p.linkage)
		if(is.matrix(x2i))
			linked.x2i<- matrix(linked.x2i, nrow=nrow(x2i), ncol=ncol(x2i), dimnames= list(rownames(x2i),colnames(x2i)))
	}
	linked.x2i
}
###############################################################################
#' compute the number of incident sequences that can be linked to a prevalent case with a perfect phylogenetic method
#' @param x					data frame of trial sites with randomized arms and population numbers
#' @param method			sampling scheme. Options are "PC after yr 1 and HCC", "PC only incident and HCC", "PC and HCC", "only HCC"
#' @param rtn.int			indicator if integer values are returned
#' @param p.vhcc.prev.AB 	proportion of initially +ve individuals who visit HCC at some point in 3 year period in Arms A & B
#' @param p.vhcc.prev.C  	proportion of initially +ve individuals who visit HCC at some point in 3 year period in Arm C
#' @param p.vhcc.inc.AB		proportion of incident cases who visit HCC at some point in 3 year period in Arms A and B
#' @param p.vhcc.inc.C		proportion of incident cases who visit HCC at some point in 3 year period in Arm C
#' @param consent.PC		consent rate to phylogenetics in population cohort
#' @param consent.HCC		consent rate to phylogenetics in health care centres
#' @param p.lab				proportion of isolates that make it through to full genomes
#' @param p.community		proportion of infections that come from within the community
#' @return x2i	the number of sampled transmissions
#' @examples 	sites<- popart.getdata.randomized.arm(1)
#' 				sites<- popart.getdata.randomized.n(sites, 2500, 3)
#' 				popart.get.sampled.transmissions(sites)
popart.get.sampled.transmissions<- function(x, method="PC and HCC", rtn.int=1, p.vhcc.prev.AB=0.95, p.vhcc.prev.C=0.5, p.vhcc.inc.AB=0.8, p.vhcc.inc.C=0.25, consent.PC=0.9, consent.HCC=0.5, p.lab=0.8, p.community=0.85)
{
	visit.prev	<- numeric(nrow(x))
	visit.inc	<- numeric(nrow(x))
	visit.prev[ x$arm%in% c("A","B") ]	<- p.vhcc.prev.AB
	visit.prev[ x$arm=="C" ]			<- p.vhcc.prev.C
	visit.inc[ x$arm%in% c("A","B") ]	<- p.vhcc.inc.AB
	visit.inc[ x$arm=="C" ]				<- p.vhcc.inc.C
			
	if(method=="PC and HCC")
	{
		seq.PC.prev		<- x$PC.not.art*consent.PC*p.lab
		seq.PC.inc		<- x$PC.inc*consent.PC*p.lab
	}
	else if(method=="PC after yr 1 and HCC")
	{
		seq.PC.prev		<- x$PC.not.art*visit.prev*consent.HCC*p.lab	+	x$PC.inc/3*consent.PC*p.lab
		seq.PC.inc		<- x$PC.inc*2/3*consent.PC*p.lab
	}
	else if(method=="PC only incident and HCC")
	{
		seq.PC.prev		<- x$PC.not.art*visit.prev*consent.HCC*p.lab
		seq.PC.inc		<- x$PC.inc*consent.PC*p.lab
	}
	else if(method=="only HCC")
	{
		seq.PC.prev		<- x$PC.not.art*visit.prev*consent.HCC*p.lab
		seq.PC.inc		<- x$PC.inc*visit.inc*consent.HCC*p.lab
	}
	else
		stop("popart.get.sampled.transmissions: unknown method to sample")
	
	seq.nonPC.prev	<- (x$n.not.art-x$PC.not.art)*visit.prev*consent.HCC*p.lab		
	seq.nonPC.inc	<- (x$n.inc-x$PC.inc)*visit.inc*consent.HCC*p.lab
	
	seq.prev.coverage	<- (seq.PC.prev + seq.nonPC.prev)/x$n.prev
	x2i					<- (seq.PC.inc + seq.nonPC.inc)*seq.prev.coverage*p.community
	if(rtn.int)	x2i		<- floor(x2i)
	x2i
}
###############################################################################
popart.predicted.firstCD4.131017<- function(sites, opt.design="PC12+HCC", verbose=0)
{
	#
	# 	CD4 by country, arm	from HCF	
	#	assume all FirstCD4 counts exclude PC counts
	#
	samples.CD4<- data.table(	FirstCD4.all= c(4623, 4078, 1032, 7358, 6353, 1159),
								FirstCD4.inc= c(710, 929, 106, 705, 926, 57),
								country		= c("ZA","ZA","ZA","SA","SA","SA"),
								arm			= c("A","B","C","A","B","C"),
								prediction	= "central",
								p.contam.act= 0.05	)		
	#
	# 	patients by country, arm	from PC
	#	assume all PC patients provide samples to HCF
	#
	samples.PC<- data.table(	PC0.u			= c(rep(502,3),rep(574,3)),
								PC0.t			= c(rep(166,3),rep(162,3)),
								PC0.l			= c(rep(74,3),rep(82,3)),
								PC12.u			= c(164,350,448,182,370,502),
								PC12.t			= c(408,236,154,3454,274,150),
								PC12.l			= c(134,138,140,148,152,154),
								SC.y1			= c(32,46,64,22,32,46),
								SC.y2			= c(22,40,58,16,28,42),
								SC.y3			= c(20,36,54,14,26,38),									
								country			= c("ZA","ZA","ZA","SA","SA","SA"),
								arm				= c("A","B","C","A","B","C"),
								prediction		= "central",
								p.contam.act	= 0.05		)									
	#
	#	total incidence and prevalence by community
	#
	file		<- paste(CODE.HOME,"data","popart.SummaryResultsPerCommunity.central.131016.Rdata",sep='/')
	tmp			<- load(file)
	if(verbose) cat(paste("loaded",paste(tmp,collapse=' ')))
	result		<- as.data.table(result)
	setnames(result, c("comm.name","allocatedArms"), c("comid_old","arm"))
	pred.model	<- result[, list( arm=arm, all.prev= round(adultpopsize*prevalence36m), all.inc= round(adultpopsize*sum(cumulincyear1, cumulincyear2, cumulincyear3)) ), by="comid_old" ]
	pred.model[,prediction:= "central"]
	set(pred.model, which(pred.model[,comid_old]=="Town II"), "comid_old", "TownII")
	set(pred.model, which(pred.model[,comid_old]=="Delft South"), "comid_old", "DelftSouth")
	
	#	recalculate CD4 by community from HCF	
	tmp			<- sites[,list( triplet.id=triplet.id, arm=arm, comid_old=comid_old, country=country,  adultsize=popsize*p.adults )]		
	tmp			<- merge( tmp[,	list(adultsize.arm.country=sum(adultsize)), by= c("country","arm")], tmp, by=c("country","arm"))
	samples.CD4	<- merge(samples.CD4, tmp, by=c("country","arm"))
	samples.CD4	<- samples.CD4[, list(country=country, prediction=prediction, FirstCD4.all=round(FirstCD4.all*adultsize/adultsize.arm.country), FirstCD4.inc=round(FirstCD4.inc*adultsize/adultsize.arm.country)), by="comid_old"]
	samples.CD4[, FirstCD4.prev:= FirstCD4.all-FirstCD4.inc]		
	setnames(samples.CD4, c("FirstCD4.inc","FirstCD4.prev"), c("HCF.inc", "HCF.prev"))
	samples.CD4	<- subset( samples.CD4, select=c(comid_old, country, prediction, HCF.inc, HCF.prev) )
	#	recalculate CD4 by community from PC
	samples.PC	<- merge(samples.PC, tmp, by=c("country","arm"))
	samples.PC[, scale:= adultsize/adultsize.arm.country]
	samples.PC	<- samples.PC[, list( prediction=prediction, PC0.prev= round(PC0.u*scale), PC12.prev= round((PC12.u+SC.y1)*scale), PC0.inc= round((SC.y1+SC.y2+SC.y3)*scale), PC12.inc=round((SC.y2+SC.y3)*scale) ), by="comid_old"]
	#	add design column
	samples.PC	<- rbindlist( list( 	samples.PC[, list(comid_old=comid_old, prediction=prediction, design="PC+HCC", PC.prev=PC0.prev,  PC.inc=PC0.inc)],
					samples.PC[, list(comid_old=comid_old, prediction=prediction, design="PC12+HCC", PC.prev=PC12.prev, PC.inc=PC12.inc)]	)		)
	#	add PC and HCF samples
	samples.CD4	<- merge(samples.PC, samples.CD4, by=c("comid_old","prediction"))
	samples.CD4	<- merge(samples.CD4, pred.model, by=c("comid_old","prediction"))
	#
	samples.CD4	<- subset(samples.CD4, design==opt.design)
	samples.CD4
}	
###############################################################################
popart.predicted.firstCD4<- function(version= "130627")
{
	if(version=="130627")	
	{
		#Anne: Results obtained for 6 months CHiPS rounds, assuming that individuals come for preART CD4 monitoring once a year 
		#	and  ART CD4 monitoring very 3 months. The definition of incident cases is individuals infected during the trial.
		
		#	this is assuming central targets for linkage and consenting
		first.CD4.HCC		<- matrix( c(	2397, 289, 3078, 604, 705, 57, 
											3960, 477, 2270, 445, 771, 62, 
											2231, 168, 3482, 431, 793, 32, 
											7567, 566, 4578, 568, 773, 31), nrow=12, ncol=2, byrow=1, 
										dimnames=list(	c("Ndeke","Chimwemwe","Ngungu","Maramba","Dambwa","Shampande","Luvuyo","Kuyasa","TownII","Ikhwezi","DelftSouth","Bloekombos"),
														c("HCC Total","HCC Incident")))
		first.CD4.PC		<- matrix( c(	307, 42, 320, 74, 333, 107, 
											307, 42, 320, 74, 333, 107, 
											343, 30, 352, 53, 359, 77, 
											343, 30, 352, 53, 359, 77), nrow=12, ncol=2, byrow=1, 
										dimnames=list(	c("Ndeke","Chimwemwe","Ngungu","Maramba","Dambwa","Shampande","Luvuyo","Kuyasa","TownII","Ikhwezi","DelftSouth","Bloekombos"),
														c("PC Total","PC Incident")))
		#	this is assuming no loss due to linkage and consenting
		first.CD4.All		<- matrix( c(	2410, 332, 3237, 750, 2985, 957, 
											3981, 549, 2386, 553, 3262, 1046, 
											2281, 199, 3735, 562, 4253, 917, 
											7736, 677, 4910, 738, 4144, 893), nrow=12, ncol=2, byrow=1, 
										dimnames=list(	c("Ndeke","Chimwemwe","Ngungu","Maramba","Dambwa","Shampande","Luvuyo","Kuyasa","TownII","Ikhwezi","DelftSouth","Bloekombos"),
														c("total prev","total inc")))
	}
	ans<- cbind(first.CD4.PC,first.CD4.HCC,first.CD4.All)
	#print(ans); stop()
	ans
}
###############################################################################
popart.predicted.sequences<- function(sites, first.CD4, p.consent.PC, p.consent.HCC, p.lab, p.vhcc.prev.AB, p.vhcc.inc.AB, p.vhcc.prev.C, p.vhcc.inc.C, method="PC and HCC")
{
	x					<- sites									
	if(any(x$comid_old != rownames(first.CD4))) stop("first.CD4 not ordered by sites")									
			 
	x$CD4.1st.PC.prev	<- first.CD4[,"PC Total"] - first.CD4[,"PC Incident"] 
	x$CD4.1st.PC.inc	<- first.CD4[,"PC Incident"] 
	x$CD4.1st.HCC.prev	<- first.CD4[,"HCC Total"] - first.CD4[,"HCC Incident"] 
	x$CD4.1st.HCC.inc	<- first.CD4[,"HCC Incident"]
	s<- matrix(NA,nrow=nrow(x),ncol=9,dimnames= list(x$comid_old, c("visit.prev","visit.inc","PC.prev","nonPC.prev","PC.inc","nonPC.inc","%prev","%inc","%avg")))	
	s[x$arm!="C", "visit.prev"]	<- p.vhcc.prev.AB
	s[x$arm=="C", "visit.prev"]	<- p.vhcc.prev.C
	s[x$arm!="C", "visit.inc"]	<- p.vhcc.inc.AB
	s[x$arm=="C", "visit.inc"]	<- p.vhcc.inc.C
	
	s[, "PC.prev"]		<- p.consent.PC*p.lab
	s[, "PC.inc"]		<- p.consent.PC*p.lab
	s[, "nonPC.prev"]	<- s[,"visit.prev"]*p.consent.HCC*p.lab		
	s[, "nonPC.inc"]	<- s[,"visit.inc"]*p.consent.HCC*p.lab		
	
	s[, "nonPC.inc"]	<- round(switch(	method,
									"PC and HCC"				= x$CD4.1st.HCC.inc * s[, "nonPC.inc"],
									"PC12+HCC"					= (x$CD4.1st.HCC.inc + x$CD4.1st.PC.inc/3) * s[, "nonPC.inc"],
									"PC only incident and HCC"	= x$CD4.1st.HCC.inc * s[, "nonPC.inc"],
									"only HCC"					= (x$CD4.1st.HCC.inc + x$CD4.1st.PC.inc/3) * s[, "nonPC.inc"],
									NA
									))						
	s[, "nonPC.prev"]	<- round(switch(	method,
									"PC and HCC"				= x$CD4.1st.HCC.prev * s[, "nonPC.prev"],
									"PC12+HCC"					= (x$CD4.1st.HCC.prev+x$CD4.1st.PC.prev/3) * s[, "nonPC.prev"],
									"PC only incident and HCC"	= (x$CD4.1st.HCC.prev+x$CD4.1st.PC.prev) * s[, "nonPC.prev"],
									"only HCC"					= (x$CD4.1st.HCC.prev+x$CD4.1st.PC.prev) * s[, "nonPC.prev"],
									NA
									))									
	s[, "PC.prev"]		<- round(switch(	method,
									"PC and HCC"				= x$CD4.1st.PC.prev*s[, "PC.prev"],
									"PC12+HCC"					= x$CD4.1st.PC.prev*s[, "PC.prev"],
									"PC only incident and HCC"	= 0,
									"only HCC"					= 0,
									NA
									))										
	s[, "PC.inc"]		<- round(switch(	method,
									"PC and HCC"				= x$CD4.1st.PC.inc*s[, "PC.inc"],
									"PC12+HCC"					= x$CD4.1st.PC.inc*2/3*s[, "PC.inc"],
									"PC only incident and HCC"	= x$CD4.1st.PC.inc*s[, "PC.inc"],
									"only HCC"					= 0,
									NA
									))
	s[, "%prev"]		<- apply(s[, c("PC.prev","nonPC.prev")],1,sum)	/ 	first.CD4[,"total prev"]
	s[, "%inc"]			<- apply(s[, c("PC.inc","nonPC.inc")],1,sum)	/ 	first.CD4[,"total inc"]
	s[, "%avg"]			<- (s[, "%prev"] + s[, "%inc"])/2				
	as.data.frame(s[,-c(1,2)])	
}
###############################################################################
popart.set.hypo<- function(sites, theta.EE.H0, theta.EE.H1, opt.analysis, df.hyp=NULL)
{
	if(opt.analysis=="1040")
	{
		sites[,"mu.inc.rate.H0"]	<- sites[,inc.rate]
		sites[,"mu.inc.rate.H1"]	<- 0.013
		sites[,"mu.pE2E.H0"]		<- theta.EE.H0
		sites[,"mu.pE2E.H1"]		<- theta.EE.H1
	}
	else
	{
		if(is.null(df.hyp))								stop("df.hyp is NULL")
		opt.analysis<- ifelse( grepl("central",opt.analysis), "central", ifelse( grepl("optimistic",opt.analysis), "optimistic", "pessimistic" ) )
		if(!opt.analysis%in%unique(df.hyp[,target]))	stop("Unknown target in df.hyp")
		tmp			<- subset(df.hyp,target==opt.analysis & h=="H0", c(country, arm, Inc, E2E))	
		setnames(tmp,	c("Inc","E2E"),	c("mu.inc.rate.H0","mu.pE2E.H0"))
		sites		<- merge(sites, tmp, by=c("country","arm"))
		tmp			<- subset(df.hyp,target==opt.analysis & h=="H1", c(country, arm, Inc, E2E))	
		setnames(tmp,	c("Inc","E2E"),	c("mu.inc.rate.H1","mu.pE2E.H1"))
		sites		<- merge(sites, tmp, by=c("country","arm"))				
	}	
	sites		<- sites[order(sites[,triplet.id],sites[,arm]),]
	sites
}
###############################################################################
popart.predicted.sequences.130717<- function(samples.CD4, df.nocontam, opt.analysis, p.lab)
{
	opt.prediction	<- ifelse( grepl("central",opt.analysis), "central", ifelse( grepl("optimistic",opt.analysis), "optimistic", "pessimistic" ) )
	tmp				<- subset(df.nocontam, target==opt.prediction, c(country, arm, p.nocontam))
	samples.seq		<- subset(samples.CD4, prediction==opt.prediction)
	samples.seq		<- merge(samples.seq, tmp, by=c("country","arm"))
	samples.seq		<- samples.seq[, list(country=country, arm=arm, PC.prev=round(PC.prev*p.lab*p.nocontam), PC.inc=round(PC.inc*p.lab*p.nocontam), HCF.inc=round(HCF.inc*p.lab*p.nocontam), HCF.prev=round(HCF.prev*p.lab*p.nocontam), all.prev=all.prev, all.inc=all.inc), by="comid_old"]
	samples.seq[,"%prev":=round( samples.seq[, (PC.prev+HCF.prev) / all.prev], d=3)]
	samples.seq[,"%inc":=round( samples.seq[, (PC.inc+HCF.inc) / all.inc], d=3)]
	samples.seq[,"%avg":=round( samples.seq[, (PC.inc+HCF.inc+PC.prev+HCF.prev) / (all.prev+all.inc)], d=3)]
	setnames(samples.seq, c("HCF.inc","HCF.prev"), c("nonPC.inc","nonPC.prev") )		#to match previous data table
	samples.seq
}
###############################################################################
popart.sampling.init<- function(x, p.consent.PC, p.consent.HCC, p.lab, p.vhcc.prev.AB, p.vhcc.inc.AB, p.vhcc.prev.C, p.vhcc.inc.C, method="PC and HCC")
{
	s<- matrix(NA,nrow=nrow(x),ncol=7,dimnames= list(x$comid, c("visit.prev","visit.inc","PC.prev","nonPC.prev","PC.inc","nonPC.inc","baseline")))
	
	s[x$arm!="C", "visit.prev"]	<- p.vhcc.prev.AB
	s[x$arm=="C", "visit.prev"]	<- p.vhcc.prev.C
	s[x$arm!="C", "visit.inc"]	<- p.vhcc.inc.AB
	s[x$arm=="C", "visit.inc"]	<- p.vhcc.inc.C
	
	s[, "PC.prev"]		<- p.consent.PC*p.lab
	s[, "PC.inc"]		<- switch(	method,
									"PC and HCC"				= rep( p.lab*p.consent.PC, nrow(s) ),
									"PC only incident and HCC"	= rep( p.lab*p.consent.PC, nrow(s) ),
									"PC after yr 1 and HCC"		= rep( p.lab*p.consent.PC*2/3, nrow(s) ),
									"only HCC"					= p.lab*p.consent.HCC*s[,"visit.inc"],
									NA
									)
	if(any(is.na(s[, "PC.inc"])))	stop("unknown method")		
	
	s[, "nonPC.prev"]	<- s[,"visit.prev"]*p.consent.HCC*p.lab		
	s[, "nonPC.inc"]	<- s[,"visit.inc"]*p.consent.HCC*p.lab			
	
	s[, "baseline"]		<- switch(	method,
									"PC and HCC"				= x$PC.not.art*s[,"PC.prev"]	+	(x$n.not.art-x$PC.not.art)*s[,"nonPC.prev"],
									"PC after yr 1 and HCC"		= x$PC.inc/3*s[,"PC.prev"]		+ 	x$n.not.art*s[,"nonPC.prev"],
									"PC only incident and HCC"	= x$n.not.art*s[,"nonPC.prev"],
									"only HCC"					= x$n.not.art*s[,"nonPC.prev"],
									NA
									)
	if(any(is.na(s[, "baseline"])))	stop("unknown method")								
	s[, "baseline"]		<- s[, "baseline"] / x$n.prev		#coverage at baseline is sampled HIV+ at baseline / total HIV+ at baseline	
	as.data.frame(s[,-c(1,2)])	
}
###############################################################################
popart.get.sampled.transmissions.from.tipc<- function(x, tipc.p, clu.n, theta, sampling, method="PC and HCC", rtn.int=1, mx.sampled.ntr=ncol(clu.n), exclude.O= 1, verbose=0)
{	
	if(abs(sum(tipc.p, na.rm=1)-1)>EPS)	stop("invalid tipc.p")
	
	tipc.PC		<- lapply( x$PC.inc, 			function(z){ tmp<- clu.simulate(tipc.p, z, rtn.int=rtn.int); tmp[,1]<- 0; tmp })
	tipc.nonPC	<- lapply( x$n.inc-x$PC.inc, 	function(z){ tmp<-  clu.simulate(tipc.p, z, rtn.int=rtn.int); tmp[,1]<- 0; tmp }  )
	#if(verbose)	print( lapply(tipc.nonPC, function(x)	apply(x,1,sum)) )
	print(tipc.PC[[1]])
	print("HERE")
	tipc.PC.s	<- lapply(	seq_along(tipc.PC), function(i)
						{
							tmp			<- as.numeric(sampling[i,c("baseline","PC.inc")])
							names(tmp)	<- c("Idx","E")							
							clu.sample(tipc.PC[[i]], tmp, rtn.exp=!rtn.int)
						})
	print(tipc.PC.s[[1]])
	stop()			
	tipc.nonPC.s<- lapply(	seq_along(tipc.nonPC), function(i)
						{
							tmp			<- as.numeric(sampling[i,c("baseline","nonPC.inc")])
							names(tmp)	<- c("Idx","E")
							clu.sample(tipc.nonPC[[i]], tmp, rtn.exp=!rtn.int)
						})					
	if(verbose)
	{
		NCLU<<- rbind(NCLU, matrix(range(sapply(tipc.nonPC.s, function(x)	sum(apply(x,1,sum)))),1,2) )
		cat(paste("\nrange of cluster numbers found: ", paste(range( sapply(tipc.nonPC.s, function(x)	sum(apply(x,1,sum))) ), collapse=', ',sep='') ))
	}
	ans<- sapply(seq_along(tipc.PC.s),function(i)
						{
							tmp			<- c(1,1)
							names(tmp)	<- c("Idx","E")		
														
							tipc.pooled	<- tipc.PC[[i]] + tipc.nonPC[[i]]
							#if(verbose) print(tipc.pooled)
							tr.complete	<- clu.exp.transmissions(tipc.pooled, 		clu.n, theta, tmp, mx.s.ntr=mx.sampled.ntr, exclude.O=exclude.O )
							
							tmp[]		<- as.numeric(sampling[i,c("baseline","PC.inc")])
							#print(tmp)
							tr.PC.s		<- clu.exp.transmissions(tipc.PC.s[[i]], 	clu.n, theta, tmp, mx.s.ntr=mx.sampled.ntr, exclude.O=exclude.O )
							
							tmp[]		<- as.numeric(sampling[i,c("baseline","nonPC.inc")])
							#print(tmp)
							tr.nonPC.s	<- clu.exp.transmissions(tipc.nonPC.s[[i]], clu.n, theta, tmp, mx.s.ntr=mx.sampled.ntr, exclude.O=exclude.O )
							
							c( tr.PC.s["E2E"]+tr.nonPC.s["E2E"], tr.complete["E2E"], sum(tr.PC.s)+sum(tr.nonPC.s), sum(tr.complete) )							
						})
	colnames(ans)<- x$comid			
	rownames(ans)<- c("i2i.s","i2i.c","x2i.s","x2i.c")	
	if(rtn.int)
		ans<- floor(ans)		
	ans
}
###############################################################################
#' compute the number of acute to acute transmissions from tip clusters up to order 3
#' @param x					data frame of trial sites with randomized arms and population numbers
#' @param tipc.p			expected frequencies of tip cluster counts under hypothesis
#' @param method			sampling scheme. Options are "PC after yr 1 and HCC", "PC and HCC", "only HCC"
#' @param rtn.int			indicator if integer values are returned 
#' @param p.vhcc.prev.AB 	proportion of initially +ve individuals who visit HCC at some point in 3 year period in Arms A & B
#' @param p.vhcc.prev.C  	proportion of initially +ve individuals who visit HCC at some point in 3 year period in Arm C
#' @param p.vhcc.inc.AB		proportion of incident cases who visit HCC at some point in 3 year period in Arms A and B
#' @param p.vhcc.inc.C		proportion of incident cases who visit HCC at some point in 3 year period in Arm C
#' @param consent.PC		consent rate to phylogenetics in population cohort
#' @param consent.HCC		consent rate to phylogenetics in health care centres
#' @param p.lab				proportion of isolates that make it through to full genomes
#' @return For each trial site, a vector containing
#' \item{i2i.s}{number of sampled acute to acute transmissions}
#' \item{i2i.c}{number of acute to acute transmissions}
#' \item{x2i.s}{number of sampled transmissions}
#' \item{x2i.c}{number of true transmissions}
#' @examples 	sites<- popart.getdata.randomized.arm(1)
#' sites<- popart.getdata.randomized.n(sites, 2500, 3)
#' popart.get.sampled.transmissions(sites)
popart.get.sampled.acute2acute<- function(x, tipc.p, method="PC and HCC", rtn.int= 1, p.vhcc.prev.AB=0.95, p.vhcc.prev.C=0.5, p.vhcc.inc.AB=0.8, p.vhcc.inc.C=0.25, consent.PC=0.9, consent.HCC=0.5, p.lab=0.8)
{
	#print("")
	#tmp<- c(p.vhcc.prev.AB, p.vhcc.prev.C, p.vhcc.inc.AB, p.vhcc.inc.C, consent.PC, consent.HCC, p.lab, p.community); names(tmp)<- c("p.vhcc.prev.AB","p.vhcc.prev.C","p.vhcc.inc.AB","p.vhcc.inc.C","consent.PC","consent.HCC","p.lab","p.community"); print(tmp)
	verbose<- 0
	visit.prev							<- numeric(nrow(x))
	visit.inc							<- numeric(nrow(x))
	visit.prev[ x$arm%in% c("A","B") ]	<- p.vhcc.prev.AB
	visit.prev[ x$arm=="C" ]			<- p.vhcc.prev.C
	visit.inc[ x$arm%in% c("A","B") ]	<- p.vhcc.inc.AB
	visit.inc[ x$arm=="C" ]				<- p.vhcc.inc.C
				
	tipc.p	<- tipc.p / sum(tipc.p)								#tipc.p may not be standardized
	inc.p	<- apply(tipc.p,2,sum)*seq_len(ncol(tipc.p))		#fraction of incidence in tipc column
	
	#total expected tip clusters with at least 1 infection
	#less than n.inc because some tip clusters have more than 1 incident case within them
	tipc.PC		<- x$PC.inc / sum(inc.p)
	tipc.nonPC	<- (x$n.inc-x$PC.inc) / sum(inc.p)
	if(verbose){	print(x$PC.inc); print(tipc.PC); print(x$n.inc-x$PC.inc); print(tipc.nonPC)	}
	if(length(tipc.PC)!=length(tipc.nonPC))	stop("popart.get.sampled.acute2acute: tipc.nonPC, tipc.PC length must be the same")
	
	#proportion of prev individuals 
	#whose sequences are isolated successfully in trial and
	#which lead to transmission in community
	if(method=="PC and HCC")
		seq.prev		<- (x$PC.not.art*consent.PC	+	(x$n.not.art-x$PC.not.art)*visit.prev*consent.HCC)	* p.lab
	else if(method=="PC after yr 1 and HCC")
		seq.prev		<- x$n.not.art*visit.prev*consent.HCC*p.lab		+		x$PC.inc/3*consent.PC*p.lab
	else if(method=="PC only incident and HCC")
		seq.prev		<- x$n.not.art*visit.prev*consent.HCC*p.lab
	else if(method=="only HCC")
		seq.prev		<- x$n.not.art*visit.prev*consent.HCC*p.lab
	else
		stop("unknown method to sample")	
	seq.prev.coverage	<- seq.prev / x$n.prev
	s.prev				<- seq.prev.coverage #there is no p.community here because contamination is dealt with excluding the tipc.p[3,] row below
	
	#incidence sampling matrices of tipc clusters up to order 3 for PC and nonPC
	#transpose so we re ready for matrix multiplication
	#since there are higher order terms of si.XX, need to do cohort and community separately				
	s.inc.nonPC	<- lapply(seq_along(visit.inc),function(i)
			{
				s.inc.nonPC	<- p.lab*consent.HCC*visit.inc[i]							
				t(matrix(	c(	s.inc.nonPC,		2*s.inc.nonPC*(1-s.inc.nonPC),		3*s.inc.nonPC*(1-s.inc.nonPC)^2,	
								0,			s.inc.nonPC*s.inc.nonPC,				3*s.inc.nonPC*s.inc.nonPC*(1-s.inc.nonPC),
								0,			0,							s.inc.nonPC*s.inc.nonPC*s.inc.nonPC 	),
								nrow=ncol(tipc.p),ncol=ncol(tipc.p),byrow=1))
			})
	s.inc.PC	<- lapply(seq_along(visit.inc),function(i)
			{
				s.inc.PC	<-	switch(	method,
										"PC and HCC"=				p.lab*consent.PC,
										"PC only incident and HCC"=	p.lab*consent.PC,
										"PC after yr 1 and HCC"=	p.lab*consent.PC*2/3,
										"only HCC"	=				p.lab*consent.HCC*visit.inc[i]
										)
				t(matrix(	c(	s.inc.PC,		2*s.inc.PC*(1-s.inc.PC),		3*s.inc.PC*(1-s.inc.PC)^2,	
								0,			s.inc.PC*s.inc.PC,				3*s.inc.PC*s.inc.PC*(1-s.inc.PC),
								0,			0,							s.inc.PC*s.inc.PC*s.inc.PC 	),
								nrow=ncol(tipc.p),ncol=ncol(tipc.p),byrow=1))
			})
	if(length(s.inc.PC)!=length(s.inc.nonPC))	stop("popart.get.sampled.acute2acute: s.inc.PC, s.inc.nonPC length must be the same")
	#print(s.inc.nonPC); print(s.inc.PC)
	#do not observe transitions from outside community
	tipc.p<- tipc.p[-3,]	
	#compute sampled acute2acute and any2acute transmissions from tip clusters
	ans<- sapply(seq_along(tipc.PC),function(i)
			{
				out			<- numeric(4)
				names(out)	<- c("i2i.s","i2i.c","x2i.s","x2i.c")
				#all tip cluster counts before sampling
				tipc.pooled	<- (tipc.PC[i]+tipc.nonPC[i])*tipc.p							
				if(verbose) print(tipc.pooled)
				tmp			<- c(apply(tipc.pooled,1,sum), sum(tipc.pooled[,2])+2*sum(tipc.pooled[,3]))		#all U->I, T->I, I->I before sampling
				if(verbose)	print(tmp)
				out[c("i2i.c","x2i.c")]<- c(tmp[3],sum(tmp))					#I->I and U->I + T->I + I->I	before sampling
				
				#all tip cluster counts after sampling
				tipc.pooled	<- tipc.PC[i]*s.prev[i]*tipc.p%*%s.inc.PC[[i]]	+ 	tipc.nonPC[i]*s.prev[i]*tipc.p%*%s.inc.nonPC[[i]]	
				if(verbose) print(tipc.pooled)	
				tmp			<- c(apply(tipc.pooled,1,sum), sum(tipc.pooled[,2])+2*sum(tipc.pooled[,3]))		#all U->I, T->I, I->I after sampling
				if(verbose) print(tmp)	
				out[c("i2i.s","x2i.s")]<- c(tmp[3],sum(tmp))					#I->I and U->I + T->I + I->I	after sampling
				if(rtn.int)	out<- floor(out)
				out				
			})
	if(verbose) print(ans)
#print(ans); stop()
	ans
}
###############################################################################
#' compute the fraction of linked incident sequences based on overall proportions
#' this is replaced by "popart.get.sampled.transmissions"
popart.get.sampled.transmissions.p<- function(pop.size,f.coh,f.susc,f.untreated,cohort.dur,perc.inc,p.lab,p.nocontam,p.consent.coh,p.consent.clu,p.prev.instudy.clu,p.inc.instudy.clu, opt.power="All")
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
#'	pool the number of transmissions according to a pooling option
#' 	@param sites	data frame of trial sites
#' 	@param inc		matrix of transmissions
#'  @param method 	poopling method. Supported methods are "no pooling", "pooled across country", "pooled across ZA", "pooled across SA", "pooled across trial".
#' 	@return \item{inc}{pooled matrix of acute to acute transmissions}
#' 			\item{idx.A}{row index of arm A sites}
#' 			\item{idx.B}{row index of arm B sites}
#'			\item{idx.C}{row index of arm C sites}
#' 	@examples 	
#' 	sites<- popart.getdata.randomized.arm(1)
#' 	sites<- popart.getdata.randomized.n(sites, 2500, 3)
#' 	sampled.x2i<- popart.get.sampled.transmissions(sites)
#' 	popart.pool(sites, sampled.x2i, "pooled across ZA")
popart.pool<- function(sites, transm, method="no pooling")
{
	if(method=="no pooling")
	{
		idx.A<-	which(sites$arm=="A")
		idx.B<-	which(sites$arm=="B")
		idx.C<-	which(sites$arm=="C")				
	}
	else if(method=="pooled across country")
	{
		transm<- do.call(rbind,tapply(seq_len(nrow(sites)),sites$rid, function(ris)
						{
							rbind(	apply(transm[ris[which(sites$arm[ris]=="A" & sites$country[ris]==1)],,drop=0],2,sum),
									apply(transm[ris[which(sites$arm[ris]=="B" & sites$country[ris]==1)],,drop=0],2,sum),
									apply(transm[ris[which(sites$arm[ris]=="C" & sites$country[ris]==1)],,drop=0],2,sum),
									apply(transm[ris[which(sites$arm[ris]=="A" & sites$country[ris]==2)],,drop=0],2,sum),
									apply(transm[ris[which(sites$arm[ris]=="B" & sites$country[ris]==2)],,drop=0],2,sum),
									apply(transm[ris[which(sites$arm[ris]=="C" & sites$country[ris]==2)],,drop=0],2,sum))					
						}))
		idx.A<-	seq(1,nrow(transm),by=3)
		idx.B<-	seq(2,nrow(transm),by=3)
		idx.C<-	seq(3,nrow(transm),by=3)
	}
	else if(method=="pooled across ZA")
	{
		transm<- do.call(rbind,tapply(seq_len(nrow(sites)),sites$rid, function(ris)
						{
							rbind(	apply(transm[ris[which(sites$arm[ris]=="A" & sites$country[ris]==1)],,drop=0],2,sum),
									apply(transm[ris[which(sites$arm[ris]=="B" & sites$country[ris]==1)],,drop=0],2,sum),
									apply(transm[ris[which(sites$arm[ris]=="C" & sites$country[ris]==1)],,drop=0],2,sum))					
						}))
		idx.A<-	seq(1,nrow(transm),by=3)
		idx.B<-	seq(2,nrow(transm),by=3)
		idx.C<-	seq(3,nrow(transm),by=3)
	}
	else if(method=="pooled across SA")
	{
		transm<- do.call(rbind,tapply(seq_len(nrow(sites)),sites$rid, function(ris)
						{
							rbind(	apply(transm[ris[which(sites$arm[ris]=="A" & sites$country[ris]==2)],,drop=0],2,sum),
									apply(transm[ris[which(sites$arm[ris]=="B" & sites$country[ris]==2)],,drop=0],2,sum),
									apply(transm[ris[which(sites$arm[ris]=="C" & sites$country[ris]==2)],,drop=0],2,sum))					
						}))
		idx.A<-	seq(1,nrow(transm),by=3)
		idx.B<-	seq(2,nrow(transm),by=3)
		idx.C<-	seq(3,nrow(transm),by=3)
	}
	else if(method=="pooled across trial")
	{
		transm<- do.call(rbind,tapply(seq_len(nrow(sites)),sites$rid, function(ris)
						{
							rbind(	apply(transm[ris[which(sites$arm[ris]=="A")],,drop=0],2,sum),
									apply(transm[ris[which(sites$arm[ris]=="B")],,drop=0],2,sum),
									apply(transm[ris[which(sites$arm[ris]=="C")],,drop=0],2,sum)	)					
						}))
		idx.A<-	seq(1,nrow(transm),by=3)
		idx.B<-	seq(2,nrow(transm),by=3)
		idx.C<-	seq(3,nrow(transm),by=3)				
	}
	else if(method=="pooled across T5")
	{
		transm<- do.call(rbind,tapply(seq_len(nrow(sites)),sites$rid, function(ris)
						{
							rbind(	apply(transm[ris[which(sites$arm[ris]=="A" & sites$triplet.id[ris]==5)],,drop=0],2,sum),
									apply(transm[ris[which(sites$arm[ris]=="B" & sites$triplet.id[ris]==5)],,drop=0],2,sum),
									apply(transm[ris[which(sites$arm[ris]=="C" & sites$triplet.id[ris]==5)],,drop=0],2,sum))					
						}))
		idx.A<-	seq(1,nrow(transm),by=3)
		idx.B<-	seq(2,nrow(transm),by=3)
		idx.C<-	seq(3,nrow(transm),by=3)		
	}
	else stop("popart.pool: pooling option not recognized")
	
	ans<- list(transm=transm, idx.A=idx.A, idx.B=idx.B, idx.C=idx.C)
	ans
}	
###############################################################################
#' internal function to perform power calculations
.phdes.binom.power<- function(x,test.prop0, test.prop1, test.alpha, method="homebrew")
{
	if(!method%in%c("homebrew","cloglog", "logit", "probit", "asymp", "lrt", "exact"))
		stop("phdes.binom.power: cannot find method")
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
		ans[ which(is.nan(ans)) ]<- 0
	}
	ans
}	
###############################################################################
phdes.power.ttest.cl<- function(x2i, p.H0, p.H1, var.p.H0, var.p.H1, alpha= 0.05, verbose=0)
{
	ncl<- length(x2i)
	if(length(p.H0)!=ncl)	p.H0<- rep(p.H0,ncl)
	if(length(p.H1)!=ncl)	p.H1<- rep(p.H1,ncl)
	if(length(var.p.H0)!=1)	stop("invalid var.p.H0")
	
	delta	<- mean( p.H1 - p.H0 )
	var.D	<- sum( p.H0*(1-p.H0)/x2i + p.H1*(1-p.H1)/x2i + var.p.H0 + var.p.H1) / ( ncl-1 )
	se.D	<- sqrt(var.D)
	tmp		<- (qnorm(alpha/2)*se.D - delta) / se.D
	pw		<- 1 - pnorm( tmp )
	#print(delta); print(tmp); print(pw)
	pw
}
###############################################################################
#' performs power calculations and computes binomial confidence intervals for a given number of acute 2 acute transmissions
#' @param x2i 			number of transmissions, typically for a specific arm
#' @param i2i			number of acute 2 acute transmissions, typically for a specific arm
#' @param test.prop0	null proportion of transmission during acute phase
#' @param test.prop1 	alternate proportion of transmission during acute phase
#' @param test.alpha	level of test
#' @param method.pw		method used to calculate Binomial power. Possible options are "homebrew","cloglog", "logit", "probit", "asymp", "lrt", "exact".
#' @param method.ci		method used to calculate Binomial confidence intervals. Possible options are "exact", "ac", "asymptotic", "wilson", "prop.test", "bayes", "logit", "cloglog", "probit".
#' @return \item{conf}{binomial confidence interval for each site in x2i and i2i}
#' 		   \item{is.conf}{boolean flag for each site, indicating if test.prop0 is below the confidence interval AND if the i2i exceed 5}
#' 		   \item{power}{binomial power for each site in x2i and i2i}
#' @examples 	sites<- popart.getdata.randomized.arm(1)
#' 	sites<- popart.getdata.randomized.n(sites, 2500, 3)
#' 	x2i<- popart.get.sampled.transmissions(sites)
#'	x2i<- phdes.get.linked.transmissions(x2i,0.8)			
#'	tmp<- popart.pool(sites, x2i)
#'	x2i<- tmp[["transm"]]
#'	idx.A<- tmp[["idx.A"]]
#' 	phdes.binom.power(	x2i[idx.A,,drop=0], round(x2i[idx.A,,drop=0]*test.prop0), test.prop0, test.prop1, test.alpha, verbose=0)
phdes.binom.power<- function(x2i, i2i, test.prop0, test.prop1, test.alpha=0.05, verbose=0, method.pw="homebrew", method.ci="bayes")
{		
	if(!method.ci%in%c("exact", "ac", "asymptotic", "wilson", "prop.test", "bayes", "logit", "cloglog", "probit"))
		stop("method.ci not known")
	if(!method.pw%in%c("homebrew","cloglog", "logit", "probit", "asymp", "lrt", "exact"))
		stop("method.pw not known")	
	if(length(test.prop0)!=ncol(x2i))
		test.prop0<- rep(test.prop0,ncol(x2i))
	if(length(test.prop1)!=ncol(x2i))
		test.prop1<- rep(test.prop1,ncol(x2i))
		
	x			<- apply(x2i,2,mean)		
	y			<- apply(i2i,2,mean)	
	
	if(verbose){	cat("\nphylog linked incident cases\n");	print(x)	}	
	conf<- binom.confint(as.vector(y), as.vector(x), conf.level = 0.95, methods=method.ci)[,c("lower","upper")]	
	rownames(conf)<- colnames(x2i)
	is.conf<- test.prop0<conf[,"lower"] & y>=5
	#suppose  test.prop0<test.prop1<0.5, then the alternative Binom|test.prop0 has smaller variance
	power<- .phdes.binom.power(x, test.prop1, test.prop0, test.alpha, method=method.pw)	 
	names(power)<- colnames(x2i)
	if(verbose){ cat(paste("\npower to distinguish",test.prop0, "from", test.prop1,"\n"));	print(power.hg.A)	}
	
	ans<- list(conf=conf, is.conf=is.conf, power=power)
	ans		
}
