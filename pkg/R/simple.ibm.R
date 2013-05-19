
#' @export
IBM.MODELT<- seq.int(1,1)
names(IBM.MODELT)<- c("Acute")

#' @export
IBM.LOCT<- seq.int(1,3)
names(IBM.LOCT)<- c("SA-A","ZA-A","ZA-C")

#' @export
IBM.IND.ATTRIBUTES<- seq.int(1,5)
names(IBM.IND.ATTRIBUTES)<- c("status","gender","risk","circ","age")


###############################################################################
#Gillespie SSA simulation of tip clusters
#' @export
ibm.init.attributes<- function(mt)
{
	if(!mt%in%names(IBM.MODELT))	stop("ibm.def.attributes: 1a")
	mt<- IBM.MODELT[mt]
	if(mt>1)	stop("ibm.def.attributes: 1b")
	a		<-  vector("list",length(IBM.IND.ATTRIBUTES))
	names(a)<- 	names(IBM.IND.ATTRIBUTES)
	a$status<- 	factor(c('s','i','t','u'))
	a$gender<-	factor(c('f','m'))
	a$risk	<- 	factor(c('l','m','h')) 
	a$circ	<- 	factor(c('y','n'))
	a$age	<- 	factor(c('<17','<19','<22','<25','<30','<35','<40','<45','<50','<55'))
	a
}
###############################################################################
#' @export
ibm.init.incidence<- function(loc)
{
	if(regexpr("ZA",loc)>0)	ans<- 0.0117
	else if(regexpr("SA",loc)>0) ans<- 0.0149
	else	stop("ibm.init.incidence: unknown loc")	
	if(regexpr("-A",loc)>0)	ans<- ans*0.4
	else if(regexpr("-B",loc)>0)	ans<- ans*0.7
	ans
}
###############################################################################
ibm.init.popinitdistributions.popart<- function(attr, popart.community)
{
	data(popart.triplets.130207)
	comm		<- popart.triplets.130207[ popart.community==popart.triplets.130207$comid_old, , drop=0]	
	tmp			<- lapply(seq_along(attr),function(i)
					{
						switch(	names(attr)[i],
								status=	{		
											s	<- 1-comm[,"hivcomb"]/100
											u	<- comm[,"hivcomb"] / 100 * ( 1 - comm[,"artadjust"]/100 )
											t	<- comm[,"hivcomb"] / 100 * comm[,"artadjust"]/100 
											ans	<- c(s,0,t,u)
											names(ans)<- c('s','i','t','u')
											ans
										},
								gender= rep(1/length(attr[[i]]),length(attr[[i]])),
								risk= 	c(0.2,0.4,0.4),
								circ= 	c(0.2,0.8),
								age= 	rep(1/length(attr[[i]]),length(attr[[i]]))
						)
					})
	names(tmp)	<- names(attr)
	tmp$npop	<- round( comm[,"popsize"] * comm[,"p.adults"] ) 
	tmp	
}
###############################################################################
#' @export
ibm.init.popinitdistributions.vanilla<- function(attr, loc, popsize)
{
	if(!loc%in%c("SA-A","ZA-A","ZA-C"))	stop("ibm.init.popinitdistributions: error at 1a")
	
	tmp			<- lapply(seq_along(attr),function(i)
					{
						switch(	names(attr)[i],
								status=	{
									if(loc=="ZA-A")			tmp<- c(0.135, 0.7, 0)					#c(estim prevalence, targeted treatment uptake)							
									else if(loc=="SA-A")	tmp<- c(0.178, 0.7, 0)
									else if(loc=="ZA-C")	tmp<- c(0.135, 0.3*0.6, 0)				#c(estim prevalence, estimated ART of those known HIV+ * known HIV+)
									#else if(loc%in%rownames(POPART.SITES))
									#						tmp<- c(POPART.SITES[loc,"adult HIV prev"],POPART.SITES[loc,"HIV inf on ART"]*POPART.SITES[loc,"known HIV status"],ibm.init.incidence(loc))
														
									ans<- c(1-tmp[1]-tmp[3]/4,tmp[3]/4,tmp[1]*tmp[2],tmp[1]*(1-tmp[2]) )
									names(ans)<- c('s','i','t','u')
									ans
								},
								gender= rep(1/length(attr[[i]]),length(attr[[i]])),
								risk= 	c(0.2,0.4,0.4),
								circ= 	c(0.2,0.8),
								age= 	rep(1/length(attr[[i]]),length(attr[[i]]))
						)
					})
	names(tmp)	<- names(attr)
	tmp$npop	<- popsize
	tmp
}
###############################################################################
#' @export
ibm.init.beta<- function(attr)
{
	#set up components for beta: 
	#rel infectiousness		rel susceptibility		mixing		
	ans<- vector("list",6)
	names(ans)<- c("base","missing","i","s","m","dT")	#base rate, infectiousness, susceptibility, mixing, sampling interval
	ans[["base"]]<- 1 
	ans[["missing"]]<- c(0,0)
	names(ans[["missing"]])<- c("st1","st2")
	ans[["dT"]]<- 1 
	#rel infectiousness
	ans[['i']]<- lapply(seq_along(attr),function(i)
			{
				tmp<- rep(1,length(attr[[i]])) 
				names(tmp)<- attr[[i]]
				tmp								
			})
	names(ans[['i']])<- names(attr)	
	#rel susceptibility
	ans[['s']]<- lapply(seq_along(attr),function(i)
			{
				tmp<- rep(1,length(attr[[i]])) 
				names(tmp)<- attr[[i]]
				tmp				
			})
	names(ans[['s']])<- names(attr)	
	#mixing
	ans[['m']]<- lapply(seq_along(attr),function(i)
			{
				matrix( 1, length(attr[[i]]), length(attr[[i]]), dimnames=list(attr[[i]],attr[[i]]))
			})
	names(ans[['m']])<- names(attr)
	
	ans
}
###############################################################################
#' @export
ibm.set.modelbeta<- function(mt, beta.template, theta)
{
	if(!mt%in%names(IBM.MODELT))	stop("ibm.init.modelbeta: 1a")
	mt<- IBM.MODELT[mt]
	if(mt>1)	stop("ibm.init.modelbeta: 1b")			
	
	if("base"%in%names(theta))	beta.template[["base"]]<- theta["base"]
	if("m.st1"%in%names(theta)) beta.template[["missing"]]["st1"]<- theta["m.st1"]
	if("m.st2"%in%names(theta)) beta.template[["missing"]]["st2"]<- theta["m.st2"]
	
	theta.acute<- ifelse("acute"%in%names(theta), theta["acute"], beta.template[['i']][["status"]][2])	
	beta.template[['i']][["status"]][]<- c(0,theta.acute, 0.1, 1)		#I causes 'theta.acute' times more infections
	beta.template[['s']][["status"]][]<- c(1,0,0,0)						#no superinfection
	
	beta.template
}
###############################################################################
#' @export
ibm.init.pop<- function(attr, distr)
{
	size<- distr$npop
	distr<- lapply(distr,function(x)
			{
				paste("c(",paste(x,sep='',collapse=','),")")
			})	
	tmp<- sapply(seq_along(attr),function(i)
			{
				paste(names(attr)[i],"=my.sample(attr$",names(attr)[i],",size=size, replace=TRUE, prob=",distr[[i]],")",sep='')				
			})
	tmp<- paste("data.table(id=seq_len(size),",paste(tmp,sep='',collapse=','),")")
	ans<- eval(parse(text= tmp))
	setkey(ans,"id","status")
	ans
}
###############################################################################
#' @export
ibm.init.model<- function(m.type, loc.type, m.popsize, theta, save='', resume= 1, init.pop=1)
{	
	data(popart.phylo.com)
	
	old.warn<- getOption("warn")
	options(warn=2)			#turn warnings into error
	require(data.table)
	options(warn=old.warn)
	theta.acute	<- theta[1]
	theta.base	<- theta[2]
	
	if(resume && nchar(save))
	{
		options(show.error.messages = FALSE)		
		cat(paste("\nibm.init.model: try to resume file ",save))
		readAttempt	<-try(suppressWarnings(load(save)))
		options(show.error.messages = TRUE)
		if(!inherits(readAttempt, "try-error"))
			cat(paste("\nibm.init.model: resumed file ",save))
	}
	if(!resume || inherits(readAttempt, "try-error"))
	{
		ibm						<- vector("list",4)
		names(ibm)				<- c("init.pop","init.pop.distr","curr.pop","beta")		
		ibm.att					<- ibm.init.attributes(m.type)
		ibm$beta				<- ibm.init.beta(ibm.att)					
		ibm$beta				<- ibm.set.modelbeta(m.type, ibm$beta, theta)
		if(loc.type%in%popart.phylo.com)
			ibm.distr			<- ibm.init.popinitdistributions.popart(ibm.att, loc.type)
		else
			ibm.distr			<- ibm.init.popinitdistributions.vanilla(ibm.att, loc.type,m.popsize)
		if(init.pop)
		{
			ibm$init.pop		<- ibm.init.pop(ibm.att,ibm.distr)
			ibm$curr.pop		<- ibm$init.pop
		}
		else
			ibm$init.pop.distr	<- ibm.distr
		
		if(nchar(save))
		{
			cat(paste("\nibm.init.model: save ibm to",save))
			save(ibm,file=save)
		}
	}
	ibm
}
###############################################################################
#' @export
ibm.sample<- function(ibm, include.prob)
{
	ibm[["curr.pop"]][runif(nrow(ibm[["curr.pop"]]))<=include.prob,id]	
}
###############################################################################
#' @export
ibm.as.data.table<- function(ibm)
{
	require(data.table)
	if(!is.data.table(ibm[["curr.pop"]]))
	{
		ibm[["curr.pop"]]<- as.data.table(ibm[["curr.pop"]])
		setkey(ibm[["curr.pop"]],"id","status")
	}
	if(!is.data.table(ibm[["init.pop"]]))
	{
		ibm[["init.pop"]]<- as.data.table(ibm[["init.pop"]])
		setkey(ibm[["init.pop"]],"id","status")
	}	
	ibm
}
###############################################################################
#' @export
ibm.collapse<- function(ibm)
{
	inc.i					<- sapply(ibm[["beta"]][['i']],function(x)	any(x!=1)	)
	inc.s					<- sapply(ibm[["beta"]][['s']],function(x)	any(x!=1)	)
	inc.m					<- sapply(ibm[["beta"]][['m']],function(x)	any(x!=1)	)
	collapse.to				<- which( inc.i | inc.s | inc.m )
	if(!length(collapse.to))	stop("ibm.collapse: error at 1a")
	
	ibm[["beta"]][['s']]	<- lapply(collapse.to, function(i)	ibm[["beta"]][['s']][[i]]	)
	ibm[["beta"]][['i']]	<- lapply(collapse.to, function(i)	ibm[["beta"]][['i']][[i]]	)
	ibm[["beta"]][['m']]	<- lapply(collapse.to, function(i)	ibm[["beta"]][['m']][[i]]	)
	
	collapse.to				<- c("id",names(ibm[["beta"]][['s']]))
	if(!is.null(ibm[["curr.pop"]]))
		ibm[["curr.pop"]]	<- ibm[["curr.pop"]][,collapse.to,with=FALSE]
	if(!is.null(ibm[["init.pop"]]))
		ibm[["init.pop"]]	<- ibm[["init.pop"]][,collapse.to,with=FALSE]	
	ibm
}
