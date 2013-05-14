###############################################################################
ssa.get.transition.compartmental<- function(ibm, verbose=0)
{
#print(ibm)	
	#st.lv<- levels( ibm[["init.pop"]][["status"]] )
	if(ncol(ibm[["curr.pop"]])>2)	stop("ssa.get.transprob: error at 1a")
	ans<- vector("list",4)
	names(ans)<- c("from","from.attr","to","time")
	propens<- acute.get.rates(ibm)	
	#determine transition by risk group x->y	 from: donor risk group x, to: recipient risk group y
	tmp <- 	my.sample(seq_len(length(propens)), size = 1, prob = as.vector(propens), replace= 1)
	ans[["from.attr"]]<- colnames(propens)[ (tmp-1)%%nrow(propens)+1 ]
	to<-	colnames(propens)[ (tmp-1)%/%nrow(propens)+1 ]
#print(from); print(to)
	#determine donor and recpipient individuals 
	setkey(ibm[["curr.pop"]],"status","id")
	old.warn<- getOption("warn")
	options(warn=-1)		#TODO known inefficiency: should convert ibm factors to characters
	ans[["from"]]<- my.sample(ibm[["curr.pop"]][ans[["from.attr"]]][,id], size=1)
	ans[["to"]]<- my.sample(ibm[["curr.pop"]][to][,id], size=1)
	options(warn=old.warn)
	setkey(ibm[["curr.pop"]],"id","status")
	if(verbose)	cat(paste("\n mean reaction time",1/sum(propens)))
	#determine time to transition
	ans[["time"]]<- -log(runif(1)) / sum(propens)
	ans
}
###############################################################################
ssa.run<- function(ibm, ssa.ctime= 0, ssa.etime= 1, verbose= 0, save= '', resume= 1, record.tpc= 1)
{
	#catch<-try({ 
	
	
	#library(compiler)
	#c.tpc.find.treeidx<- cmpfun(tpc.find.treeidx)
	#c.tpc.add.tree<- cmpfun(tpc.add.tree)
	if(resume && nchar(save))
	{
		options(show.error.messages = FALSE)		
		cat(paste("\nssa.run: try to resume file ",save))
		readAttempt<-try(suppressWarnings(load(save)))
		options(show.error.messages = TRUE)
		if(!inherits(readAttempt, "try-error"))
			cat(paste("\nssa.run: resumed file ",save))
	}
	if(!resume || inherits(readAttempt, "try-error"))
	{
		ssa.attackrate<- 0
		ssa.nevent<- 0
		
		if(record.tpc)
		{ 
			tpc.data<- tpc.init(ibm)
			tpc.tree.idx<- rep(NA,nrow(ibm[["curr.pop"]]))
			tpc.tree.idx[ ibm[["curr.pop"]][ ibm[["curr.pop"]][,status!='s'], id] ]<- seq.int(1,length(tpc.data[["trees"]]))			
		}
		else	
			tpc.tree.idx<- tpc.data<- NULL
		while(ssa.ctime<ssa.etime)
		{
			#determine next infection
			ssa.update<- ssa.get.transition.compartmental(ibm, verbose=0)
			#update and store
			ibm[["curr.pop"]][ssa.update$to,"status"]<- 'i'
			if(is.na(tpc.tree.idx[ssa.update$from]))	stop("ssa.run: cannot find tree index")
			tpc.tree.idx[ssa.update$to]<- tpc.tree.idx[ssa.update$from] 
			ssa.ctime<- ssa.ctime + ssa.update$time
			#print(ssa.update)
			#print(tpc.tree.idx)
			#print(tpc.data[["trees"]][[tpc.tree.idx]])
			#next line: expensive
			if(record.tpc) 
				tpc.data[["trees"]][[ tpc.tree.idx[ssa.update$from] ]]<- tpc.add.tree(tpc.data[["trees"]][[ tpc.tree.idx[ssa.update$from] ]], ssa.update$from, ibm[["curr.pop"]][ssa.update$to],ssa.ctime)				
			ssa.nevent<- ssa.nevent + 1						
			if(verbose)
				cat(paste("\nat time",ssa.ctime,"\tindividual",ssa.update$from,"infects",ssa.update$to,"nevent",ssa.nevent,"attack rate",ssa.nevent/nrow(ibm[["curr.pop"]])))
			#store: from, to, branch length		as well as data frame of node attributes of donor at time of infection				
		}
		if(record.tpc)
			tpc.data[["attack.rate"]]<- ssa.nevent / nrow(ibm[["curr.pop"]]) / ssa.etime
		
		if(nchar(save) && record.tpc)
		{
			cat(paste("\nssa.run: save ibm to",save))
			save(tpc.data,file=save)
		}
	}	
	
	#})
	#if(TPC.DEBUG && inherits(catch, "try-error"))
	#{
	#	geterrmessage()
	#	dump.frames(dumpto= "tipc.dump")
	#	cat(paste("\nssa.run: save file dump to",paste(DATA, paste("ssa.run",date(),".rda",sep=''),sep='/')))
	#	save(tipc.dump, file=paste(DATA, paste("ssa.run",date(),".rda",sep=''),sep='/'))
	#	q()
	#}
	
	if(record.tpc)
		return( tpc.data )
	else
		return( ssa.nevent/nrow(ibm[["curr.pop"]]) )
}