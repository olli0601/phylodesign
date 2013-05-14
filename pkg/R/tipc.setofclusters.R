###############################################################################
tpc.init.tree<- function(x, n.init= 5)
{
	tmp<- vector("list",2)
	names(tmp)<- c("nodes","edges")
	tmp$nodes<- rbind(x, matrix(NA, nrow = n.init, ncol = ncol(x), dimnames= list(c(),colnames(x))) )	
	tmp$edges<- matrix(rep(NA,n.init*3),ncol=3,nrow=n.init,dimnames=list(c(),c("from","to","bl")))
	tmp
}
###############################################################################
tpc.add.tree<- function(tree,from,to,ctime, n.init= 5)
{		
	tmp<- as.data.frame(c(to,ctime))
	names(tmp)<- c(names(to),"time")		
	if(!any(is.na(tree[["nodes"]])))
	{
		tree[["nodes"]]<- rbind(tree[["nodes"]], matrix(NA, nrow = n.init, ncol = ncol(tree[["nodes"]]), dimnames= list(c(),colnames(tree[["nodes"]]))) )
	}
	tree[["nodes"]][which(apply(is.na(tree[["nodes"]]),1,any) )[1], ]<- tmp		
	if(any(!is.na(tree[["edges"]])))
	{
		tree[["edges"]]<- rbind(tree[["edges"]], matrix(NA,ncol=ncol(tree[["edges"]]),nrow=n.init,dimnames=list(c(),colnames(tree[["edges"]]))) )
	}
	tree[["edges"]][ which(apply(is.na(tree[["edges"]]),1,any) )[1], 	]<- c(from, tmp[1,"id"], ctime - subset(tree[["nodes"]],id==from,select=time, drop=1))
	tree
}
###############################################################################
tpc.find.treeidx<- function(tpc,catt,cfrom)
{
	tmp<- which(sapply(	seq_along(tpc[["trees"]]),function(i)	
						nrow( subset(tpc[["trees"]][[i]][["nodes"]], id==cfrom & status==catt ) )>0	
			))
	#if there is more than one match (eg because of co-infection) then pick randomly one		
	my.sample(tmp,size=1)
}
###############################################################################
tpc.init<- function(ibm, ctime= 0)
{
	tpc<- vector("list",2)
	names(tpc)<- c("trees","attack.rate")
	tmp<- ibm[["curr.pop"]][ ibm[["curr.pop"]][,status!='s'], ] 	
	tmp<- cbind(as.data.frame(tmp),data.frame(time=rep(ctime,nrow(tmp))))		
	tpc[["trees"]]<- lapply(seq_len(nrow(tmp)),function(i)		tpc.init.tree(tmp[i,,drop=0])	)
	tpc[["attack.rate"]]<- NA
	tpc
}
###############################################################################
tpc.collapse<- function(tpc)
{
	tpc[["trees"]]<- lapply(tpc[["trees"]],function(x)
			{
				tmp<- which(apply(is.na(x[["nodes"]]),1,any) )[1]
				if(is.na(tmp))	#can happen that 'nodes' have been completely filled	
					tmp<- nrow(x[["nodes"]]) + 1
				x[["nodes"]]<- x[["nodes"]][seq_len(tmp-1),]
				x[["edges"]]<- x[["edges"]][seq_len(tmp-2),]
				x
			})
	tpc
}
###############################################################################
tpc.tabulate<- function(tpc)
{	 
	st.lv<- levels( tpc[["trees"]][[1]][["nodes"]][["status"]] )
	st.m<- sapply(tpc[["trees"]],function(x)
			{
				tabulate(x[["nodes"]][,"status"],nbins= length(st.lv))				
			})
	rownames(st.m)<- st.lv
	st.mx<- max(st.m)
	#count all combinations XI^n up to n=st.m
	roots<- st.lv[ !st.lv%in%c('s') ]
	ans<- matrix(NA,nrow=length(roots),ncol=1+st.mx,dimnames=list(roots,seq.int(1,1+st.mx)))
	for(i in seq_along(roots))
		for(j in seq.int(st.mx,0,-1))
		{			
			if(roots[i]=='i')
				tmp<- st.m['u',]==0 & st.m['t',]==0 & st.m['i',]==j
			else
				tmp<- st.m[roots[i],]==1 & st.m['i',]==j
			
			ans[roots[i],j+1]<- ncol(st.m[,tmp,drop=0])
			st.m<- st.m[,!tmp,drop=0]					
		}		
	ans
}
###############################################################################
tpc.tabulate.after.sample<- function(tpc, sid)
{
	st.lv<- levels( tpc[["trees"]][[1]][["nodes"]][["status"]] )
	st.m<- .Call("tipc_tabulate_after_sample", tpc, as.integer(sort(sid)))
	rownames(st.m)<- st.lv
	st.mx<- max(st.m)
	#count all combinations XI^n up to n=st.m
	roots<- st.lv[ !st.lv%in%c('s') ]
	ans<- matrix(NA,nrow=length(roots),ncol=1+st.mx,dimnames=list(roots,seq.int(1,1+st.mx)))
	for(i in seq_along(roots))
		for(j in seq.int(st.mx,0,-1))
		{			
			if(roots[i]=='i')
				tmp<- st.m['u',]==0 & st.m['t',]==0 & st.m['i',]==j
			else
				tmp<- st.m[roots[i],]==1 & st.m['i',]==j
			
			ans[roots[i],j+1]<- ncol(st.m[,tmp,drop=0])
			st.m<- st.m[,!tmp,drop=0]					
		}	
	ans	
}
###############################################################################
tipc.mle<- function(tpc, ibm, theta.acute, theta.base)
{
	theta<- expand.grid(acute= theta.acute, base= theta.base)
	log.lkls<- sapply(seq_len(nrow(theta)),function(i)
			{								
				ibm[["beta"]][['i']][["status"]]['i']<- theta[i,"acute"]
				ibm[["beta"]][["base"]]<- theta[i,"base"]
				lkl.rates<- acute.get.rates(ibm, per.capita.i=1)
				acute.loglkl(tpc, lkl.rates, ibm[["beta"]][["dT"]])									
			})
	theta<- cbind(theta,log.lkls)
	#theta[ which.max(log.lkls), ]
	t(theta[rev(sort(log.lkls,index.return=1)$ix),])					
}
###############################################################################
tpc.mean<- function(tpc)
{
	if(length(tpc)<2)	stop("tpc.mean: only 1 list element")
	tpc.dim<- sapply(tpc,dim)
	tpc.dim<- apply(tpc.dim,1,max)
	ans<- matrix(NA,nrow=tpc.dim[1],ncol=tpc.dim[2],dimnames=list(rownames(tpc[[1]],seq.int(1,tpc.dim[2]))))
	ans<- numeric(prod(tpc.dim))
	ans<- sapply(seq_along(tpc),function(i)
			{		
				tmp<- numeric(prod(tpc.dim))
				tmp[seq_len(prod(dim(tpc[[i]])))]<- as.vector(tpc[[i]])
				tmp
			})
	ans<- apply(ans,1,sum) / length(tpc)
	ans<- matrix(ans, nrow=tpc.dim[1],ncol=tpc.dim[2],dimnames=list(rownames(tpc[[1]]),seq.int(1,tpc.dim[2])))
	ans	
}
###############################################################################
tpc.sd<- function(tpc)
{
	if(length(tpc)<3)	stop("tpc.sd: only <3 list elements")
	tpc.dim<- sapply(tpc,dim)
	tpc.dim<- apply(tpc.dim,1,max)
	ans<- matrix(NA,nrow=tpc.dim[1],ncol=tpc.dim[2],dimnames=list(rownames(tpc[[1]],seq.int(1,tpc.dim[2]))))
	ans<- numeric(prod(tpc.dim))
	ans<- sapply(seq_along(tpc),function(i)
			{		
				tmp<- numeric(prod(tpc.dim))
				tmp[seq_len(prod(dim(tpc[[i]])))]<- as.vector(tpc[[i]])
				tmp
			})
	ans<- apply(ans,1,function(x)
			{
				tmp<- which(x>0)
				if(length(tmp)<3)	return(NA)
				sd(x[tmp])
			} )
	ans<- matrix(ans, nrow=tpc.dim[1],ncol=tpc.dim[2],dimnames=list(rownames(tpc[[1]]),seq.int(1,tpc.dim[2])))
	ans	
}
