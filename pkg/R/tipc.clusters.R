###############################################################################
#' For a given number \code{nIdx} of free slots, compute all possible combinations of subtrees whose total size is \code{nE}.
#' @export
clu.subtrees.find<- function(nIdx, nE)
{
	options(warn=2)
	if(nIdx<1)					stop("invalid nIdx")
	else if(nIdx==1)			return( matrix(nE,1,1) )
	if(nE<0)					stop("invalid nE")
	else if(nE==0)				return(matrix(rep(0,nIdx), ncol=1))
	else if(nE==1)				return(matrix(c(1,rep(0,nIdx-1)), ncol=1 ))	
	odd			<- seq.int(1,2^(nE-1)-1,2)					 
	odd.binary	<- as.matrix( sapply(odd, function(x) rev(as.integer(intToBits(as.integer(x))))) )	
	sticks		<- odd.binary[seq.int(nrow(odd.binary)-nE+1,nrow(odd.binary)),,drop=0]
	sticks		<- cbind(sticks, rep(1,nrow(sticks)))
	
	sticks.keep	<- apply(sticks, 2, function(x)
					{
						y				<- rle(x)						#exclude sticks for which there is a longer run of 0's to the right
						tmp				<- y$values==0
						y$lengths[tmp]	<- y$lengths[tmp] + 1			
						tmp				<- y$values==1 & y$lengths>1	#if there are two consecutive 1's there is a 0 inbetween
						y$values[tmp]	<- 0
						y$lengths[tmp]	<- 1
						y				<- y$lengths[ y$values==0 ]
						all( order(y, decreasing=1)==seq_along(y) )
					})
	sticks		<- t( sticks[,sticks.keep] )	
	comb		<- apply(sticks, 1, function(x)
					{
						tmp<- which(x==1)
						c( tmp[1],diff(tmp),rep(0,ncol(sticks)-length(tmp)) )
					})
	if(nIdx<nrow(comb))		
		comb	<- comb[seq_len(nIdx), apply(comb[seq.int(nIdx+1,nrow(comb)),,drop=0]==0,2,all) ]
	else if(nIdx>nrow(comb))
		comb	<- rbind(comb, matrix(0,nIdx-nrow(comb),ncol(comb)))
	comb		
}
###############################################################################
#' For a matrix of subtrees (across col), compute the number of possible subtree combinations
clu.subtrees.n<- function(subtrees)
{
	nE	<- sum(subtrees[,1])
	apply(subtrees, 2, function(x)
		{
			cnt	<- tabulate(x+1)
			df	<- length(x)-c(0,cumsum( cnt[-length(cnt)] ))
			ans	<- prod( choose(df,cnt) )						#Comb(unique(J));			all ways the different subtrees can be arranged onto the free slots
			cnt	<- x[x!=0]
			df	<- nE-c(0,cumsum(cnt[-length(cnt)]))
			ans	<- ans*prod( choose(df,cnt) )					#times Comb(J); 			all ways the individuals are distributed into trees 
			ans*prod( (x+1)^(x-1) )								#times prod_{j\inJ} N(j)	all ways how subtrees of size j could look like
		})
}
###############################################################################
#' For each tip cluster with \code{n} transmissions (col) return the number of spanning trees that correspond to \code{i} transmissions from the index case (row)
#' @export
clu.tipc.n<- function(closure)
{
	comb			<- sapply(seq.int(1,closure),function(edges)
						{
							c( 	sapply(seq.int(1,edges),function(nIdx)
											{												
												choose(edges,nIdx) * sum( clu.subtrees.n( clu.subtrees.find(nIdx,edges-nIdx) ) )	
											}),
								rep(0,closure-edges) )
						})
	comb			<- cbind(rep(0,nrow(comb)),comb )
	comb			<- rbind(c(1,rep(0,ncol(comb)-1)),comb)
	dimnames(comb)	<- list(paste('idx',seq.int(0,closure),sep=''), paste('n',seq.int(0,closure),sep=''))
	comb
}
###############################################################################
#' @export
clu.p.of.tchain<- function(clu.n, p, log=0 )
{
	if(length(p)!=2)		stop("invalid p.names")
	if(!"E2E"%in%names(p))	stop("invalid p.names")
	if(any(colnames(clu.n)!=paste('n',seq.int(0,ncol(clu.n)-1),sep='')))	stop("invalid clu.n colnames")
	p.E2E		<- p["E2E"]
	p.x2E		<- p[which(names(p)!="E2E")]	
	r.E2E		<- p.E2E/p.x2E
	max.ntransm	<- ncol(clu.n) - 1								#substract 1 to account for n0 case
	
	clu.logp	<- sapply(seq_len(max.ntransm),function(ntr)
			{
				tmp<- log(r.E2E)*seq.int(ntr-1,0,-1)
				if(r.E2E==0)
					tmp[length(tmp)]<- 0
				c(tmp,rep(NA,max.ntransm-ntr))
			})

	clu.logp			<- cbind( rep(NA,nrow(clu.logp)), clu.logp)		#no probability for singletons can be given
	dimnames(clu.logp)	<- list(paste("idx",1:max.ntransm,sep=''),paste("n",seq.int(0,max.ntransm),sep=''))				
	clu.log.x2E			<- sapply(seq.int(0,max.ntransm),function(ntr)	log(p.x2E)*rep(ntr,max.ntransm)	)
	clu.logp			<- clu.logp+clu.log.x2E+log(clu.n)
	if(!log)
		clu.logp<- exp(clu.logp)	
	clu.logp
}
###############################################################################
#' @export
clu.p.of.tchain.rnd.sampling<- function(clu.p, s, mx.s.ntr=ncol(clu.p)-1, log=0 )
{	
	if(any(s<0) | any(s>1))				stop("invalid s")
	if(!all(c("Idx","E")%in%names(s)))	stop("invalid names of s")
	if(mx.s.ntr>ncol(clu.p)-1)			stop("mx.s.ntr must be <= the max number of transmissions")
print(log(clu.p))
stop()
	s.Idx		<- s["Idx"]
	s.E			<- s["E"]
	closure		<- ncol(clu.p)-1
	clu.ps	<- sapply(seq.int(0,mx.s.ntr),function(sntr)			#for each column in new incomplete clu.p; sntr ~ sampled number transmissions
			{
				#print(paste("sntr",sntr))
				if(log)
				{
					clu.ps<- sapply(seq.int(0,sntr),function(index)		#consider the upper triangular part of the matrix including diagonal
							{
								tmp<- sapply(seq.int(0,closure-sntr),function(missed)				#sum over complete clusters with ntr= sntr+missed
												{
													#print(paste("missed",missed))
													tmp<- clu.p[seq.int(index,index+missed)+1, sntr+missed+1]	+ 
															lchoose(seq.int(index,index+missed), index)			+
															lchoose(seq.int(sntr+missed-index, sntr+missed-(index+missed)), sntr-index)
#print(missed)
#print(lchoose(missed,seq.int(0,missed)))
#print(tmp)													
													log( sum(exp(tmp)) ) + log(1-s.E)*missed							#since in model 'Acute' U and T occur only at baseline, can only miss E					
												}) 
								sum(exp(tmp))
							})
#					print( clu.ps )
					clu.ps<- log(clu.ps) + log(s.E)*sntr + log(s.Idx) 
				}
				else
				{
					clu.ps<- sapply(seq.int(0,sntr),function(index)		#consider the upper triangular part of the matrix including diagonal
							{
								#print("AAAAAAAAAAA")
								tmp<- sum( sapply(seq.int(0,closure-sntr),function(missed)				#sum over complete clusters with ntr= sntr+missed
												{
													#see if correct rows/cols seleced:
													#print("HERE")
													#print(c(sntr+1, index+1, missed)); print("clu.p entries"); print(seq.int(index+1,index+1+missed)); print(sntr+1+missed)
													#print("first choose"); print(seq.int(index,index+missed)); print(index)
													#print("second choose"); print(seq.int(sntr+missed-index, sntr+missed-(index+missed))); print(sntr-index)													
													tmp<- clu.p[seq.int(index+1,index+1+missed), sntr+1+missed]
													#print(tmp)													
													tmp<- tmp * choose(seq.int(index,index+missed), index) #this is  * (i+j, i)													
													tmp<- tmp * choose(seq.int(sntr+missed-index, sntr+missed-(index+missed)), sntr-index) #this is  * (n+miss-(i+j), n-i)
													sum(tmp)*(1-s.E)^missed							#since in model 'Acute' U and T occur only at baseline, can only miss E					
												}) )
								#if(index==1) stop()
								tmp
							})
					clu.ps<- clu.ps * ( s.E^sntr * s.Idx )
				}	
				
				c(clu.ps, rep(ifelse(log,-Inf,0),mx.s.ntr-sntr))				
			})			
	dimnames(clu.ps)<- list(paste("idx",seq.int(0,mx.s.ntr),sep=''),paste("ns",seq.int(0,mx.s.ntr),sep=''))
	clu.ps
}
###############################################################################
#' @export
clu.p.of.tipc.rnd.sampling<- function(clu.p, s, mx.s.ntr=length(clu.p)-1, log=0, rtn.only.closure.sum=0 )
{
	if(any(clu.p[!is.na(clu.p)]<0))			stop("invalid clu.p - not log(clu.p) (?)")
	if(sum(clu.p)!=1) 						stop("clu.p must be normalized")
	if(any(s<0) | any(s>1))					stop("invalid s")
	if(!all(c("Idx","E")%in%names(s)))		stop("invalid names of s")
	if(mx.s.ntr>length(clu.p)-1)			stop("mx.s.ntr must be <= the max number of transmissions")
	
	s.Idx		<- s["Idx"]
	s.E			<- s["E"]
	closure		<- length(clu.p)-1	
	clu.ps		<- sapply(seq.int(0,mx.s.ntr),function(sntr)
			{
				tmp<- sum( sapply(seq.int(0,closure-sntr),function(missed)
								{
									(1-s.E)^missed * choose(sntr+missed,sntr) * clu.p[sntr+missed+1]
								}) )				
				tmp<- tmp * s.Idx * s.E^sntr
				tmp
			})
	names(clu.ps)<- paste("ns",seq.int(0,mx.s.ntr),sep='')	
	clu.ps		<- clu.ps / sum(clu.ps)
	if(rtn.only.closure.sum)
		clu.ps	<- clu.ps / (s.Idx * s.E^seq.int(0,mx.s.ntr))
	if(log)
		clu.ps	<- log(clu.ps)		
	clu.ps 
}
###############################################################################
#' @export
clu.exp.X2E<- function(clu.p)
{
	if(any(clu.p[!is.na(clu.p)]<0))	stop("invalid clu.p - not log(clu.p) (?)")
	
	norm	<- apply(clu.p,2,function(x)  sum(x,na.rm=1))
	if(!all(norm==1))	
		clu.p	<- clu.p / matrix( rep(norm,each=nrow(clu.p)), nrow(clu.p), ncol(clu.p) )
	#print(clu.p)	
	n.Idx2E	<- apply( clu.p*seq_len(ncol(clu.p)), 2, function(x) sum(x,na.rm=1)	)	
	tmp		<- sapply(seq_len(ncol(clu.p)),function(ntr)
			{
				c(seq.int(ntr-1,0,-1),rep(0,ncol(clu.p)-ntr))
			})
	#print(tmp)
	#stop()
	n.E2E	<- apply( clu.p*tmp, 2, function(x) sum(x, na.rm=1) )
	list(n.E2E=n.E2E, n.Idx2E=n.Idx2E)	
	
}
###############################################################################
#' @export
clu.p.init<- function(p.E2E, l.U2E, l.T2E, p.O2E)
{
	tmp<- c(p.E2E, c(l.U2E, l.T2E) / c(l.U2E+l.T2E) * (1-p.O2E-p.E2E), p.O2E)
	names(tmp)<- c("E2E","U2E","T2E","O2E")
	tmp 
}
###############################################################################
#' @export
clu.probabilities<- function(clu.n, p, with.ntr.weight=0)
{
	if(length(p)!=4)	stop("invalid transmission probabilities")
	if(!all(c("E2E","U2E","T2E","O2E")%in%names(p)))	stop("missing transmission probabilties")		
	if(any(colnames(clu.n)!=paste('n',seq.int(0,ncol(clu.n)-1),sep='')))	stop("invalid clu.n colnames")
	
	clu.p			<- lapply(list( p[c("E2E","U2E")], p[c("E2E","T2E")], p[c("E2E","O2E")] ),function(x)
			{							
				clu.p.of.tchain(clu.n, x )
			})
	clu.nchain		<- apply( clu.n, 2, sum )
	clu.p			<- t( sapply(seq_along(clu.p),function(i)
					{
						tmp			<- apply(clu.p[[i]],2,function(x)	sum(x,na.rm=1)	)					
						if(with.ntr.weight)
							tmp		<- tmp / clu.nchain
						tmp['n0']	<- NA						#tip cluster probability not known for singleton from proportions  
						tmp
					}) )				
	rownames(clu.p)	<- c("U","T","O")	
	#stop()
	clu.p			<- clu.p / sum(clu.p, na.rm=1) 		#in or out, the final clu will be the same as this only changes the 'tmp' factor
	clu.p	
}
###############################################################################
#' @export
clu.simulate<- function(clu.p, inc, rtn.int=0)
{
	tmp				<- sum(apply(clu.p,2,function(x) sum(x, na.rm=1) )*seq_len(ncol(clu.p)))		#sum of incidence in clu.p
	clu				<- lapply(inc, function(x)  x/tmp*clu.p )
	if(rtn.int)
		clu			<- lapply(inc, round)
	if(length(clu)==1)
		clu			<- clu[[1]]
	clu
}
###############################################################################
#' @export
clu.exp.transmissions<- function(clu, clu.n, p, s, mx.s.ntr, exclude.O= 1)
{
	if(length(p)!=4)	stop("invalid transmission probabilities")
	if(!all(c("E2E","U2E","T2E","O2E")%in%names(p)))	stop("missing transmission probabilties")
	if(any(s<0) | any(s>1))	stop("invalid s")
	if(!all(c("Idx","E")%in%names(s)))	stop("invalid names of s")
	
	clu.ps			<- lapply(list( p[c("E2E","U2E")], p[c("E2E","T2E")], p[c("E2E","O2E")] ),function(x)
			{							
				tmp	<- clu.p.of.tchain(clu.n, x )
				#norm<- apply(tmp,2,function(x)  sum(x,na.rm=1))
				#tmp	<- tmp / matrix( rep(norm,each=nrow(tmp)), nrow(tmp), ncol(tmp) )
				#print(tmp)
				tmp	<- clu.p.of.tchain.rnd.sampling(tmp, s, mx.s.ntr )							
				norm<- apply(tmp,2,function(x)  sum(x,na.rm=1))
				tmp	<- tmp / matrix( rep(norm,each=nrow(tmp)), nrow(tmp), ncol(tmp) )									
				#stop()
				tmp
			})
	tmp				<- lapply(clu.ps, clu.exp.X2E )	
	clu.p.E2E		<- t( sapply(seq_along(tmp),function(i)	tmp[[i]][["n.E2E"]]	) )
	clu.p.x2E		<- t( sapply(seq_along(tmp),function(i)	tmp[[i]][["n.Idx2E"]]	) )
	#print(clu.p.E2E); print(clu.p.x2E); print(clu[,seq_len(mx.s.ntr)])
	#stop()
	clu.E2E			<- clu.p.E2E * clu[,seq_len(mx.s.ntr)]
	clu.x2E			<- clu.p.x2E * clu[,seq_len(mx.s.ntr)]
	#print(clu.x2E)
	if(!exclude.O)
	{
		ans			<- c( sum(clu.E2E), apply(clu.x2E,1,sum) )
		names(ans)	<- c("E2E","U2E","T2E","O2E")
	}
	else
	{
		ans			<- c( sum(clu.E2E[-3,]), apply(clu.x2E[-3,],1,sum) )
		names(ans)	<- c("E2E","U2E","T2E")
	}
	ans
}
###############################################################################
#' For a complete tip cluster contingency table \code{clu}, compute a downsampled version based on sampling probabilities that are random conditional on dependent variables 
#' @export
clu.sample<- function(clu, s, mx.s.ntr=ncol(clu)-1, rtn.exp=0)
{
	if(any(is.na(clu))) 						stop("invalid clu")
	if(any(clu<0))								stop("invalid clu")	
	if(any(s<0) | any(s>1))						stop("invalid s")
	if(!all(c("Idx","E")%in%names(s)))									stop("invalid names of s")
	if(any(colnames(clu)!=paste('n',seq.int(0,ncol(clu)-1),sep='')))	stop("invalid clu.n colnames")
	if(mx.s.ntr>ncol(clu)-1)					stop("mx.s.ntr must be smaller than the max number of transmissions")
	print(clu)
	s.Idx		<- s["Idx"]	
	#print(s.Idx)	
	s.E			<- s["E"]
	closure		<- ncol(clu) - 1	
	clu.noE		<- rownames(clu)!='i'
	
	if(!rtn.exp)
	{	
		stop("double check")
		clu.s			<- t(sapply(seq_len(nrow(clu)),function(i)
						{
							tmp<- sapply(seq.int(0,closure),function(j)
									{
										tmp<- rbinom(clu[i,j+1],j,s.E)		#includes 0's
										tabulate( tmp, nbins=ncol(clu) )	#first elt is number of tip clusters of size 1 -- 0's are automatically dropped
									})
							tmp<- apply(tmp, 1, sum)
							rbinom(length(tmp),tmp,s.Idx)				
						}))
		colnames(clu.s)	<- paste("ns",seq_len(ncol(clu)),sep='')	
		rownames(clu.s)	<- rownames(clu)
		clu.s			<- clu.s[,seq_len(mx.s.ntr)]
	}
	else
	{
		clu.s.noIdx		<- sapply(seq.int(0,mx.s.ntr),function(sntr)			#for each column; sntr ~ sampled number of transmissions
				{					
					tmp<- sapply(seq.int(0,closure-sntr),function(missed)	#sum over larger tip clusters up to closure
							{
								(1-s.E)^missed * choose(sntr+missed,sntr) * clu[,sntr+missed+1]		#for each row, the sampling probability is the same in model 'Acute'
							})					
					apply(tmp,1,sum) * s.E^sntr					
				})
		print(clu.s.noIdx)
		clu.s			<- clu.s.noIdx * s.Idx
		clu.addToE		<- c(apply(clu.s.noIdx[clu.noE,]*(1-s.Idx),2,sum)[-1],0)
		clu.s['i',]		<- clu.s['i',] + clu.addToE
		colnames(clu.s)	<- paste("ns",seq.int(0,mx.s.ntr),sep='')
	}
	clu.s 
}