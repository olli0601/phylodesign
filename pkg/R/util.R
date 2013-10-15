###############################################################################
# Generic form
'%<-%' = function(l, r, ...) UseMethod('%<-%')
# Binary Operator
'%<-%.lbunch' = function(l, r, ...) {
	Envir = as.environment(-1)
	
	if (length(r) > length(l))
		warning("RHS has more args than LHS. Only first", length(l), "used.")
	
	if (length(l) > length(r))  {
		warning("LHS has more args than RHS. RHS will be repeated.")
		r <- extendToMatch(r, l)
	}
	
	for (II in 1:length(l)) {
		do.call('<-', list(l[[II]], r[[II]]), envir=Envir)
	}
}
###############################################################################
# Used if LHS is larger than RHS
extendToMatch <- function(source, destin) {
	s <- length(source)
	d <- length(destin)
	
	# Assume that destin is a length when it is a single number and source is not
	if(d==1 && s>1 && !is.null(as.numeric(destin)))
		d <- destin
	
	dif <- d - s
	if (dif > 0) {
		source <- rep(source, ceiling(d/s))[1:d]
	}
	return (source)
}
###############################################################################
# Grouping the left hand side
g = function(...) {
	List = as.list(substitute(list(...)))[-1L]
	class(List) = 'lbunch'
	return(List)
}
###############################################################################
my.dumpframes<- function()
{
	geterrmessage()
	dump.frames()
	cat(paste("\nmy.dumpframes dump 'last.dump' to file",paste(DATA,paste("debug_",paste(strsplit(date(),' ')[[1]],collapse='_'),".rda\n",sep=''),sep='/')))
	save(last.dump, file=paste(DATA,paste("debug_",paste(strsplit(date(),' ')[[1]],collapse='_'),".rda",sep=''),sep='/'))
	q()
}
###############################################################################
my.sample <- function(x, ...) x[sample.int(length(x), ...)]
####### vector scan over length(search) many elements ########################################################################
my.which.in<- function(search,set) which( set[ findInterval(search,set) ] == search )		
###############################################################################
print.v<- function(x,cut=3,digits=4,prefix= "simu_",print.char= TRUE, as.R= TRUE)
{
	if(as.R)
	{
		tmp<- paste("c(",paste(c(x,recursive=T),collapse=',',sep=''),')',sep='')
		if(!is.null(names(x)))
			tmp<- paste("{tmp<-", tmp, "; names(tmp)<- ", paste('c("',paste(c(names(x),recursive=T),collapse='", "',sep=''),'")',sep=''), "; tmp}", sep= '', collapse= '')
	}
	else
	{
		if(!is.null(names(x)))
		{
			m<- matrix(NA,nrow=2,ncol=length(x))
			m[1,]<- substr(names(x),1,cut)
			m[2,]<- round( x, digits=digits )
			if(cut==0)		m<- m[2,]
			tmp<- gsub('.',',',paste(prefix,paste(as.vector(m), collapse='_',sep=''),sep=''),fixed=T)
		}
		else
			tmp<- gsub('.',',',paste(prefix,paste(round( x, digits=digits ), collapse='_',sep=''),sep=''),fixed=T)
	}
	if(print.char) print(tmp)
	tmp
}
###############################################################################
my.mkdir<-function(root,data.name)
{
	if(length(dir(root,pattern=paste('^',data.name,'$',sep='')))==0)
		system(paste("mkdir ",paste(root,data.name,sep='/'),sep=''))
}
###############################################################################
my.fade.col<-function(col,alpha=0.5)
{
	return(rgb(col2rgb(col)[1]/255,col2rgb(col)[2]/255,col2rgb(col)[3]/255,alpha))
}
###############################################################################
my.plot.persplocfit<- function(x, y, z, theta= 30, phi= 10, palette= "gray",...)
{	
	nbcol <- 100
	if(palette=="gray")	
		color<- head( rev(gray(seq(0,0.95,len=trunc(nbcol*1.4)))), nbcol)
	else
	{
		jet.colors <- colorRampPalette( c("blue", "green") )
		color <- jet.colors(nbcol)
	}
	# Compute the z-value at the facet centres
	nrz <- nrow(z)
	ncz <- ncol(z)
	zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
	# Recode facet z-values into color indices
	facetcol <- cut(zfacet, nbcol)		
	par(mar=c(0,2.5,0,0))
	pmat<- persp(x, y, z, col = color[facetcol], shade= 0.1, border=NA, ticktype = "detailed", ltheta = 120,theta = theta, phi = phi, expand = 0.5, box=1, ... )
	
	list(pmat=pmat, x= x, y= y, z= z)
}
###############################################################################
my.image<- function(x,y,z, palette= "rob",nlevels=5, cex.points=0, points.col="white", points.pch=20, contour.nlevels=0, contour.col="white",...)
{	
	if(palette=="topo")					image(x,y,z, col=tail( topo.colors(trunc(50*1.4)), 50 ),...)
	else if(palette=="gray")			image(x,y,z, col=head( rev(gray(seq(0.05,.98,len=trunc(50*1.4)))), 50),...)					
	else if(palette=="rob")				image(x,y,z, col=colorRampPalette(c("blue", "orange", "red"), space = "Lab")( 50 ),...)
	else								image(x,y,z, col=colorRampPalette(c("red", "orange", "blue"), space = "Lab")( 50 ),...)	
	if(cex.points>0)
		points(x, y, cex=cex.points, col=points.col, pch=points.pch)
	if(contour.nlevels>0)
		contour(x,y,z,add=TRUE, nlevels=contour.nlevels,col=contour.col)	
}
###############################################################################
my.normalize<- function(x, scale=1, log=TRUE)
{
	if(log)		
		x	<- exp( x - max(x) )
	x / sum( scale*x )	
}
###############################################################################
my.image.get.CI<- function(look, ci=0.95, log=TRUE)
{
	dx		<- c( diff( look$x[1:2] ), diff( look$y[1:2] ) )
	lkl		<- my.normalize(look$z, scale=prod(dx), log=log)
	
	df.CI	<- cbind( data.table( expand.grid(x=look$x, y=look$y) ), lkl=as.numeric(prod(dx)*lkl))
	df.CI[,dummy:=-df.CI[,lkl]]
	setkey(df.CI, dummy)
	df.CI[,cumlkl:=cumsum(df.CI[,lkl])]
	df.CI	<- subset(df.CI, cumlkl<ci)
	setkey(df.CI, x)
	df.CI	<- df.CI[,list(xl= x+c(-0.5,0.5)*dx[1], ymin=rep(min(y)-dx[2]/2,2), ymax=rep(max(y)+dx[2]/2,2)),by='x']
	df.CI	<- rbind( rbind( df.CI[1,],df.CI ), df.CI[nrow(df.CI),] )
	set(df.CI,1L,"ymax",df.CI[1,ymin])
	set(df.CI,nrow(df.CI),"ymin",df.CI[nrow(df.CI),ymax])
	df.CI	<- subset(df.CI,select=c(xl,ymin,ymax))
	setnames(df.CI,"xl","x")
	
	list(lkl=lkl, ci=df.CI, dx=dx)
}
###############################################################################
my.image.smooth<- function(x, y, z, xlab='', ylab='', nrow=50, palette="gray", ncol=50, nlevel=50, theta=.25, tol=1e-8, plot=1, cex.points=0, points.pch=20, points.col="white", contour.nlevels=0, contour.col="white")
{
	require(fields)		
	x		<- matrix(c(x,y),ncol=2,nrow=length(x) )
	out		<- as.image(z,x=x,nrow=nrow, ncol=ncol)
	dx		<- out$x[2]-out$x[1] 
	dy		<- out$y[2]-out$y[1] 
	look	<- image.smooth( out, theta=theta, tol=tol)
	if(plot)
	{
		if(palette=="topo")				col<- tail( topo.colors(trunc(nlevel*1.4)), nlevel )
		else if(palette=="gray")		col<- head( rev(gray(seq(0,.95,len=trunc(nlevel*1.4)))), nlevel)
		else if(palette=="rob")			col<- colorRampPalette(c("blue", "orange", "red"), space = "Lab")( 50 )
		else 							col<- tim.colors(nlevel)
		image.plot(look,xlab=xlab,ylab=ylab, col=col)
		if(cex.points>0)
			points( x , cex=cex.points, col=points.col, pch=points.pch)		
		if(contour.nlevels>0)
			contour(look$x, look$y, look$z, add=TRUE, nlevels=contour.nlevels, col=contour.col)
	}
	look
}
###############################################################################
my.plot.2D.dens<- function(x,y,xlab,ylab,xlim=NA,ylim=NA,width.infl=2,n.hists=5,method="gauss", palette= "rob", persp.theta= -30, persp.phi= 30, zero.abline=TRUE, ...)
{
	if(!method%in%c("gauss","ash","persp"))	stop("plot.2D.dens: exception 1a")
	if(!palette%in%c("topo","heat","gray","rob"))	stop("plot.2D.dens: exception 1b")
	switch(method,
			gauss=
					{
						require(KernSmooth)
						require(fields)
						x.bw<- width.infl*diff(summary(x)[c(2,5)])
						y.bw<- width.infl*diff(summary(y)[c(2,5)])
						if(!x.bw) x.bw<- EPS
						if(!y.bw) y.bw<- EPS
						if(any(is.na(xlim)))	xlim<- range(x)+c(-1.5,1.5)*x.bw
						if(any(is.na(ylim)))	ylim<- range(y)+c(-1.5,1.5)*y.bw
						dens <- bkde2D(cbind(x, y), range.x=list(xlim,ylim),bandwidth=c(x.bw,y.bw))
						contour(dens$x1, dens$x2, dens$fhat,xlab=xlab,ylab=ylab)
						if(zero.abline) abline(v=0,col="black",lty=3,lwd=1.5)
						if(zero.abline) abline(h=0,col="black",lty=3,lwd=1.5)
					},
			ash={
						require(ash)
						if(any(is.na(xlim))) xlim<- range(x)*1.05
						if(any(is.na(ylim))) ylim<- range(y)*1.05
						bins<- bin2(cbind(x, y), ab=rbind(xlim,ylim),nbin=2*c(nclass.Sturges(x),nclass.Sturges(y)))
						f <- ash2(bins,rep(n.hists,2))
						#image(f$x,f$y,f$z, col=rainbow(50,start= 4/6,end=0),xlab=xlab,ylab=ylab)
						if(palette=="topo")					image(f$x,f$y,f$z, col=tail( topo.colors(trunc(50*1.4)), 50 ),xlab=xlab,ylab=ylab,...)
						else if(palette=="gray")			image(f$x,f$y,f$z, col=head( rev(gray(seq(0,.95,len=trunc(50*1.4)))), 50),xlab=xlab,ylab=ylab,...)
						else if(palette=="rob")				image(f$x,f$y,f$z, col=colorRampPalette(c("blue", "orange", "red"), space = "Lab")( 50 ),xlab=xlab,ylab=ylab,...)
						else								image(f$x,f$y,f$z, col=heat.colors( 50 ),xlab=xlab,ylab=ylab,...)
						contour(f$x,f$y,f$z,add=TRUE, nlevels= 5)
						if(zero.abline) abline(v=0,col="black",lty=2,lwd=2)
						if(zero.abline) abline(h=0,col="black",lty=2,lwd=2)
			},
			persp={
				require(ash)
				require(MASS)
				if(any(is.na(xlim))) xlim<- range(x)*1.05
				if(any(is.na(ylim))) ylim<- range(y)*1.05
				bins<- bin2(cbind(x, y), ab=rbind(xlim,ylim),nbin=2*c(nclass.Sturges(x),nclass.Sturges(y)))
				f <- ash2(bins,rep(n.hists,2))
				
				nrz <- nrow(f$z)
				ncz <- ncol(f$z)
				col<- tail(topo.colors(trunc(1.4 * 50)),50)
				fcol      <- col[trunc(f$z / max(f$z)*(50-1))+1]
				dim(fcol) <- c(nrz,ncz)
				fcol      <- fcol[-nrz,-ncz]
				par(mar=c(1/2,1/2,1/2,1/2))
				persp(x=f$x,y=f$y,z=f$z,col= fcol,theta=persp.theta,phi=persp.phi,xlab=xlab,ylab=ylab,zlab='', ticktype= "detailed" )
			})
}