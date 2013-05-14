#! /Library/Frameworks/R.framework/Versions/2.15/Resources/bin/Rscript
#--DSCR------- #! /opt/apps/R-2.15.1/lib64/R/bin/Rscript
#--CX1-------- #! /apps/R/2.13.0/lib64/R/bin/Rscript
###############################################################################
#
#
# Author: oliver ratmann
# file: phdes.startme.R
#
# usage from R:
#> setwd("/Users/cfraser/git/phylodesign/pkg")
#> setwd("/Users/Oliver/git/phylodesign/pkg")
#> source("misc/phdes.startme.R")
# usage from bash:
#> misc/phdes.startme.R --help
#
#
# Installation instructions:
#	1 make sure the first line in this file points to your Rscript
#	2.1	create a directory CODE.HOME and set the CODE.HOME variable below to this path
#	2.2	create CODE.HOME/src_tipclust and copy all the R files into this directory
#	2.3 create a directory HOME and set the HOME variable below to this path
#
# How to run popart power calculations
#	1. open a terminal 
#	2.1 	cd to CODE.HOME/src_tipclust
#	2.2		for the link power calculations, always run ./tipc.startme.R -pPOPART.POWER.LINK.CONSENT
#	2.3		for the tip cluster power calculations, always run ./tipc.startme.R -pPOPART.POWER.TIPC.CONSENT
#	2.4		this will run a standard parameterization for the power calculation
#	2.5		to change the standard parameterization, command line options are available, see -help.
#
###############################################################################
args <- commandArgs()
if(!any(args=='--args'))
	args<- vector("numeric",0)
if(any(args=='--args'))
	args<- args[-(1:match("--args", args)) ]

CODE.HOME<<- "/Users/Oliver/git/phylodesign/pkg"
#CODE.HOME<<- "/Users/cfraser/git/phylodesign/pkg"
#CODE.HOME<<- "/work/or105/libs/popartlib"
#HOME<<- "/Users/cfraser/phylodesign"
HOME<<- "/Users/Oliver/workspace_sandbox/popart"
#HOME<<- "/home/koelle/or7/popart"
#HOME<<- "/work/or105/popart"
DATA<<- paste(HOME,"data",sep='/')
PHDES.DEBUG<<- 0
#default.fun	<- "prj.popart.tchain_test"
#default.fun	<- "prj.popart.powercalc_link_consenting"
#default.fun	<- "prj.popart.powercalc_tipc_test"
#default.fun		<- "prj.popart.powercalc_tipc_test_ukhivrdb"
default.fun		<- "prj.popart.powercalc_tipc_test_residual"
#default.fun	<- "prj.popart.powercalc_tipc_consenting"
#default.fun	<- "prj.popart.power_test"
#default.fun	<- "prj.popart.powercalc_cmp_link_tipc"
#default.fun	<- "prj.popart.powercalc_tipc_contam"

#default.fun<- "prj.popart.powercalc_link_consenting"
#default.fun<- "prj.plotfisherhettransm"

###############################################################################
#if(length(args) && !is.loaded("tipc_tabulate_after_sample"))
#{
#	file<- paste(CODE.HOME,"src_tipclust",paste("libtipcr",.Platform$dynlib.ext,sep=''),sep='/')
#	cat(paste("\nloading",file,'\n',sep=' '))
#	dyn.load(file)
#}
#cat(paste("is.loaded('tipcr')->",is.loaded("tipc_tabulate_after_sample"),'\n'))
###############################################################################
function.list<-list.files(path= paste(CODE.HOME,"R",sep='/'), pattern = ".R$", all.files = FALSE,
		full.names = TRUE, recursive = FALSE)
#function.list<- function.list[-which(sapply(function.list,function(n){identical(n, paste(CODE.HOME,"/misc/phdes.startme.R",sep='') )}))]
sapply(function.list,function(x) source(x,echo=FALSE,print.eval=FALSE, verbose=FALSE))
###############################################################################
my.mkdir<-function(root,data.name)
{
	if(length(dir(root,pattern=paste('^',data.name,'$',sep='')))==0)
		system(paste("mkdir ",paste(root,data.name,sep='/'),sep=''))
}
###############################################################################
my.make.documentation<- function()
{
	require(roxygen2)		
	roxygenize(CODE.HOME)
}
###############################################################################
my.fade.col<-function(col,alpha=0.5)
{
	return(rgb(col2rgb(col)[1]/255,col2rgb(col)[2]/255,col2rgb(col)[3]/255,alpha))
}
###############################################################################
my.mkdir(HOME,"data")
my.mkdir(HOME,"pdf")
my.mkdir(HOME,"script")
argv<- list()
if(length(args))
{
	tmp<- na.omit(sapply(args,function(arg)
					{
						switch(substr(arg,2,4),
								exe= return(substr(arg,5,nchar(arg))),
								NA)
					}))
	if(length(tmp)!=0)
	{
		if(length(tmp)>1) stop("phdes.startme.R: duplicate -exe")
		else default.fun<- switch(tmp[1],
					POPART.POWER.LINK.CONSENT	= "prj.popart.powercalc_link_consenting",
					POPART.POWER.TIPC.CONSENT	= "prj.popart.powercalc_tipc_consenting",
					TEST						= "prj.test",
					SIMU.DATA					= "prj.simudata",
					ACUTE.LKL					= "prj.acute.loglklsurface",
					ACUTE.ABC					= "prj.acutesampling.rejabc",
					WH.SLEEPER					= "prj.wh.sleeper",
					MAKE.DOCUMENTATION		 	= "my.make.documentation"
					)
	}
	tmp<- na.omit(sapply(args,function(arg)
					{
						switch(substr(arg,2,6),
								debug= 1,
								NA)
					}))		
	if(length(tmp)!=0)	PHDES.DEBUG<<- tmp[1]		
	tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,9),
									phsignal= return(as.numeric(substr(arg,10,nchar(arg)))),NA)	}))
	if(length(tmp)>0) argv$p.phylosignal<- tmp[1]
	tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,9),
									nocontam= return(as.numeric(substr(arg,10,nchar(arg)))),NA)	}))
	if(length(tmp)>0) argv$p.nocontam<- tmp[1]
	tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,17),
									instudy.clu.armC= return(as.numeric(substr(arg,18,nchar(arg)))),NA)	}))
	if(length(tmp)>0) argv$p.prev.instudy.clu.armC<- tmp[1]
	tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,11),
									opt.pooled= return(substr(arg,12,nchar(arg))),NA)	}))
	if(length(tmp)>0) argv$opt.pooled<- tmp[1]
	tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,10),
									opt.power= return(substr(arg,11,nchar(arg))),NA)	}))
	if(length(tmp)>0) argv$opt.power<- tmp[1]
	tmp<- na.omit(sapply(args,function(arg)
						{	switch(substr(arg,2,5),
									help= {	cat(paste("\n
												Help file for function ",default.fun,"\n												
												Command line options are\n\t
												-nocontam\tset prob of no contamination\n\t
												-instudy.clu.armC\tset prob of linkage to HCC in arm C\n\t
												-phsignal\tset prob that truly linked individual can be detected with phylogenetic methods, defaults to 1\n\t
												-opt.power\tset prower option. All: prevalent cases in HCC and incident cases in HCC \n\t
												-opt.pooled\tset pooling option. defaults to None\n\n
												Example: ./tipc.startme.R -pPOPART.POWER.LINK.CONSENT -nocontam0.55 -instudy.clu.armC0.2 -opt.pooledSA2\n\n"))
										q()
									})}))							
}
###############################################################################
if(PHDES.DEBUG)	options(error= my.dumpframes)	
cat(paste("\nphdes.startme.R: ",ifelse(PHDES.DEBUG,"debug",""),"call",default.fun))
do.call(default.fun,argv) 	
cat("\ntipc: ",ifelse(PHDES.DEBUG,"debug","")," end\n")
