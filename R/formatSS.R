#Copyright 2024 Doug Speed.

#    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

#    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.


################################################


#R code for putting GWAS results in LDAK format
#requires cors.bim file, so can work out consistent SNPs


################################################


formatSS=function(gwasfile=NULL, corstem=NULL, outfile=NULL, NameCol=NULL, A1Col=NULL, A2Col=NULL, ZCol=NULL, EffectCol=NULL, SECol=NULL, PCol=NULL, nCol=NULL, ncasesCol=NULL, ncontrolsCol=NULL, fixedn=NULL, headerRows=NULL)
{
################
#get start time
start_time=Sys.time()

end_time=Sys.time()
cat(paste0("Start at ",start_time,"\n\n"))


################
#check required arguments

if(is.null(gwasfile))
{return(paste0("Error, you must use the argument gwasfile to specify the GWAS results file"))}

if(is.null(corstem))
{return(paste0("Error, you must use the argument corstem to specify the correlation files"))}

if(!is.null(headerRows))
{
if(headerRows<1|headerRows!=round(headerRows))
{return(paste0("Error, the argument skipRows must be a positive integer (not ",headerRows,")"))}
}


################
#check required files exist

if(file.exists(gwasfile)==FALSE)
{
return(paste0("Error, unable to find the GWAS results file ",gwasfile))
}

filename=paste0(corstem,".cors.bim")
if(file.exists(filename)==FALSE)
{
return(paste0("Error, unable to find the file ",filename," (make sure you have downloaded the files from https://dougspeed.com/correlations)"))
}


################
#read the column names of gwas file

cat(paste0("Reading the top of ",gwasfile,"\n"))

if(is.null(headerRows))	#find the header row
{
headerRows=1
while(1)
{
gwas_head=read.table(gwasfile,head=FALSE,nrow=1,comment.char="",skip=headerRows-1,sep=";")
if(substr(gwas_head[1],1,2)!="##"){break}

cat(paste0("Row ", headerRows," starts ", substr(gwas_head[1],1,30)," - this will be treated as a comment and ignored\n"))
headerRows=headerRows+1
}
}

gwas_head=read.table(gwasfile,head=FALSE,nrow=1,comment.char="",skip=headerRows-1)
cat(paste0("Row ", headerRows," starts ", substr(gwas_head[1],1,30)," - this will be used as the header row\n\n"))
gwas_head=read.table(gwasfile,head=TRUE,nrow=1,comment.char="",skip=headerRows-1)

#first assume space delimited
comma_sep=0

if(ncol(gwas_head)==1)	#try comma delimited
{
comma_sep=1
gwas_head=read.table(gwasfile,head=TRUE,nrow=1,comment.char="",comment.char="",skip=headerRows-1,sep=",")

if(ncol(gwas_head)==1)
{return(paste0("Error, ", gwasfile," appears to have only one column (have tried to read using both spaces and commas as delimiters)"))}
}

if(ncol(gwas_head)<5)
{return(paste0("Error, ", gwasfile, "should have at least five columns (not ", ncol(gwas_head),")"))}

cat(paste0("The file ", gwasfile," has the following ", ncol(gwas_head), " columns:\n"))
for(j in 1:ncol(gwas_head))
{
if(j==21){cat(paste0("Warning, only the first 20 columns are shown\n"))}

cat(paste0("Column ",j," is ", colnames(gwas_head)[j],"\n"))
}
cat("\n")


################
#check for inconsistent column arguments

if(!is.null(ZCol)&!is.null(EffectCol))
{return(paste0("Error, you can only use one of ZCol and EffectCol"))}
if(!is.null(SECol)&is.null(EffectCol))
{return(paste0("Error, you can only use SECol when also using EffectCol"))}
if(!is.null(PCol)&is.null(EffectCol))
{return(paste0("Error, you can only use PCol when also using EffectCol"))}
if(!is.null(SECol)&!is.null(PCol))
{return(paste0("Error, you can only use one of SECol and PCol"))}

if((!is.null(ncasesCol)|!is.null(ncontrolsCol))&!is.null(nCol))
{return(paste0("Error, you not use ncasesCol or ncontrolsCol when using nCol"))}
if(is.null(ncasesCol)+is.null(ncontrolsCol)==1)
{return(paste0("Error, you must use both ncasesCol and ncontrolsCol"))}
if(!is.null(fixedn)&!is.null(nCol))
{return(paste0("Error, you not use fixedn when using nCol"))}


################
#see how many pieces of information provided

count=0

if(!is.null(NameCol)){count=count+1}
if(!is.null(A1Col)){count=count+1}
if(!is.null(A2Col)){count=count+1}
if(!is.null(ZCol)|!is.null(EffectCol)!=0){count=count+1}
if(!is.null(ZCol)|!is.null(SECol)|!is.null(PCol)){count=count+1}
if(!is.null(nCol)|!is.null(ncasesCol)|!is.null(fixedn)){count=count+1}

if(count==0)
{
cat(paste0("You should now rerun this function adding the following arguments:\
use NameCol to specify the column containing predictor names (look for labels such as Predictor, SNP, Marker or rsID);\
use A1Col to specify the column containing the A1 allele (look for labels such as A1 or EffectAlle);\
use A2Col to specify the column containing the A2 allele (look for labels such as A2 or OtherAllele);\
use EffectCol to specify the column containing effect sizes (look for labels such as Effect, Beta, OR or LogOR);\
use SECol to specify the column containing standard errors (look for labels such as SE or SD);\
and use nCol to specify the column containing sample sizes (look for labels such as n, N or nEff)\n\
If SEs are not provided, you can use PCol to specify the column containing p-values,\
while if Gaussian test statistics are provided, you can use ZCol instead of EffectCol and SECol\n\
If the sample size is divided into numbers of cases and controls, you can replace nCol with ncasesCol and ncontrolsCol\
(if sample size is not provided, you can instead use fixedn to provide the total sample size)\n\n"))
return(invisible())
}
if(count<6)
{
cat(paste0("Error, you must provide the following arguments:\
use NameCol to specify the column containing predictor names (look for labels such as Predictor, SNP, Marker or rsID);\
use A1Col to specify the column containing the A1 allele (look for labels such as A1 or EffectAlle);\
use A2Col to specify the column containing the A2 allele (look for labels such as A2 or OtherAllele);\
use EffectCol to specify the column containing effect sizes (look for labels such as Effect, Beta, OR or LogOR);\
use SECol to specify the column containing standard errors (look for labels such as SE or SD);\
and use nCol to specify the column containing sample sizes (look for labels such as n, N or nEff)\n\
If SEs are not provided, you can use PCol to specify the column containing p-values,\
while if Gaussian test statistics are provided, you can use ZCol instead of EffectCol and SECol\n\
If the sample size is divided into numbers of cases and controls, you can replace nCol with ncasesCol and ncontrolsCol\
(if sample size is not provided, you can instead use fixedn to provide the total sample size)\n\n"))
return(invisible())
}


################
#find indexes

Name_find=which(colnames(gwas_head)==NameCol)
if(length(Name_find)==0)
{return(paste0("Error, there is not a column called ", NameCol))}

A1_find=which(colnames(gwas_head)==A1Col)
if(length(A1_find)==0)
{return(paste0("Error, there is not a column called ", A1Col))}

A2_find=which(colnames(gwas_head)==A2Col)
if(length(A2_find)==0)
{return(paste0("Error, there is not a column called ", A2Col))}

if(!is.null(ZCol))
{
Z_find=which(colnames(gwas_head)==ZCol)
if(length(Z_find)==0)
{return(paste0("Error, there is not a column called ", ZCol))}
}
else
{
Effect_find=which(colnames(gwas_head)==EffectCol)
if(length(Effect_find)==0)
{return(paste0("Error, there is not a column called ", EffectCol))}

if(!is.null(SECol))
{
SE_find=which(colnames(gwas_head)==SECol)
if(length(SE_find)==0)
{return(paste0("Error, there is not a column called ", SECol))}
}
else
{
P_find=which(colnames(gwas_head)==PCol)
if(length(P_find)==0)
{return(paste0("Error, there is not a column called ", PCol))}
}
}

if(!is.null(nCol))
{
n_find=which(colnames(gwas_head)==nCol)
if(length(n_find)==0)
{return(paste0("Error, there is not a column called ", nCol))}
}
if(!is.null(ncasesCol))
{
ncases_find=which(colnames(gwas_head)==ncasesCol)
if(length(ncases_find)==0)
{return(paste0("Error, there is not a column called ", ncasesCol))}
ncontrols_find=which(colnames(gwas_head)==ncontrolsCol)
if(length(ncontrols_find)==0)
{return(paste0("Error, there is not a column called ", ncontrolsCol))}
}

cat(paste0("Will read predictor names from the column called ", colnames(gwas_head)[Name_find],"\n"))
cat(paste0("Will read A1 alleles from the column called ", colnames(gwas_head)[A1_find],"\n"))
cat(paste0("Will read A2 alleles from the column called ", colnames(gwas_head)[A2_find],"\n"))

if(!is.null(ZCol))
{cat(paste0("Will read Gaussian test statistics from the column called ", colnames(gwas_head)[Z_find],"\n"))}
else
{
cat(paste0("Will read effect sizes from the column called ", colnames(gwas_head)[Effect_find],"\n"))
if(!is.null(SECol))
{cat(paste0("Will read Ses from the column called ", colnames(gwas_head)[SE_find],"\n"))}
else
{cat(paste0("Will read p-values from the column called ", colnames(gwas_head)[P_find],"\n"))}
}

if(!is.null(nCol))
{cat(paste0("Will read sample sizes from the column called ", colnames(gwas_head)[n_find],"\n"))}
if(!is.null(ncasesCol))
{cat(paste0("Will read numbers of cases and controls from the columns called ", colnames(gwas_head)[ncases_find]," and ",colnames(gwas_head)[ncases_find],"\n"))}
if(!is.null(fixedn))
{cat(paste0("The sample size will be set to ", fixedn,"\n"))}
cat(paste0("\n"))

if(!is.null(EffectCol))
{
if(EffectCol=="OR")
{cat(paste0("Warning, it is assumed the column called ", EffectCol, " contains Odds Ratios, so will compute its logarithm\n\n"))}
}


################
#set some file names

rootfile=paste0(corstem,".cors.root")
bimfile=paste0(corstem,".cors.bim")

if(is.null(outfile))
{
outfile=paste0(gwasfile,".summaries")
cat(paste0("The summary statistics will be saved in the file ", outfile," (change this using the argument outfile)\n\n"))
}


################
#read rootfile and extract total number of predictors

root=read.table(rootfile,head=FALSE)
num_preds=as.numeric(root[5,2])


################
#read bimfile

cat(paste0("Reading predictor details from ",bimfile,"\n"))
bim=read.table(bimfile,head=F,nrow=num_preds)
cat(paste0("In total, there are ",nrow(bim)," predictors\n\n"))


################
#read details from gwasfile

req_columns=c(Name_find,A1_find,A2_find)
req_classes=c("character","character","character")

if(!is.null(ZCol)){req_columns=c(req_columns,Z_find);req_classes=c(req_classes,"numeric")}
else
{
if(!is.null(SECol)){req_columns=c(req_columns,Effect_find,SE_find);req_classes=c(req_classes,"numeric","numeric")}
else{req_columns=c(req_columns,Effect_find,P_find);req_classes=c(req_classes,"numeric","numeric")}
}

if(!is.null(nCol)){req_columns=c(req_columns,n_find);req_classes=c(req_classes,"numeric")}
if(!is.null(ncasesCol)){req_columns=c(req_columns,ncases_find,ncontrols_find);req_classes=c(req_classes,"numeric","numeric")}

classes_all=rep("NULL",ncol(gwas_head))
classes_all[req_columns]=req_classes

cat(paste0("Reading GWAS results from ",gwasfile,"\n"))
if(comma_sep==0){gwas_all=read.table(gwasfile,head=TRUE,colClasses=classes_all,comment.char="",skip=headerRows-1)}
else{gwas_all=read.table(gwasfile,head=TRUE,colClasses=classes_all,comment.char="",skip=headerRows-1,sep=",")}

cat(paste0("In total, there are results for ", nrow(gwas_all)," predictors\n\n"))


################
#find valid predictors

diff_alleles=which(gwas_all[,2]!=gwas_all[,3])
if(length(diff_alleles)==0)
{
return(paste0("Error, none of the predictors in ",gwasfile," have distinct alleles"))
}

common_start=intersect(gwas_all[diff_alleles,1],bim[,2])
if(length(common_start)==0)
{
return(paste0("Error, there is no overlap between the predictors in ",gwasfile," and those in ",bimfile))
}

cat(paste0("There are ",length(common_start)," predictors common to ",gwasfile," and ",bimfile,"\n"))


#use only predictors with consistent alleles

match1=match(common_start,bim[,2])
match2=match(common_start,gwas_all[,1])

find_a1= (bim[match1,5]==gwas_all[match2,2]) + (bim[match1,5]==gwas_all[match2,3])
find_a2= (bim[match1,6]==gwas_all[match2,2]) + (bim[match1,6]==gwas_all[match2,3])

common_end=common_start[which(find_a1==1&find_a2==1)]
if(length(common_end)==0)
{
return(paste0("Error, none of the these have consistent alleles"))
}

if(length(common_end)==length(common_start))
{
cat(paste0("All of these have consistent alleles\n\n"))
}
else
{
cat(paste0("Only ",length(common_end)," of these have consistent alleles\n\n"))
}


#get indexes, and compute Z and n

match_preds=match(common_end,gwas_all[,1])

if(!is.null(ZCol)){Z_stats=as.numeric(gwas_all[match_preds,4])}
else
{
effect_sizes=as.numeric(gwas_all[match_preds,4])
if(EffectCol=="OR"){effect_sizes=log(effect_sizes)}

if(!is.null(SECol)){Z_stats=effect_sizes/as.numeric(gwas_all[match_preds,5])}
else{Z_stats=sign(effect_sizes)*qnorm(as.numeric(gwas_all[match_preds,5])/2,lower=FALSE)}
}

if(!is.null(nCol)){sample_sizes=as.numeric(gwas_all[match_preds,ncol(gwas_all)])}
if(!is.null(ncasesCol)){sample_sizes=as.numeric(gwas_all[match_preds,ncol(gwas_all)-1])+as.numeric(gwas_all[match_preds,ncol(gwas_all)])}
if(!is.null(fixedn)){sample_sizes=rep(fixedn,length(common_end))}

cat(paste0("The median Z statistic is ", round(median(Z_stats),2)," (this should be close to zero), while ", round(100*mean(Z_stats>0),1),"% are positive (this should be close to 50%)\n\n"))

cat(paste0("The average sample size is ", round(mean(sample_sizes),1)," (the range is ", min(sample_sizes)," to ", max(sample_sizes),")\n\n"))


#save results
final_ss=cbind(gwas_all[match_preds,c(1,2,3)],Z_stats,sample_sizes)
colnames(final_ss)=c("Predictor","A1","A2","Z","n")

write.table(final_ss,outfile,row.names=F,col.names=TRUE,quote=FALSE)

cat(paste0("The formatted summary statistics are saved in the file ",outfile,"\n\n"))


################
#get total time
end_time=Sys.time()
cat(paste0("End at ",end_time,"\n"))
cat(paste0("Total run time was ", round(difftime(end_time,start_time,units="hours"),2), " hours\n\n"))
}

