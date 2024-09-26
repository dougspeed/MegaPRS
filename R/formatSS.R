#Copyright 2024 Doug Speed.

#    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

#    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.


################################################


#R code for putting GWAS results in LDAK format
#requires cors.bim file, so can work out consistent SNPs


################################################


#' Format GWAS results
#' @param gwasfile The name of the file containing GWAS results
#' @export
formatSS=function(gwasfile=NULL, outstem=NULL, NameCol=NULL, A1Col=NULL, A2Col=NULL, ZCol=NULL, EffectCol=NULL, SECol=NULL, PCol=NULL, nCol=NULL, ncasesCol=NULL, ncontrolsCol=NULL, fixedn=NULL, FreqCol=NULL, headerRows=NULL)
{
################
#get start time
start_time=Sys.time()

end_time=Sys.time()
cat(paste0("Start at ",start_time,"\n\n"))


################
#check arguments

if(is.null(gwasfile))
{return(paste0("Error, you must use the argument gwasfile to specify the file containing the GWAS results"))}

if(file.exists(gwasfile)==FALSE)
{
return(paste0("Error, unable to find the GWAS results file ",gwasfile))
}

if(is.null(outstem))
{return(paste0("Error, you must use the argument outstem to specify the stem for the output files"))}

if(!is.null(headerRows))
{
if(headerRows<1|headerRows!=round(headerRows))
{return(paste0("Error, the argument skipRows must be a positive integer (not ",headerRows,")"))}
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

cat(paste0("Row ", headerRows," starts ", substr(gwas_head[1],1,40)," - this will be treated as a comment and ignored\n"))
headerRows=headerRows+1
}
}

gwas_head=read.table(gwasfile,head=FALSE,nrow=1,comment.char="",skip=headerRows-1)
cat(paste0("Row ", headerRows," starts ", substr(gwas_head[1],1,40)," - this will be used as the header row\n\n"))
gwas_head=read.table(gwasfile,head=FALSE,nrow=1,comment.char="",skip=headerRows-1)

#first assume space delimited
comma_sep=0

if(length(gwas_head)==1)	#try comma delimited
{
comma_sep=1
gwas_head=read.table(gwasfile,head=FALSE,nrow=1,comment.char="",skip=headerRows-1,sep=",")

if(length(gwas_head)==1)
{return(paste0("Error, ", gwasfile," appears to have only one column (have tried to read using both spaces and commas as delimiters)"))}
}

if(length(gwas_head)<5)
{return(paste0("Error, ", gwasfile, "should have at least five columns (not ", length(gwas_head),")"))}

if(substr(gwas_head[1],1,1)=="#")	#remove # from start
{gwas_head[1]=substr(gwas_head[1],2,nchar(gwas_head[1]))}

cat(paste0("The file ", gwasfile," has the following ", length(gwas_head), " columns:\n"))
for(j in 1:length(gwas_head))
{
if(j==21){cat(paste0("Warning, only the first 20 columns are shown\n"))}

cat(paste0("Column ",j," is ", gwas_head[j],"\n"))
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
use A1Col to specify the column containing the A1 alleles (look for labels such as A1 or EffectAlle);\
use A2Col to specify the column containing the A2 alleles (look for labels such as A2 or OtherAllele);\
use EffectCol to specify the column containing effect sizes (look for labels such as Effect, Beta, OR or LogOR);\
use SECol to specify the column containing standard errors (look for labels such as SE or SD);\
and use nCol to specify the column containing sample sizes (look for labels such as n, N or nEff)\n\
We recommend you also use FreqCol to specify the column containing frequencies for the A1 alleles (look for labels such as Freq or FCON)\n\
Note 1: if SEs are not provided, you can use PCol to specify the column containing p-values, while if Z-test statistics are provided, you can use ZCol instead of EffectCol and SECol\n\
Note 2: if the sample size is divided into numbers of cases and controls, you can replace nCol with ncasesCol and ncontrolsCol (if sample size is not provided, you can instead use fixedn to provide the total sample size)\n\n"))
return(invisible())
}
if(count<6)
{
cat(paste0("Error, you must provide the following arguments:\
use NameCol to specify the column containing predictor names (look for labels such as Predictor, SNP, Marker or rsID);\
use A1Col to specify the column containing the A1 alleles (look for labels such as A1 or EffectAlle);\
use A2Col to specify the column containing the A2 alleles (look for labels such as A2 or OtherAllele);\
use EffectCol to specify the column containing effect sizes (look for labels such as Effect, Beta, OR or LogOR);\
use SECol to specify the column containing standard errors (look for labels such as SE or SD);\
and use nCol to specify the column containing sample sizes (look for labels such as n, N or nEff)\n\
We recommend you also use FreqCol to specify the column containing frequencies for the A1 alleles (look for labels such as Freq or FCON)\n\
Note 1: if SEs are not provided, you can use PCol to specify the column containing p-values, while if Z-test statistics are provided, you can use ZCol instead of EffectCol and SECol\n\
Note 2: if the sample size is divided into numbers of cases and controls, you can replace nCol with ncasesCol and ncontrolsCol (if sample size is not provided, you can instead use fixedn to provide the total sample size)\n\n"))
return(invisible())
}


################
#check column names are unique
all_cols=c(NameCol, A1Col, A2Col, ZCol, EffectCol, SECol, PCol, nCol, ncasesCol, ncontrolsCol, FreqCol)
for(one_col in all_cols)
{
if(sum(all_cols==one_col)>1)
{return(paste0("Error, the column name ", one_col, " is provided more than once"))}
}


################
#find column indexes

Name_find=which(gwas_head==NameCol)
if(length(Name_find)==0)
{return(paste0("Error, there is not a column called ", NameCol, " (make sure you provide column names and not numbers, and note that names are case-sensitive)"))}
if(length(Name_find)>1)
{return(paste0("Error, there is more than one column called ", NameCol))}

A1_find=which(gwas_head==A1Col)
if(length(A1_find)==0)
{return(paste0("Error, there is not a column called ", A1Col, " (make sure you provide column names and not numbers, and note that names are case-sensitive)"))}
if(length(Name_find)>1)
{return(paste0("Error, there is more than one column called ", A1Col))}

A2_find=which(gwas_head==A2Col)
if(length(A2_find)==0)
{return(paste0("Error, there is not a column called ", A2Col, " (make sure you provide column names and not numbers, and note that names are case-sensitive)"))}
if(length(Name_find)>1)
{return(paste0("Error, there is more than one column called ", A2Col))}

if(!is.null(ZCol))
{
Z_find=which(gwas_head==ZCol)
if(length(Z_find)==0)
{return(paste0("Error, there is not a column called ", ZCol, " (make sure you provide column names and not numbers, and note that names are case-sensitive)"))}
if(length(Name_find)>1)
{return(paste0("Error, there is more than one column called ", ZCol))}
}
else
{
Effect_find=which(gwas_head==EffectCol)
if(length(Effect_find)==0)
{return(paste0("Error, there is not a column called ", EffectCol, " (make sure you provide column names and not numbers, and note that names are case-sensitive)"))}
if(length(Name_find)>1)
{return(paste0("Error, there is more than one column called ", EffectCol))}

if(!is.null(SECol))
{
SE_find=which(gwas_head==SECol)
if(length(SE_find)==0)
{return(paste0("Error, there is not a column called ", SECol, " (make sure you provide column names and not numbers, and note that names are case-sensitive)"))}
if(length(Name_find)>1)
{return(paste0("Error, there is more than one column called ", SECol))}
}
else
{
P_find=which(gwas_head==PCol)
if(length(P_find)==0)
{return(paste0("Error, there is not a column called ", PCol, " (make sure you provide column names and not numbers, and note that names are case-sensitive)"))}
if(length(Name_find)>1)
{return(paste0("Error, there is more than one column called ", PCol))}
}
}

if(!is.null(nCol))
{
n_find=which(gwas_head==nCol)
if(length(n_find)==0)
{return(paste0("Error, there is not a column called ", nCol, " (make sure you provide column names and not numbers, and note that names are case-sensitive)"))}
if(length(Name_find)>1)
{return(paste0("Error, there is more than one column called ", nCol))}
}
if(!is.null(ncasesCol))
{
ncases_find=which(gwas_head==ncasesCol)
if(length(ncases_find)==0)
{return(paste0("Error, there is not a column called ", ncasesCol, " (make sure you provide column names and not numbers, and note that names are case-sensitive)"))}
if(length(Name_find)>1)
{return(paste0("Error, there is more than one column called ", ncasesCol))}

ncontrols_find=which(gwas_head==ncontrolsCol)
if(length(ncontrols_find)==0)
{return(paste0("Error, there is not a column called ", ncontrolsCol, " (make sure you provide column names and not numbers, and note that names are case-sensitive)"))}
if(length(Name_find)>1)
{return(paste0("Error, there is more than one column called ", ncontrolsCol))}
}

if(!is.null(FreqCol))
{
freq_find=which(gwas_head==FreqCol)
if(length(freq_find)==0)
{return(paste0("Error, there is not a column called ", FreqCol, " (make sure you provide column names and not numbers, and note that names are case-sensitive)"))}
if(length(Name_find)>1)
{return(paste0("Error, there is more than one column called ", FreqCol))}
}

cat(paste0("Will read predictor names from the column called ", gwas_head[Name_find],"\n"))
cat(paste0("Will read A1 alleles from the column called ", gwas_head[A1_find],"\n"))
cat(paste0("Will read A2 alleles from the column called ", gwas_head[A2_find],"\n"))

if(!is.null(ZCol))
{cat(paste0("Will read Z-test statistics from the column called ", gwas_head[Z_find],"\n"))}
else
{
cat(paste0("Will read effect sizes from the column called ", gwas_head[Effect_find],"\n"))
if(EffectCol=="OR")
{cat(paste0("(it is assumed this column provides Odds Ratios, so will compute its logarithm\n)"))}
if(!is.null(SECol))
{cat(paste0("Will read SEs from the column called ", gwas_head[SE_find],"\n"))}
else
{cat(paste0("Will read p-values from the column called ", gwas_head[P_find],"\n"))}
}

if(!is.null(nCol))
{cat(paste0("Will read sample sizes from the column called ", gwas_head[n_find],"\n"))}
if(!is.null(ncasesCol))
{cat(paste0("Will read numbers of cases and controls from the columns called ", gwas_head[ncases_find]," and ",gwas_head[ncases_find],"\n"))}
if(!is.null(fixedn))
{cat(paste0("The sample size will be set to ", fixedn,"\n"))}

if(!is.null(FreqCol)){cat(paste0("Will read A1 allele frequencies from the column called ", gwas_head[freq_find],"\n\n"))}
else{cat(paste0("\nNote that if ", gwasfile," provides frequencies for the A1 alleles, we recommend you use the argument FreqCol to specify the corresponding column\n\n"))}


################
#set numoff (where the sample sizes start) and  gotfreq (whether we have allele frequencies)
if(!is.null(ZCol)){numoff=5}
else{numoff=6}
if(is.null(FreqCol)){gotfreq=1}
else{gotfreq=0}


################
#load included data files

data(GENO.SNPs)
data(HAPMAP.SNPs)
data(PCA.DETAILS)


################
#read details from gwasfile

req_cols=c(Name_find,A1_find,A2_find)
req_classes=c("character","character","character")

if(!is.null(ZCol)){req_cols=c(req_cols,Z_find);req_classes=c(req_classes,"numeric")}
else
{
if(!is.null(SECol)){req_cols=c(req_cols,Effect_find,SE_find);req_classes=c(req_classes,"numeric","numeric")}
else{req_cols=c(req_cols,Effect_find,P_find);req_classes=c(req_classes,"numeric","numeric")}
}

if(!is.null(nCol)){req_cols=c(req_cols,n_find);req_classes=c(req_classes,"numeric")}
if(!is.null(ncasesCol)){req_cols=c(req_cols,ncases_find,ncontrols_find);req_classes=c(req_classes,"numeric","numeric")}

if(!is.null(FreqCol)){req_cols=c(req_cols,freq_find);req_classes=c(req_classes,"numeric")}

classes_all=rep("NULL",length(gwas_head))
classes_all[req_cols]=req_classes

back_cols=match(req_cols,sort(req_cols))

cat(paste0("Reading GWAS results from ",gwasfile,"\n"))

file_size=file.info(gwasfile)$size/2^20
if(file_size<100){cat(paste0("The file is relatively small (", round(file_size),"Mb), so this should only take seconds\n\n"))}
else{cat(paste0("The file is quite large (", round(file_size),"Mb), so this can take a few minutes\n\n"))}

if(comma_sep==0){gwas_all=read.table(gwasfile,head=FALSE,colClasses=classes_all,comment.char="",skip=headerRows,fill=TRUE)[,back_cols]}
else{gwas_all=read.table(gwasfile,head=FALSE,colClasses=classes_all,comment.char="",skip=headerRows,sep=",",fill=TRUE)[,back_cols]}
colnames(gwas_all)=gwas_head[req_cols]

cat(paste0("In total, there are results for ", nrow(gwas_all)," predictors\n\n"))


################
#check for different alleles

diff_alleles=which(gwas_all[,2]!=gwas_all[,3])
if(length(diff_alleles)==0)
{
return(paste0("Error, none of the predictors have distinct alleles (e.g., the first SNP is called ", gwas_all[1,1], ", and both its alleles are ", gwas_all[1,2],")"))
}
if(length(diff_alleles)<nrow(gwas_all))
{
cat(paste0("Warning, only ",length(diff_alleles), " of these predictors have distinct alleles (the remaining ", nrow(gwas_all)-length(diff_alleles)," will be ignored)\n\n"))
}


################
#print out summary statistics for all valid predictors

if(length(diff_alleles)==nrow(gwas_all))
{cat(paste0("Formatting summary statistics for all ", nrow(gwas_all), " predictors\n"))}
else
{cat(paste0("Formatting summary statistics for the ", length(diff_alleles), " predictors with distinct alleles\n"))}

if(!is.null(ZCol)){Z_stats=as.numeric(gwas_all[diff_alleles,4])}
else
{
effect_sizes=as.numeric(gwas_all[diff_alleles,4])
if(EffectCol=="OR"){effect_sizes=log(effect_sizes)}

if(!is.null(SECol))
{
SEs=as.numeric(gwas_all[diff_alleles,5])
if(sum(SEs<0)>0){return(paste0("Error, ", sum(SEs<0), " predictors have negative SEs"))}
if(sum(SEs==0)>0){return(paste0("Error, ", sum(SEs<0), " predictors have SE zero"))}

Z_stats=effect_sizes/SEs
}
else
{
pvalues=as.numeric(gwas_all[diff_alleles,5])
if(sum(pvalues<0)>0){return(paste0("Error, ", sum(pvalues<0), " predictors have p-values below zero"))}
if(sum(pvalues>1)>0){return(paste0("Error, ", sum(pvalues>1), " predictors have p-values above one"))}

small_pvalues=which(pvalues<1e-300)
if(length(small_pvalues)>0)
{
cat(paste0("Warning, ", length(small_pvalues), " predictors have p-values below 1e-300 (the p-values will be rounded up to 1e-300)\n"))
pvalues[small_pvalues]=1e-300
}

Z_stats=sign(effect_sizes)*qnorm(pvalues/2,lower=FALSE)
}
}

if(!is.null(nCol)){sample_sizes=as.numeric(gwas_all[diff_alleles,numoff])}
if(!is.null(ncasesCol)){sample_sizes=as.numeric(gwas_all[diff_alleles,numoff])+as.numeric(gwas_all[diff_alleles,numoff+1])}
if(!is.null(fixedn)){sample_sizes=rep(fixedn,length(common_hapmap_end))}

if(!is.null(FreqCol))
{
a1_freq=as.numeric(gwas_all[diff_alleles,ncol(gwas_all)])
if(sum(a1_freq<0)>0){return(paste0("Error, ", sum(a1_freq<0), " predictors have frequency below zero"))}
if(sum(a1_freq>1)>0){return(paste0("Error, ", sum(a1_freq>1), " predictors have frequency above one"))}
}

#which predictors have valid resutls

if(is.null(FreqCol)){valid_preds=which(!is.na(Z_stats)&!is.na(sample_sizes))}
else{valid_preds=which(!is.na(Z_stats)&!is.na(sample_sizes)&!is.na(a1_freq))}

if(length(valid_preds)==0)
{return(paste0("Error, none of the ", length(diff_alleles)," predictors have valid results"))}
if(length(valid_preds)<length(diff_alleles))
{cat(paste0("Warning, only ", length(valid_preds), " of the ", length(diff_alleles)," predictors have valid results\n\n"))}

#print out some summaries

cat(paste0("The median Z statistic is ", round(median(Z_stats[valid_preds]),2)," (this should be close to zero), while ", round(100*mean(Z_stats[valid_preds]>0),1),"% are positive (this should be close to 50%)\n"))

cat(paste0("The average sample size is ", round(mean(sample_sizes[valid_preds]),1)," (the range is ", min(sample_sizes[valid_preds])," to ", max(sample_sizes[valid_preds]),")\n"))

#save results
if(is.null(FreqCol))
{
final_ss=cbind(gwas_all[,c(1,2,3)],Z_stats,sample_sizes)[valid_preds,]
colnames(final_ss)=c("Predictor","A1","A2","Z","n")
}
else
{
final_ss=cbind(gwas_all[valid_preds,c(1,2,3)],Z_stats,sample_sizes,a1_freq)[valid_preds,]
colnames(final_ss)=c("Predictor","A1","A2","Z","n","A1Freq")
}

outfile=paste0(outstem,".ALL.summaries")
write.table(final_ss,outfile,row.names=F,col.names=TRUE,quote=FALSE)
cat(paste0("The formatted summary statistics are saved in the file ",outfile,"\n\n"))


################
#find overlap with genotyped SNPs

common_geno_start=intersect(final_ss[,1],GENO.SNPs[,1])
if(length(common_geno_start)==0)
{
return(paste0("Error, there is no overlap with the ", nrow(GENO.SNPs), " genotyped SNPs"))
}
cat(paste0("Have found ", length(common_geno_start)," of the ", nrow(GENO.SNPs), " genotyped SNPs\n"))


################
#check alleles

match1=match(common_geno_start,final_ss[,1])
match2=match(common_geno_start,GENO.SNPs[,1])
find_a1= (final_ss[match1,2]==GENO.SNPs[match2,2]) + (final_ss[match1,2]==GENO.SNPs[match2,3])
find_a2= (final_ss[match1,3]==GENO.SNPs[match2,2]) + (final_ss[match1,3]==GENO.SNPs[match2,3])
common_geno_end=common_geno_start[which(find_a1==1&find_a2==1)]

if(length(common_geno_end)==0)
{
pp=final_ss[match1[1],1]
a1=final_ss[match1[1],2]
a2=final_ss[match1[1],3]
b1=GENO.SNPs[match2[1],2]
b2=GENO.SNPs[match2[1],3]
return(paste0("Error, none of the ", length(common_geno_start), " genotyped SNPs have consistent alleles (e.g., ",pp ," should have alleles ", b1," and ", b2, ", but the alleles in ", gwasfile," are ", a1," and ", a2,")"))
}

if(length(common_geno_end)<length(common_geno_start))
{
cat(paste0("Warning, after checking alleles, the number of genotyped SNPs has been reduced to ", length(common_geno_end),"\n"))
}


################
#print out summary statistics

match_preds=match(common_geno_end,final_ss[,1])
outfile=paste0(outstem,".GENO.summaries")
write.table(final_ss[match_preds,],outfile,row.names=F,col.names=TRUE,quote=FALSE)
cat(paste0("The corresponding summary statistics are saved in the file ",outfile,"\n\n"))


################
#find overlap with HapMap SNPs

common_hapmap_start=intersect(final_ss[,1],HAPMAP.SNPs[,1])
if(length(common_hapmap_start)==0)
{
return(paste0("Error, there is no overlap with the ", nrow(HAPMAP.SNPs), " HapMap SNPs"))
}
cat(paste0("Have found ", length(common_hapmap_start)," of the ", nrow(HAPMAP.SNPs), " HapMap SNPs\n"))


################
#check alleles

match1=match(common_hapmap_start,final_ss[,1])
match2=match(common_hapmap_start,HAPMAP.SNPs[,1])
find_a1= (final_ss[match1,2]==HAPMAP.SNPs[match2,2]) + (final_ss[match1,2]==HAPMAP.SNPs[match2,3])
find_a2= (final_ss[match1,3]==HAPMAP.SNPs[match2,2]) + (final_ss[match1,3]==HAPMAP.SNPs[match2,3])
common_hapmap_end=common_hapmap_start[which(find_a1==1&find_a2==1)]

if(length(common_hapmap_end)==0)
{
pp=final_ss[match1[1],1]
a1=final_ss[match1[1],2]
a2=final_ss[match1[1],3]
b1=HAPMAP.SNPs[match2[1],2]
b2=HAPMAP.SNPs[match2[1],3]
return(paste0("Error, none of the ", length(common_hapmap_start), " HapMap SNPs have consistent alleles (e.g., ",pp ," should have alleles ", b1," and ", b2, ", but the alleles in ", gwasfile," are ", a1," and ", a2,")"))
}

if(length(common_hapmap_end)<length(common_hapmap_start))
{
cat(paste0("Warning, after checking alleles, the number of HapMap SNPs has been reduced to ", length(common_hapmap_end),"\n"))
}


################
#print out summary statistics

match_preds=match(common_hapmap_end,final_ss[,1])
outfile=paste0(outstem,".HAPMAP.summaries")
write.table(final_ss[match_preds,],outfile,row.names=F,col.names=TRUE,quote=FALSE)
cat(paste0("The corresponding summary statistics are saved in the file ",outfile,"\n\n"))


################
#give some advice
if(length(common_geno_end)/nrow(GENO.SNPs)>0.5)
{cat(paste0("When constructing PRS, we recommend that you use the summary statistics in ", outstem,".GENO.summaries\n\n"))}
else
{cat(paste0("There are relatively few genotyped SNPs, so consider whether it is better to construct PRS using the summary statistics in ", outstem,".HAPMAP.summaries\n\n"))}


################
#test ancestry (if possible)

if(!is.null(FreqCol))
{
cat(paste0("Assessing the ancestry of the GWAS summary statistics\n"))

#get ancestry snps and their indexes
common_pca=intersect(common_geno_end,PCA.DETAILS[,1])
match1=match(common_pca,gwas_all[,1])
match2=match(common_pca,PCA.DETAILS[,1])

if(length(common_pca)<10000)
{return(paste0("Error, unable to computer ancestry PCs because ", gwasfile, " contains only ", length(common_pca), " of the ", nrow(PCA.DETAILS), " ancestry SNPs"))}

#get frequencies - seeing if necessary to flip
a1_freq=as.numeric(gwas_all[match1,ncol(gwas_all)])
flip_preds=which(gwas_all[match1,2]!=PCA.DETAILS[match2,2])
if(length(flip_preds)>0){a1_freq[flip_preds]=1-a1_freq[flip_preds]}

#compute mean pcs

cents=as.numeric(PCA.DETAILS[match2,4])
w1=as.numeric(PCA.DETAILS[match2,5])
w2=as.numeric(PCA.DETAILS[match2,6])
w3=as.numeric(PCA.DETAILS[match2,7])

proj1=matrix(NA,length(common_pca),ncol(PCA.DETAILS)-6)
proj2=matrix(NA,length(common_pca),ncol(PCA.DETAILS)-6)
proj3=matrix(NA,length(common_pca),ncol(PCA.DETAILS)-6)

proj1[,1]=(2*a1_freq-cents)*w1
proj2[,1]=(2*a1_freq-cents)*w2
proj3[,1]=(2*a1_freq-cents)*w3

for(j in 8:ncol(PCA.DETAILS))
{
proj1[,j-6]=(as.numeric(PCA.DETAILS[match2,j])-cents)*w1
proj2[,j-6]=(as.numeric(PCA.DETAILS[match2,j])-cents)*w2
proj3[,j-6]=(as.numeric(PCA.DETAILS[match2,j])-cents)*w3
}

means1=apply(proj1,2,mean)
means2=apply(proj2,2,mean)
means3=apply(proj3,2,mean)

#get distances between the summary statistics and the five reference panels
distances=rep(NA,5)
for(j in 1:5){distances[j]=(means1[1]-means1[1+j])^2+(means2[1]-means2[1+j])^2+(means3[1]-means3[1+j])^2}

closest_ref=which.min(distances)
if(min(distances)<0.01)
{cat(paste0("The closest reference panel is ", colnames(PCA.DETAILS)[7+closest_ref],"; the distance is ", signif(min(distances),2)," indicating a very good match\n"))}
else
{cat(paste0("The closest reference panel is ", colnames(PCA.DETAILS)[7+closest_ref],"; the distance is ", signif(min(distances),2)," indicating the panel is sub-optimal (ideally the distance should be below 0.01)\n"))}

#visualize the results

cols=rep(NA,length(means1))
cols[1]="orange"
mark=1
for(tt in c("AFR","AMR","EAS","EUR","SAS","FIN"))
{
a=grep(tt,colnames(PCA.DETAILS))-6
cols[a]=mark
mark=mark+1
}

pchs=rep(4,length(means1))
pchs[1]=3
pchs[2:6]=1

cexs=rep(1.3,length(means1))
cexs[1]=5
cexs[2:6]=3

outfile=paste0(outstem,".ancestries.pdf")
pdf(outfile,height=6,width=8)
plot_title=paste0("Closest Ref Panel: ",colnames(PCA.DETAILS)[7+closest_ref]," (distance ",signif(min(distances),2),")")
plot(means1,means2,col=cols,pch=pchs,lwd=3,cex=cexs,xlab="Principal Component Axis 1",ylab="Principal Component Axis 2",main=plot_title)
legend("bottomright",leg=c("AFR","AMR","EAS","EUR","SAS","FIN"),title="Ancestries",fill=1:6,bty="n",cex=1.5,ncol=2)
legend("topright",leg=c("GWAS Summary Stats","MegaPRS Ref Panels","1000GP Populations"),col=c("orange","darkgrey","darkgrey"),cex=1.5,pch=c(3,1,4),bty="n",lwd=3,lty=NA)
dev.off()

plot_title=paste0("Closest Ref Panel: ",colnames(PCA.DETAILS)[7+closest_ref]," (distance ",signif(min(distances),2),")")
plot(means1,means2,col=cols,pch=pchs,lwd=3,cex=cexs,xlab="Principal Component Axis 1",ylab="Principal Component Axis 2",main=plot_title)
legend("bottomright",leg=c("AFR","AMR","EAS","EUR","SAS","FIN"),title="Ancestries",fill=1:6,bty="n",cex=1.5,ncol=2)
legend("topright",leg=c("GWAS Summary Stats","MegaPRS Ref Panels","1000GP Populations"),col=c("orange","darkgrey","darkgrey"),cex=1.5,pch=c(3,1,4),bty="n",lwd=3,lty=NA)

cat(paste0("The PCA plot has been saved in ", outfile,"\n\n"))
}


################
#get total time
end_time=Sys.time()
cat(paste0("End at ",end_time,"\n"))
cat(paste0("Total run time was ", round(difftime(end_time,start_time,units="hours"),2), " hours\n\n"))
}

