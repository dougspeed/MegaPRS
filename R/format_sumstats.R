#Copyright 2024 Doug Speed.

#    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

#    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.


################################################


#R code for putting GWAS results in LDAK format
#requires cors.bim file, so can work out consistent SNPs


################################################


#' @title Tool for formatting GWAS results
#
#' @description This function converts GWAS results into the format required by LDAK, infers the ancestry of the GWAS and determines a suitable reference panel
#
#' @details We recommend you first run the function using only the arguments gwasfile and outstem, then follow the on-screen instructions
#
#' @param gwasfile The name of the file containing GWAS results
#' @param outstem The desired prefix for the output files
#' @export
#' @examples
#' #These examples use the gwas results that come with the MegaPRS package
#' ex_gwas_file=system.file("extdata", "ex_gwas.txt.gz", package="MegaPRS")
#' 
#' #Example 1 - Provide only the two required arguments
#' format_sumstats(gwasfile=ex_gwas_file, outstem="ex_out")
#' 
#' #Example 2 - Run the function based on the on-screen instructions from Example 1
#' format_sumstats(gwasfile=ex_gwas_file, outstem="ex_out", NameCol="Predictor", A1Col="A1", A2Col="A2", EffectCol="Beta", SECol="SE", nCol="n", FreqCol="A1Freq")
format_sumstats=function(gwasfile=NULL, outstem=NULL, NameCol=NULL, A1Col=NULL, A2Col=NULL, EffectCol=NULL, SECol=NULL, PCol=NULL, ZCol=NULL, nCol=NULL, ncasesCol=NULL, ncontrolsCol=NULL, fixedn=NULL, FreqCol=NULL, headerRows=NULL)
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

gwas_head=read.table(gwasfile,head=FALSE,nrow=1,comment.char="",skip=headerRows-1,sep=";")
cat(paste0("Row ", headerRows," starts ", substr(gwas_head[1],1,40)," - this will be used as the header row\n\n"))

#by default, assume gwasfile is tab delimited
space_sep=0
comma_sep=0

gwas_head=read.table(gwasfile,head=FALSE,nrow=1,comment.char="",skip=headerRows-1,sep="\t")

if(length(gwas_head)==1)	#try space delimited
{
gwas_head=read.table(gwasfile,head=FALSE,nrow=1,comment.char="",skip=headerRows-1,sep=" ")
if(length(gwas_head)>1){space_sep=1}
}

if(length(gwas_head)==1)	#try comma delimited
{
gwas_head=read.table(gwasfile,head=FALSE,nrow=1,comment.char="",skip=headerRows-1,sep=",")
if(length(gwas_head)>1){comma_sep=1}
}

if(length(gwas_head)==1)
{return(paste0("Error, ", gwasfile," appears to have only one column (have tried to read using tabs, spaces and commas as delimiters)"))}

if(space_sep==0&comma_sep==0){cat(paste0(gwasfile," appears to be tab delimited\n\n"))}
if(space_sep==1){cat(paste0(gwasfile," appears to be space delimited\n\n"))}
if(comma_sep==1){cat(paste0(gwasfile," appears to be comma delimited\n\n"))}

if(length(gwas_head)<5)
{return(paste0("Error, ", gwasfile, "should have at least five columns (not ", length(gwas_head),")"))}

if(substr(gwas_head[1],1,1)=="#")	#remove # from start
{gwas_head[1]=substr(gwas_head[1],2,nchar(gwas_head[1]))}

cat(paste0(gwasfile," has ", length(gwas_head), " columns, labelled as follows:\n"))
for(j in 1:length(gwas_head))
{
if(j==31){cat(paste0("Warning, only the first 30 columns are shown\n"))}

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
use NameCol to specify the predictor names (look for column labels such as Predictor, SNP, Marker or rsID)\
use A1Col to specify the A1 alleles (look for column labels such as A1 or EffectAlle)\
use A2Col to specify the A2 alleles (look for column labels such as A2 or OtherAllele)\
use EffectCol to specify the effect sizes (look for column labels such as Effect, Beta, OR or LogOR)\
use SECol to specify the standard errors (look for column labels such as SE or SD)\
and use nCol to specify the sample sizes (look for column labels such as n, N or nEff)\n\
In addition, we recommend you use FreqCol to specify the A1 allele frequencies (look for labels such as Freq or FCON)\n\
For example, if the predictor names are stored in the column labelled ", gwas_head[1], ", you should add NameCol=\"", gwas_head[1],"\"\n\
Note 1: if p-values are provided, you can use PCol instead of SECol, while if Z-test statistics are provided, you can use ZCol instead of EffectCol and SECol\n\
Note 2: if the sample size is divided into numbers of cases and controls, you can replace nCol with ncasesCol and ncontrolsCol (if sample size is not provided, you can instead use fixedn to provide the total sample size)\n\n"))
return(invisible())
}
if(count<6)
{
cat(paste0("Error, there are insufficient column arguments; please ensure you include the following:\
use NameCol to specify the predictor names (look for column labels such as Predictor, SNP, Marker or rsID)\
use A1Col to specify the A1 alleles (look for column labels such as A1 or EffectAlle)\
use A2Col to specify the A2 alleles (look for column labels such as A2 or OtherAllele)\
use EffectCol to specify the effect sizes (look for column labels such as Effect, Beta, OR or LogOR)\
use SECol to specify the standard errors (look for column labels such as SE or SD)\
and use nCol to specify the sample sizes (look for column labels such as n, N or nEff)\n\
In addition, we recommend you use FreqCol to specify the A1 allele frequencies (look for labels such as Freq or FCON)\n\
For example, if the predictor names are stored in the column labelled ", gwas_head[1], ", you should add NameCol=\"", gwas_head[1],"\"\n\
Note 1: if p-values are provided, you can use PCol instead of SECol, while if Z-test statistics are provided, you can use ZCol instead of EffectCol and SECol\n\
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


################
#set num_offset (where the sample sizes start), got_freq (whether we have allele frequencies), and got_OR (whether we have OR)

if(!is.null(ZCol)){num_offset=5}
else{num_offset=6}
if(is.null(FreqCol)){got_freq=1}
else{got_freq=0}
if(!is.null(EffectCol))
{
if(EffectCol=="OR"|EffectCol=="Odds"|EffectCol=="ODDS"){got_OR=1}
else{got_OR=0}
}
else{got_OR=0}


################
#explain what will happen

cat(paste0("Will read predictor names from the column called ", gwas_head[Name_find],"\n"))
cat(paste0("Will read A1 alleles from the column called ", gwas_head[A1_find],"\n"))
cat(paste0("Will read A2 alleles from the column called ", gwas_head[A2_find],"\n"))

if(!is.null(ZCol))
{cat(paste0("Will read Z-test statistics from the column called ", gwas_head[Z_find],"\n"))}
else
{
cat(paste0("Will read effect sizes from the column called ", gwas_head[Effect_find]))
if(got_OR){cat(paste0(" (it is assumed this column provides Odds Ratios, so will compute its logarithm)"))}
cat("\n")
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
else{cat(paste0("\nNote that if ", gwasfile," provides A1 allele frequencies, we recommend you use the argument FreqCol to specify the corresponding column\n\n"))}


################
#load included data files

data(GENO.SNPs)
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
else{cat(paste0("The file is relatively large (", round(file_size),"Mb), so this might take a few minutes\n\n"))}

if(space_sep==0&comma_sep==0)
{gwas_all=read.table(gwasfile,colClasses=classes_all,comment.char="",skip=headerRows,fill=TRUE,sep="\t")[,back_cols]}
if(space_sep==1)
{gwas_all=read.table(gwasfile,colClasses=classes_all,comment.char="",skip=headerRows,fill=TRUE,sep=" ")[,back_cols]}
if(comma_sep==1)
{gwas_all=read.table(gwasfile,colClasses=classes_all,comment.char="",skip=headerRows,fill=TRUE,sep=",")[,back_cols]}

colnames(gwas_all)=gwas_head[sort(req_cols)]

cat(paste0("In total, there are results for ", nrow(gwas_all)," predictors (here are the first six rows)\n\n"))
print(head(gwas_all))
cat("\n")


################
#see if any predictors are missing names

missing_names=which(gwas_all[,1]=="")
if(length(missing_names)==nrow(gwas_all))
{
return(paste0("Error, all predictor names are missing; please check the argument NameCol"))
}
if(length(missing_names)>0)
{
cat(paste0("Warning, ", length(missing_names)," predictors are missing names, and will be ignored\n\n"))
}


################
#check for different alleles

diff_alleles=which(gwas_all[,2]!=gwas_all[,3])
if(length(diff_alleles)==0)
{
return(paste0("Error, none of the predictors have distinct alleles (e.g., the first predictor is called ", gwas_all[got_names[1],1], ", and both its alleles are ", gwas_all[got_names[1],2],"); please check the arguments A1Col and A2Col"))
}
if(length(diff_alleles)<nrow(gwas_all))
{
cat(paste0("Warning, only ",length(diff_alleles), " of the predictors have distinct alleles (the remaining ", nrow(gwas_all)-length(diff_alleles)," will be ignored)\n\n"))
}


################
#print out summary statistics for all valid predictors

start_preds=setdiff(diff_alleles,missing_names)

if(!is.null(ZCol)){Z_stats=as.numeric(gwas_all[start_preds,4])}
else
{
effect_sizes=as.numeric(gwas_all[start_preds,4])
if(got_OR){effect_sizes=log(effect_sizes)}

if(!is.null(SECol))
{
SEs=as.numeric(gwas_all[start_preds,5])
if(sum(SEs<0,na.rm=TRUE)>0){return(paste0("Error, ", sum(SEs<0,na.rm=TRUE), " predictors have negative SEs"))}
if(sum(SEs==0,na.rm=TRUE)>0){return(paste0("Error, ", sum(SEs==0,na.rm=TRUE), " predictors have SE zero"))}

Z_stats=effect_sizes/SEs
}
else
{
pvalues=as.numeric(gwas_all[start_preds,5])
if(sum(pvalues<0,na.rm=TRUE)>0){return(paste0("Error, ", sum(pvalues<0,na.rm=TRUE), " predictors have p-values below zero"))}
if(sum(pvalues>1,na.rm=TRUE)>0){return(paste0("Error, ", sum(pvalues>1,na.rm=TRUE), " predictors have p-values above one"))}

small_pvalues=which(pvalues<1e-300)
if(length(small_pvalues)>0)
{
cat(paste0("Warning, ", length(small_pvalues), " predictors have p-values below 1e-300 (the p-values will be rounded up to 1e-300)\n"))
pvalues[small_pvalues]=1e-300
}

Z_stats=sign(effect_sizes)*qnorm(pvalues/2,lower=FALSE)
}
}

if(!is.null(nCol)){sample_sizes=as.numeric(gwas_all[start_preds,num_offset])}
if(!is.null(ncasesCol)){sample_sizes=as.numeric(gwas_all[start_preds,num_offset])+as.numeric(gwas_all[start_preds,num_offset+1])}
if(!is.null(fixedn)){sample_sizes=rep(fixedn,length(start_preds))}

if(!is.null(FreqCol))
{
a1_freq=as.numeric(gwas_all[start_preds,ncol(gwas_all)])
if(sum(a1_freq<0,na.rm=TRUE)>0){return(paste0("Error, ", sum(a1_freq<0,na.rm=TRUE), " predictors have frequency below zero"))}
if(sum(a1_freq>1,na.rm=TRUE)>0){return(paste0("Error, ", sum(a1_freq>1,na.rm=TRUE), " predictors have frequency above one"))}
}

#which predictors have valid resutls

if(is.null(FreqCol)){valid_preds=which(!is.na(Z_stats)&!is.na(sample_sizes))}
else{valid_preds=which(!is.na(Z_stats)&!is.na(sample_sizes)&!is.na(a1_freq))}

if(length(valid_preds)==0)
{return(paste0("Error, none of the ", nrow(gwas_all)," predictors have valid results"))}
if(length(valid_preds)<nrow(gwas_all))
{cat(paste0("Warning, only ", length(valid_preds), " of the ", nrow(gwas_all)," predictors have valid results\n\n"))}

#print out some summaries

cat(paste0("The Z statistic have medium ", round(median(Z_stats[valid_preds]),4)," (this should be close to zero), while ", round(100*mean(Z_stats[valid_preds]>0),2),"% are positive (this should be close to 50%)\n"))

cat(paste0("The average sample size is ", round(mean(sample_sizes[valid_preds]),1)," (the range is ", min(sample_sizes[valid_preds])," to ", max(sample_sizes[valid_preds]),")\n"))

#save results

if(is.null(FreqCol))
{
final_ss=cbind(gwas_all[start_preds,c(1,2,3)],Z_stats,sample_sizes)[valid_preds,]
colnames(final_ss)=c("Predictor","A1","A2","Z","n")
}
else
{
final_ss=cbind(gwas_all[start_preds,c(1,2,3)],Z_stats,sample_sizes,a1_freq)[valid_preds,]
colnames(final_ss)=c("Predictor","A1","A2","Z","n","A1Freq")
}

outfile=paste0(outstem,".ALL.summaries")
write.table(final_ss,outfile,row.names=F,col.names=TRUE,quote=FALSE)
cat(paste0("The formatted summary statistics are saved in the file ",outfile,"\n\n"))


################
#find overlap with genotyped SNPs - try to match based on name, position, and genetic

cat(paste0("Searching for overlap with the ", nrow(GENO.SNPs), " genotyped SNPs (here are the first six of these)\n"))
print(head(GENO.SNPs))
cat("\n")

common_name=intersect(final_ss[,1],GENO.SNPs[,1])
common_position=intersect(final_ss[,1],GENO.SNPs[,4])

generic_names=c(paste0(GENO.SNPs[,4],"_",GENO.SNPs[,2],"_",GENO.SNPs[,3]),paste0(GENO.SNPs[,4],"_",GENO.SNPs[,3],"_",GENO.SNPs[,2]))
common_generic=intersect(final_ss[,1],generic_names)

num_overlap=max(c(length(common_name),length(common_position),length(common_generic)))
if(num_overlap)
{
return(paste0("Error, there are no genotyped SNPs (have searched based on SNP names and positions)"))
}
cat(paste0("Have found ", num_overlap," of the ", nrow(GENO.SNPs), " genotyped SNPs\n"))


################
#find indexes of overlapping SNPs

best_search=which.max(c(length(common_name),length(common_position),length(common_generic)))

if(best_search==1)	#matching based on name
{
match1=match(common_name,final_ss[,1])
match2=match(common_name,GENO.SNPs[,1])
}
if(best_search==2)	#matching based on position
{
match1=match(common_position,final_ss[,1])
match2=match(common_position,GENO.SNPs[,4])
}
if(best_search==3)	#matching based on generic format
{
match1=match(common_generic,final_ss[,1])
match2=(match(common_position,generic_names)-1)%%nrow(GENO.SNPs)+1
}


################
#check alleles

find_a1= (final_ss[match1,2]==GENO.SNPs[match2,2]) + (final_ss[match1,2]==GENO.SNPs[match2,3])
find_a2= (final_ss[match1,3]==GENO.SNPs[match2,2]) + (final_ss[match1,3]==GENO.SNPs[match2,3])
match_alleles=which(find_a1==1&find_a2==1)

if(length(match_alleles)==0)
{
pp=final_ss[match1[1],1]
a1=final_ss[match1[1],2]
a2=final_ss[match1[1],3]
b1=GENO.SNPs[match2[1],2]
b2=GENO.SNPs[match2[1],3]
return(paste0("Error, none of the ", num_found, " genotyped SNPs have consistent alleles (e.g., ",pp ," should have alleles ", b1," and ", b2, ", but the alleles in ", gwasfile," are ", a1," and ", a2,")"))
}

if(length(match_alleles)<num_overlap)
{
cat(paste0("Warning, after checking alleles, the number of genotyped SNPs has been reduced to ", length(match_alleles),"\n"))
}


################
#print out summary statistics for overlapping genotyped SNPs

subset_ss=cbind(GENO.SNPs[match2,1],final_ss[match1,-1])[match_alleles,]
outfile=paste0(outstem,".GENO.summaries")
write.table(subset_ss,outfile,row.names=F,col.names=TRUE,quote=FALSE)
cat(paste0("The corresponding summary statistics are saved in the file ",outfile,"\n\n"))


################
#test ancestry (if possible)

if(!is.null(FreqCol))
{
cat(paste0("Assessing the ancestry of the GWAS summary statistics\n"))

#get ancestry snps and their indexes, then gwas frequencies (seeing if necessary to flip)

common_pca=intersect(subset_ss[,1],PCA.DETAILS[,1])
match1=match(common_pca,subset_ss[,1])
match2=match(common_pca,PCA.DETAILS[,1])

if(length(common_pca)<10000)
{return(paste0("Error, unable to computer ancestry PCs because ", gwasfile, " contains only ", length(common_pca), " of the ", nrow(PCA.DETAILS), " ancestry SNPs"))}

a1_freq=as.numeric(subset_ss[match1,ncol(subset_ss)])
flip_preds=which(subset_ss[match1,2]!=PCA.DETAILS[match2,2])
if(length(flip_preds)>0){a1_freq[flip_preds]=1-a1_freq[flip_preds]}

#compute mean pcs, then get distances between the gwas and the five reference panels

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

distances=rep(NA,5)
for(j in 1:5){distances[j]=((means1[1]-means1[1+j])^2+(means2[1]-means2[1+j])^2+(means3[1]-means3[1+j])^2)^0.5}

closest_ref=which.min(distances)
best_panel=colnames(PCA.DETAILS)[7+closest_ref]

if(min(distances)<1e-8)
{cat(paste0("The closest reference panel is ", best_panel, "; the distance is ", signif(min(distances),2)," indicating a very good match\n"))}
else
{cat(paste0("The closest reference panel is ", best_panel, "; the distance is ", signif(min(distances),2)," indicating the panel is sub-optimal (ideally the distance should be below 1e-8)\n"))}

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

plot_title=paste0("The closest Ref Panel is ", best_panel, " (distance ",signif(min(distances),2),")")

outfile=paste0(outstem,".ancestries.pdf")
pdf(outfile,height=6,width=8)
plot(means1,means2,col=cols,pch=pchs,lwd=3,cex=cexs,xlab="Principal Component Axis 1",ylab="Principal Component Axis 2",main=plot_title)
legend("bottomright",leg=c("AFR","AMR","EAS","EUR","SAS","FIN"),title="Ancestries",fill=1:6,bty="n",cex=1.5,ncol=2)
legend("topright",leg=c("GWAS Summary Stats","MegaPRS Ref Panels","1000GP Populations"),col=c("orange","darkgrey","darkgrey"),cex=1.5,pch=c(3,1,4),bty="n",lwd=3,lty=NA)
dev.off()

plot(means1,means2,col=cols,pch=pchs,lwd=3,cex=cexs,xlab="Principal Component Axis 1",ylab="Principal Component Axis 2",main=plot_title)
legend("bottomright",leg=c("AFR","AMR","EAS","EUR","SAS","FIN"),title="Ancestries",fill=1:6,bty="n",cex=1.5,ncol=2)
legend("topright",leg=c("GWAS Summary Stats","MegaPRS Ref Panels","1000GP Populations"),col=c("orange","darkgrey","darkgrey"),cex=1.5,pch=c(3,1,4),bty="n",lwd=3,lty=NA)

cat(paste0("The PCA plot has been saved in ", outfile,"\n\n"))

#give some advice

cat(paste0("When constructing PRS, we recommend you download the correlations with prefix ", best_panel,".GENO, and use the summary statistics in ", outstem,".GENO.summaries\n"))

if(nrow(subset_ss)/nrow(GENO.SNPs)<0.8)
{cat(paste0("However, please be aware that summary statistics are missing for a relatively large proportion of the genotyped predictors (",100-round(nrow(subset_ss)/nrow(GENO.SNPs)*100),"%), which can result in low-accuracy PRS\n"))}
cat("\n")
}


################
#get total time
end_time=Sys.time()
cat(paste0("End at ",end_time,"\n"))
cat(paste0("Total run time was ", round(difftime(end_time,start_time,units="hours"),2), " hours\n\n"))
}

