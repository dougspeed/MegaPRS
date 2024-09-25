#Copyright 2024 Doug Speed.

#    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

#    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.


################################################


#R code for running a simplified version of MegaPRS
#This file contains functions called by megaprs.R


################################################

get_postmean=function(XjTR, neff, h2j, prob_null, prob_small, prob_med, prob_big)
{
#this function returns the posterior mean of the effect size beta_j

#assumes prior distribution beta_j sim prob_null N(0,0) + prob_small N(0,V/100) + prob_med N(0,V/10) + prob_big N(0,V)
#where V is set so that E(beta_j^2)=h2j

#t(Xj)R is the dot product of Xj and current residuals; n is the effective sample size

#posterior distribution will be a point mass plus three normal distributions


#compute V
V=h2j/(prob_small/100+prob_med/10+prob_big)


#compute posterior means and variances for the normal components of the posterior
mean_small=XjTR/(neff+100/V);var_small=1/(neff+100/V)
mean_med=XjTR/(neff+10/V);var_med=1/(neff+10/V)
mean_big=XjTR/(neff+1/V);var_big=1/(neff+1/V)


#maxval is used just below to avoid very large numbers
maxval=max(0, mean_small^2/var_small, mean_med^2/var_med, mean_big^2/var_big)


#compute relative posterior probabilities
rel_null=prob_null*exp(-0.5*maxval)
rel_small=prob_small*sqrt(var_small*100/V)*exp(0.5*mean_small^2/var_small-0.5*maxval)
rel_med=prob_med*sqrt(var_med*10/V)*exp(0.5*mean_med^2/var_med-0.5*maxval)
rel_big=prob_big*sqrt(var_big*1/V)*exp(0.5*mean_big^2/var_big-0.5*maxval)


#compute absolute posterior probabilies
abs_small=rel_small/(rel_null+rel_small+rel_med+rel_big)
abs_med=rel_med/(rel_null+rel_small+rel_med+rel_big)
abs_big=rel_big/(rel_null+rel_small+rel_med+rel_big)


#return posterior mean
return(mean_small*abs_small + mean_med*abs_med + mean_big*abs_big)
}


################################################


get_postmean_multi=function(XjTR, neff, h2j, prob_null, prob_small, prob_med, prob_big)
{
#same as above, except handles multiple values


#compute Vs
V=h2j/(prob_small/100+prob_med/10+prob_big)


#compute posterior means and variances for the normal components of the posteriors
mean_small=XjTR/(neff+100/V);var_small=1/(neff+100/V)
mean_med=XjTR/(neff+10/V);var_med=1/(neff+10/V)
mean_big=XjTR/(neff+1/V);var_big=1/(neff+1/V)


#maxval is used just below to avoid very large numbers
maxval=apply(rbind(0, mean_small^2/var_small, mean_med^2/var_med, mean_big^2/var_big),2,max)


#compute relative posterior probabilities
rel_null=prob_null*exp(-0.5*maxval)
rel_small=prob_small*sqrt(var_small*100/V)*exp(0.5*mean_small^2/var_small-0.5*maxval)
rel_med=prob_med*sqrt(var_med*10/V)*exp(0.5*mean_med^2/var_med-0.5*maxval)
rel_big=prob_big*sqrt(var_big*1/V)*exp(0.5*mean_big^2/var_big-0.5*maxval)


#compute absolute posterior probabilies
abs_small=rel_small/(rel_null+rel_small+rel_med+rel_big)
abs_med=rel_med/(rel_null+rel_small+rel_med+rel_big)
abs_big=rel_big/(rel_null+rel_small+rel_med+rel_big)


#return posterior means
return(mean_small*abs_small + mean_med*abs_med + mean_big*abs_big)
}



################################################


read_summary_stats=function(ssfile, bim, bimfile, exclude_preds)
{
#this function reads the summary statistics, works out the overlap with bimfile and returns key details


#read first row of summary statistics and check have required columns

ss_cols=read.table(ssfile,nrow=1)

if(length(intersect(c("Predictor","A1","A2","n"),ss_cols))!=4)
{
return(paste0("Error, ",ssfile," should have columns labelled Predictor, A1, A2 and n (see https://dougspeed.com/summary-statistics for further details)"))
}
col_pred=which(ss_cols=="Predictor")
col_A1=which(ss_cols=="A1")
col_A2=which(ss_cols=="A2")
col_n=which(ss_cols=="n")

type=0
if(length(intersect("Z",ss_cols))==1)	#using Z
{type=1;}
if(type==0&length(intersect(c("Direction","Stat"),ss_cols))==2)	#using Direction and Stat
{type=2;}
if(type==0&length(intersect(c("Direction","P"),ss_cols))==3)	#using Direction and P
{type=3;}

if(type==0)
{
return(paste0("Error, ",ssfile," should either have a column labelled Z, or a column labelled Direction and a column labelled either Stat or P (see https://dougspeed.com/summary-statistics for further details)"))
}


#read in all summary statistics and see which have different alleles and are in the bimfile

cat(paste0("Reading summary statistics from ",ssfile,"\n"))
ss_all=read.table(ssfile,head=TRUE)
cat(paste0("In total, there are summary statistics for ",nrow(ss_all)," predictors\n\n"))

diff_alleles=which(ss_all[,col_A1]!=ss_all[,col_A2])
if(length(diff_alleles)==0)
{
return(paste0("Error, none of the predictors in ",ssfile," have distinct alleles"))
}

common_start=intersect(ss_all[diff_alleles,col_pred],bim[,2])
if(length(common_start)==0)
{
return(paste0("Error, there is no overlap between the predictors in ",ssfile," and those in ",bimfile))
}

cat(paste0("There are ",length(common_start)," predictors common to ",ssfile," and ",bimfile,"\n"))


#use only predictors with consistent alleles

match1=match(common_start,bim[,2])
match2=match(common_start,ss_all[,col_pred])

find_a1= (bim[match1,5]==ss_all[match2,col_A1]) + (bim[match1,5]==ss_all[match2,col_A2])
find_a2= (bim[match1,6]==ss_all[match2,col_A1]) + (bim[match1,6]==ss_all[match2,col_A2])

common_mid=common_start[which(find_a1==1&find_a2==1)]
if(length(common_mid)==0)
{
return(paste0("Error, none of the these have consistent alleles"))
}

if(length(common_mid)==length(common_start))
{
cat(paste0("All of these have consistent alleles\n\n"))
}
else
{
cat(paste0("Only ",length(common_mid)," of these have consistent alleles\n\n"))
}


#exclude predictors that failed qc
common_end=setdiff(common_mid, exclude_preds)

if(length(common_end)==0)
{return(paste0("Error, after filtering predictors, none remain"))}

if(length(common_end)<length(common_mid))
{
cat(paste0("After filtering predictors ",length(common_end)," predictors remain\n\n"))
}


#load preds (predictors we will use, in bimfile order), get sample sizes, signs, chis and rhos

num_preds_use=length(common_end)
match1=sort(match(common_end,bim[,2]))
preds=bim[match1,2]
match2=match(preds,ss_all[,col_pred])

nss=as.numeric(ss_all[match2,col_n])

if(type==1)	#using Z
{
col_Z=which(ss_cols=="Z")
chis=as.numeric(ss_all[match2,col_Z])^2
signs=sign(as.numeric(ss_all[match2,col_Z]))
}

if(type==2)	#using Direction and Stat
{
col_dir=which(ss_cols=="Direction")
col_stat=which(ss_cols=="Stat")
signs=sign(as.numeric(ss_all[match2,col_dir]))
chis=as.numeric(ss_all[match2,col_Stat])
}

if(type==3)	#using Direction and P
{
col_dir=which(ss_cols=="Direction")
col_p=which(ss_cols=="P")
signs=sign(as.numeric(ss_all[match2,col_dir]))
pvals=as.numeric(ss_all[match2,col_P])

small_p=which(pvals<1e-300)
if(length(small_p)>0)
{
pvals[small_p]=1e-300
cat(paste0("Warning, ",length(small_p)," predictors have p-values below 1e-300; these have been replaced with 1e-300\n\n"))
}

chis=qchisq(pvals,1,low=F)
}

#do we need to flip signs?
match1=match(preds,bim[,2])
match2=match(preds,ss_all[,col_pred])
flip_preds=which(bim[match1,5]==ss_all[match2,col_A2])
if(length(flip_preds)>0){signs[flip_preds]=-signs[flip_preds]}

#also need correlations and effective sample size 
rhos=signs*sqrt(chis/(chis+nss))
neff=mean(nss)


#return the things we require
ret_list=list("preds"=preds, "rhos"=rhos, "nss"=nss, "chis"=chis, "num_preds_use"=num_preds_use, "neff"=neff)
return(ret_list)
}


################################################


get_pseudo_stats=function(noisefile, rhos, nss, match_preds, cvprop)
{
#this function makes the pseudo summary statistics (specifically, rhos2, rhos3 and neff2)


datarands=as.numeric(read.table(noisefile)[match_preds,1])
value=sqrt(cvprop/(1-cvprop))
rhos2=rhos+datarands*value/sqrt(nss)

rhos3=(rhos-(1-cvprop)*rhos2)/cvprop
neff2=(1-cvprop)*mean(nss)


#return the things we require
ret_list=list("rhos2"=rhos2, "rhos3"=rhos3, "neff2"=neff2)
return(ret_list)
}


################################################


