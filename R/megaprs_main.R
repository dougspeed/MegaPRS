#Copyright 2024 Doug Speed.

#    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

#    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.
8

################################################


#R code for running a simplified version of MegaPRS
#This file uses functions defined in megaprs_functions.R


################################################

#' Estimate predictor effect sizes
#' @param corstem The stem of the correlations files (e.g., if you have files called example.cors.root, example.cors.windows, etc, then use "corstem=example")
#' @param ssfile The name of file containing GWAS summary statistics in LDAK format (which you either created yourself, or using the function formatSS)
#' @param  outfile The name of the file where effect sizes will be saved
#' @param power Defines the relationship between per-predictor heritability and MAF (Default value is -0.25) 
#' @param shrink How much to shrink predictor-predictor correlations (Default value is 0.8)
#' @param jackknife Whether to estimate standard deviations via jackknifing (Default value FALSE)
#' @param perpredhers (Optional) Used to provide per-predictor heritabilities (the file should have two columns, providing predictor names and heritabilities)
#' @param onechr (Optional) Used to restrict to predictors on a single chromosome
#' @param extractfile (Optional) Used to restrict to predictors in the specified file#' @export
#' @export
#' @examples
#' megaprs(corstem="example", ssfile="gwas.summaries", outfile="output.effects")
megaprs=function(corstem=NULL, ssfile=NULL, outfile=NULL, power=-0.25, shrink=0.9, jackknife=FALSE, perpredhers=NULL, onechr=NULL, extractfile=NULL)
{

################
#some basic checks and a welcome message

cat(paste0("\n"))

if(is.null(corstem))
{return(paste0("Error, you must use the argument corstem to specify the correlation files"))}

if(is.null(ssfile))
{return(paste0("Error, you must use the argument ssfile to specify the summary statistics file (in LDAK format)"))}

if(is.null(outfile))
{return(paste0("Error, you must use the argument outfile to specify the name of the output file"))}

if(power<(-1.5)|power>0.5)
{cat(paste0("Warning, we normally recommend the power is set between -1 and 0 (not ",power,")\n\n"))}

if(shrink<0|shrink>1)
{return(paste0("Error, shrink must be between 0 and 1 (not ",shrink,")"))}

if(!is.null(onechr))
{
if(onechr!=round(onechr)|onechr<1)
{return(paste0("Error, onechr must be a positive integer (not ",onechr,")"))}
}

cat(paste0("Running MegaPRS (assuming the BayesR model) using the correlations with stem ",corstem," and the summary statistics in the file ",ssfile,"; will save the estimated effects to the file ", outfile,"\n\n"))


################
#get start time
start_time=Sys.time()

end_time=Sys.time()
cat(paste0("Start at ",start_time,"\n\n"))


################
#set some values

cvprop=0.1
jackprop=0.2
max_her=0.8
max_iter=50
tol=1e-5


################
#check files exist

for(suffix in c("root","windows","bim","bin","highld","noise"))
{
filename=paste0(corstem,".cors.",suffix)
if(file.exists(filename)==FALSE)
{
return(paste0("Error, unable to find the file ",filename," (make sure you have downloaded the files from https://dougspeed.com/correlations)"))
}
}

if(file.exists(ssfile)==FALSE)
{
return(paste0("Error, unable to find the summary statistics file ",ssfile))
}

if(jackknife==TRUE)
{
filename=paste0(corstem,".cors.jackknifes")
if(file.exists(filename)==FALSE)
{
return(paste0("Error, unable to find the file ",filename," (make sure you have downloaded the files from https://dougspeed.com/correlations)"))
}
}

if(!is.null(perpredhers))
{
if(file.exists(perpredhers)==FALSE)
{
return(paste0("Error, unable to find the per-predictor heritabilities file ",perpredhers))
}
}


################
#set some file names (and test outfile valid)

rootfile=paste0(corstem,".cors.root")
bimfile=paste0(corstem,".cors.bim")
windowfile=paste0(corstem,".cors.windows")
binfile=paste0(corstem,".cors.bin")
highldfile=paste0(corstem,".cors.highld")
noisefile=paste0(corstem,".cors.noise")
jackfile=paste0(corstem,".cors.jackknifes")

write.table("If_this_message_shows_MegaPRS_started_but_did_not_complete",outfile,row.names=FALSE,col.names=FALSE,quote=FALSE)


################
#read rootfile and extract numbers of predictors, windows, pairs and jackknifes (might not use latter)

root=read.table(rootfile,head=FALSE)
num_preds=as.numeric(root[5,2])
num_windows=as.numeric(root[6,2])
num_pairs=as.numeric(root[7,2])
num_jacks=as.numeric(root[8,2])


################
#read bimfile

cat(paste0("Reading predictor details from ",bimfile,"\n"))
bim=read.table(bimfile,head=FALSE,nrow=num_preds)
cat(paste0("In total, there are ",nrow(bim)," predictors\n\n"))


################
#take care of any filtering

exclude_preds=NULL

if(!is.null(onechr))
{
exclude_preds=bim[which(as.numeric(bim[,1])!=onechr),2]

if(length(exclude_preds)>0)
{cat(paste0("Warning, ignoring ", length(exclude_preds)," predictors that are not on Chromosome ", onechr, "\n\n"))}
}

if(!is.null(extractfile))
{
extract_preds=read.table(extractfile)
exclude_extract=setdiff(setdiff(bim[,2],extract_preds[,1]),exclude_preds)

if(length(exclude_extract)>0)
{cat(paste0("Warning, ignoring ", length(exclude_extract)," predictors that are not in ",extractfile , "\n\n"))}

exclude_preds=c(exclude_preds,exclude_extract)

rm(extract_preds)
rm(exclude_extract)
}


################
#read and process summary statistics (also decides which predictors we will use)

ret_list=read_summary_stats(ssfile, bim, bimfile, exclude_preds)
preds=ret_list$preds
rhos=ret_list$rhos
nss=ret_list$nss
chis=ret_list$chis
num_preds_use=ret_list$num_preds_use
neff=ret_list$neff
rm(ret_list)


################
#get indexes of predictors we will use
match_preds=match(preds, bim[,2])


################
#read windows and make indexes, then adjust for filtering

windows=read.table(windowfile,head=T)
window_starts=as.numeric(windows[,2])
window_ends=as.numeric(windows[,3])
window_lengths=window_ends-window_starts+1
window_indexes=cumsum(window_lengths*(window_lengths-1)/2)-window_lengths*(window_lengths-1)/2
rm(windows)

window_starts_orig=window_starts
window_ends_orig=window_ends
window_lengths_orig=window_lengths

for(window in 1:num_windows)
{
within_window=which(match_preds>=window_starts_orig[window]&match_preds<=window_ends_orig[window])

if(length(within_window)>0)
{
window_starts[window]=min(within_window)
window_ends[window]=max(within_window)
}
else
{window_starts[window]=-1;window_ends[window]=-1}
}
window_lengths=window_ends-window_starts+1

#can remove empty windows
keep_windows=which(window_starts!=-1)
num_windows=length(keep_windows)
window_starts=window_starts[keep_windows]
window_ends=window_ends[keep_windows]
window_lengths=window_lengths[keep_windows]
window_starts_orig=window_starts_orig[keep_windows]
window_ends_orig=window_ends_orig[keep_windows]
window_lengths_orig=window_lengths_orig[keep_windows]


################################################
#read top of bin file

filecon=file(binfile, "rb")
centres=readBin(filecon, numeric(), num_preds, 8)[match_preds]
mults=readBin(filecon, numeric(), num_preds, 8)[match_preds]
sqdevs=readBin(filecon, numeric(), num_preds, 8)[match_preds]
rjksums=readBin(filecon, numeric(), num_preds, 8)[match_preds]
close(filecon)


################
#get pseudo summary statistics

ret_list=get_pseudo_stats(noisefile, rhos, nss, match_preds, cvprop)
rhos2=ret_list$rhos2
rhos3=ret_list$rhos3
neff2=ret_list$neff2


################
#get high-ld predictors

temp_ld=read.table(highldfile)
{
if(temp_ld[1,1]=="None")
{
cat(paste0("Warning, there are no predictors are in ",highldfile, "\n\n"))
highld_preds=NULL
}
else
{
highld_common=intersect(preds, temp_ld[,1])
highld_preds=sort(match(highld_common,preds))

if(length(highld_common)>0)
{
cat(paste0(length(highld_common)," of the ",num_preds_use," predictors are in ", highldfile, "\n\n"))
}
else
{
cat(paste0("Warning, none of the ",num_preds_use," predictors are in ", highldfile, "\n\n"))
}
}
}
rm(temp_ld)


################
#get per-predictor heritabilities

if(is.null(perpredhers))	#estimate per-predictor heritabilities
{
#estimate total heritability using ldscores

top=sum(nss/neff*(chis-1))
bot=sum(rjksums*(nss/neff)^2)
her=top/bot*num_preds_use/neff

cat(paste0("Estimated heritability is ",round(her,4),"\n\n"))
if(her<0.01)
{
cat(paste0("Warning, this is very low, so has been increased to 0.01\n\n"))
her=0.01
}
if(her>max_her)
{
cat(paste0("Warning, this is very high, so has been reduced to ",max_her,"\n\n"))
her=max_her
}

#get per-predictor heritabilties (ensure they sum to her)
exps_temp=(centres*(1-centres/2))^(1+power)
exps=exps_temp/sum(exps_temp)*her
rm(exps_temp)
}else	#read per-predictor heritabilities from file
{
exps_all=read.table(perpredhers,head=FALSE)

overlap=intersect(exps_all,preds)
if(length(overlap)<length(preds)){return(paste0("Error, ",perpredhers," contains per-predictor heritabilities for only ", length(overlap), " of the ", length(preds)," predictors"))}

match1=match(preds,exps_all[,1])
exps=as.numeric(exps_all[match1,2])
rm(exps_all)

cat(paste0("Have read per-predictor heritabilities for all ", length(preds), "predictors from ", perpredhers,"\n\n"))

if(sum(is.na(exps))>0){return(paste0("Error, ",sum(is.na(exps))," of the values are NA"))}

exps_neg=which(exps<0)
if(length(exps_neg)>0)
{
cat(paste0("Warning, ",length(exps_neg)," of the values are negative, so will be replaced by zero"))
exps[exps_neg]=0
}

cat(paste0("The sum of the per-predictor heritabilities is ",sum(exps),"\n\n"))
}


################
#set parameter options

num_try=35
tryp1s=rep(NA,num_try)
tryp2s=rep(NA,num_try)
tryp3s=rep(NA,num_try)
tryp4s=rep(NA,num_try)

loads=c(0,0.01,0.05,0.1,0.2)

count=1
for(j4 in 1:5)
{
for(j3 in j4:5)
{
for(j2 in j3:5)
{
tryp4s[count]=loads[j4]
tryp3s[count]=loads[j3]
tryp2s[count]=loads[j2]
tryp1s[count]=1-tryp2s[count]-tryp3s[count]-tryp4s[count]
count=count+1
}}}

#make sure p4s term always non-zero
for(j in 1:num_try)
{
if(tryp4s[j]==0)
{
if(tryp3s[j]>0){tryp4s[j]=tryp3s[j];tryp3s[j]=tryp2s[j];tryp2s[j]=0;}
else	#so p4 and p3 are zero
{
if(tryp2s[j]>0){tryp4s[j]=tryp2s[j];tryp2s[j]=0;}
else	#p4, p3 and p2 are all zero (ridge model)
{tryp4s[j]=tryp1s[j];tryp1s[j]=0;}
}
}
}


################
#get jackknife summary statistics (if used)

if(jackknife==TRUE)
{
YTdata2=matrix(numeric(),num_preds_use,num_jacks)
neff3=neff*(1-jackprop)
jack_scalar=(jackprop/(1-jackprop))^0.5
jack_scalar=0
neff3=neff

filecon=file(jackfile, "rb")
for(j in 1:num_jacks)
{
data_rands=readBin(filecon, numeric(), num_preds, 4)[match_preds]
YTdata2[,j]=neff3*(rhos+data_rands*jack_scalar*nss^0.5)
}
close(filecon)
}


################
#might be helpful to clear memory before we perform the heavy computation
gc()


################
#construct and test training models

cat(paste0("Constructing ",num_try," models using pseudo training summary statistics\n\n"))

#some preparation
YTdata=rhos2*neff2
effs_all=matrix(numeric(),num_preds_use,num_try)
failed_windows=0
failed_preds=0
est_cov=matrix(0,num_try,1)
est_var=matrix(0,num_try,1)


#set starting effect sizes
for(k in 1:num_try)
{
effs_all[,k]=get_postmean_multi(YTdata, neff2, exps, tryp1s[k], tryp2s[k], tryp3s[k], tryp4s[k])/rjksums
gc()
}


#will ignore highld predictors for training models
effs_all[highld_preds,]=0


#construct models
for(window in 1:num_windows)
{
if(window%%10==1){cat(paste0("Estimating effect sizes for Window ",window," of ",num_windows,"\n"))}

#get start and end predictors of window
start=window_starts[window]
end=window_ends[window]
window_preds=setdiff(start:end,highld_preds)

if(length(window_preds)==1)	#simple case, no need to read correlations nor iterate
{
j=start

#get t(Xj) residuals and new effect
XjTR_all=rep(YTdata[j],num_try)
effs_all[j,]=get_postmean_multi(XjTR_all, neff2, exps[j], tryp1s, tryp2s, tryp3s, tryp4s)

est_cov=est_cov+effs_all[j,]*rhos3[j]
est_var=est_var+effs_all[j,]^2
}

if(length(window_preds)>1)	#hard case, must read correlations and iterate
{
#save starting effects for window, in case convergence fails
effs_all_save=effs_all[window_preds,]

#get correlations for full block
all_cors=diag(window_lengths_orig[window])/2
lower_tri_cors=which(lower.tri(all_cors))

filecon=file(binfile, "rb")
seek(filecon,4*num_preds*8+4*window_indexes[window])
all_cors[lower_tri_cors]=readBin(filecon, numeric(), window_lengths_orig[window]*(window_lengths_orig[window]-1)/2, 4)
close(filecon)

all_cors=all_cors+t(all_cors)

#reduce to predictors we are using
match_cors=match_preds[window_preds]-window_starts_orig[window]+1
use_cors=all_cors[match_cors,match_cors]
rm(all_cors)
rm(match_cors)

#get explained sum of squares
ess_all=2*t(YTdata[window_preds])%*%effs_all[window_preds,]/neff2-apply(effs_all[window_preds,]*(use_cors%*%effs_all[window_preds,]),2,sum)

for(iter in 1:max_iter)
{
for(j in window_preds)	#update effect size for predictor j
{
#get t(Xj) residuals
XjTR_all=rep(YTdata[j],num_try)-shrink*neff2*(t(effs_all[window_preds,])%*%use_cors[,j-start+1]-effs_all[j,])

#then new effect sizes
effs_all[j,]=get_postmean_multi(XjTR_all, neff2, exps[j], tryp1s, tryp2s, tryp3s, tryp4s)

#check for NAs
if(sum(is.na(effs_all[j,]))>0){return(paste0("Error, algorithm failed for Window ",window," of ",num_windows," - length ",length(window_preds)," - predictor ",j," - ",sum(is.na(effs_all[j,])),"\n\n"))}
}	#end of j loop

#save previous sum of squares, then get current
ess_all_save=ess_all
ess_all=2*t(YTdata[window_preds])%*%effs_all[window_preds,]/neff2-apply(effs_all[window_preds,]*(use_cors%*%effs_all[window_preds,]),2,sum)

#see if converged
if(sum(abs(ess_all-ess_all_save)<tol)==num_try){break}
}	#end of iter loop

if(sum(ess_all>=1)>0)	#estimates are suspect - revert to saved estimates
{
effs_all[window_preds,]=effs_all_save
failed_windows=failed_windows+1
failed_preds=failed_preds+window_lengths[window]
}

#add on contribution to variances (first setting small correlations to zero)
small_cors=which(abs(use_cors)<0.01)
if(length(small_cors)>0){use_cors[small_cors]=0}

est_cov=est_cov+t(effs_all[window_preds,])%*%rhos3[window_preds]
est_var=est_var+apply(effs_all[window_preds,]*(use_cors%*%effs_all[window_preds,]),2,sum)
}	#end of hard case

gc()
}	#end of window loop

if(failed_windows==0){cat(paste0("\nCompleted: all windows converged\n\n"))}else{cat(paste0("\nCompleted: ",failed_windows," windows failed to converge (in total, ",failed_preds," predictors)\n\n"))}

model_cors=as.numeric(est_cov)/as.numeric(est_var)^.5

for(k in 1:num_try)
{
cat(paste(k,est_cov[k], est_var[k], model_cors[k],"\n"))
}

if(sum(!is.na(model_cors))==0)
{
return(paste0("Error, it was not possible to compute the correlation for any of the models"))
}

cat(paste0("The correlations range from ", min(model_cors,na.rm=TRUE), " to ", max(model_cors,na.rm=TRUE),"\n"))

model_best=which.max(model_cors)
cat(paste0("The best model has has the following parameters: prob_big=",tryp4s[model_best],", prob_med=",tryp3s[model_best],", prob_small=",tryp2s[model_best],"\n\n"))

rm(effs_all)


################
#helpful to clear memory
gc()


################
#construct final model

cat(paste0("Constructing the final model using the original summary statistics\n\n"))

#some preparation
YTdata=rhos*neff
failed_windows=0
failed_preds=0

#set starting values
effs=get_postmean_multi(YTdata, neff, exps, tryp1s[model_best], tryp2s[model_best], tryp3s[model_best], tryp4s[model_best])/rjksums


#construct models
for(window in 1:num_windows)
{
if(window%%100==1){cat(paste0("Estimating effect sizes for Window ",window," of ",num_windows,"\n"))}

#get start and end predictors of window
start=window_starts[window]
end=window_ends[window]
window_preds=start:end

if(length(window_preds)==1)	#simple case, no need to iterate
{
j=start

#get t(Xj) residuals and new effect
XjTR=YTdata[j]
effs[j]=get_postmean(XjTR, neff, exps[j], tryp1s[model_best], tryp2s[model_best], tryp3s[model_best], tryp4s[model_best])
}

if(length(window_preds)>1)	#hard case, must read correlations and iterate
{
#save starting effects for window, in case convergence fails
effs_save=effs[window_preds]

#get correlations for full block
all_cors=diag(window_lengths_orig[window])/2
lower_tri_cors=which(lower.tri(all_cors))

filecon=file(binfile, "rb")
seek(filecon,4*num_preds*8+4*window_indexes[window])
all_cors[lower_tri_cors]=readBin(filecon, numeric(), window_lengths_orig[window]*(window_lengths_orig[window]-1)/2, 4)
close(filecon)

all_cors=all_cors+t(all_cors)

#reduce to predictors we are using
match_cors=match_preds[window_preds]-window_starts_orig[window]+1
use_cors=all_cors[match_cors,match_cors]

#get explained sum of squares
ess=2*sum(YTdata[window_preds]*effs[window_preds])/neff-sum(effs[window_preds]*(use_cors%*%effs[window_preds]))

for(iter in 1:max_iter)
{
for(j in window_preds)	#update effect size for predictor j
{
#get t(Xj) residuals
XjTR=YTdata[j]-shrink*neff*(sum(effs[window_preds]*use_cors[,j-start+1])-effs[j])

#then new effect
effs[j]=get_postmean(XjTR, neff, exps[j], tryp1s[model_best], tryp2s[model_best], tryp3s[model_best], tryp4s[model_best])

#check for NA
if(is.na(effs[j])){return(paste0("Error, algorithm failed for Window ",window," of ",num_windows," - length ",length(window_preds)," - predictor ",j))}
}	#end of j loop

#save previous sum of squares, then get current
ess_save=ess
ess=2*sum(YTdata[window_preds]*effs[window_preds])/neff-sum(effs[window_preds]*(use_cors%*%effs[window_preds]))

#see if converged
if(abs(ess-ess_save)<tol){break}
}	#end of iter loop

if(ess>=1)	#estimates are suspect - revert to saved estimates
{
effs[window_preds]=effs_save
failed_windows=failed_windows+1
failed_preds=failed_preds+window_lengths[window]
}
}	#end of hard case

gc()
}	#end of window loop

if(failed_windows==0){cat(paste0("\nCompleted: all windows converged\n\n"))}
else{cat(paste0("\nCompleted: ",failed_windows," windows failed to converge (in total, ",failed_preds," predictors)\n\n"))}


################
#construct jackknife models (if used)

if(jackknife==TRUE)
{
cat(paste0("Constructing ",num_jacks," jackknife models\n\n"))

#some preparation
effs_all=matrix(numeric(),num_preds_use,num_jacks)
failed_windows=0
failed_preds=0

#set starting values
for(k in 1:num_jacks)
{
effs_all[,k]=get_postmean_multi(YTdata2[,k], neff3, exps, tryp1s[model_best], tryp2s[model_best], tryp3s[model_best], tryp4s[model_best])/rjksums
gc()
}

#construct models
for(window in 1:num_windows)
{
if(window%%10==1){cat(paste0("Estimating effect sizes for Window ",window," of ",num_windows,"\n"))}

#get start and end predictors of window
start=window_starts[window]
end=window_ends[window]
window_preds=start:end

if(length(window_preds)==1)	#simple case, no need to read correlations nor iterate
{
j=start

#get t(Xj) residuals and new effect
XjTR_all=YTdata2[j,]
effs_all[j,]=get_postmean_multi(XjTR_all, neff3, exps[j], tryp1s[model_best], tryp2s[model_best], tryp3s[model_best], tryp4s[model_best])
}

if(length(window_preds)>1)	#hard case, must read correlations and iterate
{
#save starting effects for window, in case convergence fails
effs_all_save=effs_all[window_preds,]

#get correlations for full block
all_cors=diag(window_lengths_orig[window])/2
lower_tri_cors=which(lower.tri(all_cors))

filecon=file(binfile, "rb")
seek(filecon,4*num_preds*8+4*window_indexes[window])
all_cors[lower_tri_cors]=readBin(filecon, numeric(), window_lengths_orig[window]*(window_lengths_orig[window]-1)/2, 4)
close(filecon)

all_cors=all_cors+t(all_cors)

#reduce to predictors we are using
match_cors=match_preds[window_preds]-window_starts_orig[window]+1
use_cors=all_cors[match_cors,match_cors]
rm(all_cors)
rm(match_cors)

#get explained sum of squares
ess_all=2*t(YTdata[window_preds])%*%effs_all[window_preds,]/neff3-apply(effs_all[window_preds,]*(use_cors%*%effs_all[window_preds,]),2,sum)

for(iter in 1:max_iter)
{
for(j in window_preds)	#update effect size for predictor j
{
#get t(Xj) residuals
XjTR_all=rep(YTdata[j],num_try)-shrink*neff3*(t(effs_all[window_preds,])%*%use_cors[,j-start+1]-effs_all[j,])

#then new effect sizes
effs_all[j,]=get_postmean_multi(XjTR_all, neff3, exps[j], tryp1s[model_best], tryp2s[model_best], tryp3s[model_best], tryp4s[model_best])

#check for NAs
if(sum(is.na(effs_all[j,]))>0){return(paste0("Error, algorithm failed for Window ",window," of ",num_windows," - length ",length(window_preds)," - predictor ",j," - fail ",sum(is.na(effs_all[j,])),"\n\n"))}
}	#end of j loop

#save previous sum of squares, then get current
ess_all_save=ess_all
ess_all=2*t(YTdata[window_preds])%*%effs_all[window_preds,]/neff3-apply(effs_all[window_preds,]*(use_cors%*%effs_all[window_preds,]),2,sum)

#see if converged
if(sum(abs(ess_all-ess_all_save)<tol)==num_try){break}
}	#end of iter loop

if(sum(ess_all>=1)>0)	#estimates are suspect - revert to saved estimates
{
effs_all[window_preds,]=effs_all_save
failed_windows=failed_windows+1
failed_preds=failed_preds+window_lengths[window]
}
}	#end of hard case

gc()
}	#end of window loop

if(failed_windows==0){cat(paste0("\nCompleted: all windows converged\n\n"))}
else{cat(paste0("\nCompleted: ",failed_windows," windows failed to converge (in total, ",failed_preds," predictors)\n\n"))}
}

################
#save effects

if(jackknife==FALSE)
{
final_model=cbind(bim[match_preds,c(2,5,6)], centres, signif(effs*mults,4))
colnames(final_model)=c("Predictor","A1","A2","Centre","Effect")
}
else
{
final_model=cbind(bim[match_preds,c(2,5,6)], centres, signif(effs*mults,4))
colnames(final_model)=c("Predictor","A1","A2","Centre","Effect")
for(j in 1:num_jacks)
{
final_model=cbind(final_model, signif(effs_all[,j]*mults,4))
colnames(final_model)[5+j]=paste0("Jackknife_",jackprob)
}
}

write.table(final_model,outfile,row.names=FALSE,col.names=TRUE,quote=FALSE)

cat(paste0("The estimated effect sizes are saved in the file ",outfile,"\n\n"))


################
#get total time
end_time=Sys.time()
cat(paste0("End at ",end_time,"\n"))
cat(paste0("Total run time was ", round(difftime(end_time,start_time,units="hours"),2), " hours\n\n"))

}

