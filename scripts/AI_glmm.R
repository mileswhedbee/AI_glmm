rm(list = ls())
library(dplyr)
library(stringr)

read.csv(file = "data/Aripo_dp_mat.csv") -> aripo_dp
read.csv(file = "data/aripo.mat") -> aripo_sr
read.csv(file = "data/Quare_dp_mat.csv") -> quare_dp
read.csv(file = "data/quare.mat") -> quare_sr
read.csv(file = "data/SraRunTable.csv") -> sra_table

 # make sample table from sra run table
sra_table %>% filter(Strain == "Aripo") -> sra_aripo
#sra_table %>% filter(Strain == "Quare") -> sra_quare
strsplit(sra_aripo$Library.Name, c("_")) -> tmp
#strsplit(sra_quare$Library.Name, "_") -> tmp.quare
as.data.frame(tmp) -> tmp
#as.data.frame(tmp.quare) -> tmp.quare
t(tmp) -> tmp
#t(tmp.quare) -> tmp.quare
strsplit(tmp[1:40,2], split = "-") -> tmp2
#strsplit(tmp.quare[1:60,2], split = "-") -> tmp.quare2
as.data.frame(tmp2) -> tmp2
#as.data.frame(tmp.quare2) -> tmp.quare2
t(tmp2) -> tmp2
#t(tmp.quare2) -> tmp.quare2
cbind(tmp[1:40,1], tmp2[1:40,1:2]) -> tmp3
#cbind(tmp.quare[1:60,1], tmp.quare2[1:60,1:2]) -> tmp.quare3
as.data.frame(tmp3) -> samples
#as.data.frame(tmp.quare3) -> samples.quare
colnames(samples) = c("sample","Pop","Rearing")
#colnames(samples.quare) = c("sample","Pop","Rearing")
rownames(samples) = 1:40
#rownames(samples.quare) = 1:60

Family= c( 208,211,208,212,207,209,207,207,210,210,210,210,214,204,204,205,214,211,214,219,206,218,217,219,217,219,201,215,201,215,216,203,216,203,216,215,201,204,202,205)
samples[order(samples$sample),] -> samples
samples$Family = Family

read.csv(file = "data/BehO2_Male_decode.csv") -> samples.quare
samples.quare %>% dplyr::select(-X, -Sequenced) -> samples.quare

colnames(aripo_sr)[1] = "Id"
colnames(quare_sr)[1] = "Id"
aripo_dp$Id = paste(aripo_dp$Chrom, aripo_dp$Location, "aripo", sep = "_")
quare_dp$Id = paste(quare_dp$Chrom, quare_dp$Location, "quare", sep = "_")
# match rows between SR and DP matrices
tmp.out = data.frame()
for (id in aripo_sr$Id) {
  aripo_dp %>% filter(Id == id) -> tmp
  tmp.out = rbind(tmp.out, tmp)
}
aripo_dp = tmp.out
# match rows between SR and DP, quare
tmp.out = data.frame()
for (id in quare_sr$Id) {
  quare_dp %>% filter(Id == id) -> tmp
  tmp.out = rbind(tmp.out, tmp)
}
quare_dp = tmp.out
# convert sra to sample id
tmp4 = ""

for (col in colnames(aripo_dp)[3:38]) {
  strsplit(col, split = ".tmp") -> tmp
  tmp2 = tmp[[1]]
  sra_table %>% filter(Run == tmp2) -> tmp3
  tmp4 = append(tmp4, tmp3$Library.Name)
}
aripo_dp_ids = tmp4[2:37]
paste("X",aripo_dp_ids, sep = "") -> aripo_dp_ids
strsplit(aripo_dp_ids, split = "-") -> tmp
as.data.frame(tmp) -> tmp
t(tmp) -> tmp
as.data.frame(tmp) -> tmp
paste(tmp$V1, tmp$V2, sep = ".") -> aripo_dp_ids

# select same set of cols as dp matrix
aripo_sr %>% dplyr::select(aripo_dp_ids) -> aripo_sr

# convert sra to sample id - quare
tmp4 = ""
for (col in colnames(quare_dp)[3:62]) {
  strsplit(col, split = ".tmp") -> tmp
  tmp2 = tmp[[1]]
  sra_table %>% filter(Run == tmp2) -> tmp3
  tmp4 = append(tmp4, tmp3$Library.Name)
}
quare_dp_ids = tmp4[2:60]
strsplit(quare_dp_ids, split = "-") -> tmp
as.data.frame(tmp) -> tmp
t(tmp) -> tmp
as.data.frame(tmp) -> tmp
paste(tmp$V1, tmp$V2, sep = ".") -> quare_dp_ids

colnames(quare_sr)[2:59] -> quare_sr_ids
colnames(quare_dp)[3:61] = quare_dp_ids
# select same set of cols as dp matrix
aripo_sr %>% dplyr::select(aripo_dp_ids) -> aripo_sr
quare_dp %>% dplyr::select(quare_sr_ids) -> quare_dp

write.csv(quare_sr, file = "output/quare_sr.mat")
write.csv(quare_dp, file = "output/quare_dp.mat")
# create row names from id col
rownames(quare_dp) = quare_dp$Id
rownames(quare_sr) = quare_sr$Id
quare_sr %>% dplyr::select(-Id) -> quare_sr
# clean up dp matrix
rownames(aripo_dp) = aripo_dp$Id
rownames(aripo_sr) = aripo_dp$Id
aripo_dp %>% dplyr::select(-Chrom, -Location, -Id) -> aripo_dp
# change from srr colnames to actual Ids
colnames(aripo_dp) = aripo_dp_ids

strsplit(colnames(aripo_sr), split = "X") -> tmp
as.data.frame(tmp) -> tmp
t(tmp) -> tmp
strsplit(tmp[1:36,2], split = "_") -> tmp
as.data.frame(tmp) -> tmp
t(tmp) -> tmp
as.data.frame(tmp) -> tmp

tmp.out = data.frame()
for (id in tmp$V1) {
  samples %>% filter(sample == id) -> tmp2
  tmp.out = rbind(tmp.out, tmp2)
}
samples = tmp.out

colnames(aripo_sr) = samples$sample
colnames(aripo_dp) = samples$sample

as.factor(samples$Family) -> samples$Family

as.matrix(aripo_sr) -> aripo_sr
as.matrix(aripo_dp) -> aripo_dp
# write matrices to file
write.csv(aripo_sr, file = "output/aripo_sr.mat")
write.csv(aripo_dp, file = "output/aripo_dp.mat")

# change quare colnames to match samples.quare
strsplit(colnames(quare_sr), split = "_") -> tmp 
as.data.frame(tmp) -> tmp
t(tmp) -> tmp
colnames(quare_sr) = tmp[1:58,1] # add col names
colnames(quare_dp) = tmp[1:58,1]
colnames(samples.quare)[1] = "sample"
# set samples.quare equal to quare matrices
tmp.out = data.frame()
for (id in colnames(quare_sr)) {
  samples.quare %>% filter(sample == id) -> tmp
  tmp.out = rbind(tmp.out, tmp)
}
samples.quare = tmp.out
as.matrix(quare_dp) -> quare_dp
as.matrix(quare_sr) -> quare_sr


####################################
# run glmm
####################################
library(dplyr)
library(tidyr)
library(lme4)
#install.packages(c("lme4", "MASS", "parallel","doSNOW","foreach","optimx","lsmeans"))
library(MASS)
library(parallel)
library(doSNOW)
library(foreach)
library(optimx)
library(lsmeans)

# test first N models
#aripo_dp[1:10,] -> aripo_dp
#aripo_sr[1:10,] -> aripo_sr

#quare_dp[1:10,] -> quare_dp
#quare_sr[1:10,] -> quare_sr

i = 1
p = nrow(quare_sr)
#LRT results + DE results
DE_re = foreach(i=1:p, .combine = rbind,.packages = c("lme4","MASS","lsmeans"),
                .errorhandling = "remove",.verbose=TRUE) %dopar%{
                  #If fail to converge by Negative binomial, use poisson.
                  #If fail to converge by Negative binomial, change tolPwrss. use poisson.
                  warn_flag= FALSE
                  #If any convergence warning, set warn_flag = TRUE
                  #If any errors occur, use different optimization method for the second stage optimization.
                  #Fit by default algorithm, if fails, set error_flag=TRUE
                  err_flag=FALSE
                  tryCatch(withCallingHandlers({fit_glmer<-glmer.nb(quare_sr[i,] ~ Pop*Rearing + (1|Family) + offset(log(quare_dp[i,])),
                                                                    data = samples.quare)},
                                               warning = function(w){ warn_flag<<- TRUE
                                               invokeRestart("muffleWarning")
                                               })
                           ,error = function(e) err_flag<<-TRUE,finally=(if(!err_flag)alg_optim <<- c("bobyqa", "Nelder_Mead")))
                  #Try lme4 optimizers
                  optim_alg = c("bobyqa","Nelder_Mead","nloptwrap","nlminbwrap")
                  j=1
                  while(err_flag & (j<=length(optim_alg))){
                    err_flag = FALSE
                    tryCatch(withCallingHandlers({fit_glmer<-glmer.nb(quare_sr[i,] ~Pop*Rearing + (1|Family),
                    data = samples.quare, 
                    control=glmerControl(optimizer=optim_alg[j]))},
                    warning = function(w){
                      warn_flag<<- TRUE
                      invokeRestart("muffleWarning")
                    })
                    ,error = function(e) err_flag<<-TRUE, finally = {if(!err_flag) alg_optim<<-optim_alg[j]})
j=j+1
                  }
                  
                  #If no error
                  if(err_flag) optimx_flag=TRUE
                  else optimx_flag=FALSE
                  #Try optimx optimizers if still errors
                  j=1
                  optimx_alg = c('Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', 'nlm',
                                 'nlminb', 'spg', 'ucminf', 'newuoa', 'bobyqa',
                                 'nmkb', 'hjkb', 'Rcgmin', 'Rvmmin')
                  
                  while(err_flag & (j<=length(optimx_alg))){
                    err_flag = FALSE
                    tryCatch(withCallingHandlers({fit_glmer<-glmer.nb(quare_sr[i,] ~Pop*Rearing + (1|Family),
                    data = samples.quare, 
                    control=glmerControl(optimizer="optimx",
                                         optCtrl=list(method=optimx_alg[j])))},
                    warning = function(w){
                      warn_flag<<- TRUE
                      invokeRestart("muffleWarning")
                    })
                    ,error = function(e) err_flag<<-TRUE, finally = {alg_optim<<-optimx_alg[j]})
j=j+1
                  }
                  #Update to remove warnings: up to 10 times
                  #If error occures in updates, ignore it.
                  j=0
                  while(warn_flag & j<10){
                    warn_flag = FALSE
                    if(!optimx_flag){
                      tryCatch(withCallingHandlers(fit_glmer<-update(fit_glmer,start=getME(fit_glmer,c("theta","fixef")),
                                                                     control=glmerControl(optimizer=alg_optim)),
                                                   warning=function(w){warn_flag<<-TRUE;invokeRestart("muffleWarning")})
                               ,error = function(e) warn_flag=FALSE)
                    }else{
                      tryCatch(withCallingHandlers(fit_glmer<-update(fit_glmer,start=getME(fit_glmer,c("theta","fixef")),
                                                                     control=glmerControl(optimizer="optimx",
                                                                                          optCtrl=list(method=optimx_alg[j]))),
                                                   warning=function(w){warn_flag<<-TRUE;invokeRestart("muffleWarning")})
                               ,error = function(e) warn_flag=FALSE)
                    }
                    j=j+1
                    if( j > 8) {
                      singular_fit = "TRUE"
                    }else{
                      singular_fit = "FALSE"
                    }
                  }
                  
                  #Fit
                  sum_glmer = summary(fit_glmer)
                  
                  ## Model comparison test between glmer.nb and glm.nb to check random effects
                  fit_glm <- glm.nb(quare_sr[i,] ~Pop*Rearing, data = samples.quare)
                  
                  model_comp = anova(fit_glmer, fit_glm)[2,6:8]
                  #test contrast to check population and rearing effect
                  # main_eff=test(lsmeans(fit_glmer,~group,contr=list(pop=c(0.5,0.5,-0.5,-0.5),rear=c(0.5,-0.5,0.5,-0.5),
                  #                                                   interaction=c(1,-1,-1,1)))$contrasts)
                  main_eff=test(lsmeans(fit_glmer,~Pop*Rearing,contr=list(Pop=c(-0.5,0.5,-0.5,0.5),Rearing=c(0.5,0.5,-0.5,-0.5),
                                                                          inter=c(1,-1,-1,1),
                                                                          simple_HPP_LPP= c(-1,1,0,0),
                                                                          simple_HPNP_LPNP= c(0,0,-1,1),
                                                                          simple_LPP_LPNP= c(1,0,-1,0),
                                                                          simple_HPP_HPNP= c(0,1,0,-1)
                  ))$contrasts)
                  
                  
                  return(c(i, id = rownames(quare_sr)[i],
                           coef_glmer = getME(fit_glmer,"beta"), theta_glmer = getME(fit_glmer,"glmer.nb.theta"),
                           var_re = getME(fit_glmer,"theta")^2,
                           estim_pop = main_eff[1,2], se_pop = main_eff[1,3], 
                           stat_pop_wald=main_eff[1,5], pval_pop_wald = main_eff[1,6],
                           estim_rear = main_eff[2,2], se_rear = main_eff[2,3], 
                           stat_rear_wald=main_eff[2,5], pval_rear_wald = main_eff[2,6],
                           estim_int = main_eff[3,2], se_int = main_eff[3,3], stat_int_wald=main_eff[3,5], 
                           pval_int_wald = main_eff[3,6],
                           estim_simple_HPP_LPP = main_eff[4,2], se_simple_HPP_LPP = main_eff[4,3], 
                           stat_simple_HPP_LPP_wald=main_eff[4,5], pval_simple_HPP_LPP_wald = main_eff[4,6],
                           estim_simple_HPNP_LPNP = main_eff[5,2], se_simple_HPNP_LPNP= main_eff[5,3],
                           stat_simple_HPNP_LPNP_wald=main_eff[5,5], pval_simple_HPNP_LPNP_wald = main_eff[5,6],
                           estim_simple_LPP_LPNP = main_eff[6,2], se_simple_LPP_LPNP = main_eff[6,3],
                           stat_simple_LPP_LPNP_wald=main_eff[6,5], pval_simple_LPP_LPNP_wald = main_eff[6,6],
                           estim_simple_HPP_HPNP= main_eff[7,2], se_simple_HPP_HPNP = main_eff[7,3],
                           stat_simple_HPP_HPNP_wald=main_eff[7,5], pval_simple_HPP_HPNP_wald = main_eff[7,6],
                           model_comp_stat = model_comp[1],model_comp_df = model_comp[2],  
                           model_comp_pval = model_comp[3],
                           conv_opt = sum_glmer$optinfo$conv$opt,
                           conv_lme4 = ifelse(is.null(sum_glmer$optinfo$conv$lme4$code),0,-1),
                           warn_flag=warn_flag,
                           singular_fit = singular_fit,
                           resid = resid(fit_glmer)
                           
                  ))
                  
                }

write.csv(DE_re,"output/DE_re_quare2.csv")

message("End of Part1: Interaction model, Aripo")


