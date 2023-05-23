rm(list = ls())
library(dplyr)

read.csv(file = "output/DE_re_aripo2.csv") -> de_res
read.csv(file = "output/aripo_sr.mat") -> aripo_sr
read.csv(file = "output/DE_re_quare2.csv") -> de_res_quare
read.csv(file = "output/quare_sr.mat") -> quare_sr
load(file = "data/res_gene_overlap")

# N fish filter
N_fish_filt = 8
gene_overlap.res %>% filter(N_fish >= N_fish_filt) -> res_filt
res_filt %>% filter(Drainage == "Aripo") -> aripo_de
res_filt %>% filter(Drainage == "Quare") -> quare_de
# Sr sum filt
sr_filt = 10
aripo_de %>% filter(SR_sum >= sr_filt) -> aripo_de
quare_de %>% filter(SR_sum >= sr_filt) -> quare_de
# match de_res to filtered aripo_de
tmp.out = data.frame()
for (Id in unique(aripo_de$Id)){
  de_res %>% filter(id == Id) -> tmp
  tmp.out = rbind(tmp.out, tmp)
} 
de_res = tmp.out

tmp.out = data.frame()
for (Id in unique(quare_de$Id)){
  de_res_quare %>% filter(id == Id) -> tmp
  tmp.out = rbind(tmp.out, tmp)
} 
de_res_quare = tmp.out

i = 1
j = 1
resid_sq = 0
resid_sum = 0
for (resid in i:nrow(de_res)) {
  de_res[i,44:79] -> tmp
  t(tmp) -> tmp
  colnames(tmp) = "Resid"
  i = i + 1
  for (resid2 in tmp[1:36,]) {
    resid2^2 -> resid_sq[j]
    j = j + 1
  }
  append(resid_sum, sum(resid_sq)) -> resid_sum
  j = 1
}
resid_sum[2:length(resid_sum)] -> resid_sum
de_res$Resid_sqSum = resid_sum

hist(de_res$Resid_sqSum, breaks = 50)

de_res %>% filter(singular_fit == TRUE) -> de_res_sing
de_res %>% filter(singular_fit == FALSE) -> de_res_nonSing

hist(de_res_nonSing$Resid_sqSum, breaks = 50)
hist(de_res_sing$Resid_sqSum, breaks = 50)

hist(de_res$pval_pop_wald)
hist(de_res$pval_rear_wald)
hist(de_res$pval_int_wald)

# check high pvals
de_res %>% filter(pval_int_wald >= 0.9) -> tmp
tmp.out = data.frame()
for (id in tmp$id) {
  gene_overlap.res %>% filter(Id == id) -> tmp2
  tmp.out = rbind(tmp.out, tmp2[1,])
}
hist(tmp.out$N_fish, main = "Aripo pvals >= 0.9")
hist(tmp.out$SR_sum, breaks = 30, main = "Aripo pvals >= 0.9")
hist(tmp.out$DP_sum, breaks = 30, main = "Aripo pvals >= 0.9")

library(fdrtool)

# calc fdr as qvals - aripo
fdrtool(de_res$pval_int_wald, statistic = "pvalue", plot = F) -> tmp
tmp$qval -> de_res$qval_int
fdrtool(de_res$pval_pop_wald, statistic = "pvalue", plot = F) -> tmp
tmp$qval -> de_res$qval_pop
fdrtool(de_res$pval_rear_wald, statistic = "pvalue", plot = F) -> tmp
tmp$qval -> de_res$qval_rear
# calc fdr as qvals - quare
fdrtool(de_res_quare$pval_int_wald, statistic = "pvalue", plot = F) -> tmp
tmp$qval -> de_res_quare$qval_int
fdrtool(de_res_quare$pval_pop_wald, statistic = "pvalue", plot = F) -> tmp
tmp$qval -> de_res_quare$qval_pop
fdrtool(de_res_quare$pval_rear_wald, statistic = "pvalue", plot = F) -> tmp
tmp$qval -> de_res_quare$qval_rear

hist(de_res$Resid_sqSum, breaks = 50, main = "Sum of squares residuals - aripo")

library(ggplot2)

# filter by fdr <= 5% - aripo
de_res %>% filter(qval_pop <= 0.05) -> de_res_popSig
de_res %>% filter(qval_rear <= 0.05) -> de_res_rearSig
de_res %>% filter(qval_int <= 0.05) -> de_res_intSig

# filter by fdr <= 5% - quare 
de_res_quare %>% filter(qval_pop <= 0.05) -> de_res_popSig_quare
de_res_quare %>% filter(qval_rear <= 0.05) -> de_res_rearSig_quare
de_res_quare %>% filter(qval_int <= 0.05) -> de_res_intSig_quare

tmp.out = data.frame()
for (id in de_res_intSig$id) {
  aripo_de %>% filter(Id == id) -> tmp
  tmp.out = rbind(tmp.out, tmp)
}

for (id in unique(tmp.out$Id)) {
  
  tmp.out %>% filter(Id == id) -> tmp
  
  print(ggplot(tmp) +
    geom_col(aes(reorder(Group, Supporting_reads), Supporting_reads, fill = Pop)) +
    facet_grid(Pop ~ Rearing, margins = F, scales = "free") +
    ggtitle(id)) +
    xlab("Individual sample") +
    theme(axis.text.x = element_blank())
}

write.csv(de_res, file = "output/aripo_de_res.csv", row.names = F)
write.csv(de_res_quare, file = "output/quare_de_res.csv", row.names = F)
