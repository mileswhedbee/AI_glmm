rm(list = ls())
library(tidyverse)
load(file = "data/res_gene_overlap")
read.csv(file = "output/aripo_de_res.csv") -> de_res_aripo
read.csv(file = "output/quare_de_res.csv") -> de_res_quare

de_res_aripo %>% filter(pval_int_wald <= 0.05) -> de_int

# filter for RES data from sig list
gene_overlap.res[which(is.element(gene_overlap.res$Id, de_int$id)),] -> res_int

hpp.out = 0
hpnp.out = 0
lpp.out = 0
lpnp.out = 0
for (id in sort(unique(res_int$Id))) {
  res_int %>% filter(Id == id & Pop == "high-predation" & Rearing == "P") -> tmp.hpp
  res_int %>% filter(Id == id & Pop == "high-predation" & Rearing == "NP") -> tmp.hpnp
  res_int %>% filter(Id == id & Pop == "low-predation" & Rearing == "P") -> tmp.lpp
  res_int %>% filter(Id == id & Pop == "low-predation" & Rearing == "NP") -> tmp.lpnp

  if(nrow(tmp.hpp) > 0) {
    hpp.out = append(hpp.out, mean(tmp.hpp$Edit_level))
  } else {hpp.out = append(hpp.out, 0)}
  
  if(nrow(tmp.hpnp) > 0) {
    hpnp.out = append(hpnp.out, mean(tmp.hpnp$Edit_level))
  } else {hpnp.out = append(hpnp.out, 0)}
  
  if(nrow(tmp.lpp) > 0) {
    lpp.out = append(lpp.out, mean(tmp.lpp$Edit_level))
  } else {lpp.out = append(lpp.out, 0)}
  
  if(nrow(tmp.lpnp) > 0) {
    lpnp.out = append(lpnp.out, mean(tmp.lpnp$Edit_level))
  } else {lpnp.out = append(lpnp.out, 0)}
}

plot.df = data.frame(Id = sort(rep(unique(res_int$Id), times = 4)), Pop = rep(c("HP","HP","LP","LP"), times = 20),
                     Rearing = rep(c("P","NP","P","NP"), times = 20))


plot.df %>% filter(Pop == "HP" & Rearing == "P") -> tmp
tmp$Mean = hpp.out[2:21]

plot.df %>% filter(Pop == "HP" & Rearing == "NP") -> tmp2
tmp2$Mean = hpnp.out[2:21]

plot.df %>% filter(Pop == "LP" & Rearing == "P") -> tmp3
tmp3$Mean = lpp.out[2:21]

plot.df %>% filter(Pop == "LP" & Rearing == "NP") -> tmp4
tmp4$Mean = lpnp.out[2:21]

rbind(tmp, tmp2, tmp3, tmp4) -> plot.df

plot.df$Group = paste(plot.df$Pop, plot.df$Rearing, sep = "")

de_int %>% filter(pval_simple_HPP_HPNP_wald <= 0.05) -> plot.df.hp
de_int %>% filter(pval_simple_LPP_LPNP_wald <= 0.05) -> plot.df.lp

plot.df[which(is.element(plot.df$Id, plot.df.hp$id)),] -> plot.df.hp
plot.df.hp %>% filter(Group == "HPP" | Group == "HPNP") -> plot.df.hp
plot.df[which(is.element(plot.df$Id, plot.df.lp$id)),] -> plot.df.lp 
plot.df.lp %>% filter(Group == "LPP" | Group == "LPNP") -> plot.df.lp

ggplot(plot.df.hp) +
  geom_point(aes(x = Rearing, y = Mean, group = Id, color = Pop)) + 
  geom_line(aes(x = Rearing, y = Mean, group = Id, color = Pop)) +
  geom_point(data = plot.df.lp, aes(x = Rearing, y = Mean, group = Id, color = Pop)) + 
  geom_line(data = plot.df.lp, aes(x = Rearing, y = Mean, group = Id, color = Pop)) +
  theme_minimal()

ggplot(plot.df.hp) +
  geom_point(aes(x = Rearing, y = Mean, group = Id, color = Id)) + 
  geom_line(aes(x = Rearing, y = Mean, group = Id, color = Id)) +
  geom_point(data = plot.df.lp, aes(x = Rearing, y = Mean, group = Id, color = Id)) + 
  geom_line(data = plot.df.lp, aes(x = Rearing, y = Mean, group = Id, color = Id)) 

plot.df.hp[which(is.element(plot.df.hp$Id, plot.df.lp$Id)), ] -> plot.df.hp.share
plot.df.lp[which(is.element(plot.df.lp$Id, plot.df.hp$Id)), ] -> plot.df.lp.share
 
ggplot(plot.df.hp.share) +
  geom_point(aes(x = Rearing, y = Mean, group = Id, color = Pop)) + 
  geom_line(aes(x = Rearing, y = Mean, group = Id, color = Pop)) +
  geom_point(data = plot.df.lp.share, aes(x = Rearing, y = Mean, group = Id, color = Pop)) + 
  geom_line(data = plot.df.lp.share, aes(x = Rearing, y = Mean, group = Id, color = Pop)) 
