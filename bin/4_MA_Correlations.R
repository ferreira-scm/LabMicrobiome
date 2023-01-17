
library(relaimpo)
library(lme4)
library(MASS)
library(DescTools)
library(cocor)
library(reshape2)
library(microbiome)
library(phyloseq)

source("bin/4_MA_SA_filtering.R")

############################################
##############################################

#### sensitivity, specificity, positive predicted value and negative predictive value
############# sensitivity and specificity

#### for MA
sensit(Eim)

MA.asv1 <- prune_taxa(rownames(tax_table(Eim2))[1], Eim)
#MA.asv12 <- prune_taxa(rownames(tax_table(Eim2))[c(1,2)], Eim)
#MA.asv13 <- prune_taxa(rownames(tax_table(Eim2))[c(1,3)], Eim)
MA.asv2 <- prune_taxa(rownames(tax_table(Eim2))[2], Eim)
#MA.asv3 <- prune_taxa(rownames(tax_table(Eim2))[3], Eim)

sensit(MA.asv1)

sensit(MA.asv2)

#### for SA
sensit(Eim2)

SA.asv1 <- prune_taxa(rownames(tax_table(Eim2))[1], Eim2)
SA.asv12 <- prune_taxa(rownames(tax_table(Eim2))[c(1,2)], Eim2)
SA.asv123 <- prune_taxa(rownames(tax_table(Eim2))[c(1,2,3)], Eim2)
SA.asv13 <- prune_taxa(rownames(tax_table(Eim2))[c(1,3)], Eim2)
SA.asv2 <- prune_taxa(rownames(tax_table(Eim2))[2], Eim2)
SA.asv3 <- prune_taxa(rownames(tax_table(Eim2))[3], Eim2)
SA.asv4 <- prune_taxa(rownames(tax_table(Eim2))[4], Eim2)

sensit(SA.asv1)

sensit(SA.asv2)

sensit(SA.asv12)

sensit(SA.asv13)

sensit(SA.asv123)

sensit(SA.asv3)

sensit(SA.asv4)

# now plotting and doing correlation analysis and comparisons
## for "pooled" MA
#Plotting_cor(ps=all.PS, "MA", dir="fig/MA/")
# Seqf, TSS, RLE, CLR, ACS
#p <- c(0.0208, 0.0134, 0.1037, 0.00001, 0.00001)
#p.adjust(p, method="BH")

### for single amplicon
Plotting_cor_MA.l(ps=sin.PS18S.slv, f.sin18.slv, "SA", dir="fig/SA/")

# tss,rle,clr,acs, rare
p <- c(0.3243, 0.00001, 0.00001, 0.9665, 0.3759)
p.adjust(p, method="BH")

### for MA but individually filtered
Plotting_cor_MA.l(ps=all.PS.slv, ps.f=f.all.lp.slv, "MA_individually_filtered", dir="fig/MA/")

# tss,rle,clr,acs, rare
p <- c(0.0049, 0.00001, 0.0037, 0.0007, 0.0043)
p.adjust(p, method="BH")

### for MA but wang only

Plotting_cor_MA.l(ps=all.PS.l.slv[[37]], ps.f=f.all.l.slv[[37]], "MA_wang_TSS", dir="fig/MA/")

# tss,rle,clr,acs, rare
p <- c(0.0046, 0.00001, 0.7266, 0.0002, 0.0009)

round(p.adjust(p, method="BH"),3)

fCor.l

