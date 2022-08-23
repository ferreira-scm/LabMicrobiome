######

fil <- function(ps){
x = phyloseq::taxa_sums(ps)
keepTaxa = (x / sum(x) > 0.00001)
summary(keepTaxa)
ps = phyloseq::prune_taxa(keepTaxa, ps)
# plus prevalnce filter at 1%
KeepTaxap <- microbiome::prevalence(ps)>0.01
ps <- phyloseq::prune_taxa(KeepTaxap, ps)
# subset samples based on total read count (500 reads)
#ps <- phyloseq::subset_samples(ps, phyloseq::sample_sums(ps) > 500)
ps <- phyloseq::prune_samples(sample_sums(ps)>100, ps)
ps
}

subPS <- function(ps) {
ps <- transform_sample_counts(ps, function(x) x / sum(x))
ps <- subset_taxa(ps, family%in%"Eimeriidae")
ps <-aggregate_taxa(ps, level="family")
ps <- psmelt(ps)
ps$logA <- log(1+ps$Abundance)
ps$logGC <- log(1+ps$Genome_copies_gFaeces)
return(ps)
}

subPSslv <- function(ps) {
ps <- transform_sample_counts(ps, function(x) x / sum(x))
ps <- subset_taxa(ps, Family%in%"Eimeriorina")
ps <-aggregate_taxa(ps, level="Family")
ps <- psmelt(ps)
ps$logA <- log(1+ps$Abundance)
ps$logGC <- log(1+ps$Genome_copies_gFaeces)
return(ps)
}


p_tss <- function(df, lb, name){
ggplot(df, aes(x=logA, y=logGC))+
    geom_jitter(shape=21, position=position_jitter(0.002), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    scale_fill_brewer(palette="Spectral")+
    geom_smooth(method = "lm", se=FALSE, na.rm=TRUE) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                                          parse = TRUE) +  
    ylab("Genome copies gFaeces log(1+)")+
    xlab(paste(name, "Eimeriidae (log1+)", sep=" "))+
    ggtitle(name)+
    labs(tag=lb)+
                                        #        annotate(geom="text", x=13, y=0.5, label="Spearman rho=0.94, p<0.001")+
        theme_bw()+
    theme(text = element_text(size=12),
#          axis.title.x = element_blank(),
          legend.position= "bottom")+
    guides(fill=guide_legend(nrow=2, byrow=TRUE))
}


Plotting_cor <- function(ps, name, dir){
#No filters
library("ggpmisc")

set.seed(500)
    
sam <- data.frame(sample_data(ps))
PSeimf <- subset_taxa(ps, family%in%"Eimeriidae")

#create total sums and Eimeria sums data frame
df <- data.frame(sample_sums(otu_table(ps)))
df$labels <- rownames(df)
eimf <-as.data.frame(sample_sums(PSeimf))
eimf$labels <- rownames(eimf)
names(eimf) <- c("EimeriaSums", "labels")
names(df) <- c("TotalSums", "labels")

#merge
df <- merge(df,eimf, by="labels", all=FALSE) 
sdt <- merge(df,sam, by="labels")

#correlation tests
sdt$logOPG <- log(1+sdt$OPG)
sdt$logGC <- log(1+sdt$Genome_copies_gFaeces)
sdt$logEimeriaSums <- log(1+sdt$EimeriaSums)
sdt$logTotalSums <- log(1+sdt$TotalSums)

print(cor.test(sdt$logGC, sdt$logEimeriaSums, method="pearson"))
#print(cor.test(sdt$logOPG, sdt$logEimeriaSums, method="pearson"))
#print(cor.test(sdt$logTotalSums, sdt$logEimeriaSums, method="pearson"))
# Linear models
Alm <- lm(logGC ~ logEimeriaSums , sdt)
summary(Alm)

a <-ggplot(sdt, aes(y=logGC, x=logEimeriaSums))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    scale_fill_brewer(palette="Spectral")+
    geom_smooth(method = "lm", se=FALSE) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                                          parse = TRUE) +  
    ylab("Genome copies gFaeces(log 1+)")+
    xlab("Eimeriidae (log1+)")+
    ggtitle("Raw counts unfiltered")+
    labs(tag= "a)")+
#    annotate(geom="text", x=12, y=7, label="Spearman rho=0.93, p<0.001")+
    theme_bw()+
    theme(text = element_text(size=12),
#          axis.title.x = element_blank(),
          legend.position= "bottom")+
    guides(fill=guide_legend(nrow=2, byrow=TRUE))

##############now we filter
ppPS <- fil(ps)

#now we make the data frame
bPSeimf <- subset_taxa(ppPS, family%in%"Eimeriidae")

#create total sums and Eimeria sums data frame
bdf <- data.frame(sample_sums(otu_table(ppPS)))
bdf$labels <- rownames(bdf)
bdf$eim <-(sample_sums(otu_table(bPSeimf)))
names(bdf) <- c("Fil_TotalSums", "labels", "FilEimeriaSums")

#merge
sdt <- merge(bdf,sdt, by="labels", all=TRUE) 

#correlation tests
sdt$logFilEimeriaSums <- log(1+sdt$FilEimeriaSums)

print(cor.test(sdt$logGC, sdt$logFilEimeriaSums, method="pearson"))
#print(cor.test(sdt$logOPG, sdt$logFilEimeriaSums, method="pearson"))

# Linear models
Blm <- lm(logGC ~ logFilEimeriaSums, sdt)

# plotting with abundance and prevalence filter correlation

b <-ggplot(sdt, aes(y=logGC, x=logFilEimeriaSums))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    geom_smooth(method = "lm", se=FALSE) +
    scale_fill_brewer(palette="Spectral")+
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                                          parse = TRUE) +  
    ylab("Genome copies gFaeces(log 1+)")+
    xlab("Eimeriidae (log1+)")+
    ggtitle("Raw counts filtered")+
    labs(tag= "b)")+
#    annotate(geom="text", x=12, y=7, label="Spearman rho=0.93, p<0.001")+
    theme_bw()+
    theme(text = element_text(size=12),
#          axis.title.x = element_blank(),
          legend.position= "bottom")+
    guides(fill=guide_legend(nrow=2, byrow=TRUE))

################################################################
#### using relative abundance
PSTSS = transform_sample_counts(ppPS, function(x) x / sum(x))
cPSeimf <- subset_taxa(PSTSS, family%in%"Eimeriidae")
cPSeimf <-aggregate_taxa(cPSeimf, level="family")
    
#create total sums and Eimeria sums data frame
df <-data.frame(sample_sums(otu_table(cPSeimf)))
df$labels <- rownames(df)
names(df) <- c("TSS_Eim", "labels")

#merge
sdt <- merge(df, sdt, by="labels", all=TRUE) 

#correlation tests
sdt$logTSS_Eim <- log(1+sdt$TSS_Eim)

print(cor.test(sdt$logGC, sdt$logTSS_Eim, method="pearson"))
#print(cor.test(sdt$logOPG, sdt$logTSS_Eim, method="pearson"))
# Linear models
Clm <- lm(logGC ~ logTSS_Eim, sdt)

# plotting TSS correlation
c <-ggplot(sdt, aes(x=logTSS_Eim, y=logGC))+
    geom_jitter(shape=21, position=position_jitter(0.002),
                size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    scale_fill_brewer(palette="Spectral")+
    geom_smooth(method = "lm", se=FALSE, na.rm=TRUE) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                                          parse = TRUE) +  
    ylab("Genome copies gFaeces log(1+)")+
    xlab("Eimeriidae (log1+)")+
    ggtitle("Total sum scaling")+
        labs(tag= "c)")+
#        annotate(geom="text", x=13, y=0.5, label="Spearman rho=0.94, p<0.001")+
    theme_bw()+
    theme(text = element_text(size=12),
#          axis.title.x = element_blank(),
          legend.position= "bottom")+
    guides(fill=guide_legend(nrow=2, byrow=TRUE))

###############################################################
#### using Relative log expression
library(edgeR)
source("bin/edgeR_phyloseq.R")
edgePS <- phyloseq_to_edgeR(bPSeimf)

edgePS$samples$labels <- rownames(edgePS$samples)
df <- (edgePS$samples)

df$group <- NULL
df$lib.size <- NULL
names(df) <- c("REL_Eim", "labels")

#merge
sdt <- merge(df, sdt, by="labels", all=TRUE) 

#correlation tests
sdt$logREL_Eim <- log(1+sdt$REL_Eim)

#print(cor.test(sdt$logGC, sdt$logREL_Eim, method="pearson"))
print(cor.test(sdt$logGC, sdt$REL_Eim, method="pearson"))

# Linear models
Dlm <- lm(logGC ~ logREL_Eim, sdt)
#summary(Dlm)

# plotting REL correlation
d <-ggplot(sdt, aes(x=REL_Eim, y=logGC))+
    geom_jitter(shape=21, position=position_jitter(0.002), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    scale_fill_brewer(palette="Spectral")+
    geom_smooth(method = "lm", se=FALSE, na.rm=TRUE) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                                          parse = TRUE) +  
    ylab("Genome copies gFaeces log(1+)")+
    xlab("Eimeriidae")+
    ggtitle("Relative log expression")+
        labs(tag= "d)")+
    theme_bw()+
    theme(text = element_text(size=12),
#          axis.title.x = element_blank(),
          legend.position= "bottom")+
    guides(fill=guide_legend(nrow=2, byrow=TRUE))

########### centered log Ration
PSclr = microbiome::transform(ppPS, transform="clr")
ePSeimf <- subset_taxa(PSclr, family%in%"Eimeriidae")

df <-data.frame(sample_sums(otu_table(ePSeimf)))
df$labels <- rownames(df)
names(df) <- c("clr_Eim", "labels")

#merge
sdt <- merge(df, sdt, by="labels", all=TRUE) 

#correlation tests
print(cor.test(sdt$logGC, sdt$clr_Eim, method="pearson"))
# Linear models
Elm <- lm(logGC ~ clr_Eim, sdt)

# plotting CLR correlation
e <-ggplot(sdt, aes(x=clr_Eim, y=logGC))+
    geom_jitter(shape=21, position=position_jitter(0.002),
                size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    scale_fill_brewer(palette="Spectral")+
    geom_smooth(method = "lm", se=FALSE, na.rm=TRUE) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                                          parse = TRUE) +  
    ylab("Genome copies gFaeces log(1+)")+
    xlab("Eimeriidae")+
    ggtitle("Centered log-ratio")+
        labs(tag= "e)")+
#        annotate(geom="text", x=13, y=0.5, label="Spearman rho=0.94, p<0.001")+
    theme_bw()+
    theme(text = element_text(size=12),
#          axis.title.x = element_blank(),
          legend.position= "bottom")+
    guides(fill=guide_legend(nrow=2, byrow=TRUE))


#### experimental quantification
#ABsolute Count Scaling: scaled to DNA/g/faeces

sdt$ACS_Eim <- sdt$TSS_Eim*sdt$DNA_g_feces
sdt$logACS_Eim <- log(1+sdt$ACS_Eim)

#correlation tests
print(cor.test(sdt$logGC, sdt$logACS_Eim, method="pearson"))
# Linear models
Flm <- lm(logGC ~ logACS_Eim, sdt)
#summary(Flm)

# plotting REL correlation
f <-ggplot(sdt, aes(x=logACS_Eim, y=logGC))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    scale_fill_brewer(palette="Spectral")+
    geom_smooth(method = "lm", se=FALSE) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE) +  
    ylab("Genome copies gFaeces log(1+)")+
    xlab("Eimeriidae (log1+)")+
    ggtitle("Absolute count scaling")+
        labs(tag= "f)")+
    theme_bw()+
    theme(text = element_text(size=12),
#          axis.title.x = element_blank(),
          legend.position= "bottom")+
    guides(fill=guide_legend(nrow=2, byrow=TRUE))

# save plots of what we have so far
plot_grid(a,b,c,d,e,f) -> fCor

ggplot2::ggsave(file=paste(dir, name, "COR.pdf", sep=""), fCor, width = 12, height = 10, dpi = 600)

ggplot2::ggsave(file=paste(dir, name, "COR.png", sep=""), fCor, width = 12, height = 10, dpi = 600)

ggplot2::ggsave(file=paste(dir, name, "ACS_COR.png", sep=""), e, width = 5, height = 5, dpi = 600)

#names(sdt)
cor.df <- sdt[,c(69, 68, 70,72, 73, 3, 2, 76)]
#quick fix on RLE, to force it into a positive correlation
cor.df$REL_Eim <- cor.df$REL_Eim*-1
    
print(cocor(~logGC + logEimeriaSums | logGC + logFilEimeriaSums, data = cor.df,
            test = c("hittner2003", "zou2007")))

print(cocor(~logGC + logEimeriaSums | logGC + logTSS_Eim, data = cor.df,
            test = c("hittner2003", "zou2007")))

print(cocor(~logGC + logEimeriaSums | logGC + REL_Eim, data = cor.df,
            test = c("hittner2003", "zou2007")))

print(cocor(~logGC + logEimeriaSums | logGC + clr_Eim, data = cor.df,
            test = c("hittner2003", "zou2007")))

print(cocor(~logGC + logEimeriaSums | logGC + logACS_Eim, data = cor.df,
            test = c("hittner2003", "zou2007")))
    
#png(filename=paste(dir, name, "COR.png", sep=""),
#    width =14, height = 14, units = "in", res= 300)
#fCor
#dev.off()
}


                                        # Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
      }


NoFilPlotting_cor <- function(ps, name, dir){
#No filters
library("ggpmisc")

sam <- data.frame(sample_data(ps))
PSeimf <- subset_taxa(ps, family%in%"Eimeriidae")

#create total sums and Eimeria sums data frame
df <- data.frame(sample_sums(otu_table(ps)))
df$labels <- rownames(df)
eimf <-as.data.frame(sample_sums(PSeimf))
eimf$labels <- rownames(eimf)
names(eimf) <- c("EimeriaSums", "labels")
names(df) <- c("TotalSums", "labels")

#merge
df <- merge(df,eimf, by="labels", all=FALSE) 
sdt <- merge(df,sam, by="labels")

#correlation tests
sdt$logOPG <- log(1+sdt$OPG)
sdt$logGC <- log(1+sdt$Genome_copies_gFaeces)
sdt$logEimeriaSums <- log(1+sdt$EimeriaSums)
sdt$logTotalSums <- log(1+sdt$TotalSums)

print(cor.test(sdt$logGC, sdt$logEimeriaSums, method="pearson"))
print(cor.test(sdt$logOPG, sdt$logEimeriaSums, method="pearson"))
print(cor.test(sdt$logTotalSums, sdt$logEimeriaSums, method="pearson"))

# Linear models
Alm <- lm(logGC ~ logEimeriaSums , sdt)
summary(Alm)

a <-ggplot(sdt, aes(y=logGC, x=logEimeriaSums))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    geom_smooth(method = "lm", se=FALSE) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                                          parse = TRUE) +  
    ylab("Genome copies gFaeces(log 1+)")+
    xlab(paste(name, "Eimeriidae (log1+)", sep=" "))+
    ggtitle("Raw counts unfiltered")+
    labs(tag= "a)")+
#    annotate(geom="text", x=12, y=7, label="Spearman rho=0.93, p<0.001")+
    theme_bw()+
    theme(text = element_text(size=16))

#lazy fix
ppPS <- ps

#### using relative abundance
PSTSS = transform_sample_counts(ppPS, function(x) x / sum(x))
cPSeimf <- subset_taxa(PSTSS, family%in%"Eimeriidae")

#create total sums and Eimeria sums data frame
bdf <- data.frame(sample_sums(otu_table(ppPS)))
bdf$labels <- rownames(bdf)
bdf$eim <-(sample_sums(otu_table(bPSeimf)))
names(bdf) <- c("Fil_TotalSums", "labels", "FilEimeriaSums")


df <-data.frame(sample_sums(otu_table(cPSeimf)))
df$labels <- rownames(df)
names(df) <- c("TSS_Eim", "labels")

#merge
sdt <- merge(df, sdt, by="labels", all=TRUE) 

#correlation tests
sdt$logTSS_Eim <- log(1+sdt$TSS_Eim)

print(cor.test(sdt$logGC, sdt$logTSS_Eim, method="pearson"))
print(cor.test(sdt$logOPG, sdt$logTSS_Eim, method="pearson"))

# Linear models
Clm <- lm(logGC ~ logTSS_Eim, sdt)

# plotting TSS correlation
c <-ggplot(sdt, aes(x=logTSS_Eim, y=logGC))+
    geom_jitter(shape=21, position=position_jitter(0.002), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    geom_smooth(method = "lm", se=FALSE, na.rm=TRUE) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                                          parse = TRUE) +  
    ylab("Genome copies gFaeces log(1+)")+
    xlab(paste(name, "Eimeriidae (log1+)", sep=" "))+
    ggtitle("TSS")+
        labs(tag= "c)")+
#        annotate(geom="text", x=13, y=0.5, label="Spearman rho=0.94, p<0.001")+
        theme_bw()+
    theme(text = element_text(size=16))

#### using Relative log expression
library(edgeR)
source("bin/edgeR_phyloseq.R")
edgePS <- phyloseq_to_edgeR(bPSeimf)

edgePS$samples$labels <- rownames(edgePS$samples)
df <- (edgePS$samples)

df$group <- NULL
df$lib.size <- NULL
names(df) <- c("REL_Eim", "labels")

#merge
sdt <- merge(df, sdt, by="labels", all=TRUE) 

#correlation tests
sdt$logREL_Eim <- log(1+sdt$REL_Eim)

print(cor.test(sdt$logGC, sdt$logREL_Eim, method="pearson"))
print(cor.test(sdt$logOPG, sdt$logREL_Eim, method="pearson"))

# Linear models
Dlm <- lm(logGC ~ logREL_Eim, sdt)
summary(Dlm)

# plotting REL correlation
d <-ggplot(sdt, aes(x=logREL_Eim, y=logGC))+
    geom_jitter(shape=21, position=position_jitter(0.002), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    geom_smooth(method = "lm", se=FALSE, na.rm=TRUE) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                                          parse = TRUE) +  
    ylab("Genome copies gFaeces log(1+)")+
    xlab("SA Eimeriidae")+
    ggtitle("REL")+
        labs(tag= "d)")+
        theme_bw()+
    theme(text = element_text(size=16))

#### experimental quantification
#ABsolute Count Scaling: scaled to DNA/g/faeces

sdt$ACS_Eim <- sdt$TSS_Eim*sdt$DNA_g_feces
sdt$logACS_Eim <- log(1+sdt$ACS_Eim)

#correlation tests
print(cor.test(sdt$logGC, sdt$logACS_Eim, method="pearson"))
print(cor.test(sdt$logOPG, sdt$logACS_Eim, method="pearson"))

# Linear models
Elm <- lm(logGC ~ logACS_Eim, sdt)
summary(Elm)

# plotting REL correlation
e <-ggplot(sdt, aes(x=logACS_Eim, y=logGC))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    geom_smooth(method = "lm", se=FALSE) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                                          parse = TRUE) +  
    ylab("Genome copies gFaeces log(1+)")+
    xlab("SA Eimeriidae log(1+)")+
    ggtitle("ACS")+
        labs(tag= "e)")+
        theme_bw()+
    theme(text = element_text(size=16))

# save plots of what we have so far
plot_grid(a,b,c,d,e) -> fCor

ggplot2::ggsave(file=paste(dir, name, "COR.pdf", sep=""), fCor, width = 14, height = 14, dpi = 600)

ggplot2::ggsave(file=paste(dir, name, "COR.png", sep=""), fCor, width = 14, height = 14, dpi = 600)

#png(filename=paste(dir, name, "COR.png", sep=""),
#    width =14, height = 14, units = "in", res= 300)
#fCor
#dev.off()
}
