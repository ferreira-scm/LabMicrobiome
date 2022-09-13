######

fil <- function(ps){
    x = phyloseq::taxa_sums(ps)
    # abundance filtering at 0.001%
    keepTaxa = (x / sum(x) > 0.00001)
#    keepTaxa = (x / sum(x) > 0.0005)
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
ps <- subset_taxa(ps, genus%in%"Eimeria")
ps <-aggregate_taxa(ps, level="genus")
ps <- psmelt(ps)
ps <- ps[ps$Abundance>0,]
ps <- ps[ps$Genome_copies_ngDNA>0,]
ps$logA <- log(ps$Abundance)
ps$logGC <- log(ps$Genome_copies_ngDNA)
return(ps)
}

subPSslv <- function(ps) {
ps <- transform_sample_counts(ps, function(x) x / sum(x))
ps <- subset_taxa(ps, Family%in%"Eimeriorina")
ps <-aggregate_taxa(ps, level="Family")
ps <- psmelt(ps)
ps$logA <- log(1+ps$Abundance)
ps$logGC <- log(1+ps$Genome_copies_ngDNA)
return(ps)
}


p_tss <- function(df, lb, name){
ggplot(df, aes(x=logA, y=logGC))+
    geom_jitter(shape=21, position=position_jitter(0.002), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    scale_fill_brewer(palette="Spectral")+
    geom_smooth(method = "lm", se=FALSE, na.rm=TRUE) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                                          parse = TRUE) +  
    ylab("Genome copies/ng DNA (log)")+
    xlab(paste(name, "Eimeria (log)", sep=" "))+
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

sam <- data.frame(sample_data(ps))
PSeimf <- subset_taxa(ps, genus%in%"Eimeria")

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
#sdt$logGC <- log(sdt$Genome_copies_gFaeces)
sdt$logGC <- log(sdt$Genome_copies_ngDNA)
sdt$logEimeriaSums <- log(sdt$EimeriaSums)
sdt$logTotalSums <- log(sdt$TotalSums)

sdta <- sdt[sdt$EimeriaSums>0,]
#sdta <- sdta[sdta$Genome_copies_gFaeces>0,]
sdta <- sdta[sdta$Genome_copies_ngDNA>0,]

print(cor.test(sdta$logGC, sdta$logEimeriaSums, method="pearson"))
#print(cor.test(sdt$logOPG, sdt$logEimeriaSums, method="pearson"))
#print(cor.test(sdt$logTotalSums, sdt$logEimeriaSums, method="pearson"))
# Linear models

Alm <- lm(logGC ~ logEimeriaSums+logTotalSums, sdta)
Alm.0 <- lm(logGC ~ logEimeriaSums, sdta)
summary(Alm)

anova(Alm, Alm.0)

a <- ggplot(sdta, aes(y=logGC, x=logEimeriaSums))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    scale_fill_brewer(palette="Spectral")+
    geom_smooth(aes(y=logGC, x=logEimeriaSums),method = "lm", se=TRUE, colour="black") +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),parse = TRUE) +  
    ylab("Genome copies ngDNA(log)")+
    xlab("Eimeria (log)")+
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
bPSeimf <- subset_taxa(ppPS, genus%in%"Eimeria")

#create total sums and Eimeria sums data frame
bdf <- data.frame(sample_sums(otu_table(ppPS)))
bdf$labels <- rownames(bdf)
bdf$eim <-(sample_sums(otu_table(bPSeimf)))
names(bdf) <- c("Fil_TotalSums", "labels", "FilEimeriaSums")

#merge
sdt <- merge(bdf,sdt, by="labels", all=TRUE) 
#correlation tests
sdt$logFilEimeriaSums <- log(sdt$FilEimeriaSums)

sdtb <- sdt[sdt$FilEimeriaSums>0,]
#sdtb <- sdtb[sdtb$Genome_copies_gFaeces>0,]
sdtb <- sdtb[sdtb$Genome_copies_ngDNA>0,]

print(cor.test(sdtb$logGC, sdtb$logFilEimeriaSums, method="pearson"))
#print(cor.test(sdt$logOPG, sdt$logFilEimeriaSums, method="pearson"))

# Linear models
Blm <- lm(logGC ~ logFilEimeriaSums+logTotalSums, sdtb)
Blm0 <- lm(logGC ~ logFilEimeriaSums, sdtb)
summary(Blm)
anova(Blm, Blm0)

# plotting with abundance and prevalence filter correlation

b <- ggplot(sdtb, aes(y=logGC, x=logFilEimeriaSums))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    geom_smooth(method = "lm", se=TRUE, colour="black") +
    scale_fill_brewer(palette="Spectral")+
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                                          parse = TRUE) +  
    ylab("Genome copies/ngDNA (log)")+
    xlab("Eimeria (log)")+
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
cPSeimf <- subset_taxa(PSTSS, genus%in%"Eimeria")
cPSeimf <-aggregate_taxa(cPSeimf, level="genus")
    
#create total sums and Eimeria sums data frame
df <-data.frame(sample_sums(otu_table(cPSeimf)))
df$labels <- rownames(df)
names(df) <- c("TSS_Eim", "labels")

#merge
sdt <- merge(df, sdt, by="labels", all=TRUE) 

#correlation tests
sdt$logTSS_Eim <- log(sdt$TSS_Eim)

sdtc <- sdt[sdt$TSS_Eim>0,]
#sdtc <- sdtc[sdtc$Genome_copies_gFaeces>0,]
sdtc <- sdtc[sdtc$Genome_copies_ngDNA>0,]

print(cor.test(sdtc$logGC, sdtc$logTSS_Eim, method="pearson"))
#print(cor.test(sdt$logOPG, sdt$logTSS_Eim, method="pearson"))
# Linear models
Clm <- lm(logGC ~ logTSS_Eim, sdtc)

summary(Clm)

# plotting TSS correlation
c <-ggplot(sdtc, aes(x=logTSS_Eim, y=logGC))+
    geom_jitter(shape=21, position=position_jitter(0.002),
                size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    scale_fill_brewer(palette="Spectral")+
    geom_smooth(method = "lm", se=TRUE, na.rm=TRUE, colour="black") +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                                          parse = TRUE) +  
    ylab("Genome copies/ngDNA (log)")+
    xlab("Eimeria (log)")+
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
sdtd <- sdt[sdt$FilEimeriaSums>0,]
sdtd <- sdtd[sdtd$Genome_copies_gFaeces>0,]
sdtd <- sdtd[sdtd$Genome_copies_ngDNA>0,]
#print(cor.test(sdt$logGC, sdt$logREL_Eim, method="pearson"))
print(cor.test(sdtd$logGC, sdtd$REL_Eim, method="pearson"))

# Linear models
Dlm <- lm(logGC ~ REL_Eim, sdtd)
summary(Dlm)

# plotting REL correlation
d <-ggplot(sdtd, aes(x=REL_Eim, y=logGC))+
    geom_jitter(shape=21, position=position_jitter(0.002), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    scale_fill_brewer(palette="Spectral")+
    geom_smooth(method = "lm", se=TRUE, na.rm=TRUE, colour="black") +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                                          parse = TRUE) +  
    ylab("Genome copies/ngDNA (log)")+
    xlab("Eimeria")+
    ggtitle("Relative log expression")+
        labs(tag= "d)")+
    theme_bw()+
    theme(text = element_text(size=12),
#          axis.title.x = element_blank(),
          legend.position= "bottom")+
    guides(fill=guide_legend(nrow=2, byrow=TRUE))

########### centered log ratio
PSclr = microbiome::transform(ppPS, transform="clr")
ePSeimf <- subset_taxa(PSclr, genus%in%"Eimeria")

df <-data.frame(sample_sums(otu_table(ePSeimf)))
df$labels <- rownames(df)
names(df) <- c("clr_Eim", "labels")

#merge
sdt <- merge(df, sdt, by="labels", all=TRUE) 

sdte <- sdt[sdt$TSS_Eim>0,] #removing zeros, this is correct
sdte <- sdte[sdte$Genome_copies_gFaeces>0,]

#correlation tests
print(cor.test(sdte$logGC, sdte$clr_Eim, method="pearson"))
# Linear models
Elm <- lm(logGC ~ clr_Eim, sdte)

summary(Elm)

# plotting CLR correlation
e <-ggplot(sdte, aes(x=clr_Eim, y=logGC))+
    geom_jitter(shape=21, position=position_jitter(0.002),
                size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    scale_fill_brewer(palette="Spectral")+
    geom_smooth(method = "lm", se=TRUE, na.rm=TRUE, colour="black") +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                                          parse = TRUE) +  
    ylab("Genome copies gFaeces (log)")+
    xlab("Eimeria")+
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

#sdt$ACS_Eim <- sdt$TSS_Eim*sdt$DNA_g_feces
sdt$ACS_Eim <- sdt$TSS_Eim*sdt$Total_DNA
sdt$logACS_Eim <- log(sdt$ACS_Eim)

sdtf <- sdt[sdt$ACS_Eim>0,]
sdtf <- sdtf[sdtf$Genome_copies_gFaeces>0,]
sdtf <- sdtf[sdtf$Genome_copies_ngDNA>0,]

#correlation tests
print(cor.test(sdtf$logGC, sdtf$logACS_Eim, method="pearson"))

# Linear models
Flm <- lm(logGC ~ logACS_Eim, sdtf)
summary(Flm)

# plotting ACS correlation
f <-ggplot(sdt, aes(x=logACS_Eim, y=logGC))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    scale_fill_brewer(palette="Spectral")+
    geom_smooth(method = "lm", se=TRUE, na.rm=TRUE, color="black") +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE) +  
    ylab("Genome copies/ngDNA (log)")+
    xlab("Eimeria (log)")+
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

### terrible coding here now...
names(sdt)

cor.df <- sdt[,c("logGC", "logEimeriaSums", "logFilEimeriaSums", "logTSS_Eim", "logACS_Eim", "clr_Eim", "REL_Eim")]

cor.df$REL_Eim <- cor.df$REL_Eim*-1

cor.df1 <- na.omit(cor.df[!is.infinite(cor.df$logEimeriaSums),])

#cor.df1 <- na.omit(cor.df1[cor.df1$logFilEimeriaSums>0,])
cor.df1 <- cor.df1[!is.infinite(cor.df1$logGC),]
cor.df1$logGC
print(cocor(~logGC + logEimeriaSums | logGC + logFilEimeriaSums, data = cor.df1,
            test = c("hittner2003", "zou2007")))

cor.df2 <- na.omit(cor.df[!is.infinite(cor.df$logTSS_Eim),])
cor.df2 <- cor.df2[!is.infinite(cor.df2$logGC),]
print(cocor(~logGC + logEimeriaSums | logGC + logTSS_Eim, data = cor.df2, test = c("hittner2003", "zou2007")))

cor.df3 <- na.omit(cor.df[!is.infinite(cor.df$logFilEimeriaSums),])
cor.df3 <- cor.df3[!is.infinite(cor.df3$logGC),]
print(cocor(~logGC + logEimeriaSums | logGC + REL_Eim, data = cor.df3,test = c("hittner2003", "zou2007")))

cor.df4 <- na.omit(cor.df[!is.infinite(cor.df$logTSS_Eim),])
cor.df4 <- cor.df4[!is.infinite(cor.df4$logGC),]
print(cocor(~logGC + logEimeriaSums | logGC + clr_Eim, data = cor.df4,
            test = c("hittner2003", "zou2007")))

cor.df5 <- na.omit(cor.df[!is.infinite(cor.df$logACS_Eim),])
cor.df5 <- cor.df5[!is.infinite(cor.df5$logGC),]
print(cocor(~logGC + logEimeriaSums | logGC + logACS_Eim, data = cor.df5,
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
PSeimf <- subset_taxa(ps, genus%in%"Eimeria")

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
    xlab(paste(name, "Eimeria (log1+)", sep=" "))+
    ggtitle("Raw counts unfiltered")+
    labs(tag= "a)")+
#    annotate(geom="text", x=12, y=7, label="Spearman rho=0.93, p<0.001")+
    theme_bw()+
    theme(text = element_text(size=16))

#lazy fix
ppPS <- ps

#### using relative abundance
PSTSS = transform_sample_counts(ppPS, function(x) x / sum(x))
cPSeimf <- subset_taxa(PSTSS, genus%in%"Eimeria")

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
    xlab(paste(name, "Eimeria (log1+)", sep=" "))+
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
#sdt$logREL_Eim <- log(1+sdt$REL_Eim)

print(cor.test(sdt$logGC, sdt$REL_Eim, method="pearson"))
print(cor.test(sdt$logOPG, sdt$REL_Eim, method="pearson"))

# Linear models
Dlm <- lm(logGC ~ REL_Eim, sdt)
summary(Dlm)

# plotting REL correlation
d <-ggplot(sdt, aes(x=REL_Eim, y=logGC))+
    geom_jitter(shape=21, position=position_jitter(0.002), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    geom_smooth(method = "lm", se=FALSE, na.rm=TRUE) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                                          parse = TRUE) +  
    ylab("Genome copies gFaeces log(1+)")+
    xlab("SA Eimeria")+
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
    xlab("SA Eimeria log(1+)")+
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


Plotting_cor_MA.l <- function(ps, ps.f, name, dir){
#No filters
library("ggpmisc")

sam <- data.frame(sample_data(ps))
PSeimf <- subset_taxa(ps, genus%in%"Eimeria")

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
#sdt$logGC <- log(sdt$Genome_copies_gFaeces)
sdt$logGC <- log(sdt$Genome_copies_ngDNA)
sdt$logEimeriaSums <- log(sdt$EimeriaSums)
sdt$logTotalSums <- log(sdt$TotalSums)

sdta <- sdt[sdt$EimeriaSums>0,]
#sdta <- sdta[sdta$Genome_copies_gFaeces>0,]
sdta <- sdta[sdta$Genome_copies_ngDNA>0,]

print(cor.test(sdta$logGC, sdta$logEimeriaSums, method="pearson"))
#print(cor.test(sdt$logOPG, sdt$logEimeriaSums, method="pearson"))
#print(cor.test(sdt$logTotalSums, sdt$logEimeriaSums, method="pearson"))
# Linear models

Alm <- lm(logGC ~ logEimeriaSums+logTotalSums, sdta)
Alm.0 <- lm(logGC ~ logEimeriaSums, sdta)
summary(Alm)

anova(Alm, Alm.0)

a <- ggplot(sdta, aes(y=logGC, x=logEimeriaSums))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    scale_fill_brewer(palette="Spectral")+
    geom_smooth(aes(y=logGC, x=logEimeriaSums),method = "lm", se=TRUE, colour="black") +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),parse = TRUE) +  
    ylab("Genome copies ngDNA(log)")+
    xlab("Eimeria (log)")+
    ggtitle("Raw counts unfiltered")+
    labs(tag= "a)")+
#    annotate(geom="text", x=12, y=7, label="Spearman rho=0.93, p<0.001")+
    theme_bw()+
    theme(text = element_text(size=12),
#          axis.title.x = element_blank(),
          legend.position= "bottom")+
    guides(fill=guide_legend(nrow=2, byrow=TRUE))

##############now we filter
ppPS <- ps.f

#now we make the data frame
bPSeimf <- subset_taxa(ppPS, genus%in%"Eimeria")

#create total sums and Eimeria sums data frame
bdf <- data.frame(sample_sums(otu_table(ppPS)))
bdf$labels <- rownames(bdf)
bdf$eim <-(sample_sums(otu_table(bPSeimf)))
names(bdf) <- c("Fil_TotalSums", "labels", "FilEimeriaSums")

#merge
sdt <- merge(bdf,sdt, by="labels", all=TRUE) 
#correlation tests
sdt$logFilEimeriaSums <- log(sdt$FilEimeriaSums)

sdtb <- sdt[sdt$FilEimeriaSums>0,]
#sdtb <- sdtb[sdtb$Genome_copies_gFaeces>0,]
sdtb <- sdtb[sdtb$Genome_copies_ngDNA>0,]

print(cor.test(sdtb$logGC, sdtb$logFilEimeriaSums, method="pearson"))
#print(cor.test(sdt$logOPG, sdt$logFilEimeriaSums, method="pearson"))

# Linear models
Blm <- lm(logGC ~ logFilEimeriaSums+logTotalSums, sdtb)
Blm0 <- lm(logGC ~ logFilEimeriaSums, sdtb)
summary(Blm)
anova(Blm, Blm0)

# plotting with abundance and prevalence filter correlation

b <- ggplot(sdtb, aes(y=logGC, x=logFilEimeriaSums))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    geom_smooth(method = "lm", se=TRUE, colour="black") +
    scale_fill_brewer(palette="Spectral")+
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                                          parse = TRUE) +  
    ylab("Genome copies/ngDNA (log)")+
    xlab("Eimeria (log)")+
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
cPSeimf <- subset_taxa(PSTSS, genus%in%"Eimeria")
cPSeimf <-aggregate_taxa(cPSeimf, level="genus")
    
#create total sums and Eimeria sums data frame
df <-data.frame(sample_sums(otu_table(cPSeimf)))
df$labels <- rownames(df)
names(df) <- c("TSS_Eim", "labels")

#merge
sdt <- merge(df, sdt, by="labels", all=TRUE) 

#correlation tests
sdt$logTSS_Eim <- log(sdt$TSS_Eim)

sdtc <- sdt[sdt$TSS_Eim>0,]
#sdtc <- sdtc[sdtc$Genome_copies_gFaeces>0,]
sdtc <- sdtc[sdtc$Genome_copies_ngDNA>0,]

print(cor.test(sdtc$logGC, sdtc$logTSS_Eim, method="pearson"))
#print(cor.test(sdt$logOPG, sdt$logTSS_Eim, method="pearson"))
# Linear models
Clm <- lm(logGC ~ logTSS_Eim, sdtc)

summary(Clm)

# plotting TSS correlation
c <-ggplot(sdtc, aes(x=logTSS_Eim, y=logGC))+
    geom_jitter(shape=21, position=position_jitter(0.002),
                size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    scale_fill_brewer(palette="Spectral")+
    geom_smooth(method = "lm", se=TRUE, na.rm=TRUE, colour="black") +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                                          parse = TRUE) +  
    ylab("Genome copies/ngDNA (log)")+
    xlab("Eimeria (log)")+
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
sdtd <- sdt[sdt$FilEimeriaSums>0,]
sdtd <- sdtd[sdtd$Genome_copies_gFaeces>0,]
sdtd <- sdtd[sdtd$Genome_copies_ngDNA>0,]
#print(cor.test(sdt$logGC, sdt$logREL_Eim, method="pearson"))
print(cor.test(sdtd$logGC, sdtd$REL_Eim, method="pearson"))

# Linear models
Dlm <- lm(logGC ~ REL_Eim, sdtd)
summary(Dlm)

# plotting REL correlation
d <-ggplot(sdtd, aes(x=REL_Eim, y=logGC))+
    geom_jitter(shape=21, position=position_jitter(0.002), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    scale_fill_brewer(palette="Spectral")+
    geom_smooth(method = "lm", se=TRUE, na.rm=TRUE, colour="black") +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                                          parse = TRUE) +  
    ylab("Genome copies/ngDNA (log)")+
    xlab("Eimeria")+
    ggtitle("Relative log expression")+
        labs(tag= "d)")+
    theme_bw()+
    theme(text = element_text(size=12),
#          axis.title.x = element_blank(),
          legend.position= "bottom")+
    guides(fill=guide_legend(nrow=2, byrow=TRUE))

########### centered log ratio
PSclr = microbiome::transform(ppPS, transform="clr")
ePSeimf <- subset_taxa(PSclr, genus%in%"Eimeria")

df <-data.frame(sample_sums(otu_table(ePSeimf)))
df$labels <- rownames(df)
names(df) <- c("clr_Eim", "labels")

#merge
sdt <- merge(df, sdt, by="labels", all=TRUE) 

sdte <- sdt[sdt$TSS_Eim>0,] #removing zeros, this is correct
sdte <- sdte[sdte$Genome_copies_gFaeces>0,]

#correlation tests
print(cor.test(sdte$logGC, sdte$clr_Eim, method="pearson"))
# Linear models
Elm <- lm(logGC ~ clr_Eim, sdte)

summary(Elm)

# plotting CLR correlation
e <-ggplot(sdte, aes(x=clr_Eim, y=logGC))+
    geom_jitter(shape=21, position=position_jitter(0.002),
                size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    scale_fill_brewer(palette="Spectral")+
    geom_smooth(method = "lm", se=TRUE, na.rm=TRUE, colour="black") +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                                          parse = TRUE) +  
    ylab("Genome copies gFaeces (log)")+
    xlab("Eimeria")+
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

#sdt$ACS_Eim <- sdt$TSS_Eim*sdt$DNA_g_feces
sdt$ACS_Eim <- sdt$TSS_Eim*sdt$Total_DNA
sdt$logACS_Eim <- log(sdt$ACS_Eim)

sdtf <- sdt[sdt$ACS_Eim>0,]
sdtf <- sdtf[sdtf$Genome_copies_gFaeces>0,]
sdtf <- sdtf[sdtf$Genome_copies_ngDNA>0,]

#correlation tests
print(cor.test(sdtf$logGC, sdtf$logACS_Eim, method="pearson"))

# Linear models
Flm <- lm(logGC ~ logACS_Eim, sdtf)
summary(Flm)

# plotting ACS correlation
f <-ggplot(sdt, aes(x=logACS_Eim, y=logGC))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    scale_fill_brewer(palette="Spectral")+
    geom_smooth(method = "lm", se=TRUE, na.rm=TRUE, color="black") +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE) +  
    ylab("Genome copies/ngDNA (log)")+
    xlab("Eimeria (log)")+
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

### terrible coding here now...
names(sdt)

cor.df <- sdt[,c("logGC", "logEimeriaSums", "logFilEimeriaSums", "logTSS_Eim", "logACS_Eim", "clr_Eim", "REL_Eim")]

cor.df$REL_Eim <- cor.df$REL_Eim*-1

cor.df1 <- na.omit(cor.df[!is.infinite(cor.df$logEimeriaSums),])

#cor.df1 <- na.omit(cor.df1[cor.df1$logFilEimeriaSums>0,])
cor.df1 <- cor.df1[!is.infinite(cor.df1$logGC),]
cor.df1$logGC
print(cocor(~logGC + logEimeriaSums | logGC + logFilEimeriaSums, data = cor.df1,
            test = c("hittner2003", "zou2007")))

cor.df2 <- na.omit(cor.df[!is.infinite(cor.df$logTSS_Eim),])
cor.df2 <- cor.df2[!is.infinite(cor.df2$logGC),]
print(cocor(~logGC + logEimeriaSums | logGC + logTSS_Eim, data = cor.df2, test = c("hittner2003", "zou2007")))

cor.df3 <- na.omit(cor.df[!is.infinite(cor.df$logFilEimeriaSums),])
cor.df3 <- cor.df3[!is.infinite(cor.df3$logGC),]
print(cocor(~logGC + logEimeriaSums | logGC + REL_Eim, data = cor.df3,test = c("hittner2003", "zou2007")))

cor.df4 <- na.omit(cor.df[!is.infinite(cor.df$logTSS_Eim),])
cor.df4 <- cor.df4[!is.infinite(cor.df4$logGC),]
print(cocor(~logGC + logEimeriaSums | logGC + clr_Eim, data = cor.df4,
            test = c("hittner2003", "zou2007")))

cor.df5 <- na.omit(cor.df[!is.infinite(cor.df$logACS_Eim),])
cor.df5 <- cor.df5[!is.infinite(cor.df5$logGC),]
print(cocor(~logGC + logEimeriaSums | logGC + logACS_Eim, data = cor.df5,
            test = c("hittner2003", "zou2007")))
    
#png(filename=paste(dir, name, "COR.png", sep=""),
#    width =14, height = 14, units = "in", res= 300)
#fCor
#dev.off()
}

