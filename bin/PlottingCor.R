######

fil <- function(ps){
    x = phyloseq::taxa_sums(ps)
    # abundance filtering at 0.005%
    keepTaxa = (x / sum(x) > 0.00005)
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
    otu_table(ps) <- otu_table(ps)*sample_data(ps)$Total_DNA
    ps <- subset_taxa(ps, Genus%in%"g__Eimeria")
    ps <-aggregate_taxa(ps, level="Genus")
    ps <- psmelt(ps)
    ps <- ps[ps$Abundance>0,]
    ps <- ps[ps$Genome_copies_ngDNA>0,]
    ps$logA <- log(ps$Abundance)
    ps$logGC <- log(ps$Genome_copies_ngDNA)
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
    theme_bw()+
    theme(text = element_text(size=12),
#          axis.title.x = element_blank(),
          legend.position= "bottom")+
    guides(fill=guide_legend(nrow=2, byrow=TRUE))
}


sensit <- function(Eim_nf){
# GC- Eim -
    tneg <- summary(sample_sums(Eim_nf@otu_table)==0&Eim_nf@sam_data$Genome_copies_gFaeces==0)[[3]]

# GC+ Eim -
    fneg <- summary(sample_sums(Eim_nf@otu_table)==0&Eim_nf@sam_data$Genome_copies_gFaeces>0)[[3]]

# GC+ Eim+
    tpos <- summary(sample_sums(Eim_nf@otu_table)>0&Eim_nf@sam_data$Genome_copies_gFaeces>0)[[3]]

# GC- Eim+
    fpos <- try(summary(sample_sums(Eim_nf@otu_table)>0&Eim_nf@sam_data$Genome_copies_gFaeces==0)[[3]])

    if (class(fpos)=="try-error"){
        fpos <- "0"
    }

    print(fpos)
# GC+ GC-
gcpos <- summary(Eim_nf@sam_data$Genome_copies_gFaeces==0)[[2]]
gcneg <- summary(Eim_nf@sam_data$Genome_copies_gFaeces==0)[[3]]

# Eim + Eim -
epos <- summary(sample_sums(Eim_nf@otu_table)==0)[[2]]
eneg <- summary(sample_sums(Eim_nf@otu_table)==0)[[3]]

e.s <- matrix(as.numeric(c(tneg, fpos, fneg, tpos)), ncol=2, byrow=TRUE)
colnames(e.s) <- c("Eim-","Eim+")
rownames(e.s) <- c("GC-", "GC+")
margin1 <- margin.table(e.s, margin=1)
margin2 <- margin.table(e.s, margin=2)

#sensitivity
print("Sensitivity:")
print(e.s[2,2]/margin1[2]*100)
#specificity
print("Specificity:")
print(e.s[1,1]/margin1[1]*100)
#ppv
print("positive predictive value")
print(e.s[2,2]/margin2[2]*100)
#npv
print("negative predictive value")
print(e.s[1,1]/margin2[1]*100)
}


Plotting_cor_MA.l <- function(ps, ps.f, name, dir){
#No filters
library("ggpmisc")

sam <- data.frame(sample_data(ps))
PSeimf <- subset_taxa(ps, Genus%in%"g__Eimeria")

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

#print(cor.test(sdta$logGC, sdta$logEimeriaSums, method="pearson"))
##############now we filter
ppPS <- ps.f

#now we make the data frame
bPSeimf <- subset_taxa(ppPS, Genus%in%"g__Eimeria")

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

sdtb <- sdtb[!is.na(sdtb$Genome_copies_ngDNA),]

b.c <- cor.test(sdtb$logGC, sdtb$logFilEimeriaSums, method="pearson")
print(b.c)

coul <- c("#F6E8C3","#DFC27D","#BF812D","#F1B6DA", "#C51B7D", "#DE77AE","#01665E", "#35978F", "#80CDC1","#C7EAE5")

# plotting with abundance and prevalence filter correlation
b <- ggplot(sdtb, aes(y=logGC, x=logFilEimeriaSums))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    scale_fill_manual(values=coul)+
    ylab("Eimeria genome copies (log)")+
    xlab("Eimeria ASV abundance (log)")+
    ggtitle("No normalization")+
    annotate(geom="text", x=min(sdtb$logFilEimeriaSums), y=max(sdtb$logGC), hjust=0.05, label=paste("Spearman rho=", round(b.c$estimate, 2), ", p<0.001", sep=""))+
    theme_bw(base_size=10)+
    theme(legend.position= "none")

################################################################
#### using relative abundance
PSTSS = transform_sample_counts(ppPS, function(x) x / sum(x))
cPSeimf <- subset_taxa(PSTSS, Genus%in%"g__Eimeria")
cPSeimf <-aggregate_taxa(cPSeimf, level="Genus")
    
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
sdtc <- sdtc[!is.na(sdtc$Genome_copies_ngDNA),]

c.c <- cor.test(sdtc$logGC, sdtc$logTSS_Eim, method="pearson")
print(c.c)

# plotting TSS correlation
c <-ggplot(sdtc, aes(x=logTSS_Eim, y=logGC))+
    geom_jitter(shape=21, position=position_jitter(0.002),
                size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    scale_fill_manual(values=coul)+
    ylab("Eimeria genome copies (log)")+
    xlab("Eimeria relative abundance (log)")+
    ggtitle("Total sum scaling")+
    annotate(geom="text", x=min(sdtc$logTSS_Eim), y=max(sdtc$logGC), hjust=0.05, label=paste("Spearman rho=", round(c.c$estimate, 2), ", p<0.001", sep=""))+
    theme_bw(base_size=10)+
    theme(legend.position= "none")

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
sdtd <- sdtd[!is.na(sdtd$Genome_copies_ngDNA),]

d.c <- cor.test(sdtd$logGC, sdtd$REL_Eim, method="pearson")

print(d.c)

# plotting REL correlation
d <-ggplot(sdtd, aes(x=REL_Eim, y=logGC))+
    geom_jitter(shape=21, position=position_jitter(0.002), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    scale_fill_manual(values=coul)+
    ylab("Eimeria genome copies (log)")+
    xlab("Eimeria ASV abundance")+
    annotate(geom="text", x=max(sdtd$REL_Eim), y=max(sdtd$logGC), hjust=0.9, label=paste("Spearman rho=", round(d.c$estimate, 2), ", p<0.001", sep="")) +
    ggtitle("Relative log expression")+
    theme_bw(base_size=10)+
    theme(legend.position= "none")

########### centered log ratio
PSclr = microbiome::transform(ppPS, transform="clr")
ePSeimf <- subset_taxa(PSclr, Genus%in%"g__Eimeria")

df <-data.frame(sample_sums(otu_table(ePSeimf)))
df$labels <- rownames(df)
names(df) <- c("clr_Eim", "labels")

#merge
sdt <- merge(df, sdt, by="labels", all=TRUE) 
sdte <- sdt[sdt$TSS_Eim>0,] #removing zeros, this is correc
sdte <- sdte[sdte$Genome_copies_gFaeces>0,]
sdte <- sdte[!is.na(sdte$Genome_copies_gFaeces),]

#correlation tests
e.c <- (cor.test(sdte$logGC, sdte$clr_Eim, method="pearson"))
print(e.c)

# plotting CLR correlation
e <-ggplot(sdte, aes(x=clr_Eim, y=logGC))+
    geom_jitter(shape=21, position=position_jitter(0.002),
                size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    scale_fill_manual(values=coul)+
    ylab("Eimeria genome copies (log)")+
    xlab("Eimeria ASV abundance")+
    ggtitle("Centered log-ratio")+
    annotate(geom="text", x=min(sdte$clr_Eim), y=max(sdte$logGC), hjust=0.05, label=paste("Spearman rho=", round(e.c$estimate, 2), ", p<0.001", sep=""))+
    theme_bw(base_size=10)+
    theme(legend.position= "none")

#### experimental quantification
#ABsolute Count Scaling: scaled to DNA/g/faeces

#sdt$ACS_Eim <- sdt$TSS_Eim*sdt$DNA_g_feces
sdt$ACS_Eim <- sdt$TSS_Eim*sdt$Total_DNA
sdt$logACS_Eim <- log(sdt$ACS_Eim)

sdtf <- sdt[sdt$ACS_Eim>0,]
sdtf <- sdtf[sdtf$Genome_copies_gFaeces>0,]
sdtf <- sdtf[sdtf$Genome_copies_ngDNA>0,]
sdtf <- sdtf[!is.na(sdtf$Genome_copies_ngDNA),]

#correlation tests
f.c <- (cor.test(sdtf$logGC, sdtf$logACS_Eim, method="pearson"))
print(f.c)

# plotting ACS correlation
f <-ggplot(sdtf, aes(x=logACS_Eim, y=logGC))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    scale_fill_manual(values=coul)+
    ylab("Eimeria genome copies (log)")+
    xlab("Eimeria ASV abundance (log)")+
    annotate(geom="text", x=min(sdtf$logACS_Eim), y=max(sdte$logGC), hjust=0.05, label=paste("Spearman rho=", round(f.c$estimate, 2), ", p<0.001", sep=""))+
    ggtitle("Sample DNA concentration scaling")+
    theme_bw(base_size=10)+
    theme(legend.position= "none")

############### rarefaction
rare <- rarefy_even_depth(ppPS, rngseed=1234)
gPSeimf <- subset_taxa(rare, Genus%in%"g__Eimeria")

#create total sums and Eimeria sums data frame
gdf <- data.frame(sample_sums(otu_table(rare)))
gdf$labels <- rownames(gdf)
gdf$eimf <-(sample_sums(otu_table(gPSeimf)))
names(gdf) <- c("Fil_TotalSums_rare", "labels", "Eim_rare")

#merge
sdt <- merge(gdf,sdt, by="labels", all=TRUE) 
#correlation tests
sdt$logEim_rare <- log(sdt$Eim_rare)

sdtg <- sdt[sdt$Eim_rare>0,]
sdtg <- sdtg[sdtg$Genome_copies_ngDNA>0,]
sdtg <- sdtg[!is.na(sdtg$Genome_copies_ngDNA),]
g.c <- (cor.test(sdtg$logGC, sdtg$logEim_rare, method="pearson"))
print(g.c)

# plotting with rarefied abundance
g <- ggplot(sdtg, aes(y=logGC, x=logEim_rare))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
#    geom_smooth(method = "lm", se=TRUE, colour="black") +
    scale_fill_manual(values=coul)+
#    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),parse = TRUE) +  
    ylab("Eimeria genome copies (log)")+
    xlab("Eimeria ASV abundance (log)")+
    ggtitle("Rarefaction")+
#    labs(tag= "f)")+
    annotate(geom="text", x=min(sdtg$logEim_rare), y=max(sdtg$logGC), hjust=0.05, label=paste("Spearman rho=", round(g.c$estimate, 2), ", p<0.001", sep=""))+
    theme_bw(base_size=10)+
theme(legend.position= "none")

g
    
# save plots of what we have so far
fCor <-plot_grid(b, c, f, d, e, g, 
              align="vh",
              labels="auto",
              nrow=3)

legend <- get_legend(g+
          guides(fill=guide_legend(nrow=1, byrow=TRUE))+
          theme(legend.position="top"))

fCor.l <- plot_grid(legend,fCor, nrow=2, rel_heights=c(0.1, 1))

fCor.l  <<- fCor.l
    
ggplot2::ggsave(file=paste(dir, name, "COR.pdf", sep=""), fCor.l, width = 8, height = 12, dpi = 600)
ggplot2::ggsave(file=paste(dir, name, "COR.png", sep=""), fCor.l, width = 8, height = 12, dpi = 600)

ggplot2::ggsave(file=paste(dir, name, "SCS_COR.png", sep=""), f, width = 5, height = 5, dpi = 600)

### terrible coding here now...
names(sdt)

cor.df <- sdt[,c("logGC", "logEimeriaSums", "logFilEimeriaSums", "logTSS_Eim", "logACS_Eim", "clr_Eim", "REL_Eim", "logEim_rare")]

cor.df$REL_Eim <- cor.df$REL_Eim*-1

cor.df1 <- na.omit(cor.df[!is.infinite(cor.df$logEimeriaSums),])

#cor.df1 <- na.omit(cor.df1[cor.df1$logFilEimeriaSums>0,])
cor.df2 <- na.omit(cor.df[!is.infinite(cor.df$logTSS_Eim),])
cor.df2 <- cor.df2[!is.infinite(cor.df2$logGC),]
print(cocor(~logGC + logFilEimeriaSums | logGC + logTSS_Eim, data = cor.df2, test = c("hittner2003", "zou2007")))

cor.df3 <- na.omit(cor.df[!is.infinite(cor.df$logFilEimeriaSums),])
cor.df3 <- cor.df3[!is.infinite(cor.df3$logGC),]
print(cocor(~logGC + logFilEimeriaSums | logGC + REL_Eim, data = cor.df3,test = c("hittner2003", "zou2007")))

cor.df4 <- na.omit(cor.df[!is.infinite(cor.df$logTSS_Eim),])
cor.df4 <- cor.df4[!is.infinite(cor.df4$logGC),]
print(cocor(~logGC + logFilEimeriaSums | logGC + clr_Eim, data = cor.df4,
            test = c("hittner2003", "zou2007")))

cor.df5 <- na.omit(cor.df[!is.infinite(cor.df$logACS_Eim),])
cor.df5 <- cor.df5[!is.infinite(cor.df5$logGC),]
print(cocor(~logGC + logFilEimeriaSums | logGC + logACS_Eim, data = cor.df5,
            test = c("hittner2003", "zou2007")))

cor.df6 <- na.omit(cor.df[!is.infinite(cor.df$logEim_rare),])
cor.df6 <- cor.df6[!is.infinite(cor.df6$logGC),]
print(cocor(~logGC + logFilEimeriaSums | logGC + logEim_rare, data = cor.df6,
            test = c("hittner2003", "zou2007")))

#png(filename=paste(dir, name, "COR.png", sep=""),
#    width =14, height = 14, units = "in", res= 300)
#fCor
#dev.off()

}


