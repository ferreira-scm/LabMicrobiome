## Code for 
## 1) Standard curve of qPCR Eimeria

### We should have ONE model able to predict for different species and
### cyclers, IF they are significantly different

## 2) Determine Eimeria amount for infection experiment samples

### then make the predictions based on the data (including cycler and
### species) and the models (with those factors if significant)

## only load packages that are needed!!!
library(ggpubr)
library(rcompanion)
library(dplyr)
library(gridExtra)
library(rstatix)

##Load data
if(!exists("sample.data")){
    source("bin/1_Data_preparation.R")
}
##Standard curves
data.std<- read.csv("data/Eimeria_quantification_Std_Curve_data.csv")
data.std%>%
    dplyr::mutate(Genome_copies= Oocyst_count*8)-> data.std

##Define numeric and factor variables 
num.vars <- c("Ct", "Ct_mean", "Sd_Ct", "Qty", "Qty_mean", "Sd_Qty", "Oocyst_count", "Feces_weight", "Qubit", "NanoDrop", "Beads_weight", "Tm", "Genome_copies")
fac.vars <- c("Well", "Sample.Name", "Detector", "Task",  "Std_series","Date", "Operator", "Cycler", "Parasite", "Sample_type", "Extraction")  

## as.numeric alone will likely fail if stringsAsfactors is TRUE! 
data.std[, num.vars] <- apply(data.std[, num.vars], 2,
                                 function (x) as.numeric(as.character(x)))
data.std[, fac.vars] <- apply(data.std[, fac.vars], 2, as.factor)

##Correct zero in NTC with not-detected results 
data.std$Ct[data.std$Ct == 0] <- NA
data.std$Ct_mean[data.std$Ct_mean == 0 & data.std$Task == "NTC"] <- NA
data.std$Sd_Ct[data.std$Sd_Ct == 0] <- NA

##Correct labels to have them homogeneous 
data.std$Sample.Name<- gsub(pattern = " ", replacement = "_", x = data.std$Sample.Name)

## Select just standards data
## Estimate the number of genome copies per ng of gDNA
data.std.lm<- subset(data.std, Task== "Standard") ## Select just data from standards 
data.std.lm %>% 
  select(Sample.Name, Task, Ct, Cycler, Oocyst_count, Parasite, Genome_copies)%>%
  dplyr::mutate(Oocyst_DNA= Oocyst_count*(3.8E-4))%>% ##Estimation of DNA (ng) derived from Oocyst
  dplyr::mutate(DNA_PCR= Oocyst_DNA/30)%>% ##DNA (ng) in PCR considering 1uL from a stock of 30uL
  ##Considering that 1 ng of Eimeria gDNA is equivalent to 2.11E4 genome copies
  dplyr::mutate(Genome_copies_ngDNA= (2.11E4)*DNA_PCR)-> data.std.lm  

##Inter-sample variation and spiked samples
data.unk<-read.csv("data/Eimeria_quantification_Sample_data.csv")
  
##Define numeric and factor variables 
num.vars2 <- c("Ct", "Ct_mean", "Sd_Ct", "Qty", "Qty_mean", "Sd_Qty", "Oocyst_count", "Feces_weight", "Qubit", "NanoDrop", "Beads_weight", "Tm", 
                "Oocyst_1", "Oocyst_2", "Oocyst_3", "Oocyst_4", "Oocyst_5", "Oocyst_6", "Oocyst_7", "Oocyst_8", "Dilution_factor", "Volume", "Sporulated")
  fac.vars2 <- c("Well", "Sample.Name", "Detector", "Task", "Date", "Operator", "Cycler", "Parasite", "Sample_type", "Extraction", "Strain")  

  data.unk[, num.vars2] <- apply(data.unk[, num.vars2], 2,
                               function (x) as.numeric(as.character(x)))
  data.unk[, fac.vars2] <- apply(data.unk[, fac.vars2], 2, as.factor)

##Information from intersample variation experiment 
data.unk%>%
  dplyr::select(Sample.Name, Task, Ct, Cycler, Parasite, Sample_type, Extraction, Tm, Qubit,
                Oocyst_1, Oocyst_2, Oocyst_3, Oocyst_4, Oocyst_5,Oocyst_6, Oocyst_7, Oocyst_8, 
                Sporulated, Dilution_factor, Volume, Strain)%>%
  filter(Sample_type=="Oocysts" & Task=="Unknown")%>%
  dplyr::group_by(Sample.Name)%>%
  dplyr::mutate(N= n())%>%
  dplyr::mutate(Total_oocysts= (sum(Oocyst_1, Oocyst_2, Oocyst_3, Oocyst_4, Oocyst_5,Oocyst_6,
                                    Oocyst_7, Oocyst_8))/N)%>%
  ##Concentration of Oocyst in the solution
  dplyr::mutate(Oocyst_count= (((Total_oocysts*10000)/8)*Dilution_factor))%>%
  ##Concentration of sporulated oocyst in the solution
  dplyr::mutate(Sporulated_count= (((Sporulated*10000)/8)*Dilution_factor))%>%
  dplyr::mutate(Sporulation_rate= (Sporulated_count/Oocyst_count)*100)%>%
  dplyr::mutate(Sporulation_rate= as.numeric(Sporulation_rate))-> data.unk.lm

##Spiked data 
data.unk%>%
  select(Sample.Name, Task, Ct, Cycler,Parasite, Sample_type, Feces_weight, Extraction, Oocyst_count, Qubit, NanoDrop)%>%
  filter(Sample_type=="Feces" & Task=="Unknown")-> data.spk
  
##Infection experiment
data.inf.exp<-read.csv("data/Eimeria_quantification_Inf_exp_data.csv")
data.inf.exp%>%
    select(Content, Sample, Plate_number, Cq, Melt_Temperature)%>%
    dplyr::rename(Ct= Cq, labels= Sample, Task= Content, Tm= Melt_Temperature)%>%
    dplyr::mutate(Cycler= "BioRad")-> data.inf.exp
  
##Define numeric and factor variables 
num.vars3 <- c("Ct", "Tm")
fac.vars3 <- c("labels", "Task", "Plate_number", "Cycler")  
data.inf.exp[, num.vars3] <- apply(data.inf.exp[, num.vars3], 2,
                               function (x) as.numeric(as.character(x)))
data.inf.exp[, fac.vars3] <- apply(data.inf.exp[, fac.vars3], 2, as.factor)
  

rm(fac.vars, num.vars, fac.vars2, num.vars2, fac.vars3, num.vars3)

####### Standard curves #######
set.seed(2020)
## Compute simple linear models from standards
## "Genome copies modeled by Ct"

##Model 8: Genome copies modeled by Ct and cycler as predictors
lm.SCCyc<- lm(log10(Genome_copies)~Ct+Cycler, data.std.lm, na.action = na.exclude)

######### Infection experiment data############
## Define real positive and negatives based on Tm 
data.inf.exp %>% 
    dplyr::mutate(Infection = case_when(is.na(Tm)  ~ "Negative",
                                        Tm >= 80   ~ "Negative",
                                        Tm < 80 ~ "Positive")) -> data.inf.exp 

##Estimate number of genome copies with qPCR Ct value (Model 8)
data.inf.exp$Genome_copies<- 10^predict(lm.SCCyc, data.inf.exp)

##Summarise genome copies by sample  
data.inf.exp %>%
    select(Genome_copies,labels) %>% # select variables to summarise
    na.omit()%>%
    dplyr::group_by(labels)%>%
    dplyr::summarise_each(funs(Genome_copies_min = min, Genome_copies_q25 = quantile(., 0.25),
                               Genome_copies_median = median, Genome_copies_q75 = quantile(., 0.75), 
                               Genome_copies_max = max, Genome_copies_mean = mean, Genome_copies_sd = sd)) -> Sum.inf
##qPCR data from 236 samples
## Tm values were not sumarised to avoid problems in the function 

##Join summirised data
data.inf.exp<- inner_join(data.inf.exp, Sum.inf, by= "labels")

##Eliminate an unprocessed sample and controls
data.inf.exp%>%
    select(labels, Genome_copies_mean, Infection)%>%
    filter(!labels%in%c("Pos_Ctrl","Neg_Ctrl","FML"))%>% ## Replace NAs in real negative samples to 0 
    dplyr::mutate(Genome_copies_mean= replace_na(Genome_copies_mean, 0))%>%
    ##Get unique labels from qPCR data
    distinct(labels, .keep_all = TRUE)-> data.inf.exp

##Merging Infection experiment oocyst and weight loss data with qPCR data
##Check differences between two dataframes
setdiff(sample.data$labels, data.inf.exp$labels)

##Samples MTU, FJL, EHM, DRT, CEY were not taken 
##Sample CPY was taken but not extracted (Faeces not found in boxes)
##Sample FLM was collected and extracted but not processed for qPCR (DNA not found)
##We end with 235 samples out of the original 242

###Join all the data in the same dataframe
sdt<- left_join(sample.data, data.inf.exp, by="labels") ## Add qPCR data
###Tiny adjustment  
sdt$dpi<- as.factor(sdt$dpi)

sdt%>%
  dplyr::mutate(Genome_copies_ngDNA= Genome_copies_mean/50, ## copies by ng of fecal DNA considering 1uL from 50 ng/uL DNA
                DNA_sample= Conc_DNA*30, ## Estimate total gDNA of sample considering 30uL of elution buffer
                DNA_g_feces= DNA_sample/fecweight_DNA,
                ## Transform it to ng fecal DNA by g of faeces
                Genome_copies_gFaeces= Genome_copies_ngDNA*DNA_g_feces) -> sdt ## Estimate genome copies by g of faeces

##Transform to zero OPGs for DPI 1 and 2 
sdt$OPG[sdt$dpi==1] <- 0
sdt$OPG[sdt$dpi==2] <- 0                 
##Remove dataframes with data not related to the infection experiment data that won't be used in the following scripts
rm(data.std, data.std.lm, data.unk, data.unk.lm, data.spk, data.spk.lm, Sum.inf)