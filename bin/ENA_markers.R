
require(devtools)

devtools::load_all("/SAN/Susanas_den/taxAnyEukAmp/")

setwd("/SAN/Susanas_den/gitProj/LabMicrobiome/")

source("/SAN/Susanas_den/taxAnyEukAmp/R/download.R")

library(RCurl)

X18SDownloads <- getENAdownloads("COI", "/SAN/Susanas_den/AmpMarkers/ENA_Marker/")

X18SDownloads

X18Seq <- Biostrings::readDNAStringSet(getFiles(X18SDownloads))

X18Seq
