# Script for AffyArray Normalization

#install.packages("http://mbni.org/customcdf/25.0.0/ensg.download/hgu133plus2hsensgcdf_25.0.0.tar.gz", type = "source", repos = NULL)
library(affy)
library(oligo)
library("hgu133plus2hsensgcdf")

setwd("~/Research/Data_Analysis/ASPIS_OMICS/Data/celfiles/")

# Unzip the files and store in unzipped directory
#zip_file_paths <- list.files("TG_GATES_HUMAN_IN_VITRO_FULLDATA/", pattern = ".zip", full.names = TRUE)

#for (zip_file_path in zip_file_paths) {
 # unzip(zip_file_path, exdir = "Unzipped_files/")
#}

rawData <- ReadAffy(cdfname="HGU133Plus2_Hs_ENSG")
normExp <- expresso(rawData, bgcorrect.method="rma", normalize.method="constant", 
                    normalize.param=list(), pmcorrect.method="pmonly", summary.method="medianpolish")

#add column of IDs and normalized data to normDataTable
normDataTable <- exprs(normExp)
normDataTable <- normDataTable[-grep("AFFX",rownames(normDataTable)),]
rownames(normDataTable) <- gsub(pattern="_at", replacement="",rownames(normDataTable))

#output normalised expression data to file
normFileName <- paste("NormData_","TG_Gates_Full_Hum_invitro_ENS",".txt",sep="")
write.table(normDataTable, normFileName, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

save.image("Full_TG_GATES.RData")

