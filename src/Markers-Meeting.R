### Andrew Trister
### Sage Bionetworks
### Seattle, WA
### 20120927

### Look at the relationship between D, rho and gene expression in TCGA
###

require(synapseClient)
require(ggplot2)
require(survival)
require(reshape)
require(scales)
require(HDclassif)
require(limma)
source("./src/ggkm.R")

synapseLogin()


#load in the D/rho data
#swanson.data <- read.delim("./swanson_input_modified.txt", skip=1,header = TRUE,stringsAsFactors = FALSE)

swanson.data <- read.delim("./data/swanson_PINHA.txt",header = TRUE,stringsAsFactors = FALSE)

metadataLoad <- loadEntity('syn673127') #all of the coherent metadata
metadata <- metadataLoad$objects$metadata #extract the R object

#here get rid of the extraneous rows
metadata <- metadata[!duplicated(metadata$bcr_patient_barcode),]
rownames(metadata) <- as.character(metadata$bcr_patient_barcode)

rownames(swanson.data) <- paste("TCGA-",rownames(swanson.data),sep="")

#and now let's just get the metadata for patients for whom we have imaging
temp <- intersect(rownames(swanson.data),rownames(metadata))

swanson.data.culled <- swanson.data[temp,]
metadata.culled <- metadata[temp,]

swanson.data.culled$velocity <- 2*sqrt(swanson.data.culled$PIHNA.D..mm.*swanson.data.culled$PIHNA.rho..mm.)


#make a survival object

metadata.culled$vital_status[metadata.culled$vital_status=="[Not Available]"] <- NA
gbmPat <- c()
gbmPat$survTime <- as.numeric(metadata.culled$days_to_death)
gbmPat$surv <- ifelse(metadata.culled$vital_status=="DECEASED", 1,0)
gbmPat$survTime[ which(metadata.culled$vital_status=="LIVING")] <- as.numeric(metadata.culled$days_to_last_followup)[ which(metadata.culled$vital_status=="LIVING")]

tmpSurv <- Surv(gbmPat$survTime,gbmPat$surv)



ggkm(survfit(tmpSurv ~ swanson.data.culled$velocity>median(swanson.data.culled$velocity)), 
     ystratalabs = (c("Slow", "Fast")), 
     timeby = 365,
     main = "GBM K-M Plot By Growth Velocity")


ggkm(survfit(tmpSurv ~ swanson.data.culled$PIHNA.rho..mm.>median(swanson.data.culled$PIHNA.rho..mm.)), 
     ystratalabs = (c("Low p", "High p")), 
     timeby = 365,
     main = "GBM K-M Plot By p")

ggkm(survfit(tmpSurv ~ swanson.data.culled$PIHNA.D..mm.>median(swanson.data.culled$PIHNA.D..mm.)), 
     ystratalabs = (c("Low D", "High D")), 
     timeby = 365,
     main = "GBM K-M Plot By D")

#now make a nice scatter plot with ggplot

p <- ggplot(swanson.data.culled)

p + geom_point(aes(PIHNA.D..mm.,PIHNA.rho..mm.),size=4) +
  scale_y_log10(breaks=trans_breaks("log10",function(x) 10^x),
                labels=trans_format("log10",math_format(10^.x)),name="rho (1/years)") + 
                  scale_x_log10(breaks=trans_breaks("log10",function(x) 10^x),
                                labels=trans_format("log10",math_format(10^.x)),name="D (mm2/years)")  +
                                  theme_grey(30)                









#download the REMBRANDT set:

rembrandtDataReturn <- loadEntity('syn376920') #REMBRANDT CEL SNM Normalized
rembrandtEset <- rembrandtDataReturn$objects$rembrandt.matrix.matched
rembrandtPat <- rembrandtDataReturn$objects$phenotype.data.matched





#get the genes that are important for GBM



design <- model.matrix(~0+factor(rembrandtPat$Disease))
colnames(design) <- substring(levels(factor(rembrandtPat$Disease)),2,20)
fit <- lmFit(rembrandtEset,design)
fit <- eBayes(fit)
toptable(fit, coef=2)


contrast.matrix <- makeContrasts(GBM-ASTROCYTOMA,GBM-MIXED,GBM-OLIGODENDROGLIOMA, levels =design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
toptable(fit2, coef=2)



results <- decideTests(fit2)
vennDiagram(results)


genes <- rownames(rembrandtEset[which(results[,1]!=0 & results[,2]!=0 & results[,3]!=0),])
genes <- unlist(lapply(genes,function(x){
  sub("_mt","_eg",x)
}))

#now grab the full normalized matrix from TCGA


#'syn313583', 'syn372761' - these are the SNM normalized data
dataReturn <- loadEntity('syn313583')
eset <- exprs(dataReturn$objects$eset)

colnames(eset) <- substring(colnames(eset),1,12)

temp <- intersect(colnames(eset),rownames(metadata.culled))
temp <- temp[-grep("06-0173",temp)] #remove the cystic patient!

#reduce the genes to only the ones that are important in GBM!
eset.culled <- eset[,temp]
metadata.culled <- metadata.culled[temp,]
swanson.data.culled <- swanson.data.culled[temp,]


temp <- intersect(rownames(eset.culled),genes)
eset.culled <- eset.culled[temp,]


design <- model.matrix(~swanson.data.culled$PIHNA.rho..mm.)

fit <- lmFit(eset.culled,design)
fit <- eBayes(fit)
toptable(fit,coef=2)


results <- decideTests(fit,method="global")
table(results[,2])

esetSel <- eset.culled[which(results[,2]!=0),]

heatmap(esetSel[,order(swanson.data.culled$velocity)],scale="none",Colv=NA)