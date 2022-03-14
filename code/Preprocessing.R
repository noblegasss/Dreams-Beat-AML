# Preprocessing

library(factoextra)
library(tidyverse)
library(glmnet)
library(broom)
library(purrr)
library(WGCNA)

aucs = read.csv("aucs.csv")
clinical_categ_legnd = read.csv("clinical_categorical_legend.csv")
clinical_categ = read.csv("clinical_categorical.csv")
clinical_num = read.csv("clinical_numerical.csv")
dnaseq = read.csv("dnaseq.csv")
response = read.csv("response.csv")
rnaseq = read.csv("rnaseq.csv")

### Clean RNA data

dim(rnaseq)

rnaExpr = as.data.frame(t(rnaseq[,-c(1:2)]))
names(rnaExpr) = rnaseq$Gene
rownames(rnaExpr) = names(rnaseq[,-c(1:2)])

gsg = goodSamplesGenes(rnaExpr,verbose = 3)
gsg$allOK

####16387 is missing####

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(rnaExpr)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(rnaExpr)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  rnaExpr = rnaExpr[gsg$goodSamples, gsg$goodGenes]
}

### Clean Clinical Data

trait = as.data.frame(merge(clinical_categ,clinical_num,by = "lab_id"))

trait[,1] = rownames(rnaExpr)

Sample = data.frame(lab_id=rownames(rnaExpr))

traitRow = merge(index,trait,by="lab_id")
traitAll0 = traitRow[,-1]
rownames(traitAll0) = traitRow[,1]
traitAll = traitAll_reduced1 = traitAll0[,-c(2,4,6,8,11,12,21)]
traitAll_reduced2 = trait[,-c(2,4,6,8,11,12,21)]

## Random Forest method

clinical_all0 = missForest::missForest(traitAll_reduced2[,-1])[[1]]

write.csv(clinical_all0,"clinical_all2.csv")
write.csv(rnaExpr,"RNA.csv")


#Detect Module

rnaExprN = read.csv("RNA.csv")
traitAll = read.csv("clinical_all2.csv")

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(rnaExprN, powerVector = powers, verbose = 5)

sizeGrWindow(9,5)
par(mfrow = c(1,2))
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

## we choose power = 7 based on the plot

net = blockwiseModules(rnaExprN,maxBlockSize = 5000, power = 7,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, 
                       saveTOMs = TRUE,
                       saveTOMFileBase = "TOM")


table(net$colors)

sizeGrWindow(12,9)

mergedColors = labels2colors(net$colors)
for(i in c(1:10)){
  plotDendroAndColors(net$dendrograms[[i]],mergedColors[net$blockGenes[[i]]],
                      "Module colors",dendroLabels = FALSE,hang = 0.03,addGuide = TRUE,
                      guideHang = 0.05)
}

moduleLabels = net$colors

moduleColors = labels2colors(net$colors)

MEs0 = net$MEs
blockwiseMEs = moduleEigengenes(rnaExprN, moduleColors)$eigengenes



nGenes = ncol(rnaExprN);
nSamples = nrow(rnaExprN);
# Recalculate MEs with color labels

MEs = orderMEs(blockwiseMEs)
moduleTraitCor = cor(MEs, traitAll, use = "p",method = "spearman");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(100,50)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traitAll),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               setStdMargins = FALSE,
               cex.text = 0.01,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))



SpecificDx = as.data.frame(traitAll$specificDxAtAcquisition);
names(SpecificDx) = "SpecificDx"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(rnaExprN, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(rnaExprN, SpecificDx, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(SpecificDx), sep="");
names(GSPvalue) = paste("p.GS.", names(SpecificDx), sep="");

moduleTraitCor_0 = as.data.frame(moduleTraitCor[,7])
abs(moduleTraitCor_0)>=0.3
#MEblue

moduleTraitPvalue_0 = as.data.frame(moduleTraitPvalue[,7])
moduleTraitPvalue_0["MEblue",] #1.951396e-06 significant

module = "blue"
column = match(module,modNames)
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for SpecificDx",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

new_rna = rnaExprN[moduleColors=="blue"]

write.csv(new_rna,"rna_blue.csv")

new_rna = read.csv("rna_blue.csv")

rownames(new_rna) = new_rna[,1]

new_rna = new_rna[,-1]
### PCA

rna_pca0 = prcomp(new_rna)
plot(rna_pca0)

summary(rna_pca0)

##pc100

rna_pca = as.data.frame(rna_pca0$x[,c(1:100)])
rna_pca = cbind(lab_id=Sample,rna_pca)

write.csv(rna_pca,"rpca.csv")

