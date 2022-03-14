#Detect Module

rnaE = read.csv("First/rnaseq.csv")
rnaExpr = rnaE[,-c(1:2)]

rownames(rnaExpr) = rnaE$Gene
cv = apply(rnaExpr,1,sd)/rowMeans(rnaExpr)

rnaExpr.use = rnaExpr[order(cv,decreasing = T),]
rnaExpr.use = scale(t(rnaExpr[1:1500,]))

geneuse <- data.frame(Gene = colnames(rnaExpr.use))

rnaExprN = rnaExpr.use
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

net = blockwiseModules(rnaExprN,maxBlockSize = 5000, power = 8,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE)


table(net$colors)

sizeGrWindow(12,9)

mergedColors = labels2colors(net$colors)
for(i in c(1:6)){
  plotDendroAndColors(net$dendrograms[[i]],mergedColors[net$blockGenes[[i]]],
                    "Module colors",dendroLabels = FALSE,hang = 0.03,addGuide = TRUE,
                    guideHang = 0.05)
}

moduleLabels = net$colors

moduleColors = labels2colors(net$colors)

MEs0 = net$MEs
blockwiseMEs = moduleEigengenes(rnaExprN, moduleColors)$eigengenes

write.csv(moduleColors,"colorcluster.csv",row.names = F)



net2 = blockwiseModules(new_rna,maxBlockSize = 1000, power = 8,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE)


table(net2$colors)

sizeGrWindow(12,9)

mergedColors = labels2colors(net2$colors)
for(i in c(1:6)){
  plotDendroAndColors(net2$dendrograms[[i]],mergedColors[net2$blockGenes[[i]]],
                      "Module colors",dendroLabels = FALSE,hang = 0.03,addGuide = TRUE,
                      guideHang = 0.05)
}
