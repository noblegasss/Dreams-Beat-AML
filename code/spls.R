library(mixOmics)

gene_response = data.frame(lab_id = index,rnaExpr.use)
gene_response = merge(response,gene_response,by = "lab_id")

gene_response = gene_response %>%
  mutate(survival = ifelse(survival == 1,0,1))

pca.gene = pca(gene_response[,-c(1,2,3)],ncomp = 20)

pca.gene
plot(pca.gene)

plotIndiv(pca.gene,comp=c(1,2),group = gene_response$survival)

spls = splsda(X=gene_response[,-c(1,2,3)],Y = gene_response$survival ,ncomp = 3,keepX = c(100,50))

tune.spls = perf(spls,validation = "Mfold", folds = 3, progressBar = FALSE, nrepeat = 5)
plot(tune.spls, overlay = 'measure', sd=TRUE)

grid.keepX = c(seq(50,500, 5))  
# if you dont understand what this means, type
# grid.keepX  # adjust this grid as necessary for your own data

# this chunk takes ~2 min to run
set.seed(33)  # for reproducible results for this code, remove for your own code
tune.splsda = tune.splsda(X = gene_response[,-c(1,2,3)],
                                Y = gene_response$survival,
                                ncomp = 3,
                                test.keepX = grid.keepX,
                                validation = c('Mfold'),
                                folds = 5,
                                dist = 'mahalanobis.dist', # prediction distance can be chosen according to tune.plsda results
                                nrepeat = 10,
                                progressBar = FALSE)

plot(tune.splsda)

tune.splsda$choice.keepX 

spls = splsda(X=gene_response[,-c(1,2,3)],Y = gene_response$survival ,ncomp = 3,keepX = c(185,210,445))

plotIndiv(spls, ind.names = gene_response$survival, ellipse = TRUE, legend = TRUE)
plotVar(spls, comp =2:3)

color.edge <- color.GreenRed(50)  
# to save as a pdf
network(spls, comp = 1, shape.node = c("rectangle", "rectangle"),
        color.node = c("white", "pink"), color.edge = color.edge)
cim(spls,comp = 1)

plotLoadings(spls, comp = 3, method = 'mean', contrib = 'max',
             size.title = 1, ndisplay = 85, size.name = 0.5, size.legend = 0.3)

spls_gene1 = selectVar(spls)
spls_gene2 = selectVar(spls,comp=2)
spls_gene3 = selectVar(spls,comp=3)

p.spls = predict(spls,newdata = test_rna.final)

features = c("ageAtDiagnosis","FLT3.ITD1","priorMDS1",
             "NPM11" ,"specificDxAtAcquisition12","auccluster","dxAtSpecimenAcquisition2",
             "DNMT3A","SRSF2","WT1","IDH1","STAG2","RUNX1")
features.gene = c(features,spls_gene)

train_rna = data.frame(clin,auccluster = auc.cluster,rnacluster = rna.cluster,m01_de,rna)
train_rna = train_rna[,features.gene]
train_rna = data.frame(index,train_rna)

train_rna.final = merge(train_rna,response,by="lab_id")

#ggplot(data = train_rna.final,aes(x = index, y = ENSG00000079691)) + geom_density()

test_rna.final = data.frame(lab_id = sample,clin.test,auccluster =auc.cluster2,
                            rnacluster = rna.cluster2,m01_de.test,rna2)

train_rna.final = train_rna.final[,-1]

train_rna.final =train_rna.final%>%
  mutate(survival = ifelse(survival==1,0,1))
test_rna.final =test_rna.final%>%
  mutate(survival = ifelse(survival==1,0,1))

splstrain = train_rna.final[,spls_gene1$name]
splstrain = data.frame(train_rna.final[,c(1:3)],splstrain)

rf.sur = rfsrc(Surv(overallSurvival,survival)~.,
               data = splstrain[,-1],ntree = 6000,block.size = 1)

rf.sur = rfsrc(Surv(overallSurvival,survival)~ageAtDiagnosis+consensus_sex1+timeOfSampleCollectionRelativeToInclusion,
               nsplit =10,tree.err=TRUE,na.impute = T,
               data = train_rna.final[,-1],ntree = 100,block.size = 1)

ageAtDiagnosis+FLT3.ITD1+priorMDS1+
  NPM11 +specificDxAtAcquisition12+auccluster+dxAtSpecimenAcquisition2+
  ENSG00000163606+ENSG00000160191+ENSG00000181754+ENSG00000230295+ENSG00000176927+
  ENSG00000226174+
  DNMT3A+SRSF2+WT1+IDH1+STAG2+RUNX1

p = predict(rf.sur,newdata = test_rna.final,na.impute = T)
survival = p$predicted

aucs_final = data.frame(lab_id = sample,survival = survival)

write.csv(aucs_final,"prediction/predictions.csv",row.names = FALSE)


