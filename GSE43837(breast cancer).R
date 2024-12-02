#WEEK 3

```{r}
#BiocManager::install("useful")
library("GEOquery")
# download the GSE26922 expression data series
gse <- getGEO("GSE43837")
              sep = "\t")
# getGEread.delim()# getGEO returns a list of expression objects, but ...
length(gse)
# ... shows us there is only one object in it. We assign it to
# the same variable.
gse <- gse[[1]]
# The actual expression data are accessible in the "exprs" ... an
# Expression Set, the generic data class that BioConductor uses
# for expression data
head(exprs(gse))
length(exprs(gse)) 
# exprs(gse) is a matrix that we can assign to its own variable, to
# conveniently access the data rows and columns
ex <- exprs(gse)
dim(ex)      #61259 rows and 38 columns
colnames(ex)
#ex[1:5,]
# Analyze value distributions
#boxplot(ex)
ex2 <- log2(ex)
boxplot(ex2)

gse$title   #19 brain_met and 19 breast primary
gse@assayData

exprs(gse)   #access the data matrix

#scale normalization in R
normalized.ex <- scale(ex2)   #
boxplot(normalized.ex)
```


#week 4 

```{r}
# PCA

ex2 <- na.omit(as.matrix(ex2))
pca <- prcomp(t(ex2))
summary(pca)
screeplot(pca)

# draw PCA plot
grpcol <- c(rep("blue",19), rep("red",19))
plot(pca$x[,1], pca$x[,2], xlab="PCA1", ylab="PCA2", main="PCA for components 1&2", type="p", pch=10, col=grpcol)
text(pca$x[,1], pca$x[,2], rownames(pca$x), cex=0.75) 


#WEEK 4
#k-means   
set.seed(99)
library("useful")

k <- 2  
kmeans_result <- kmeans(t(ex2), k)
table(kmeans_result$cluster)
#Note that the cluster visualization is not trivial (space of 7815 dimensions!). We use the
# plot function in the ‘useful’ package that performs a dimensionality reduction using PCA
plot(kmeans_result, data=t(ex2)) + geom_text(aes(label=colnames(ex2)),hjust=0,vjust=0)

#hierarchical clustering          
dist_matrix <- dist(t(ex2))
hc_result <- hclust(dist_matrix, method = "complete")  #ave = average linkage clustering --> change it
k <- 2                               #k = number of group
groups <- cutree(hc_result, k=k)
table(groups)
#plot(hc_result, hang <- -1, labels=groups)
plot(hc_result, hang <- -1, labels=groups)
rect.hclust(hc_result, k = 2, which = NULL, x = NULL, h = NULL, border = 2, cluster = NULL) # red boxes to show groups
```



#week 5
```{r}
#WEEK 5-6
#simpler method if you don’t need the complexity of genefilter
small.eset <- log2(na.omit(ex))
dim(small.eset) 
group <- c(rep('B',19),rep('T',19)) # classification, in order
#B --> brain metastais T --> breast primary

# build RF
library(randomForest)
set.seed(99)
rf <- randomForest(x=t(small.eset), y=as.factor(group), ntree=1000)
# a trivial test
predict(rf, t(small.eset[, 1:5]))
# graph of sorted importance values
plot(sort(rf$importance, decreasing=TRUE)) # can also use: varImpPlot(rf)
#extract the most 'important' genes
probe.names <- rownames(rf$importance)
top200 <- probe.names[order(rf$importance, decreasing=TRUE)[1:200]]
write.csv(top200, file = "probes-top200-breastdataset.txt", quote=FALSE, row.names = FALSE, col.names=FALSE)
          
#Draw a heatmap
# Look at variable importance
imp.temp <- abs(rf$importance[,])
t <- order(imp.temp,decreasing=TRUE)
plot(c(1:nrow(small.eset)),imp.temp[t],log='x',cex.main=1.5,
     xlab='gene rank',ylab='variable importance',cex.lab=1.5,
     pch=16,main='ALL subset results')
# Get subset of expression values for 25 most 'important' genes
gn.imp <- names(imp.temp)[t]
gn.25 <- gn.imp[1:25] # vector of top 25 genes, in order
t <- is.element(rownames(small.eset),gn.25)
sig.eset <- small.eset[t,] # matrix of expression values, not necessarily in order

## Make a heatmap, with group differences obvious on plot
library(RColorBrewer)
hmcol <- colorRampPalette(brewer.pal(11,"PuOr"))(256)
colnames(sig.eset) <- group # This will label the heatmap columns
csc <- rep(hmcol[50],30)
csc[group=='T'] <- hmcol[200]
# column side color will be purple for T and orange for B
heatmap(sig.eset, scale="row", col=hmcol, ColSideColors=csc)

```


#WEEK 7

```{r}
#LDA + ROC
# install.packages("BiocManager")
# BiocManager::install("genefilter")
# BiocManager::install("GEOquery")
# BiocManager::install("pROC")

ex2 <- na.omit(as.matrix(ex))
f <- factor(c(rep("affected",19), rep("control",19)))
library("genefilter")
tt38 <- rowttests(ex2,f)
# keepers <- which(p.adjust(tt40$p.value)<0.1)
keepers <- which(tt38$p.value<0.01)
ex3 <- ex2[keepers,]
tex3 <- t(ex3)
dat <- cbind(as.data.frame(tex3),f)
colnames(dat)[ncol(dat)] <- "AFFECTED"
n.controls <- 19
n.affected <- 19

train <- sample(1:(n.controls), (n.controls-5))
test <- setdiff(1:(n.controls),train)
test<- c(test, test+20)
train <- c(train, train+20)

library("MASS")
mod <- lda(AFFECTED ~ ., data=dat, prior = c(0.5,0.5),
           subset = train)
plot(mod)

mod.values <- predict(mod, dat[train,])
mod.values$class
plot(mod.values$x[,1], ylab=c("LDA Axis"))
text(mod.values$x[,1],
     col=c(as.numeric(dat[train,"AFFECTED"])+10))

preds<-predict(mod, dat[test,])
preds$class
table(as.numeric(preds$class),
      as.numeric(dat[test, "AFFECTED"]) )

library("pROC")
roc_lda <- plot.roc(as.numeric(preds$class),
                    as.numeric(dat[test, "AFFECTED"]) )

#CARET: another way to perform LDA, more efficient and sophisticate way
#ex2 <- na.omit(as.matrix(ex))

f <- factor(c(rep("affected",19), rep("control",19)))
library("genefilter")
tt38 <- rowttests(ex2,f)
# keepers <- which(p.adjust(tt40$p.value)<0.1)
keepers <- which(tt38$p.value<0.01)
ex3 <- ex2[keepers,]
tex3 <- t(ex3)
dat <- cbind(as.data.frame(tex3),f)
colnames(dat)[ncol(dat)] <- "AFFECTED"

library("caret")
# Run algorithms using 10-fold cross validation
#
control <- trainControl(method="cv", number=10)
metric <- "Accuracy"
fit.lda <- train(AFFECTED~., data=dat, method="lda",
                 metric=metric, trControl=control)
fit.rf <- train(AFFECTED~., data=dat, method="rf",
                metric=metric, trControl=control)
results <- resamples(list(LDA=fit.lda, RF=fit.rf))
summary(results)
ggplot(results) + labs(y = "Accuracy")

# Run algorithms using 10-fold cross validation, 10 times
#I can not do more than 10 because of the small dataset but i can do 
#control <- trainControl(method = "repeatedcv", number = 10, repeats = 10)
#fit.lda.2 <- train(AFFECTED~., data=dat, method="lda",
#                  metric=metric, trControl=control)
#fit.rf.2 <- train(AFFECTED~., data=dat, method="rf",
  #                metric=metric, trControl=control)
#results <- resamples(list(LDA=fit.lda.2, RF=fit.rf.2))  
#ggplot(results) + labs(y = "Accuracy")


```


#WEEK 8

```{r}
#LASSO
y <- c(rep(0,19),rep(1,19))
f <- factor(y)
# IMPORTANT: the type of factor changes the
# way the function glmnet works:
# - f with numerical values: regression mode
# - f with categorical values: classification mode
#install.packages("glmnet")

library("glmnet")
dat <- t(ex2)
fit=glmnet(dat,y,standardize=FALSE,family="binomial")
plot(fit, xvar = "lambda", label=TRUE)
cfit=cv.glmnet(dat,y,standardize=FALSE,family="binomial")
plot(cfit)
coef(cfit, s=cfit$lambda.min)

# repeat analysis but by using train + test sample subsets
n.controls<-19
n.affected<-19
train <- sample(1:(n.controls), (n.controls-5))
test <- setdiff(1:(n.controls),train)
test<- c(test, test+19)
train <- c(train, train+19)
fit=glmnet(dat[train,],y[train],standardize=FALSE,family="binomial")
plot(fit)
cfit=cv.glmnet(dat[train,],y[train],standardize=FALSE,family="binomial")
plot(cfit)
predict(fit,dat[test,], type="class", s= cfit$lambda.min)

# plot ROCR curve
library("ROCR")
pred2 <- predict(fit,dat[test,], type="response", s=cfit$lambda.min)
plot(performance(prediction(pred2, y[test]), 'tpr', 'fpr'))
# compute Area Under the Curve (AUC)
auc.tmp <- performance(prediction(pred2, y[test]),"auc")
auc <- as.numeric(auc.tmp@y.values)

#CARET + LASSO
# factor vector of categorical type
f <- factor(c(rep("affected",19), rep("control",19)))

# optional feature selection step
library("genefilter")
tt38 <- rowttests(ex2,f)

# keepers <- which(p.adjust(tt40$p.value)<0.1)
keepers <- which(tt38$p.value<0.01)
ex3 <- ex2[keepers,]
tex3 <- t(ex3)
#tex3 <- t(ex2) # disable feature selection
dat <- cbind(as.data.frame(tex3),f)
colnames(dat)[ncol(dat)] <- "AFFECTED"

library("caret")
control <- trainControl(method="cv", number=10)
metric <- "Accuracy"
library("glmnet")
fit.lasso <- train(AFFECTED~., data=dat, method="glmnet",
                   family = "binomial",
                   tuneGrid = expand.grid(alpha = 1,
                                          lambda = seq(0,1,by=0.05)),
                   trControl = control,
                   metric = metric)

plot(fit.lasso)

# comparison with other classification methods
fit.lda <- train(AFFECTED~., data=dat, method="lda",
                 metric=metric, trControl=control)
fit.rf <- train(AFFECTED~., data=dat, method="rf",
                metric=metric, trControl=control)
results <- resamples(list(RF=fit.rf, LDA=fit.lda, Lasso=fit.lasso))
summary(results)
ggplot(results) + labs(y = "Accuracy") 


#parte 2
#RSCUDO 
#install.packages("igraph")
#BiocManager::install("rScudo")

y <- c(rep(0,19),rep(1,19))
f <- factor(y, labels = c("Affected",
                          "Control"))
library("caret")
set.seed(123)
inTrain <- createDataPartition(f, list = FALSE)
#trainData <- dat[, inTrain]
#testData <- dat[, -inTrain]
trainData <- ex3[, inTrain]
testData <- ex3[, -inTrain]

# analyze training set
library("rScudo")
trainRes <- scudoTrain(trainData, groups = f[inTrain],
                       nTop = 10, nBottom = 10, alpha = 0.05)
trainRes

# inspect signatures
upSignatures(trainRes)[1:5,1:5]
consensusUpSignatures(trainRes)[1:5, ]
# generate and plot map of training samples
trainNet <- scudoNetwork(trainRes, N = 0.2)

scudoPlot(trainNet, vertex.label = NA)
# perform validation using testing samples
testRes <- scudoTest(trainRes, testData, f[-inTrain],
                     nTop = 25, nBottom = 25)
testNet <- scudoNetwork(testRes, N = 0.2)
scudoPlot(testNet, vertex.label = NA)
# identify clusters on map
library("igraph")
testClust <- igraph::cluster_spinglass(testNet, spins = 2)
plot(testClust, testNet, vertex.label = NA)
# perform classification
classRes <- scudoClassify(trainData, testData, N = 0.25,
                          nTop = 12, nBottom = 12,
                          trainGroups = f[inTrain], alpha = 0.5)
caret::confusionMatrix(classRes$predicted, f[-inTrain])



#RSCUDO + CARET
y <- c(rep(0,19),rep(1,19))
f <- factor(y, labels = c("Affected",
                          "Control"))
library("caret")
set.seed(123)
inTrain <- createDataPartition(f, list = FALSE)
#trainData <- dat[, inTrain]
#testData <- dat[, -inTrain]

trainData <- ex3[, inTrain]
testData <- ex3[, -inTrain]
# use caret to test a grid a values for nTop & nBottom
# using cross validation
model <- scudoModel(nTop = (2:6)*5, nBottom = (2:6)*5,
                    N = 0.25)
control <- caret::trainControl(method = "cv", number = 5,
                               summaryFunction =
                                 caret::multiClassSummary)
cvRes <- caret::train(x = t(trainData), y = f[inTrain],
                      method = model,
                      trControl = control)
# plot map of testing samples using best nTop & nBottom
# values
testRes <- scudoTest(trainRes, testData, f[-inTrain],
                     cvRes$bestTune$nTop,
                     cvRes$bestTune$nBottom5)
testNet <- scudoNetwork(testRes, N = 0.2)
scudoPlot(testNet, vertex.label = NA)
# perform classification of testing samples using best
# nTop & nBottom values
classRes <- scudoClassify(ex3[, inTrain], ex3[, -inTrain],
                          0.25,
                          cvRes$bestTune$nTop,
                          cvRes$bestTune$nBottom, f[inTrain], alpha = 0.05)
caret::confusionMatrix(classRes$predicted, f[-inTrain])

```



#WEEK 9
Functional enrichment analysis:
- DAVID 
- amiGO 
- GSEA 
- g profiler

```{r}
#install.packages("gprofiler2")
library(gprofiler2)
# see vignette at https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html
# Caricare i dati dal file probes-top200-breastdataset
genes <- read.table("C:/Users/eleon/Desktop/probes-top200-breastdataset.txt", header = TRUE, sep = "\t")
gene_list <- as.vector(genes)
gostres <- gost(query = gene_list,   
                organism = "hsapiens", ordered_query = FALSE,
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                measure_underrepresentation = FALSE, evcodes = FALSE,
                user_threshold = 0.05, correction_method = "g_SCS",
                domain_scope = "annotated", custom_bg = NULL,
                numeric_ns = "", sources = NULL, as_short_link = FALSE)
names(gostres)
head(gostres$result)
# visualize results using a Manhattan plot
gostplot(gostres, capped = TRUE, interactive = TRUE)
# when ready, create publication quality (static) plot + table of interesting terms/pathways
p <- gostplot(gostres, capped = TRUE, interactive = FALSE)
publish_gostplot(p, highlight_terms = c("?????"),           # ????????geni di interesse????
                 width = NA, height = NA, filename = NULL )


```

#week 10

```{r}
common_genes <- intersect(top200, rownames(ex2))
filtered_expression <- ex2[common_genes, , drop=FALSE]
f <- factor(c(rep("affected",19), rep("control",19)))
input_pathfindR <- rowttests(filtered_expression,f)
gene_pvalue_df <- data.frame(
  Gene = rownames(input_pathfindR),      # Nomi dei geni
  P_Value = input_pathfindR$p.value     # P-value calcolati
)

gene_pvalue_df <- gene_pvalue_df[order(gene_pvalue_df$Gene), ]
# Scrivi il data frame in un file CSV

write.csv(gene_pvalue_df, file = "rf_genes_pvalue.csv", quote = FALSE, row.names = FALSE)

rf_genes_pvalue = read.csv("C:/Users/eleon/Desktop/rf_genes_pvalue.csv", sep = ",")
rf_genes_pvalue$P_Value <- as.numeric(rf_genes_pvalue$P_Value)
rf_genes_pvalue <- lda_genes_pvalue[!is.na(lda_genes_pvalue$P_Value), ]

write.csv(gene_pvalue_df, file = "rf_genes_pvalue.csv", quote = FALSE, row.names = FALSE)


```





```{r}

#BiocManager::install("u133x3p.db")
library("u133x3p.db")

rf_genes_pvalue <- read.csv("C:/Users/eleon/Desktop/rf_genes_pvalue.csv", stringsAsFactors = FALSE)

probe_ids <- rf_genes_pvalue$Gene
my.map <- u133x3pENTREZID
mapped.probes <- mappedkeys(my.map)
gene_ids <- sapply(probe_ids, function(probe) {
  if (probe %in% mapped.probes) {
    return(my.map[[probe]])
  } else {
    return(NA)  
  }
})

rf_genes_pvalue$Gene_ID <- gene_ids
head(lda_genes_pvalue)

write.csv(rf_genes_pvalue, "rf_genes_pvalue_with_gene_ids.csv", row.names = FALSE)

# Load necessary packages
library(AnnotationDbi)
library(u133x3p.db)

# Load your data
lda_genes_pvalue <- read.csv("C:/Users/eleon/Desktop/rf_genes_pvalue.csv", stringsAsFactors = FALSE)

# Assuming the 'Gene' column contains the probe IDs
probe_ids <- lda_genes_pvalue$Gene

# Use the mapIds function to map the probe IDs to gene symbols
gene_symbols <- mapIds(u133x3p.db, keys = probe_ids, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")

# Add the gene symbols to your data frame
lda_genes_pvalue$Gene_Symbol <- gene_symbols

# Check the updated dataframe
head(lda_genes_pvalue)

# Assuming lda_genes_pvalue already contains the Gene_Symbol and p-value columns
# Ensure that you have the Gene_Symbol and p-value columns
gene_pvalue_data <- lda_genes_pvalue[, c("Gene_Symbol", "P_Value")]
gene_pvalue_data_clean <- gene_pvalue_data[!is.na(gene_pvalue_data$Gene), ]
gene_pvalue_data_clean <- na.omit(gene_pvalue_data)

# Check if the columns are correct
head(gene_pvalue_data)

# Save this new data frame to a CSV file
write.csv(gene_pvalue_data_clean, "gene_symbols_pvalues.csv", row.names = FALSE)
input = read.csv("C:/Users/eleon/Desktop/gene_symbols_pvalues.csv")


```





```{r}

library("KEGGREST")
library("KEGGgraph")
library("AnnotationDbi")
library("org.Hs.eg.db")
#install.packages("pathfindR")
#BiocManager::install("GeneMania")

library("pathfindR")

## pathway enrichment
RA_demo <- run_pathfindR(input, iterations = 1) # keeps execution time short - default is 10
head(RA_demo)
## clutser enriched terms
RA_demo_clu <- cluster_enriched_terms(RA_demo)
## term-gene graph of top 10 terms
term_gene_graph(RA_demo)

## pathway enrichment can use different networks, different gene sets
RA_demo2 <- run_pathfindR(input,
 gene_sets = "KEGG",
 pin_name_path = "GeneMania")
head(RA_demo2)
```
