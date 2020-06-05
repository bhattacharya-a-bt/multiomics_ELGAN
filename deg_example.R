### EXAMPLE OF DEG ANALYSIS
### CAN BE MODIFIED DEMIR ANALYSIS
require(DESeq2)
require(RUVSeq)
require(qvalue)
require(data.table)

### Create a SeqExpressionSet using:
## 1. cts - counts of mRNA expression
## 2. fullCovar - covariates as the phenodata
set <- newSeqExpressionSet(as.matrix(round(cts)),
                           phenoData = data.frame(fullCovar, 
                                                  row.names=colnames(cts)))

#### FIND GENES NOT ASSOCIATED WITH OUTCOME OF INTEREST FOR RUV
design <- model.matrix(~srstt,data = pData(set))
y = DGEList(counts = counts(set), group = pData(set)$srstt)
y <- calcNormFactors(y, method="TMM")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)

### USE THE TOP 1000 UNASSOCIATED GENES AS EMPIRICAL HOUSEKEEPERS
top <- topTags(lrt, n=nrow(set))$table
empirical <- rownames(top)[1:1000]
empirical = empirical[-which(empirical %in% 
                                 genelocs$geneid[genelocs$chr == 'Y'])]


#### REMOVE TECHNICAL VARIATION BASED ON NEGATIVE CONTROLS
set2 <- RUVg(set, empirical, k=1)


dds <- DESeqDataSetFromMatrix(countData = counts(set2),
                              colData = pData(set2),
                              design = ~ W_1 + as.factor(race1) + 
                                  gadays + 
                                  as.factor(medu1) + as.factor(sex) + srstt)


# Estimate size factors by the median ratio method (Anders 2010, Genome Biology)
dds <- estimateSizeFactors(dds)

# Check for overdispersion
dds <- estimateDispersionsGeneEst(dds)
cts <- counts(dds, normalized=TRUE)

# Estimate dispersion per gene
disp <- pmax((rowVars(cts) - rowMeans(cts)),0)/rowMeans(cts)^2

# Assign dispersion to metadata of dds
mcols(dds)$dispGeneEst <- disp

# Correct for overdispersion
dds <- estimateDispersionsFit(dds, fitType="mean")

# Induce homoskedasticity based on fitted dispersion-mean relations
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
mat <- assay(vsd)

# Remove UV
covars <- as.matrix(colData(dds)[,grep("W",colnames(colData(dds))),drop=FALSE])
mat <- removeBatchEffect(mat, covariates=covars)

dds <- DESeq(dds, parallel=TRUE, BPPARAM=MulticoreParam(8))
res <- as.data.frame(results(dds))
res$qvalue = qvalue(res$pvalue)$qvalue
fwrite(res,'srstt_res_deg_nozcat.csv',quote=F,
       col.names=T,row.names=T)
