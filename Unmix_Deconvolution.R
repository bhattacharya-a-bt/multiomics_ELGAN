require(data.table)
require(RUVSeq)
require(DESeq2)

### cts is raw counts of mRNA expression from ELGAN
### genelocs is MatrixEQTL format gene locations

cts = cts[,-which(colnames(cts) %in% exclude),with=F]
genes = cts$Gene
cts = as.matrix(cts[,-1])
rownames(cts) = genes
colnames(cts) = sapply(strsplit(colnames(cts),'[.]'),function(vec) vec[1])

### FILTER OUT NON-EXPRESSED GENES
filter <- apply(cts, 1, function(x) length(as.numeric(x)[as.numeric(x)>5])>=2)
cts <- cts[filter,]


### CREATE AN EXPERIMENT
### fullCovar is the data frame for covariates from ELGAN
set <- newSeqExpressionSet(as.matrix(round(cts)),
                           phenoData = data.frame(fullCovar, row.names=colnames(cts)))

#### FIND GENES NOT ASSOCIATED WITH OUTCOME OF INTEREST
design <- model.matrix(~asd,data = pData(set))
y = DGEList(counts = counts(set), group = pData(set)$asd)
y <- calcNormFactors(y, method="TMM")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)

top <- topTags(lrt, n=nrow(set))$table
empirical <- rownames(top)[1:1000]
empirical = empirical[-which(empirical %in% genelocs$geneid[genelocs$chr == 'Y'])]


#### REMOVE TECHNICAL VARIATION BASED ON NEGATIVE CONTROLS
set2 <- RUVg(set, empirical, k=1)
counts = counts(set2)

sc = fread('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE89nnn/GSE89497/suppl/GSE89497_Human_Placenta_TMP_V2.txt.gz')
sc = subset(sc,V1 %in% rownames(counts))
counts.dc = counts[match(sc$V1,
                      rownames(counts)),]
hash = data.frame(cn = colnames(sc)[-1])
hash$cellType = sapply(strsplit(hash$cn,'_'),
                       function(x) x[2])
hash$cellType = as.factor(hash$cellType)
levels(hash$cellType) = c('CTB','CTB','CTB','CTB',
                          'EVT','STB','STR','STR')
pure.mat = matrix(ncol=4,nrow=nrow(sc))
colnames(pure.mat) = c('CTB','EVT','STB','STR')
pure.mat[,1] = rowMeans(sc[,hash$cn[hash$cellType == 'CTB'],with=F])
pure.mat[,2] = rowMeans(sc[,hash$cn[hash$cellType == 'EVT'],with=F])
pure.mat[,3] = rowMeans(sc[,hash$cn[hash$cellType == 'STB'],with=F])
pure.mat[,4] = rowMeans(sc[,hash$cn[hash$cellType == 'STR'],with=F])
rownames(pure.mat) = sc$V1

unmix.props = unmix(x = counts.dc,
                    pure = pure.mat,
                    shift = 1e2)
require(DeconRNASeq)
decon.props = as.matrix(DeconRNASeq(as.data.frame(counts.dc),
                                    as.data.frame(pure.mat))$out.all)
rownames(decon.props) = rownames(unmix.props)
unmix.props[which(rowSums(is.na(unmix.props)) == 4),] = 
    decon.props[which(rowSums(is.na(unmix.props)) == 4),]
unmix.props[which(rowSums(unmix.props == 1) == 1),] = 
    decon.props[which(rowSums(unmix.props == 1) == 1),]

props = unmix.props
library(DESeq2)
require(qvalue)

pData(set2)$CTB = props[,1]
pData(set2)$EVT = props[,2]
pData(set2)$STB = props[,3]
pData(set2)$STR = props[,4]

df = pData(set2)
results = data.frame(Gene = rownames(counts),
                     Beta_CTB = NA,
                     P_CTB = NA,
                     Beta_EVT = NA,
                     P_EVT = NA,
                     Beta_STB = NA,
                     P_STB = NA,
                     Beta_STR = NA,
                     P_STR = NA)
counts = vst(counts)

for (i in 1:nrow(counts)){
    
    print(i)
    df.new = df
    df.new$Gene = log(counts[i,]+1)
    df.new$Gene = resid(lm(Gene~W_1,data=df.new))
    lm1 = lm(srstt ~ Gene*CTB + Gene + as.factor(race1) + 
                 gadays + 
                 as.factor(medu1) + as.factor(sex) +
                 CTB, data = df.new)
    results[i,2:3] = as.numeric(coef(summary(lm1))[10,c(1,4)])
    
    lm1 = lm(srstt ~ Gene*EVT + Gene + as.factor(race1) + 
                 gadays + 
                 as.factor(medu1) + as.factor(sex) +
                 EVT, data = df.new)
    results[i,4:5] = as.numeric(coef(summary(lm1))[10,c(1,4)])
    
    lm1 = lm(srstt ~ Gene*STB + Gene + as.factor(race1) + 
                 gadays + 
                 as.factor(medu1) + as.factor(sex) +
                 STB, data = df.new)
    results[i,6:7] = as.numeric(coef(summary(lm1))[10,c(1,4)])
    
    lm1 = lm(srstt ~ Gene*STR + Gene + as.factor(race1) + 
                 gadays + 
                 as.factor(medu1) + as.factor(sex) +
                 STR, data = df.new)
    results[i,8:9] = as.numeric(coef(summary(lm1))[10,c(1,4)])
    
}
results$FDR_CTB = p.adjust(results$P_CTB,method = 'BH')
results$FDR_EVT = p.adjust(results$P_EVT,method = 'BH')
results$FDR_STB = p.adjust(results$P_STB,method = 'BH')
results$FDR_STR = p.adjust(results$P_STR,method = 'BH')

fwrite(results,'results_srstt_int_umix.tsv',sep='\t',quote=F,col.names=T,row.names=F)

for (i in 1:nrow(counts)){
    
    print(i)
    df.new = df
    df.new$Gene = log(counts[i,]+1)
    df.new$Gene = resid(lm(Gene~W_1,data=df.new))
    lm1 = lm(diq ~ Gene*CTB + Gene + as.factor(race1) + 
                 gadays + 
                 as.factor(medu1) + as.factor(sex) +
                 CTB, data = df.new)
    results[i,2:3] = as.numeric(coef(summary(lm1))[10,c(1,4)])
    
    lm1 = lm(diq ~ Gene*EVT + Gene + as.factor(race1) + 
                 gadays + 
                 as.factor(medu1) + as.factor(sex) +
                 EVT, data = df.new)
    results[i,4:5] = as.numeric(coef(summary(lm1))[10,c(1,4)])
    
    lm1 = lm(diq ~ Gene*STB + Gene + as.factor(race1) + 
                 gadays + 
                 as.factor(medu1) + as.factor(sex) +
                 STB, data = df.new)
    results[i,6:7] = as.numeric(coef(summary(lm1))[10,c(1,4)])
    
    lm1 = lm(diq ~ Gene*STR + Gene + as.factor(race1) + 
                 gadays + 
                 as.factor(medu1) + as.factor(sex) +
                 STR, data = df.new)
    results[i,8:9] = as.numeric(coef(summary(lm1))[10,c(1,4)])
    
}
results$FDR_CTB = p.adjust(results$P_CTB,method = 'BH')
results$FDR_EVT = p.adjust(results$P_EVT,method = 'BH')
results$FDR_STB = p.adjust(results$P_STB,method = 'BH')
results$FDR_STR = p.adjust(results$P_STR,method = 'BH')

fwrite(results,'results_diq_int_umix.tsv',sep='\t',quote=F,col.names=T,row.names=F)

