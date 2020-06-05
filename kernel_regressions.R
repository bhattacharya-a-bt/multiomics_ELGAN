library(data.table)

###### FUNCTIONS USED IN SCRIPT

### Transform a matrix into a linear kernel matrix
matToKernel <- function(mat){
    
    return(mat %*% t(mat) / ncol(mat))
    
}


### Transform a matrix into a Gaussian kernel matrix 
### with width = no. of features
matToGaussKernel <- function(mat){
    
    require(KRLS)
    mod <- gausskernel(X = mat,sigma = ncol(mat))
    return(mod)
    
}

### Measure similarity between two kernels k1 and k2
kernelSim = function(k1, k2){
    
    test.stat = sum(k1 * k2)/sqrt(sum(k1 * k1) * sum(k2 * k2))
    return(test.stat)
    
    
}

### FUNCTION TO RUN KERNEL REGRESSION AND OUTPUT PREDICTIVE R2
### Inputs: 
###         y, outcome
###         kernel, kernel matrix
###         partition, training partition (use same seed)
###         clinical, logical to use clinical fixed effects
###         clinMat, clinical model matrix
linearKernelRegression = function(y,
                                  kernel,
                                  partition,
                                  clinical,
                                  clinMat){
    
    if (!clinical){
        require(rrBLUP)
        mod = mixed.solve(y = y[partition],
                          Z = kernel[partition,partition],
                          method='REML',
                          SE=T)
        pred = as.numeric(mod$beta) + 
            as.numeric(mod$u %*% 
                           solve(kernel[partition,partition]) %*% 
                           kernel[partition,-partition])
        return(summary(lm(pred ~ y[-partition]))$adj.r.squared)
    }
    
    if (clinical){
        require(rrBLUP)
        mod = mixed.solve(y = y[partition],
                          Z = kernel[partition,partition],
                          X = clinMat[partition,],
                          method='REML',
                          SE=T)
        pred = as.numeric( mod$beta %*% 
                               t(clinMat[-partition,])) + 
            as.numeric(mod$u %*% 
                           solve(kernel[partition,partition]) %*% 
                           kernel[partition,-partition])
        return(summary(lm(pred ~ y[-partition]))$adj.r.squared)
    }
}	

### FUNCTION TO RUN ELASTIC NET REGRESSION AND OUTPUT PREDICTIVE R2
### Inputs: 
###         y, outcome
###         kernel, kernel matrix
###         partition, training partition (use same seed)
###         clinical, logical to use clinical fixed effects
###         clinMat, clinical model matrix
enetRegression = function(y,
                          fullMat,
                          partition,
                          method,
                          clinical,
                          clinMat){
    
    if (!clinical){
        require(glmnet)
        mod = cv.glmnet(y = y[partition],
                        x = fullMat[partition,],
                        nfolds = 5,
                        alpha = .5)
        pred = as.numeric(predict(mod,
                                  newx = fullMat[-partition,],
                                  s = 'lambda.min'))
        return(summary(lm(pred ~ y[-partition]))$adj.r.squared)
    }
    
    if (clinical){
        require(glmnet)
        mod = cv.glmnet(y = y[partition],
                        x = cbind(clinMat[partition,],
                                  fullMat[partition,]),
                        nfolds = 5, 
                        alpha = .5)
        pred = as.numeric(predict(mod,
                                  newx = cbind(clinMat[-partition,],
                                               fullMat[-partition,]),
                                  s = 'lambda.min'))
        return(summary(lm(pred ~ y[-partition]))$adj.r.squared)
    }
    
}	

### FUNCTION TO RUN ELASTIC NET REGRESSION AND OUTPUT PREDICTIVE R2
### Inputs: 
###         y, outcome
###         fullMat, full omic matrix (not a kernel)
###         partition, training partition (use same seed)
gaussKernelRegression = function(y,
                                 fullMat,
                                 partition){
    require(KRLS)
    mod = krls(y = y[partition],
               X = fullMat[partition,])
    pred = as.numeric(predict(mod,
                              newdata = fullMat[-partition,])$fit)
    return(summary(lm(pred ~ y[-partition]))$adj.r.squared)
}	


###############################################################################

## DEFINE PHENOTYPE AND P-VALUE CUTOFF
pheno = 'SRS'

## METHYLATION/mRNA EXPRESSION/miRNA EXPRESSION
## MATRICES ARE REDUCED TO THE CURRENT TEST AND TRAINING FOLDS
## FEATURES ARE SELECTED WITH THE FOLD-WISE EWAS/DEG/DEMIR ANALYSIS

## REFER TO DESEQ2 DOCUMENTATION (Love 2014, Genome Biology)
## FOR DEG AND DEMIR ANALYSIS CODE

## REFER TO EWAS_example.R FOR EWAS CODE

### REGRESSION OF SRS and IQ AGAINST ASD
asd_cog_reg = glm(ASD ~ SRS + IQ, data = covs, family = 'binomial')
summary(asd_cog_reg)
null.mod = glm(ASD ~ 1, data = covs, family = 'binomial')
mcfad.r2 = 1-logLik(asd_cog_reg)/logLik(null.mod)

### CREATE COVARIATES MATRIX
covs = as.data.frame(covs)
Y = covs[,pheno]
metadata = covs[,c('Race','Sex',
                   'Maternal Age',
                   'Smoker','Insurance','Gestational Days')]

### REMOVE DEMOGRAPHICS FROM OMICS
require(limma)
clinical.design = model.matrix(~.-1,data=metadata)
mRNA = t(removeBatchEffect(x = t(exp), covariates = clinical.design))
miRNA = t(removeBatchEffect(x = t(mirna), covariates = clinical.design))
meth = t(removeBatchEffect(x = t(meth), covariates = clinical.design))

### MAKE LISTS OUT OF OMIC MATRICES TO FACILITATE KERNEL PROJECTION
require(pbapply)
full.Y = Y
full.X = list(mRNA = pbapply(mRNA ,2,scale),
              miRNA = pbapply(miRNA ,2,scale),
              methylation = pbapply(meth ,2,scale),
              mRNA.miRNA = cbind(pbapply(mRNA ,2,scale),
                                 pbapply(miRNA ,2,scale)),
              mRNA.methylation = cbind(pbapply(mRNA ,2,scale),
                                       pbapply(meth,2,scale)),
              miRNA.methylation = cbind(pbapply(miRNA ,2,scale),
                                        pbapply(meth,2,scale)),
              mRNA.miRNA.methylation = cbind(pbapply(mRNA ,2,scale),
                                             pbapply(miRNA ,2,scale),
                                             pbapply(meth,2,scale)))


### TRANSFORM INTO KERNELS
full.X.kernel = lapply(full.X,matToKernel)
full.X.gauss = lapply(full.X[c(1,2,3,4)],
                      matToKernel)

require(rlist)
full.X.gauss = list.append(full.X.gauss, 
                           full.X.kernel$mRNA + 
                               matToGaussKernel(full.X$methylation))
full.X.gauss = list.append(full.X.gauss, 
                           full.X.kernel$miRNA + 
                               matToGaussKernel(full.X$methylation))
full.X.gauss = list.append(full.X.gauss, 
                           full.X.kernel$mRNA.miRNA + 
                               matToGaussKernel(full.X$methylation))
names(full.X.gauss) = names(full.X)


clinical.kernel = matToKernel(clinical.design)


### REGRESS OUTCOME AGAINST DEMOGRAPHICS AND RESIDUALIZE
totalData = data.frame(Y = full.Y, metadata)
meta.reg = lm(Y ~ ., data = totalData)
print(summary(meta.reg))
r2.metadata = summary(meta.reg)$adj.r.squared
metaDesign = model.matrix(~.-1,data=metadata)
clinical.design = metaDesign
rY = as.numeric(resid(meta.reg))

set.seed(seed)
pblapply(full.X,
         gaussKernelRegression,
         y = rY,
         partition = partition)
pblapply(full.X.kernel,
         linearKernelRegression,
         y = rY,
         partition = partition,
         clinical = F,
         clinMat = clinMat)
pblapply(full.X,
         enetRegression,
         y = rY,
         partition = partition,
         clinical = F,
         clinMat = clinMat)

