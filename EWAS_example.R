vectorClass <- function(vec,class){
    if (class == 'numeric'){
        return(as.numeric(vec))
    }
    if (class == 'factor'){
        return(as.factor(round(vec,0)))
    }
}

##############################################
# svaEWAS() runs an EWAS with including the top correlated SVs from SVA analysis
#
# Inputs:
# fileBeta -> the file path to the beta values
# fileCovaasr -> the file path to the metadata
# stratified -> set to 1 or 0 depending on if you want a sex-stratified analysis
#               default is 'No', meaning an unstratified analysis
# pheno -> the variable name for your phenotype of interest
# class -> the class of your phenotype (either 'factor' or 'numeric')
# numSV -> the number of SVs you want in your analysis
#          default is 10


svaEWAS <- function(fileBeta,
                    fileCovar,
                    impute = T,
                    stratified = 'No',
                    phenoName,
                    phenoClass,
                    numSV = 10,
                    covarNames,
                    covarClasses,
                    smart,
                    fileQQ,
                    fdrcut,
                    cut){
    
    vectorClass <- function(vec,class){
        if (class == 'numeric'){
            return(as.numeric(vec))
        }
        if (class == 'factor'){
            return(as.factor(round(vec,0)))
        }
    }
    
    
    
    #################################################
    ## REQUIRED PACKAGES FOR THIS FUNCTION        ##
    #################################################
    require(data.table)
    require(MASS)
    require(limma)
    require(sva)
    require(SmartSVA)
    require(missForest)
    require(qqman)
    require(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    require(statmod)
    
    #################################################
    ## READ IN BETA VALUES AND CONVERT TO M-VALUES ##
    #################################################
    y = fread(fileBeta)
    colnames(y)[which(grepl('cg',y[1,]))] = 'V1'
    y1 = as.data.frame(y[,-which(grepl('cg',y[1,])),with=F])
    rownames(y1) = as.character(as.vector(y$V1))
    edata = as.matrix(log2(y1/(1-y1)))
    edata = edata[,cut]
    rm(y,y1)
    
    #################################################
    ## DEAL WITH COVARIATES                        ##
    #################################################
    
    if (stratified != 'No'){
        covarNames = covarNames[covarNames != 'sex']
        covarClasses = covarClasses[covarNames != 'sex']
        edata =  edata[rownames(edata) %in% 
                           rownames(Locations)[!(Locations$chr %in% c('chrX','chrY'))],]
    }
    
    ## 1. READ IN COVARIATES AND ARRANGE AS EDATA
    covar <- as.data.frame(fread(fileCovar))
    covar <- covar[match(colnames(edata),covar$Array_ID),]
    
    if (stratified == 0){
        covar <- covar[covar$sex==0,]
        edata <- edata[,which(colnames(edata) %in% covar$Array_ID)]
        covar <- subset(covar,Array_ID %in% colnames(edata))}
    if (stratified == 1){
        covar <- covar[covar$sex==1,]
        edata <- edata[,which(colnames(edata) %in% covar$Array_ID)]
        covar <- subset(covar,Array_ID %in% colnames(edata))}
    
    covar <- covar[,c(phenoName,covarNames)]
    
    
    ## 2. IMPUTE MISSING VALUES
    set.seed(1218)
    if (impute == T){
        covar.imp <- missForest(covar)
        fullCovar <- covar.imp$ximp
        rm(covar,covar.imp)}
    if (impute == F){
        fullCovar = covar
        rm(covar)
    }
    
    
    ## 3. SEPARATE OUT OUTCOME AND COVARIATES
    pheno = vectorClass(fullCovar[,phenoName],phenoClass)
    covs = paste0('x',1:length(covarNames))
    
    
    for (i in 1:length(covarNames)){
        assign(covs[i],vectorClass(fullCovar[,covarNames[i]],covarClasses[i]))
    }
    
    
    formula1 = paste('~pheno',paste(covs,collapse = '+'),sep='+')
    formula01 = paste0('~',paste(covs,collapse='+'))
    
    #################################################
    ## SVA FOR EWAS                                ##
    #################################################
    
    ## 1. MODELS WITH AND WITHOUT
    mod1 <- model.matrix(as.formula(formula1))
    mod01 <- model.matrix(as.formula(formula01))
    
    ## 2. FIT TWO-STEP SVA WITHOUT SPECIFYING NUMBER OF SURROGATES
    if (smart == F){
        svobj1 = sva(edata,mod1,mod01,n.sv=numSV,method="two-step")
        modSv1 = cbind(mod1,svobj1$sv[,1:numSV])}
    
    ## 2B. SMARTSVA
    if (smart == T){
        svobj1 = smartsva.cpp(dat = edata, mod = mod1, 
                              mod0 = mod01, n.sv = numSV, 
                              B = 100, alpha = 0.25,
                              epsilon = 0.001, VERBOSE = F)
        modSv1 = cbind(mod1,svobj1$sv[,1:numSV])}
    
    
    ## 3. FIT ROBUST REGRESSION WITH SURROGATES AS COVARIATES
    fit1 = lmFit(edata,modSv1,method="robust")
    
    ## 4. EXTRACT EMPIRICAL BAYES STATISTICS
    fite1 = eBayes(fit1,trend=T,robust=T)
    SE = sqrt(fite1$s2.post) * fite1$stdev.unscaled
    pdf(fileQQ,height=6,width=6)
    qqt(fite1$t[,2], df=fite1$df.total,
        main="Q-Q plot of EWAS Results",
        pch = 18, col = "blue4", cex = .5, las = 1)
    abline(0,1,col='red')
    dev.off()
    
    
    ## 5. P-VALUE ADJUSTMENT
    if (phenoClass == 'numeric'){
        tab1 = topTable(fite1, coef = "pheno",number=nrow(edata),
                        p.val=1,adjust.method = "fdr", 
                        genelist = row.names(edata))
        tab1$SE = as.numeric(SE[rownames(tab1),'pheno'])
    }
    
    if (phenoClass == 'factor'){
        tab1 = topTable(fite1, coef = "pheno1",number=nrow(edata),
                        p.val=1,adjust.method = "fdr", 
                        genelist = row.names(edata))
        tab1$SE = as.numeric(SE[rownames(tab1),'pheno1'])
    }
    
    tab1$Upper_unadjusted = tab1$logFC + tab1$SE * 1.96
    tab1$Lower_unadjusted = tab1$logFC - tab1$SE * 1.96
    if (nrow(subset(tab1,adj.P.Val < fdrcut)) == 0){
        print('No results meet FDR cutoff!')}
    if (nrow(subset(tab1,adj.P.Val < fdrcut)) > 0){
        tab1$Upper_adjusted = tab1$logFC + 
            tab1$SE * qnorm(1 - (nrow(subset(tab1,
                                             adj.P.Val < fdrcut))/dim(fite1)[[1]] * fdrcut))
        tab1$Lower_adjusted = tab1$logFC - 
            tab1$SE * qnorm(1 - (nrow(subset(tab1,
                                             adj.P.Val < fdrcut))/dim(fite1)[[1]] * fdrcut))}
    tab1$AveM = 2^(tab1$AveExpr)
    tab1$betaDiff = (2^(tab1$AveM*2^tab1$logFC))/(1+
                                                      2^(tab1$AveM*2^tab1$logFC)) - (2^tab1$AveM)/(1+2^tab1$AveM)
    return(tab1)
}
