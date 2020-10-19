lm_test_interaction <- function(N, 
                                bRace1,
                                bRace2,
                                bSex,
                                bAge,
                                bSmoke,
                                bMedicaid,
                                bGest,
                                bZ,
                                Binflam) {
    
    threeway = function(N,prob){
        u = runif(N,0,1)
        prob = cumsum(prob)
        a = ifelse(u < prob[1],'Race0',
                   ifelse(u > prob[2],'Race2','Race1'))
        return(a)
    }
    
    race = as.factor(threeway(N,prob = c(.615,.295,.09)))
    sex = rbinom(N,1,prob = 1-.475)
    age = rnorm(N,29.6,6.61)
    smoke = (rbinom(N,1,prob = .103))
    insurance = (rbinom(N,1,prob=.338))
    gest = rnorm(N,182.5,9.17)
    zscore = rnorm(N)
    inflame = (rbinom(N,1,prob=.336))
    
    
    yvar <- sqrt(1 - bRace1^2 - bRace2^2 + bSex^2 +
                     bAge^2 + bSmoke^2 + bMedicaid^2 + 
                     bGest^2 + bZ^2 + Binflam) # residual variance
    y <- rnorm(N, bRace1*(ifelse(race == 'Race1',1,0)) + 
                   bRace2*(ifelse(race == 'Race2',1,0)) + 
                   bSex*sex +
                   bAge*age +
                   bSmoke*(smoke) +
                   bMedicaid*(insurance) +
                   bGest*gest +
                   bZ*zscore +
                   Binflam*(inflame), 
               yvar)
    model <- lm(y ~ race + sex + age + smoke + insurance + gest + zscore + inflame)
    
    # pull output from model (two main effects and interaction)
    sig = coef(summary(model))[-1,4] < .05
    
    return(sig)
}

race1_coef = c(.2,.3,.5)
race2_ceof = c(.3,.5,.7)
sex_coef = c(.3,.5,.7)
age_coef = c(-0.01,-0.05,-.1)
smoke_coef = c(.2,.3,.5)
insur_coef = c(.4,.5,.6)
gest_coef = c(-0.01,-0.05,-0.1)
z_coef = c(-0.05,-.1,-.2)
inflam_coef = c(-0.03,-0.05,-0.07)
a = expand.grid(race1_coef,
                race2_ceof,
                sex_coef,
                age_coef,
                smoke_coef,
                insur_coef,
                gest_coef,
                z_coef,
                inflam_coef)
colnames(a) = c('Race1_Beta',
                'Race2_Beta',
                'Sex_Beta',
                'Age_Beta',
                'Smoker_Beta',
                'Insurace_Beta',
                'Gestational_Beta',
                'Z_Beta',
                'Inflammation_Beta')

df = as.data.frame(matrix(nrow = nrow(a),
                         ncol = 9))
colnames(df) =  c('raceRace1',
                  'raceRace2',
                  'sex',
                  'age',
                  'smoke',
                  'insurance',
                  'gest',
                  'zscore',
                  'inflame')
a = cbind(a,df)
for (i in seq(1,nrow(a),by=1)){
    print(i)
    ppp = a[i,]
    qqq = pbapply::pbreplicate(1e3,
              lm_test_interaction(N = 379, 
                              bRace1 = as.numeric(ppp[1]),
                              bRace2 = as.numeric(ppp[2]),
                              bSex = as.numeric(ppp[3]),
                              bAge = as.numeric(ppp[4]),
                              bSmoke = as.numeric(ppp[5]),
                              bMedicaid = as.numeric(ppp[6]),
                              bGest = as.numeric(ppp[7]),
                              bZ = as.numeric(ppp[8]),
                              Binflam = as.numeric(ppp[9])))
    a[i,10:ncol(a)] = as.numeric(rowMeans(qqq))
}

colnames(a.tot)[10:18] = paste0('Power_',
                                c('raceRace1',
                                  'raceRace2',
                                  'sex',
                                  'age',
                                  'smoke',
                                  'insurance',
                                  'gest',
                                  'zscore',
                                  'inflame'))
data.table::fwrite(a.tot,'C:/Users/Arjun/OneDrive - University of North Carolina at Chapel Hill/Multiomics - autism/power_linear_regression.tsv',
       sep='\t',quote=F,col.names=T,row.names=F)


##### RNA-seq power
require(RnaSeqSampleSizeData)
require(RnaSeqSampleSize)
fold_change = seq(1,2,by = .1)
lambda = c(5,10,20)
phi0 = c(.5,3)
eee = expand.grid(fold_change,
                  lambda,
                  phi0)
for (i in 1:nrow(eee)){
    print(i)
    eee$Power[i] = est_power(n=379,
                             rho=eee[i,1],
                             lambda0=eee[i,2],
                             phi0=eee[i,3],
                             alpha = .05,
                             f=.01,
                             m = 12000,
                             m1 = .1*12000)
       }

colnames(eee)[1:3] = c('Fold change',
                       'Mean',
                       'Dispersion')
eee = eee[order(eee$Mean,eee$Dispersion),]
require(ggplot2)
eee$Group = paste0('Mean = ',eee$Mean,'; ',
                   'Dispersion = ',eee$Dispersion)

eee$Mean = paste0('Mean = ',eee$Mean)
eee$Dispersion = paste0('Dispersion = ',eee$Dispersion)
eee = eee[order(eee$Group),]
mrna_powerplot = ggplot(data = eee,
                        aes(x = `Fold change`,
                            y = Power,
                            color = Group)) +
    geom_line() +
    geom_point(size = 4) +
    facet_wrap(~Mean*Dispersion) +
    theme_minimal() +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=24),
          plot.title = element_text(size = 30),
          legend.title=element_text(size=24),
          legend.text=element_text(size=24),
          strip.text = element_text(size=16),
          panel.spacing=unit(1, "lines"),
          panel.border = element_rect(color = "grey",
                                      fill = NA, size = .1),
          legend.position="bottom", legend.box = "horizontal") +
    guides(color = F) +
    scale_color_brewer(palette = 'Dark2') +
    xlab('Fold change (mRNA DEG analysis)')

##### miRNA-seq power
require(RnaSeqSampleSizeData)
require(RnaSeqSampleSize)
fold_change = seq(1,2,by = .1)
lambda = c(5,10,20)
phi0 = c(.5,3)
eee = expand.grid(fold_change,
                  lambda,
                  phi0)
for (i in 1:nrow(eee)){
    print(i)
    eee$Power[i] = est_power(n=379,
                             rho=eee[i,1],
                             lambda0=eee[i,2],
                             phi0=eee[i,3],
                             alpha = .05,
                             f=.01,
                             m = 2000,
                             m1 = .1*2000)
}

colnames(eee)[1:3] = c('Fold change',
                       'Mean',
                       'Dispersion')
eee = eee[order(eee$Mean,eee$Dispersion),]
require(ggplot2)
eee$Group = paste0('Mean = ',eee$Mean,'; ',
                   'Dispersion = ',eee$Dispersion)

eee$Mean = paste0('Mean = ',eee$Mean)
eee$Dispersion = paste0('Dispersion = ',eee$Dispersion)
eee = eee[order(eee$Group),]
mirna_powerplot = ggplot(data = eee,
                        aes(x = `Fold change`,
                            y = Power,
                            color = Group)) +
    geom_line() +
    geom_point(size = 4) +
    facet_wrap(~Mean*Dispersion) +
    theme_minimal() +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=24),
          plot.title = element_text(size = 30),
          legend.title=element_text(size=24),
          legend.text=element_text(size=24),
          strip.text = element_text(size=16),
          panel.spacing=unit(1, "lines"),
          panel.border = element_rect(color = "grey",
                                      fill = NA, size = .1),
          legend.position="bottom", legend.box = "horizontal") +
    guides(color = F) +
    scale_color_brewer(palette = 'Dark2') +
    xlab('Fold change (miRNA DEG analysis)')

## ewas power - collected from https://epigenetics.essex.ac.uk/shiny/EPICDNAmPowerCalcs/ shiny app
power = data.table::fread('C:/Users/Arjun/OneDrive - University of North Carolina at Chapel Hill/Multiomics - autism/ewas_power.csv')
power$`Mean Difference` = as.factor(power$`Mean Difference`)
power$`Mean Difference` = paste0(power$`Mean Difference`,'%')
ewas_powerplot = ggplot(data = power,
                         aes(x = `Proportion of sites`,
                             y = Power,
                             color = `Mean Difference`)) +
    geom_line() +
    geom_point(size = 4) +
    theme_minimal() +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=24),
          plot.title = element_text(size = 30),
          legend.title=element_text(size=24),
          legend.text=element_text(size=24),
          strip.text = element_text(size=16),
          panel.spacing=unit(1, "lines"),
          panel.border = element_rect(color = "grey",
                                      fill = NA, size = .1),
          legend.position="bottom", legend.box = "horizontal") +
    scale_color_brewer(palette = 'Dark2')

setwd('C:/Users/Arjun/OneDrive - University of North Carolina at Chapel Hill/Multiomics - autism/')
require(cowplot)
a = plot_grid(mrna_powerplot, mirna_powerplot, ewas_powerplot,
              ncol = 3,
              labels = c('A','B','C'),
              label_size = 32)
ggsave(plot = a,
       filename = 'power_plot.png',
       width = 24,
       height = 10)
