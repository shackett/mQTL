setwd("~/Desktop/BrueingNMR/2013_01_17_IRA2_switch")
library(reshape)
library(ggplot2)

#import an ROI summary table
#maximum height
NMRdata <- read.table("roiSummary_sean.txt", header = T, stringsAsFactors = F)
#peak area
NMRarea <- read.table("roiArea.txt", header = TRUE, stringsAsFactors = F)
#peak area and height are proportional across all concentrations -> use peak height

#normalize signal against DSS to account for smolarity differences
nNMRdata <- NMRdata[,-1]
nNMRdata <- nNMRdata/NMRdata$DSS

sampleTable <- data.frame(t(sapply(NMRdata$File, function(x){strsplit(x, split = '_')[[1]][4:5]})), stringsAsFactors = F)
colnames(sampleTable) <- c("genotype", "time")
sampleTable$time <- sub('.ucsf', '', sampleTable$time)

#input actual sampled time-points
timeCorr <- data.frame(timePoint = c("t0", "t1", "t2", "t3", "t4"), time = c(0, 97, 244, 395, 558))

#input OD                                            
sampleTable$OD <- c(0.05, 0.077, 0.188, 0.26, 0.514, 0.05, 0.075, 0.128, 0.246, 0.482, 0.05, 0.089, 0.133, 0.434, 0.887, 0.05, 0.094, 0.207, 0.476, 0.923, rep(NA, 12))
sampleTable$time <- as.numeric(sapply(sampleTable$time, function(x){if(x %in% timeCorr$timePoint){timeCorr$time[timeCorr$timePoint == x]}else{x}}))
sampleTable$is_sample <- T; sampleTable$is_sample[grep('std', sampleTable$genotype)] <- F

od_plot <- ggplot(sampleTable[sampleTable$is_sample,], aes(x = time, y = log(OD), group = genotype, col = genotype))
od_plot + geom_line()

#determine concentration of EtOH/Glucose as a function of NMR peak height using standards
standards <- sampleTable[!sampleTable$is_sample,]

#t0 molarity = 2g/100mL glucose * 1 mol /180.16g * 10 = 0.11M glucose, 0.22M EtOH if glucose is fully fermented
stdConc = data.frame(glucose = rep(c(0.11, 0.11 * (1/2), 0.11 * (1/2)^2, 0.11 * (1/2)^2*(1/5), 0.11 * (1/2)^2*(1/5)^2, 0.11 * (1/2)^2*(1/5)^3), each = 2), ethanol = rep(c(0.22, 0.22 * (1/2), 0.22 * (1/2)^2, 0.22 * (1/2)^2*(1/5), 0.22 * (1/2)^2*(1/5)^2, 0.22 * (1/2)^2*(1/5)^3), each = 2))

#where do actual samples fall relative to standards
plot(nNMRdata$Gluc.1, col = ifelse(sampleTable$is_sample, "RED", "BLUE"))
plot(nNMRdata$EtOH, col = ifelse(sampleTable$is_sample, "RED", "BLUE"))

#the two glucose peaks are essentially perfectly correlated
plot(nNMRdata$Gluc.1 ~ nNMRdata$Gluc.2)



#peak height scales linearly with abundance 
plot(nNMRdata[grep('std', sampleTable$genotype),]$Gluc.1 ~ stdConc$glucose)
plot(nNMRdata$EtOH[!sampleTable$is_sample] ~ stdConc$ethanol)
#one of the highest concentrations standards was obviously inappropriatley prepared and will be excluded


ethanol_fit <- lm(stdConc$ethanol[-2] ~ nNMRdata[!sampleTable$is_sample,]$EtOH[-2])$coef
#gorgeous linear fit
plot(stdConc$ethanol[-2] ~ nNMRdata[!sampleTable$is_sample,]$EtOH[-2], xlab = "peak height", ylab = "concentration", main =  "Predicting ethanol concentration from NMR peak heigt", pch = 16, cex = 2)
abline(a = ethanol_fit[1], b = ethanol_fit[2], lwd = 3, col = "RED")

#determine the growth rate of each sample - linear fit through log(OD) ~ time
samples <- data.frame(sampleTable[sampleTable$is_sample,], nNMRdata[sampleTable$is_sample,colnames(nNMRdata) %in% c("EtOH", "Gluc.1", "Gluc.2")])
#remove a clear outlier (probably contaminated)
samples <- samples[rownames(samples) != "2013_01_16_BY_t3.ucsf",]

samples$background <- "RM"; samples$background[grep('^BY', samples$genotype)] <- "BY"
samples$ira2 <- ifelse(samples$genotype %in% c("BY", "RMira2BY"), "BY", "RM")

  
#create a design matrix for regression - OD ~ intercept + G*slope + residual
Design <- cbind(1, ifelse(samples$genotype == "BY", 1, 0)*samples$time, ifelse(samples$genotype == "BYira2RM", 1, 0)*samples$time, ifelse(samples$genotype == "RM", 1, 0)*samples$time, ifelse(samples$genotype == "RMira2BY", 1, 0)*samples$time)  

D <- model.matrix(log(samples$OD) ~ samples$background*samples$time + samples$ira2*samples$time)[,-c(2,4)]
summary(lm(log(samples$OD) ~ D + 0))
GRcoef <- lm(log(samples$OD) ~ Design + 0)$coef; names(GRcoef) <- c("offset", "BY", "BYira2RM", "RM", "RMira2BY")
anova(lm(log(samples$OD) ~ samples$genotype*samples$time))
anova(lm(log(samples$OD) ~ samples$genotype + samples$background*samples$time + samples$ira2*samples$time))

fittedAbund <- as.data.frame(sapply(unique(samples$background), function(geno){
  exp(GRcoef[1] + GRcoef[names(GRcoef) == geno]*seq(0, max(samples$time)))
  }))
fittedAbund$time <- seq(0, max(samples$time))


mfittedAbund <- melt(fittedAbund, id.vars = "time")
colnames(mfittedAbund) <- c("time", "genotype", "OD")

growth_plot <- ggplot(samples, aes(x = time, y = OD, group = genotype, col = genotype))
growth_plot + geom_point(size = 4) + geom_line(data = mfittedAbund, size = 1.5)

#the amount of ethanol excreted, glucose consumed at any interval should be proprotional to the density of cells at this point
#if the exponent of ethanol excretion differs from the growth-rate, this is evidence of differential rates of EtOH excretion per cell

growth_fit <- Design %*% GRcoef
plot(log(samples$OD) ~ growth_fit, col = factor(samples$genotype))

samples$ethanolConc <- ethanol_fit[1] + ethanol_fit[2]*samples$EtOH
plot(log(samples$ethanolConc) ~ growth_fit, col = factor(samples$genotype))
summary(lm(log(samples$ethanolConc) ~ Design2))

summary(lm(log(ethanolConc) ~ background*time + ira2*time, data = samples))
anova(lm(log(ethanolConc) ~ factor(genotype) + factor(background)*time + factor(ira2)*time, data = samples[rownames(samples) != "2013_01_16_RMira2BY_t1.ucsf",]))

D <- model.matrix(log(ethanolConc) ~ background*time + ira2*time, data = samples)[,-c(2,4)]
#not looking at background/ira2 main effects (i.e. the effect on the intercept) because these should be equal by virtue of study design
summary(lm(log(samples$ethanolConc) ~ D + 0))
#remove the first time point, since this is meaningless becasuse of transfer
summary(lm(log(samples$ethanolConc)[as.data.frame(D)$time != 0] ~ D[as.data.frame(D)$time != 0,] + 0))
ethanolExCoef <- lm(log(samples$ethanolConc)[as.data.frame(D)$time != 0] ~ D[as.data.frame(D)$time != 0,] + 0)$coef

plot(log(samples$ethanolConc) ~ samples$time, col = )
D[as.data.frame(D)$time != 0,] %*% lm(log(samples$ethanolConc)[as.data.frame(D)$time != 0] ~ D[as.data.frame(D)$time != 0,] + 0)$coef

fittedAbund <- as.data.frame(sapply(unique(samples$genotype), function(geno){
  exp(ethanolExCoef[1] + ethanolExCoef[2]*seq(0, max(samples$time)) + 
   ethanolExCoef[3]*ifelse(geno %in% c("RM", "RMira2BY"), 1, 0)*seq(0, max(samples$time)) + 
   ethanolExCoef[4]*ifelse(geno %in% c("RM", "BYira2RM"), 1, 0)*seq(0, max(samples$time))
  )
  }))
fittedAbund$time <- seq(0, max(samples$time))

mfittedAbund <- melt(fittedAbund, id.vars = "time")
colnames(mfittedAbund) <- c("time", "genotype", "ethanolConc")

growth_plot <- ggplot(samples, aes(x = time, y = ethanolConc, group = genotype, col = genotype))
growth_plot + geom_point(size = 3, alpha = 0.8 ) + geom_line(data = mfittedAbund, size = 1, alpha = 0.8)



