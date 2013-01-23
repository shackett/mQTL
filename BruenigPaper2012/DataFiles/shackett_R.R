setwd("~/Desktop/BruenigPaper2012/DataFiles")

library(stringr)
library(qtl)

colinfo <- scan("segAbund.txt", what = "character", nlines = 1)
seginfo <- data.frame(date = sapply(colinfo, function(entry){unlist(strsplit(entry, "-"))[1]}), seg = sapply(colinfo, function(entry){str_extract(str_extract(entry, '-[a-z0-9]{2,10}-'), '[a-z0-9]{2,10}')}), stringsAsFactors = FALSE)
seginfo$date[grep("071013-2", colinfo)] <- "071013-2"
seginfo$date[grep("070926-2", colinfo)] <- "070926-2"
seginfo$date[grep("071029-2", colinfo)] <- "071029-2"

rminfo <- scan("rmAbund.txt", what = "character", nlines = 1)
rminfo <- data.frame(date = sapply(rminfo, function(entry){unlist(strsplit(entry, "-"))[1]}), stringsAsFactors = FALSE)
rminfo$date[grep("071013-2", rownames(rminfo))] <- "071013-2"
rminfo$date[grep("070926-2", rownames(rminfo))] <- "070926-2"
rminfo$date[grep("071029-2", rownames(rminfo))] <- "071029-2"

rmdata <-  read.delim("rmAbund.txt", header = TRUE)
rownames(rmdata) <- rmdata[,1]; rmdata <- rmdata[,-1]
rmdata <- log2(as.matrix(rmdata))

rmcol_list <- list()
for(num in 1:length(seginfo$date)){
	rmcol_list[[num]] <- c(1:length(rminfo[,1]))[rminfo$date == seginfo$date[num]]
	}

#sig associations


GWAShits <- read.delim("figure 2.csv", sep = ",", stringsAsFactors = FALSE)


abund_dataMD <- read.delim("sampleDat32.txt", header = TRUE, stringsAsFactors = FALSE)
abund_data <- read.delim("segAbund.txt", header = TRUE, stringsAsFactors = FALSE)

baseline <- 32
abund_dataMD[abund_dataMD < baseline] <- baseline
mv_cutoff <- 0.25#only use metabolites found in >=25% of samples

abundMat <- as.matrix(abund_data[,-1])
rownames(abundMat) <- abund_data[,1]
abund_dataMD <- abund_dataMD[,-1]

abund_dataMD[apply(((is.na(log_abundMatr) * 1) %*% sampleToSeg) != 0, 1, sum) == 96,]


sampleToSeg <- matrix(0, nrow = length(seginfo[,1]), ncol = length(unique(seginfo$seg)))
colnames(sampleToSeg) <- unique(seginfo$seg)
for(samp in 1:length(seginfo[,1])){
	sampleToSeg[samp,] <- ifelse(colnames(sampleToSeg) == seginfo$seg[samp], 1, 0)
	}

good_compounds <- rowSums(((abund_dataMD != baseline)*1) %*% sampleToSeg != 0) >= length(sampleToSeg[1,])*mv_cutoff
missingValFrac <- 1 - (rowSums(((abund_dataMD != baseline)*1) %*% sampleToSeg != 0)/length(sampleToSeg[1,]))[good_compounds]
	
abundMatr <- abundMat[good_compounds,]
abundMatr[abund_dataMD[good_compounds,] == baseline] <- NA
log_abundMatr <- log(abundMatr, base = 2)	


rm_normMat <- matrix(NA, ncol = length(abundMatr[1,]), nrow = length(abundMatr[,1]))
for(num in 1:length(rm_normMat[1,])){
	rm_normMat[,num] <- apply(matrix(rmdata[good_compounds, rmcol_list[[num]]], ncol = length(rmcol_list[[num]])), 1, mean) 
	}
	
normMat <- log_abundMatr - rm_normMat
	
rowMean <- apply(normMat, 1, mean, na.rm = TRUE)	
rowSS <- rowSums((normMat - t(t(rowMean)) %*% matrix(1, ncol = length(normMat[1,]), nrow = 1))^2, na.rm =TRUE)

#recalculate heritability
herit <- rep(NA, times = length(normMat[,1]))
for(met in 1:length(herit)){
  #reduce data to metabolites with only 
  segN <- t(!is.na(normMat[met,])*1) %*% sampleToSeg
  matchMAT <- (sampleToSeg * normMat[met,])[,segN >= 2]
  matchMAT <- matchMAT[!is.na(rowSums(matchMAT)) & rowSums(matchMAT) != 0,]
  matMean <- mean(matchMAT[matchMAT != 0])
  
  TSS <- sum((rowSums(matchMAT) - matMean)^2)
  
  #residuals calculated as nV - where V was adjusted with the small sample correction (n/n-1)
  RSS <- sum(segN[segN >= 2]*apply(matchMAT, 2, function(var_calc){var(var_calc[var_calc != 0])}))
  
  herit[met] <- (TSS - RSS)/TSS
  
  }

library(lme4)

for(met in 1:length(herit)){
  #reduce data to metabolites with only 
  segN <- t(!is.na(normMat[met,])*1) %*% sampleToSeg
  matchMAT <- sampleToSeg[,segN >= 2]
  anova(lm(normMat[met,rowSums(matchMAT) == 1] ~ matchMAT[rowSums(matchMAT) == 1,]))
  
  
}



  

envSSmat <- matrix(NA, ncol = length(sampleToSeg[1,]), nrow = length(normMat[,1]))
#if there is no estimate of the unconfounded SS for a segregant (i.e. no biological replicates exist), remove that fraction of the SS from the rowSS
unnestedSS <- rep(0, times = length(rowSS))

for(seg in 1:length(unique(seginfo$seg))){
	
	subData <- normMat[,c(1:length(normMat[1,]))[sampleToSeg[,seg] == 1]]
	can_calc_Venv <- rowSums(!is.na(subData)) >= 2
	if(sum(!can_calc_Venv) != 0){
		envSSmat[!can_calc_Venv,seg] <- NA
		if(sum(rowSums(!is.na(subData)) == 1) >= 2){
			unnestedSS[rowSums(!is.na(subData)) == 1] <- unnestedSS[rowSums(!is.na(subData)) == 1] + (rowSums(subData[rowSums(!is.na(subData)) == 1,], na.rm = TRUE) - rowMean[rowSums(!is.na(subData)) == 1])^2
			}
		if(sum(rowSums(!is.na(subData)) == 1) == 1){
			unnestedSS[rowSums(!is.na(subData)) == 1] <- unnestedSS[rowSums(!is.na(subData)) == 1] + (subData[rowSums(!is.na(subData)) == 1,][!is.na(subData[rowSums(!is.na(subData)) == 1,])] - rowMean[rowSums(!is.na(subData)) == 1])^2
			}	
		}
	
	envSSmat[can_calc_Venv,seg] <- apply((subData[can_calc_Venv,] - apply(subData[can_calc_Venv,], 1, mean, na.rm = TRUE) %*% matrix(1, ncol = length(subData[1,])))^2 * t(matrix(1, nrow = length(subData[1,]), ncol = 1) %*% apply(!is.na(subData[can_calc_Venv,]), 1, sum)/(apply(!is.na(subData[can_calc_Venv,]), 1, sum) - 1)), 1, sum)
		
	}
	
#met_herit <- data.frame(compound = rownames(normMat), heritability = ((rowSS - unnestedSS) - rowSums(envSSmat, 1, na.rm = TRUE))/(rowSS - unnestedSS), nassoc = 0, frac_missingVals = missingValFrac, stringsAsFactors = FALSE)
met_herit <- data.frame(compound = rownames(normMat), heritability = herit, nassoc = 0, frac_missingVals = missingValFrac, stringsAsFactors = FALSE)


GWAShitnum <- table(GWAShits$name)
#sapply(met_herit$compound, function(cmpd){unlist(strsplit(cmpd, '-[0-9]'))[1]})

for(assoc in names(GWAShitnum)){
	met_herit$nassoc[met_herit$compound == assoc] <- GWAShitnum[names(GWAShitnum) == assoc]
	}


met_herit$nassoc[met_herit$compound == "glutathione-0"] <- met_herit$nassoc[met_herit$compound == "glutathione-0"] + met_herit$nassoc[met_herit$compound == "glutathione disulfide-posi-1"]
met_herit$heritability[met_herit$compound == "glutathione-0"] <- mean(met_herit$heritability[met_herit$compound %in% c("glutathione-0", "glutathione disulfide-posi-1")])
met_herit$frac_missingVals[met_herit$compound == "glutathione-0"] <- min(met_herit$frac_missingVals[met_herit$compound %in% c("glutathione-0", "glutathione disulfide-posi-1")])
met_herit <- met_herit[!met_herit$compound == "glutathione disulfide-posi-1",]

met_herit$nassoc[met_herit$compound == "S-adenosyl-L-homoCysteine-posi-0"] <- met_herit$nassoc[met_herit$compound == "S-adenosyl-L-homoCysteine-posi-0"] + met_herit$nassoc[met_herit$compound == "S-adenosyl-L-homocysteine-nega-1"]
met_herit$heritability[met_herit$compound == "S-adenosyl-L-homoCysteine-posi-0"] <- mean(met_herit$heritability[met_herit$compound %in% c("S-adenosyl-L-homoCysteine-posi-0", "S-adenosyl-L-homocysteine-nega-1")])
met_herit$frac_missingVals[met_herit$compound == "S-adenosyl-L-homoCysteine-posi-0"] <- min(met_herit$frac_missingVals[met_herit$compound %in% c("S-adenosyl-L-homoCysteine-posi-0", "S-adenosyl-L-homocysteine-nega-1")])
met_herit <- met_herit[!met_herit$compound == "S-adenosyl-L-homocysteine-nega-1",]

met_herit$fraction_missing_lines <- as.factor(sapply(floor(met_herit$frac_missingVals * 10)*10, function(val){
	paste(val, "-", val+10, " %", sep = "")
	}))

met_herit$nassoc <- as.factor(met_herit$nassoc)

library(ggplot2)

pdf("QTLnum.pdf")
herit_plot <- ggplot(met_herit, aes(x = heritability, fill = fraction_missing_lines)) + facet_grid(nassoc ~ .)
herit_plot + geom_histogram()


herit_plot <- ggplot(met_herit, aes(x = heritability, fill = fraction_missing_lines)) + facet_grid(. ~ nassoc)
herit_plot + geom_histogram(binwidth = 0.05) + coord_flip() + scale_fill_discrete(name = 'Fraction of \nMissing Lines', h = c(200, 360)) + scale_x_continuous(name = "Heritability", limits = c(0,1), expand = c(0,0)) + scale_y_continuous(name = "Number of Metabolites") + 
  theme(axis.text.x = element_text(size = 12, face = "bold"), axis.text.y = element_text(size = 12, face = "bold"), strip.text.x = element_text(size = 16, colour = "RED", face = "bold")) +
  theme(axis.title.x = element_text(size = 12, face = "bold"), axis.title.y = element_text(size = 12, face = "bold"), legend.position = "top", legend.title = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size = 10))
  
herit_plot <- ggplot(met_herit, aes(x = 1, y = heritability, fill = fraction_missing_lines)) + facet_wrap( ~ nassoc, ncol = 2)
herit_plot + geom_dotplot(binaxis = "y", stackgroups = TRUE, binwidth = 0.05, method = "histodot")


glm(as.numeric(met_herit$nassoc)-1 ~ met_herit$heritability, family = poisson)
glm((as.numeric(met_herit$nassoc)-1)/max(as.numeric(met_herit$nassoc)-1) ~ met_herit$heritability, family = binomial)

disc_plot <- ggplot(met_herit, aes(x = fraction_missing_lines, fill = nassoc))
disc_plot + geom_bar()

disc_plot <- ggplot(met_herit, aes(x = fraction_missing_lines)) + facet_grid(nassoc ~ .)
disc_plot + geom_bar()
dev.off()


###### read in supplementary metabolomics data #########
library(qvalue); library(gplots)

suppMSdata <- read.delim("~/Desktop/krugMetabData/srh_bruenigDat.txt", stringsAsFactors = FALSE)
suppMSlist <- list()
suppMSlist$mets <- suppMSdata[,1:5]; colnames(suppMSlist$mets) <- c("Instrument", "groupRank", "Ngood", "medianRT", "Compound")
suppMSlist$abund <- suppMSdata[,-c(1:5)]
rownames(suppMSlist$mets) <- rownames(suppMSlist$abund) <- apply(suppMSlist$mets, 1, function(x){paste(x[5], x[2], sep = "-")})
suppMSlist$samples <- data.frame(sampleName = colnames(suppMSlist$abund), cond = c(rep("BL", times = 4), rep("BY", times = 2), rep("BYira2RM", times = 2), rep("BYura3del", times = 2), rep("NULL", times = 2), rep("RM", times = 2), rep("RMira2BY", times = 2)), stringsAsFactors = FALSE)

ODdf <- data.frame(conditions = c("BL", "BY", "BYira2RM", "BYura3del", "NULL", "RM", "RMira2BY"), od1 = c(NA, 0.321, 0.497, 0.493, NA, 0.649, 0.749), od2 = c(NA, 0.336, 0.568, 0.532, NA, 0.693, 0.759), stringsAsFactors = FALSE)
ODdf$average <- (ODdf$od1 + ODdf$od2)/2
suppMSlist$samples$OD <- sapply(suppMSlist$samples$cond, function(cond){
	ODdf$average[ODdf$conditions == cond]
	})

#scale biological samples by OD
suppMSlist$abund[suppMSlist$abund < baseline] <- NA

suppMSlist$log_abund <- log2(suppMSlist$abund)

average_sf <- mean(suppMSlist$samples$OD, na.rm = TRUE)
norm_mat <- suppMSlist$log_abund[,!is.na(suppMSlist$samples$OD)] - t(t(rep(1, times = length(suppMSlist$abund[,1])))) %*% t(log2(suppMSlist$samples$OD/average_sf)[!is.na(suppMSlist$samples$OD)])

suppMSlist$log_abund_norm <- suppMSlist$log_abund
suppMSlist$log_abund_norm[,!is.na(suppMSlist$samples$OD)] <- norm_mat; suppMSlist$log_abund_norm <- as.matrix(suppMSlist$log_abund_norm)
suppMSlist$log_abund_norm[is.na(suppMSlist$log_abund_norm)] <- log2(baseline)
suppMSlist$log_abund_norm[suppMSlist$log_abund_norm < log2(baseline)] <- log2(baseline)

#test to see how many features have null signal (filter-culture using just media) comparable to sample signal
#determine sd(log_val) from all pairs

residual_mat <- matrix(NA, ncol = length(unique(suppMSlist$samples$cond[!(suppMSlist$samples$cond == "BL")])), nrow = length(suppMSlist$abund[,1]))
rownames(residual_mat) <- rownames(suppMSlist$abund); colnames(residual_mat) <- unique(suppMSlist$samples$cond[!(suppMSlist$samples$cond == "BL")])

for(cond in unique(suppMSlist$samples$cond[!(suppMSlist$samples$cond == "BL")])){
	residual_mat[,colnames(residual_mat) == cond] <- apply(suppMSlist$log_abund_norm[,suppMSlist$samples$cond == cond], 1, diff)
	}

#fit the metabolite-specific MSE from row residuals 
met_sd <- apply(abs(residual_mat) * 2, 1, function(a_row){median(a_row[a_row != 0])})

#average biological signal < #t

n1 <- sum((suppMSlist$samples$cond %in% c("NULL")))
n2 <- sum(!(suppMSlist$samples$cond %in% c("NULL", "BL")))

tstat_null <- (apply(suppMSlist$log_abund_norm[,!(suppMSlist$samples$cond %in% c("NULL", "BL"))], 1, mean) - apply(suppMSlist$log_abund_norm[,suppMSlist$samples$cond %in% c("NULL")], 1, mean))/(((n1 + n2 -2)*met_sd^2)/(n1 + n2 -2)*sqrt(1/n1 + 1/n2))
pt_null <- 1 - pt((tstat_null), (n1 + n2 - 2))


#
names(qvalue(pt_null)$qvalues)[qvalue(pt_null)$qvalues > 0.05]

#4-Pyridoxic acid-1: B6 provided
#leucine/isoleucine-0: provided
#Acetylcarnitine DL-1 ?
#nicotinate: provided
#deoxyribose-phosphate-1: uracil breakdown?
#pantothenate-1: B5 provided
#D-glucono-delta-lactone-6-phosphate-0

heatmap.2(suppMSlist$log_abund_norm, trace = "none")

#reduced data to determine if ura3 deletiom differs from BY

met_lm <- matrix(NA, ncol = 4, nrow = length(suppMSlist$log_abund_norm[,1]))
rownames(met_lm) <- rownames(suppMSlist$log_abund_norm)

for(met in 1:length(suppMSlist$log_abund_norm[,1])){
	met_lm[met,1:4] <- summary(lm(suppMSlist$log_abund_norm[,suppMSlist$samples$cond %in% c("BY", "BYura3del")][met,] ~ factor(suppMSlist$samples[suppMSlist$samples$cond %in% c("BY", "BYura3del"),]$cond)))$coef[2,]
	
	}

qvalue(met_lm[,4][!is.nan(met_lm[,4])])$qvalue[qvalue(met_lm[,4][!is.nan(met_lm[,4])])$qvalue < 0.05]
met_lm[,3][rownames(met_lm) %in% names(qvalue(met_lm[,4][!is.nan(met_lm[,4])])$qvalue[qvalue(met_lm[,4][!is.nan(met_lm[,4])])$qvalue < 0.05]
)]

#orotate
#N-carbamoyl-aspartate
#orotidine-5'phosphate
#IDP

met_lm[,3][rownames(met_lm) %in% names(qvalue(met_lm[,4][!is.nan(met_lm[,4])])$qvalue[qvalue(met_lm[,4][!is.nan(met_lm[,4])])$qvalue < 0.05])]

met_herit$spec_compound = sapply(met_herit$compound, function(cmpd){
	paste(strsplit(cmpd, "-")[[1]][1:(length(strsplit(cmpd, "-")[[1]])-1)], collapse = "-")
	})

suppMSlist$mets$Compound %in% met_herit$spec_compound
met_herit$spec_compound %in% suppMSlist$mets$Compound

suppMSlist$mets$Compound[!(suppMSlist$mets$Compound %in% met_herit$spec_compound)]
met_herit$spec_compound[!(met_herit$spec_compound %in% suppMSlist$mets$Compound )]





	