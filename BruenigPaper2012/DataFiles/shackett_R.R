setwd("~/Desktop/mQTL/BruenigPaper2012/DataFiles")

library(stringr)
library(qtl)
library(ggplot2)

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

#abund_dataMD[apply(((is.na(log_abundMatr) * 1) %*% sampleToSeg) != 0, 1, sum) == 96,]


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

met_herit <- data.frame(compound = rownames(normMat), heritability = herit, nassoc = 0, heritExplained = NA, frac_missingVals = missingValFrac, stringsAsFactors = FALSE)



#overwrite incorrect annotations
GWAShits$name[GWAShits$name == "Unk_131@10-+149-0"] <- "Unk_131@10--149-0"
GWAShits$name[GWAShits$name == "Unk_135@10-+225.001-0"] <- "Unk_135@10--225.001-0"
GWAShits$name[GWAShits$name == "Unk_239@10-+337-0"] <- "Unk_239@10--337-0"
GWAShits$name[GWAShits$name == "Unk_267@20-+387-0"] <- "Unk_267@20--387-0"
GWAShits$name[GWAShits$name == "thiamine-1"] <- "thiamine-0"
GWAShits$name[GWAShits$name == "lysine-G"] <- "lysine-0"
GWAShits$name[GWAShits$name == "hexose-phosphate-2"] <- "hexose-phosphate-2-0"
GWAShits$name[GWAShits$name == "ADP-nega-0"] <- "ADP-posi-0"

GWAShitnum <- table(GWAShits$name)

#table of QTLs

table(table(GWAShits$name[-(grep('Unk', GWAShits$name))]))
#glut are combined so -2 ones and +1 2
table(table(GWAShits$name[(grep('Unk', GWAShits$name))]))

colinfo <- scan("josh_formatted_genotypes.txt", what = "character", nlines = 1)
segGeno <- read.delim("josh_formatted_genotypes.txt", header = TRUE, sep = "\t", stringsAsFactors = F)

#determine genotype corresponding to QTL peak
#marker is indicated by "X" column in GWAShits, matches rowname of genotypeMat
segName <- colinfo[grep('[0-9]+_[0-9]', colinfo)]
rownames(segGeno) <- segGeno[,1]

segMeta <- segGeno[,!(colinfo %in% segName)]
segGeno <- segGeno[,colinfo %in% segName]

collapsedSegName <- sapply(segName, function(segname){
  paste(strsplit(segname, split = '_')[[1]], collapse = "")
  })

#missing genotypes for 2 segregants
#table(colnames(sampleToSeg) %in% collapsedSegName)

segGenoQTL <- matrix(NA, ncol = sum(colnames(sampleToSeg) %in% collapsedSegName), nrow = length(GWAShits[,1]))
rownames(segGenoQTL) <- GWAShits$X; colnames(segGenoQTL) <- colnames(sampleToSeg)[colnames(sampleToSeg) %in% collapsedSegName]
  
for(a_row in 1:length(segGenoQTL[,1])){
  matchMeta <- segMeta[segMeta$RQTL_name %in% GWAShits$X[a_row],]
  orderedGeno <- segGeno[segMeta$chromosome == matchMeta$chromosome,][order(abs(matchMeta$position - segMeta$position[segMeta$chromosome == matchMeta$chromosome]), decreasing = FALSE),]
  calledGeno <- apply(orderedGeno, 2, function(genoCall){
    genoCall[genoCall != 2][1]
    })
  
  segMatches <- unlist(lapply(colnames(segGenoQTL), function(segMatch){
    c(1:length(collapsedSegName))[collapsedSegName == segMatch]
    }))
  
  segGenoQTL[a_row,] <- calledGeno[segMatches]
  }



GWAShits$effectSize <- NA
GWAShits$varExplained <- NA

#reduce sampleSeg to align it to genotyped segregants

sampleToSegGen <- sampleToSeg[,sapply(colnames(segGenoQTL), function(a_seg){
  c(1:length(sampleToSeg[1,]))[colnames(sampleToSeg) == a_seg]
  })]
  

#unique compounds 
unqSigCmpds <- unique(GWAShits$name)[unique(GWAShits$name) %in% rownames(normMat)]
warnbox_plotList <- list()

for(cmpd in unqSigCmpds){
  cmpdAbund <- normMat[rownames(normMat) == cmpd,]
  if(cmpd %in% c("glutathione-0", "glutathione disulfide-posi-1")){
    cmpd <- c("glutathione-0", "glutathione disulfide-posi-1")
    }
  
  QTLgenotypes <- t(segGenoQTL[GWAShits$name %in% cmpd,] %*% t(sampleToSegGen))
  #ls regression for effect size
  GWAShits$effectSize[GWAShits$name %in% cmpd] <- summary(lm(cmpdAbund ~ QTLgenotypes))$coef[-1,1]
  #anova to determine the fraction of variance explained by all QTLs
  #GWAShits$varExplained[GWAShits$name %in% cmpd] <- anova(lm(cmpdAbund ~ QTLgenotypes))[1,2]/sum(anova(lm(cmpdAbund ~ QTLgenotypes))[,2])
  genSS <- anova(lm(as.formula(paste("cmpdAbund", paste(paste("QTLgenotypes[,", c(1:length(QTLgenotypes[1,])), "]", sep = ""), collapse = " + "), sep = ' ~ '))))[,2]
  GWAShits$varExplained[GWAShits$name %in% cmpd] <- genSS[-length(genSS)]/sum(genSS)
  
  
  #generate a boxplot of
  
  boxplotGenoPrep <- QTLgenotypes[,order(GWAShits$varExplained[GWAShits$name %in% cmpd], decreasing = TRUE)]
  boxplotGenoPrep[boxplotGenoPrep == 0] <- "BY"; boxplotGenoPrep[boxplotGenoPrep == 1] <- "RM"

  boxplot_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background =  element_blank(), legend.position = "none", 
  panel.grid.minor =  element_blank(), panel.grid.major.y =  element_blank(), axis.text.y = element_text(size = 20, face = "bold"), axis.ticks.y = element_line(size = 1), axis.ticks.x = element_blank(), axis.line = element_blank()) 
  
  if(is.vector(boxplotGenoPrep)){
      boxplotDF <- data.frame(abundance = cmpdAbund, genotype = boxplotGenoPrep)  
  }else{
      boxplotDF <- data.frame(abundance = cmpdAbund, genotype = apply(boxplotGenoPrep, 1, paste, collapse = "/"))
    }
  
  #bootstrap to determine a 95% CI for each genotypes' median
  B <- 10000
  unq_geno <- unique(boxplotDF$genotype)
  genoCI <- data.frame(genotype = unq_geno, ymin = NA, ymax = NA)
  for(gen in unq_geno){
    genAbunds <- boxplotDF$abundance[boxplotDF$genotype == gen][!is.na(boxplotDF$abundance[boxplotDF$genotype == gen])]
    genMedians <- sapply(1:B, function(b){
      median(sample(genAbunds, length(genAbunds), replace = T))
      })
    genoCI[genoCI$genotype == gen, c(2,3)] <- sort(genMedians)[c(B*0.025, B*0.975)]
    }
  
  boxplotPlot <- ggplot(boxplotDF, aes(x = genotype)) + boxplot_theme
  box_plotList[[cmpd[1]]] <- boxplotPlot + geom_violin(aes(y = abundance, fill = "RED")) + geom_vline(aes(xintercept = 0), size = 1) + scale_x_discrete(name = "Segregant Genotype", expand = c(0,0)) + scale_y_continuous(name = "Relative abundance (log2)") + geom_errorbar(data = genoCI, aes(x = genotype, ymin = ymin, ymax = ymax), size = 1)
  
} 

plot(GWAShits$varExplained ~ GWAShits$effectSize)




#bin ahead of time - determine how many compounds are in each bin and the max bin



#bin heritability by 0.05 on [0,1] - i.e. 20 bins
herit_frame <- data.frame(herit = c(0, 0, 0, 0.2, 0.2, 0.2, 0.2, 0.2), frac = c(0.2, 0.4, 0.4, 0.3, 0.1, 0.6, 0.5, 0.5), metStack = c(1,1,1,1,1,1,2,2), color = c("A", "B", "R", "A", "B", "R", "A", "R"))

herit_plot <- ggplot(herit_frame, aes(x = factor(1), fill = color, y = frac))
herit_plot + geom_bar(width = 1) + facet_grid(metStack ~ herit) + coord_polar(theta = "y")


#"S-adenosyl-L-homoCysteine-posi-0" tossed because of too few segregants

for(assoc in names(GWAShitnum)){
  met_herit$nassoc[met_herit$compound == assoc] <- GWAShitnum[names(GWAShitnum) == assoc]
  met_herit$heritExplained[met_herit$compound == assoc] <- paste(GWAShits$varExplained[GWAShits$name == assoc], collapse ='_')
  if(sum(met_herit$compound == assoc) == 0){print(paste(assoc, "missing"))}
}



qplot(GWAShits$effectSize)
qplot(met_herit$heritExplained/met_herit$heritability)
qplot(GWAShits$varExplained)



length(met_herit[grep('Unk', met_herit$compound),][,1])





#### look at overlap of Zhu/Schadt NMR data and ours ####
met_herit$Zhu <- "Not-Measured"

met_herit$Zhu[met_herit$compound %in% c("ATP-nega-0", "ADP-posi-0", "orotate-0", "dihydroorotate-0", "acetyl-CoA-posi-1", "alanine-0",
                                       "arginine-0", "aspartate-0", "phenylpyruvate-0", "phenylalanine-0", "phosphoenolpyruvate-0",
                                       "pyruvate-0", "S-adenosyl-L-methionine-0", "UDP-D-glucose-0", "S-adenosyl-L-homocysteine-nega-1",
                                        "glutamate-0", "glutamine-0", "glutathione-0", "inosine-0", "leucine/isoleucine-0", "lysine-0",
                                        "N-acetyl-glutamate-0", "NAD+_posi-0", "serine-0", "threonine-0", "tryptophan-0", "valine-0"
                                        )] <- "Measured"

met_herit$Zhu[met_herit$compound %in% c("phenylpyruvate-0", "alanine-0", "arginine-0", "N-acetyl-glutamate-0", "orotate-0", "dihydroorotate-0", "S-adenosyl-L-methionine-0",
"S-adenosyl-L-homocysteine-nega-1", "leucine/isoleucine-0", "threonine-0", "valine-0", "lysine-0")] <- "mQTL Found"

table(met_herit$Zhu)

#our study quantified 27/X metabolites found in Zhu et al., with 12/16 mQTL-associated metabolites measureable on our platform


met_herit$nassoc[met_herit$compound == "glutathione-0"] <- met_herit$nassoc[met_herit$compound == "glutathione-0"] + met_herit$nassoc[met_herit$compound == "glutathione disulfide-posi-1"]
met_herit$heritability[met_herit$compound == "glutathione-0"] <- mean(met_herit$heritability[met_herit$compound %in% c("glutathione-0", "glutathione disulfide-posi-1")])
met_herit$frac_missingVals[met_herit$compound == "glutathione-0"] <- min(met_herit$frac_missingVals[met_herit$compound %in% c("glutathione-0", "glutathione disulfide-posi-1")])
met_herit <- met_herit[!met_herit$compound == "glutathione disulfide-posi-1",]

met_herit$fraction_missing_lines <- as.factor(sapply(floor(met_herit$frac_missingVals * 10)*10, function(val){
	paste(val, "-", val+10, " %", sep = "")
	}))

met_herit$nassoc <- as.factor(met_herit$nassoc)


pdf("QTLnum.pdf")
herit_plot <- ggplot(met_herit, aes(x = heritability, y = fill = fraction_missing_lines)) + facet_grid(nassoc ~ .)
herit_plot + geom_histogram()

#reduced code
hotCols <- colorRampPalette(c("lightgoldenrodyellow", "firebrick1"))
assocCols <- hotCols(max(as.numeric(met_herit$nassoc)-1)+1)


histodot_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_blank(), legend.position = "top", 
                        panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line = element_blank()) 

herit_plot <- ggplot(met_herit, aes(x = heritability, fill = nassoc)) 
herit_plot <- herit_plot + geom_dotplot(binwidth = 0.05, method = "histodot", binpositions = "all") + histodot_theme + geom_hline(aes(yintercept = 0), size = 0.5) + geom_vline(aes(xintercept = 0), size = 0.5) +
  scale_x_continuous(name = "Heritability", expand = c(0,0)) + scale_y_discrete(name = "Number of Metabolites", expand = c(0,0)) +
  scale_fill_manual(name = 'Number of QTLs', values = c("0" = assocCols[1], "1" = assocCols[2], "2" = assocCols[3], "3" = assocCols[4], "4" = assocCols[5]))


ggsave(plot = herit_plot, file = "heritPlot.pdf", width = 10, height = 8, dpi = 500)  

#### comapre explained heritability to total broad for metabolites with QTLs
herit_explained_df <- met_herit[!is.na(met_herit$heritExplained),]
herit_explained_df <- herit_explained_df[order(herit_explained_df$heritability),]
herit_explained_df_unpacked <- NULL
for(metN in 1:length(herit_explained_df[,1])){
  herit_row <- herit_explained_df[metN,]
  QTLvar <- sort(as.numeric(strsplit(herit_row$heritExplained, split = "_")[[1]]), decreasing = T)  
  
  herit_explained_df_unpacked <- rbind(herit_explained_df_unpacked, data.frame(x = metN, compound = herit_row$compound, QTLnumber = c(c(1:length(QTLvar)), "Residual Heritability"),  varExplained = c(QTLvar, herit_row$heritability - sum(QTLvar))))
  }

herit_explained_df_unpacked$QTLnumber <- factor(herit_explained_df_unpacked$QTLnumber, levels = sort(levels(herit_explained_df_unpacked$QTLnumber)))


barCols <- colorRampPalette(c("darkslategray", "gold"))
assocCols <- barCols(max(as.numeric(met_herit$nassoc)-1))
                                        
barplot_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_blank(), legend.position = "none", 
  panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.line = element_blank()) 
                                      
herit_explained_plot <- ggplot(herit_explained_df_unpacked, aes(x = factor(x), y = varExplained, fill = as.factor(QTLnumber)))
herit_explained_plot <- herit_explained_plot + geom_bar(binwidth = 1) + barplot_theme + geom_vline(aes(xintercept = 0), size = 0.5) + geom_hline(aes(yintercept = 0), size = 0.5) + 
  scale_x_discrete(name = "Associated Metabolites", expand = c(0,0)) + scale_y_continuous(name = "Heritability", expand = c(0,0), limits = c(0,1)) +
  scale_fill_manual(values = c("Residual Heritability" = "gray80", "1" = assocCols[1], "2" = assocCols[2], "3" = assocCols[3], "4" = assocCols[4]))
ggsave(plot = herit_explained_plot, file = "heritExplained.pdf", width = 14, height = 7) 


#pick out interesting examples and save boxplots
ggsave(plot = box_plotList[["Unk_111@20--287.003-0"]], file = "Unknown.pdf", width = 10, height = 8) 
ggsave(plot = box_plotList[["valine-0"]], file = "Valine.pdf", width = 15, height = 8) 
ggsave(plot = box_plotList[["aspartate-0"]], file = "Aspartate.pdf", width = 10, height = 8) 





herit_plot <- ggplot(met_herit, aes(x = heritability, fill = nassoc)) 
herit_plot + geom_dotplot(binwidth = 0.05, method = "histodot", binpositions = "all") + histodot_theme + geom_hline(aes(yintercept = 0), size = 0.5) + geom_vline(aes(xintercept = 0), size = 0.5) +
  scale_x_continuous(name = "Heritability", expand = c(0,0)) + scale_y_discrete(name = "Number of Metabolites", expand = c(0,0)) +
  scale_fill_manual(name = 'Number of QTLs', values = c("0" = assocCols[1], "1" = assocCols[2], "2" = assocCols[3], "3" = assocCols[4], "4" = assocCols[5]))

#library(colorRamps)
#hotCols <- colorRampPalette(c("black", "firebrick1"))

>
herit_plot <- ggplot(met_herit, aes(x = heritability, fill = fraction_missing_lines)) + facet_grid(. ~ nassoc)
herit_plot + geom_histogram(binwidth = 0.05) + coord_flip() + scale_fill_discrete(name = 'Fraction of \nMissing Lines', h = c(200, 360)) + scale_x_continuous(name = "Heritability", limits = c(0,1), expand = c(0,0)) + scale_y_continuous(name = "Number of Metabolites") + 
  theme(axis.text.x = element_text(size = 12, face = "bold"), axis.text.y = element_text(size = 12, face = "bold"), strip.text.x = element_text(size = 16, colour = "RED", face = "bold")) +
  theme(axis.title.x = element_text(size = 12, face = "bold"), axis.title.y = element_text(size = 12, face = "bold"), legend.position = "top", legend.title = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size = 10))
 
#color by zhu status
herit_plot <- ggplot(met_herit, aes(x = heritability, fill = Zhu)) + facet_grid(. ~ nassoc)
herit_plot + geom_histogram(binwidth = 0.05) + coord_flip() + scale_x_continuous(name = "Heritability", limits = c(0,1), expand = c(0,0)) + scale_y_continuous(name = "Number of Metabolites") + 
  theme(axis.text.x = element_text(size = 12, face = "bold"), axis.text.y = element_text(size = 12, face = "bold"), strip.text.x = element_text(size = 16, colour = "RED", face = "bold")) +
  theme(axis.title.x = element_text(size = 12, face = "bold"), axis.title.y = element_text(size = 12, face = "bold"), legend.position = "top", legend.title = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size = 10)) + scale_fill_discrete(name = 'Comparison with \nZhu et. al. 2012')




prob_matchMat <- met_herit[met_herit$Zhu != "Not-Measured",]
prob_matchMat$Zhu <- ifelse(prob_matchMat$Zhu == "mQTL Found", 1, 0)
#positive relationship between the probability of a Zhu QTL and the heritability of a trait
summary(glm(prob_matchMat$Zhu ~ prob_matchMat$heritability, family = binomial))

table(prob_matchMat$Zhu, (prob_matchMat$nassoc != 0)*1)



herit_plot <- ggplot(met_herit, aes(x = 1, y = heritability, fill = fraction_missing_lines)) + facet_wrap( ~ nassoc, ncol = 2)
herit_plot + geom_dotplot(binaxis = "y", stackgroups = TRUE, binwidth = 0.05, method = "histodot")

summary(glm(as.numeric(met_herit$nassoc)-1 ~ met_herit$heritability, family = poisson))
glm((as.numeric(met_herit$nassoc)-1)/max(as.numeric(met_herit$nassoc)-1) ~ met_herit$heritability, family = binomial)

disc_plot <- ggplot(met_herit, aes(x = fraction_missing_lines, fill = nassoc))
disc_plot + geom_bar()

disc_plot <- ggplot(met_herit, aes(x = fraction_missing_lines)) + facet_grid(nassoc ~ .)
disc_plot + geom_bar()



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





	