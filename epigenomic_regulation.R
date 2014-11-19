
##This will calculate all relevant methylation features: avg M val, avg logFC, and number of hyper and hypo methylated probes
##Two user files are needed as Rdata files: one containing the differential methylation information, and another containing the transcripts with respective upregulation/downregulation information
##See default input files for examples

load("reannotate_cpg_withTSS.Rdata") #provided
load("reannotate_specialmarks.Rdata") #provided
load("cpg_locations.Rdata") #provided
load("user_input/lung_cancer_meth450.Rdata") #USER INPUT FILE 1
load("user_input/transcript_regulation.Rdata") #USER INPUT FILE 2

relevant_transcripts <- transcript_regulation[,1]

cpg_locations <- cpg_locations[rownames(cpg_locations) %in% rownames(cpgs_limma),]
cpgs_limma <- cpgs_limma[rownames(cpgs_limma) %in% rownames(cpg_locations),]
cpg_locations <- cpg_locations[order(rownames(cpg_locations)),]
cpgs_limma <- cpgs_limma[order(rownames(cpgs_limma)),]

cpgs_limma <- cbind(cpg_locations, cpgs_limma)

reannotate_introns <- reannotate_introns[reannotate_introns[,7] %in% relevant_transcripts,]
transcript_regulation <- transcript_regulation[transcript_regulation[,1] %in% reannotate_introns[,7],]

for(i in 3:17) {
	meth450[,i] <- as.character(meth450[,i])
}

meth450 <- meth450[meth450$transcript != "" | meth450$mRNA != "" | meth450$TSS1500 != "" | meth450$TSS200 != "",] #get rid of unmarked methylation sites

transcript_id_vec <- unique(as.vector(reannotate_introns[,"transcript_id"])) #all transcripts to search (dictates loop flow)
logFC <- vector(mode = "numeric", length = length(transcript_id_vec)) #empty vector to fill (sum of avg expression / number of hits)
logFC_unmatched <- as.vector(cpgs_limma[,"logFC"]) #vector containing rows of logFC data
avgExpr <- vector(mode = "numeric", length = length(transcript_id_vec)) #empty vector to fill (sum of avg expression / number of hits)
avgExpr_unmatched <- as.vector(cpgs_limma[,"AveExpr"]) #vector containing rows of avgExpr data
nHyper <- vector(mode = "numeric", length = length(transcript_id_vec)) #empty vector to fill 
nHypo <- vector(mode = "numeric", length = length(transcript_id_vec)) #empty vector to fill 



meth450_transcripts <- as.vector(meth450[,"transcript"])
meth450_mRNAs <- as.vector(meth450[,"mRNA"])
meth450_aggregate <- paste(meth450_transcripts, meth450_mRNAs, sep = "") 
meth450 <- cbind(meth450, meth450_aggregate)
meth450[,18] <- as.character(meth450[,18])

nrowUnique <- length(transcript_id_vec)

logFC_avgs_matrix <- matrix(,nrow = nrowUnique, ncol = 15)
avgMval_avgs_matrix <- matrix(,nrow = nrowUnique, ncol = 15)
nHyper_matrix <- matrix(,nrow = nrowUnique, ncol = 15)
nHypo_matrix <- matrix(,nrow = nrowUnique, ncol = 15)

print("start methylation feature extraction (ignore UTR regions)")
for (j in 5:18) {
	meth450_current <- as.vector(meth450[,j])
	for (i in 1:nrowUnique) {
		if (i %% 300 == 0) {
			progress <- (j - 5)*nrowUnique + i
			progress <- progress / (14*nrowUnique)
      progress <- progress*100
			print(paste(progress, "%", sep = ""))
			}
		rowHits_logFC <- grep(paste(transcript_id_vec[i],":",sep = ""), meth450_current) #find how many transcripts have a cpg site on their "first exon"
	
		if (length(rowHits_logFC) == 0) { #if no hits, no need to go further
			logFC[i] <- "?"
      avgExpr[i] <- "?"
      nHyper[i] <- "?" 
      nHypo[i] <- "?"
		}
		else {
			cpgNamesMeth450 <- rownames(meth450)[rowHits_logFC] #get the cpg ID's because the two data frames do not have cpg sites organized in the same way
			matchingRowsCpgs_limma_index <- which(rownames(cpgs_limma) %in% cpgNamesMeth450) #don't care about the order, we just want to find the rows with matching hits
				if(length(matchingRowsCpgs_limma_index) == 0) { #because not all cpg sites have matching statistical data, possible that the probes exist but have no relevant data.  this avoids dividing by 0.
					logFC[i] <- "?"
					avgExpr[i] <- "?"
          nHyper[i] <- "?"
          nHypo[i] <- "?"
					next
				}
			aveLogFCSum <- sum(logFC_unmatched[matchingRowsCpgs_limma_index])
			aveLogFC <- aveLogFCSum/length(matchingRowsCpgs_limma_index)
			logFC[i] <- aveLogFC
      nHyper[i] <- sum(logFC_unmatched[matchingRowsCpgs_limma_index] > 0, 0, na.rm = T)
      nHypo[i] <- sum(logFC_unmatched[matchingRowsCpgs_limma_index] < 0, 0, na.rm = T)
			aveavgExprSum <- sum(avgExpr_unmatched[matchingRowsCpgs_limma_index])
			aveavgExpr <- aveavgExprSum/length(matchingRowsCpgs_limma_index)
			avgExpr[i] <- aveavgExpr
		}
	}
	logFC_avgs_matrix[,j-3] <- logFC
  nHyper_matrix[,j-3] <- nHyper
  nHypo_matrix[,j-3] <- nHypo
	avgMval_avgs_matrix[,j-3] <- avgExpr
}
logFC_avgs_matrix[,1] <- transcript_id_vec
avgMval_avgs_matrix[,1] <- transcript_id_vec
nHyper_matrix[,1] <- transcript_id_vec
nHypo_matrix[,1] <- transcript_id_vec
colnames(logFC_avgs_matrix) <- c("transcript_id", "first_exon_avglogFC", "first_intron_avglogFC", "exon_avglogFC", "intron_avglogFC", "last_exon_avglogFC", "last_intron_avglogFC", "5'-UTR_avglogFC", "3'-UTR_avglogFC", "CDS_avglogFC", "single_exon_avglogFC", "single_intron_avglogFC", "TSS1500_avglogFC", "TSS200_avglogFC", "fullTranscript_avglogFC")
colnames(avgMval_avgs_matrix) <- c("transcript_id", "first_exon_avgMval", "first_intron_avgMval", "exon_avgMval", "intron_avgMval", "last_exon_avgMval", "last_intron_avgMval", "5'-UTR_avgMval", "3'-UTR_avgMval", "CDS_avgMval", "single_exon_avgMval", "single_intron_avgMval", "TSS1500_avgMval", "TSS200_avgMval", "fullTranscript_avgMval")
colnames(nHyper_matrix) <- c("transcript_id", "first_exon_numHyper", "first_intron_numHyper", "exon_numHyper", "intron_numHyper", "last_exon_numHyper", "last_intron_numHyper", "5'-UTR_numHyper", "3'-UTR_numHyper", "CDS_numHyper", "single_exon_numHyper", "single_intron_numHyper", "TSS1500_numHyper", "TSS200_numHyper", "fullTranscript_numHyper")
colnames(nHypo_matrix) <- c("transcript_id", "first_exon_numHypo", "first_intron_numHypo", "exon_numHypo", "intron_numHypo", "last_exon_numHypo", "last_intron_numHypo", "5'-UTR_numHypo", "3'-UTR_numHypo", "CDS_numHypo", "single_exon_numHypo", "single_intron_numHypo", "TSS1500_numHypo", "TSS200_numHypo", "fullTranscript_numHypo")

logFC_avgs_matrix <- logFC_avgs_matrix[,c(1:7, 10:15)]
avgMval_avgs_matrix <- avgMval_avgs_matrix[,c(1:7, 10:15)]
nHyper_matrix <- nHyper_matrix[,c(1:7, 10:15)]
nHypo_matrix <- nHypo_matrix[,c(1:7, 10:15)]
print("First round of methylation features done.  Starting UTR calculations now.")







load("reannotate_UTR_hotfix.Rdata") #PROVIDED
reannotate_UTR <- reannotate_UTR[reannotate_UTR[,7] %in% relevant_transcripts,]


UTR_meth_hotfix <- matrix(nrow = length(unique(reannotate_UTR$transcript_id)), ncol = 9)

colnames(UTR_meth_hotfix) <- c("transcript_id", "UTR5_avgMval", "UTR5_avglogFC", "UTR3_avgMval", "UTR3_avglogFC", "UTR5_numHyper", "UTR5_numHypo", "UTR3_numHyper", "UTR3_numHypo")
UTR_meth_hotfix[,1] <- unique(reannotate_UTR$transcript_id)
rownames(UTR_meth_hotfix) <- unique(reannotate_UTR$transcript_id)

cpgs_limma$CHR <- paste("chr", cpgs_limma$CHR, sep = "")

i <- 1


while (i < nrow(reannotate_UTR)) {
	if (i %% 300 == 0) {
		progress <- i / nrow(reannotate_UTR)
    progress <- progress*100
		print(progress)
	}
	rowhits <- NULL
	rowhits <- which(cpgs_limma$CHR == reannotate_UTR$chr[i] & cpgs_limma$MAPINFO > reannotate_UTR$start[i] & cpgs_limma$MAPINFO < reannotate_UTR$end[i])
	while (reannotate_UTR$transcript_id[i] == reannotate_UTR$transcript_id[i+1] && reannotate_UTR$element[i] == reannotate_UTR$element[i+1]) {
		i <- i+1
		rowhits <- c(rowhits, which(cpgs_limma$CHR == reannotate_UTR$chr[i] & cpgs_limma$MAPINFO > reannotate_UTR$start[i] & cpgs_limma$MAPINFO < reannotate_UTR$end[i]))
	}
	if (length(rowhits) != 0) {
		avglogFC <- mean(cpgs_limma$logFC[rowhits], na.rm = T)
		avgexpr <- mean(cpgs_limma$AveExpr[rowhits], na.rm = T)
		nhyper <- sum(cpgs_limma$logFC[rowhits] > 0, na.rm = T)
		nhypo <- sum(cpgs_limma$logFC[rowhits] < 0, na.rm = T)

		if (reannotate_UTR$element[i] == "5'-UTR") {
			UTR_meth_hotfix[reannotate_UTR$transcript_id[i], 2] <- avgexpr
			UTR_meth_hotfix[reannotate_UTR$transcript_id[i], 3] <- avglogFC
			UTR_meth_hotfix[reannotate_UTR$transcript_id[i], 6] <- nhyper
			UTR_meth_hotfix[reannotate_UTR$transcript_id[i], 7] <- nhypo
		}

		else {
			UTR_meth_hotfix[reannotate_UTR$transcript_id[i], 4] <- avgexpr
			UTR_meth_hotfix[reannotate_UTR$transcript_id[i], 5] <- avglogFC
			UTR_meth_hotfix[reannotate_UTR$transcript_id[i], 8] <- nhyper
			UTR_meth_hotfix[reannotate_UTR$transcript_id[i], 9] <- nhypo
		}
	}

	i <- i + 1
}


UTR_meth_hotfix <- replace(UTR_meth_hotfix, is.na(UTR_meth_hotfix), "?")



load("precomputed_data_lungcancermodel.Rdata") #provided
precomputed_data <- precomputed_data[precomputed_data[,1] %in% transcript_regulation[,1],]
transcript_regulation <- transcript_regulation[transcript_regulation[,1] %in% precomputed_data[,1],]
precomputed_data <- precomputed_data[order(precomputed_data[,1]),]
transcript_regulation <- transcript_regulation[order(transcript_regulation[,1]),]
logFC_avgs_matrix <- logFC_avgs_matrix[order(logFC_avgs_matrix[,1]),]
avgMval_avgs_matrix <- avgMval_avgs_matrix[order(avgMval_avgs_matrix[,1]),]
nHyper_matrix <- nHyper_matrix[order(nHyper_matrix[,1]),]
nHypo_matrix <- nHypo_matrix[order(nHypo_matrix[,1]),]
UTR_meth_hotfix <- UTR_meth_hotfix[order(rownames(UTR_meth_hotfix)),]
finalModel <- cbind(transcript_regulation[,1], avgMval_avgs_matrix[,2:13], logFC_avgs_matrix[,2:13], UTR_meth_hotfix[,2:5], precomputed_data[,2:1369], nHyper_matrix[,2:13], nHypo_matrix[,2:13], UTR_meth_hotfix[,6:9], transcript_regulation[,2])
colnames(finalModel)[1] <- "transcript_id"
colnames(finalModel)[1426] <- "regulated"

write.csv(finalModel, file = "finalModel_forWEKA.csv", quote = F, row.names = F)

print("Model calculation done!  Saved to current directory as finalModel_forWEKA.csv")
