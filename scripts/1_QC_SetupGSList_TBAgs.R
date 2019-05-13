#!/usr/bin/env Rscript
library(here)
library(flowWorkspace)
library(plyr)

projectDir <- here::here()
source(file.path(projectDir, "scripts/flowHelperFunctions.R")) # for boxplot.cell.counts and prepare.gating.set.list.4.compass

# Task: Read in FlowJo workspaces, perform QC, and prepare a GatingSetList for the COMPASS run.

flowJoXmlPath1 <- file.path(projectDir, "data/TBAgs/20170605_RSTR_OMIP14_ICS_Batch1/20170607_RSTR_ICS_Batch1.xml")
flowJoXmlPath2 <- file.path(projectDir, "data/TBAgs/20170612_RSTR_OMIP14_ICS_Batch2/20170614_RSTR_ICS_Batch2.xml")
batch1FcsDirPath <- file.path(projectDir, "data/TBAgs/20170605_RSTR_OMIP14_ICS_Batch1/20170607_RSTR_OMIP14_ICS")
batch2FcsDirPath <- file.path(projectDir, "data/TBAgs/20170612_RSTR_OMIP14_ICS_Batch2/20170614_RSTR_OMIP14_ICS_Batch2")

# Additional variables
qcOutDir <- file.path(projectDir, "out/QC/TBAgs")
patientStatusFilePath <- file.path(projectDir, "data/20170518_HiRisk_VisitA_Only_1.txt")
gatingSetListOutDir <- file.path(projectDir, "out/GatingSets/AllBatchesForCompass_TBAgs")
flowJoWorkspace_sampleID_FileMapping <- read.table(file.path(projectDir, "data/flowJoWorkspace_sampleID_FileMapping.tsv"), sep = "\t",
                                                   header = T, colClasses = c("character", "numeric", "character"))


# First read in the flowJoXmlPaths with desired keywords
ws1 <- openWorkspace(flowJoXmlPath1)
keywords2import=c("PATIENT ID", "Comp", "Antigen", "Status", "PLATE NAME", "TUBE NAME", "WELL ID")
sampleGroup <- "Samples"
ws1_filemap <- subset(flowJoWorkspace_sampleID_FileMapping, Experiment == "TBAgs_B1")[, c("sampleID", "file")]
missingFiles_b1 <- c("Specimen_001_B2_B02_012.fcs", "Specimen_001_B4_B04_014.fcs", "Specimen_001_B6_B06_016.fcs", "Specimen_001_B8_B08_018.fcs", "Specimen_001_B10_B10_020.fcs", "Specimen_001_C6_C06_066.fcs", "Specimen_001_C8_C08_068.fcs", "Specimen_001_A2_A02_112.fcs", "Specimen_001_A4_A04_114.fcs", "Specimen_001_A6_A06_116.fcs", "Specimen_001_A8_A08_118.fcs", "Specimen_001_A10_A10_120.fcs", "Specimen_001_A2_A02_192.fcs", "Specimen_001_A4_A04_194.fcs", "Specimen_001_A6_A06_196.fcs", "Specimen_001_A8_A08_198.fcs", "Specimen_001_A10_A10_200.fcs")
ws1_filemap <- ws1_filemap[which(!(ws1_filemap$file %in% missingFiles_b1)),]
ws1_filemap$file <- file.path(batch1FcsDirPath, ws1_filemap$file)
gs1 <- parseWorkspace(ws1, name=sampleGroup, keywords=keywords2import,
                      path=ws1_filemap) # batch1FcsDirPath

#########################

ws2 <- openWorkspace(flowJoXmlPath2)
keywords2import=c("PATIENT ID", "Comp", "Antigen", "Status", "PLATE NAME", "TUBE NAME", "WELL ID")
sampleGroup <- "Samples"
ws2_filemap <- subset(flowJoWorkspace_sampleID_FileMapping, Experiment == "TBAgs_B2")[, c("sampleID", "file")]
missingFiles_b2 <- c("Specimen_001_F10_F10_100.fcs")
ws2_filemap <- ws2_filemap[which(!(ws2_filemap$file %in% missingFiles_b2)),]
ws2_filemap$file <- file.path(batch2FcsDirPath, ws2_filemap$file)
gs2 <- parseWorkspace(ws2, name=sampleGroup, keywords=keywords2import,
                      path=ws2_filemap)
png(file.path(qcOutDir, "Batch2GatingTree.png"))
plot(gs2)
dev.off()

# The only inconsistency between the batches is the extra "8+/4+" gate in Batch 2. Remove it:
Rm("8+/4+", gs2)

# We want to run COMPASS on DNeg and DPos populations as well.
# Cytokine gates don't exist on DNeg and DPos, so
# copy over those gates from the 4+ node to DNeg and DPos
# Load the GatingSetList and copy gates from the 4+ subset over onto DN and DP
# getNodes(gs1[[1]], path="auto", order="tsort")
nodes2Copy <- c("4+/TNF+", "4+/IL17a+", "4+/IL4+", "4+/IL2+", "4+/IFNg+", "4+/CD154+", "4+/CD107a+")
parents <- c("DPos", "DNeg")
for (parent in parents) {
  for (nodeName in nodes2Copy) {
    add(gs1, getGate(gs1, nodeName), parent=parent)
    add(gs2, getGate(gs2, nodeName), parent=parent)
  }
  recompute(gs1, parent)
  recompute(gs2, parent)
}


# Add treatment and control columns to pData.
# Add column to pData labeling each row as Control or Treatment, based on the Antigen
getTrt <- function(x) {
  if (x == "DMSO") {
    "Control"
  } else {
    "Treatment"
  }
}
pData(gs1)$trt <- sapply(pData(gs1)$Antigen, getTrt)
pData(gs2)$trt <- sapply(pData(gs2)$Antigen, getTrt)

# Read in the patient status mappings and add to GatingSet metadata
patientStatusFile <- read.table(patientStatusFilePath, header = T, sep = "\t")[,c("RS_SUB_ACCESSION_NO", "STATUS_TST")]
colnames(patientStatusFile) <- c("PATIENT ID", "Status")
patientStatusFile$Status <- revalue(patientStatusFile$Status, c("PERSISTENT NEG." = "RSTR", "TST+ CONTROL" = "LTBI"))
pData(gs1)$Batch <- "Batch 1"
pData(gs2)$Batch <- "Batch 2"
pData(gs1) <- subset(pData(gs1), select = -c(Status))
pData(gs2) <- subset(pData(gs2), select = -c(Status))
pData(gs1)$tmpRownames <- rownames(pData(gs1))
pData(gs2)$tmpRownames <- rownames(pData(gs2))
pDataGsList1Tmp <- merge(pData(gs1), patientStatusFile, by = "PATIENT ID")
rownames(pDataGsList1Tmp) <- pDataGsList1Tmp$tmpRownames
pDataGsList2Tmp <- merge(pData(gs2), patientStatusFile, by = "PATIENT ID")
rownames(pDataGsList2Tmp) <- pDataGsList2Tmp$tmpRownames
pData(gs1) <- pDataGsList1Tmp
pData(gs2) <- pDataGsList2Tmp

# Now perform QC
# First on CD3 counts
boxplot.cell.counts(gatingSet=subset(gs1, `PLATE NAME` == "Plate 1"),
                    outdir=qcOutDir,
                    stratifyByLevel1="PATIENT ID",
                    subpopulation="/S/LV/L/3+",
                    batch="Batch 1 Plate 1",
                    threshold=10000)
boxplot.cell.counts(gatingSet=subset(gs1, `PLATE NAME` == "Plate 2"),
                    outdir=qcOutDir,
                    stratifyByLevel1="PATIENT ID",
                    subpopulation="/S/LV/L/3+",
                    batch="Batch 1 Plate 2",
                    threshold=10000)
boxplot.cell.counts(gatingSet=subset(gs1, `PLATE NAME` == "Plate 3"),
                    outdir=qcOutDir,
                    stratifyByLevel1="PATIENT ID",
                    subpopulation="/S/LV/L/3+",
                    batch="Batch 1 Plate 3",
                    threshold=10000)
boxplot.cell.counts(gatingSet=subset(gs1, `PLATE NAME` == "Plate 4"),
                    outdir=qcOutDir,
                    stratifyByLevel1="PATIENT ID",
                    subpopulation="/S/LV/L/3+",
                    batch="Batch 1 Plate 4",
                    threshold=10000)
boxplot.cell.counts(gatingSet=subset(gs2, `PLATE NAME` == "Plate 1"),
                    outdir=qcOutDir,
                    stratifyByLevel1="PATIENT ID",
                    subpopulation="/S/LV/L/3+",
                    batch="Batch 2 Plate 1",
                    threshold=10000)
boxplot.cell.counts(gatingSet=subset(gs2, `PLATE NAME` == "Plate 2"),
                    outdir=qcOutDir,
                    stratifyByLevel1="PATIENT ID",
                    subpopulation="/S/LV/L/3+",
                    batch="Batch 2 Plate 2",
                    threshold=10000)
boxplot.cell.counts(gatingSet=subset(gs2, `PLATE NAME` == "Plate 3"),
                    outdir=qcOutDir,
                    stratifyByLevel1="PATIENT ID",
                    subpopulation="/S/LV/L/3+",
                    batch="Batch 2 Plate 3",
                    threshold=10000)
boxplot.cell.counts(gatingSet=subset(gs2, `PLATE NAME` == "Plate 4"),
                    outdir=qcOutDir,
                    stratifyByLevel1="PATIENT ID",
                    subpopulation="/S/LV/L/3+",
                    batch="Batch 2 Plate 4",
                    threshold=10000)

# # Samples to remove due to low cell count:
# PATIENT IDs:
#  Batch 1:
#   RS102057
#   RS102042
#  Batch 2:
#   None
# These samples are removed by filename later in this file

# This time with CD4 counts:
boxplot.cell.counts(gatingSet=subset(gs1, `PLATE NAME` == "Plate 1"),
                    outdir=qcOutDir,
                    stratifyByLevel1="PATIENT ID",
                    subpopulation="/S/LV/L/3+/4+",
                    batch="Batch 1 Plate 1",
                    threshold=3000)
boxplot.cell.counts(gatingSet=subset(gs1, `PLATE NAME` == "Plate 2"),
                    outdir=qcOutDir,
                    stratifyByLevel1="PATIENT ID",
                    subpopulation="/S/LV/L/3+/4+",
                    batch="Batch 1 Plate 2",
                    threshold=3000)
boxplot.cell.counts(gatingSet=subset(gs1, `PLATE NAME` == "Plate 3"),
                    outdir=qcOutDir,
                    stratifyByLevel1="PATIENT ID",
                    subpopulation="/S/LV/L/3+/4+",
                    batch="Batch 1 Plate 3",
                    threshold=3000)
boxplot.cell.counts(gatingSet=subset(gs1, `PLATE NAME` == "Plate 4"),
                    outdir=qcOutDir,
                    stratifyByLevel1="PATIENT ID",
                    subpopulation="/S/LV/L/3+/4+",
                    batch="Batch 1 Plate 4",
                    threshold=3000)
boxplot.cell.counts(gatingSet=subset(gs2, `PLATE NAME` == "Plate 1"),
                    outdir=qcOutDir,
                    stratifyByLevel1="PATIENT ID",
                    subpopulation="/S/LV/L/3+/4+",
                    batch="Batch 2 Plate 1",
                    threshold=3000)
boxplot.cell.counts(gatingSet=subset(gs2, `PLATE NAME` == "Plate 2"),
                    outdir=qcOutDir,
                    stratifyByLevel1="PATIENT ID",
                    subpopulation="/S/LV/L/3+/4+",
                    batch="Batch 2 Plate 2",
                    threshold=3000)
boxplot.cell.counts(gatingSet=subset(gs2, `PLATE NAME` == "Plate 3"),
                    outdir=qcOutDir,
                    stratifyByLevel1="PATIENT ID",
                    subpopulation="/S/LV/L/3+/4+",
                    batch="Batch 2 Plate 3",
                    threshold=3000)
boxplot.cell.counts(gatingSet=subset(gs2, `PLATE NAME` == "Plate 4"),
                    outdir=qcOutDir,
                    stratifyByLevel1="PATIENT ID",
                    subpopulation="/S/LV/L/3+/4+",
                    batch="Batch 2 Plate 4",
                    threshold=3000)
cd4CountsBatch1 <- merge(getPopStats(gs1)[which(getPopStats(gs1)$Population == "4+"),], pData(gs1), by.x = "name", by.y = "tmpRownames")
cd4CountsBatch1 <- cd4CountsBatch1[order(cd4CountsBatch1$Count),]
cd4CountsBatch2 <- merge(getPopStats(gs2)[which(getPopStats(gs2)$Population == "4+"),], pData(gs2), by.x = "name", by.y = "tmpRownames")
cd4CountsBatch2 <- cd4CountsBatch2[order(cd4CountsBatch2$Count),]

# Rename columns, e.g. Count and ParentCount columns to CD4_Count and CD3_Count.
colnames(cd4CountsBatch1)[which(colnames(cd4CountsBatch1) == "Count")] <- "CD4_Count"
colnames(cd4CountsBatch2)[which(colnames(cd4CountsBatch2) == "Count")] <- "CD4_Count"
colnames(cd4CountsBatch1)[which(colnames(cd4CountsBatch1) == "ParentCount")] <- "CD3_Count"
colnames(cd4CountsBatch2)[which(colnames(cd4CountsBatch2) == "ParentCount")] <- "CD3_Count"
colnames(cd4CountsBatch1)[which(colnames(cd4CountsBatch1) == "name.y")] <- "name_short"
colnames(cd4CountsBatch2)[which(colnames(cd4CountsBatch2) == "name.y")] <- "name_short"

cd3_cd4_counts <- rbind(subset(cd4CountsBatch1, select = -c(Population, Parent)), subset(cd4CountsBatch2, select = -c(Population, Parent)))
write.csv(cd3_cd4_counts, file.path(qcOutDir, "CD3_CD4_Counts_Both_Batches.csv"), row.names = F)

# Remove samples deemed unfit:
# tmpRownames values
samples2Drop <- c("5920.fcs_17827", "5916.fcs_41828", "5922.fcs_69785", "5914.fcs_39744", "5918.fcs_34582",
                  "5694.fcs_135306", "5696.fcs_113469", "5688.fcs_119558", "5690.fcs_98464", "5692.fcs_102349",
                  "5785.fcs_4721", "5787.fcs_6723", "7182.fcs_268521", "6079.fcs_260438", "6081.fcs_299427",
                  "6075.fcs_253894", "6083.fcs_313780", "6077.fcs_265292")
gs1 <- subset(gs1, !(tmpRownames %in% samples2Drop))
gs2 <- subset(gs2, !(tmpRownames %in% samples2Drop))

# Finally, save the GatingSetList
prepare.gating.set.list.4.compass(gsList=list(gs1, gs2),
                                  outDir=gatingSetListOutDir)
