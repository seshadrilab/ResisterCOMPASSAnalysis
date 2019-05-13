#!/usr/bin/env Rscript
library(here)
library(flowWorkspace)
library(plyr)
library(data.table)

projectDir <- here::here()
source(file.path(projectDir, "scripts/flowHelperFunctions.R")) # for boxplot.cell.counts and prepare.gating.set.list.4.compass

# Task: Read in FlowJo workspaces, perform QC, and prepare a GatingSetList for the COMPASS run.

flowJoXmlPath1 <- file.path(projectDir, "data/NonTBAgs/20180207_OMIP14_Batch1/20180209_NonTB_Antgns_OMIP14.xml")
flowJoXmlPath2 <- file.path(projectDir, "data/NonTBAgs/20180214_OMIP14_Batch2/20180216_NonTB_Antgns_OMIP14.xml")
batch1FcsDirPath <- file.path(projectDir, "data/NonTBAgs/20180207_OMIP14_Batch1/20180209_RSTR_NonTB_Atgns_Omip14_FCS")
batch2FcsDirPath <- file.path(projectDir, "data/NonTBAgs/20180214_OMIP14_Batch2/20180216_RSTR_NonTB_Atgns_Omip14_FCS")

# Additional variables
qcOutDir <- file.path(projectDir, "out/QC/NonTBAgs")
patientStatusFilePath <- file.path(projectDir, "data/20170518_HiRisk_VisitA_Only_1.txt")
gatingSetListOutDir <- file.path(projectDir, "out/GatingSets/AllBatchesForCompass_NonTBAgs")
flowJoWorkspace_sampleID_FileMapping <- read.table(file.path(projectDir, "data/flowJoWorkspace_sampleID_FileMapping.tsv"), sep = "\t",
                                                   header = T, colClasses = c("character", "numeric", "character"))

# First read in the flowJoXmlPaths with desired keywords
ws1 <- openWorkspace(flowJoXmlPath1)
# getKeywords(ws1, "14521.fcs")
keywords2import=c("PATIENT ID", "Comp", "Antigen", "PLATE NAME", "TUBE NAME", "WELL ID")
sampleGroup <- "Samples"
ws1_filemap <- subset(flowJoWorkspace_sampleID_FileMapping, Experiment == "NonTBAgs_B1")[, c("sampleID", "file")]
missingFiles_b1 <- c("Specimen_001_D2_D02_122.fcs", "Specimen_001_D4_D04_124.fcs", "Specimen_001_D6_D06_126.fcs", "Specimen_001_D8_D08_128.fcs", "Specimen_001_E2_E02_130.fcs", "Specimen_001_E4_E04_132.fcs", "Specimen_001_E6_E06_134.fcs", "Specimen_001_E8_E08_136.fcs", "Specimen_001_F2_F02_138.fcs", "Specimen_001_F4_F04_140.fcs", "Specimen_001_F6_F06_142.fcs", "Specimen_001_F8_F08_144.fcs")
ws1_filemap <- ws1_filemap[which(!(ws1_filemap$file %in% missingFiles_b1)),]
ws1_filemap$file <- file.path(batch1FcsDirPath, ws1_filemap$file)
gs1 <- parseWorkspace(ws1, name=sampleGroup, keywords=keywords2import,
                      path=ws1_filemap)
png(file.path(qcOutDir, "Batch1GatingTree.png"))
plot(gs1)
dev.off()

ws2 <- openWorkspace(flowJoXmlPath2)
keywords2import=c("PATIENT ID", "Comp", "Antigen", "PLATE NAME", "TUBE NAME", "WELL ID")
sampleGroup <- "Samples"
ws2_filemap <- subset(flowJoWorkspace_sampleID_FileMapping, Experiment == "NonTBAgs_B2")[, c("sampleID", "file")]
missingFiles_b2 <- c("Specimen_001_F2_F02_138.fcs", "Specimen_001_F4_F04_140.fcs", "Specimen_001_F6_F06_142.fcs", "Specimen_001_F8_F08_144.fcs")
ws2_filemap <- ws2_filemap[which(!(ws2_filemap$file %in% missingFiles_b2)),]
ws2_filemap$file <- file.path(batch2FcsDirPath, ws2_filemap$file)
gs2 <- parseWorkspace(ws2, name=sampleGroup, keywords=keywords2import,
                      path=ws2_filemap)
png(file.path(qcOutDir, "Batch2GatingTree.png"))
plot(gs2)
dev.off()

# The batches have the same Gating tree layout. However, remove the extra 8+/4+ gates from both batches
Rm("8+/4+", gs1)
Rm("8+/4+", gs2)

# Read in the patient status mappings and add to GatingSet metadata
patientStatusFile <- read.table(patientStatusFilePath, header = T, sep = "\t")[,c("RS_SUB_ACCESSION_NO", "STATUS_TST")]
colnames(patientStatusFile) <- c("PATIENT ID", "Status")
patientStatusFile$Status <- revalue(patientStatusFile$Status, c("PERSISTENT NEG." = "RSTR", "TST+ CONTROL" = "LTBI"))
pData(gs1)$Batch <- "Batch 1"
pData(gs2)$Batch <- "Batch 2"
pData(gs1)$tmpRownames <- rownames(pData(gs1))
pData(gs2)$tmpRownames <- rownames(pData(gs2))
pDataGsList1Tmp <- merge(pData(gs1), patientStatusFile, by = "PATIENT ID")
rownames(pDataGsList1Tmp) <- pDataGsList1Tmp$tmpRownames
pDataGsList2Tmp <- merge(pData(gs2), patientStatusFile, by = "PATIENT ID")
rownames(pDataGsList2Tmp) <- pDataGsList2Tmp$tmpRownames
pData(gs1) <- pDataGsList1Tmp
pData(gs2) <- pDataGsList2Tmp

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

# # Samples to investigate due to low CD3 or CD4 cell count:
# PATIENT IDs:
#  Batch 1:
#   None
#  Batch 2:
#   RS102372 (Batch 2 Plate 3 TST+ Control), has both low CD3 (3,311 DMSO events) and CD4 (1,044 DMSO events) counts
#      "15712.fcs_751608", "15718.fcs_695046", "15716.fcs_613734", "15714.fcs_781110"
# # Samples to remove due to observed contamination:
# PATIENT IDs:
#  Batch 1
#   RS102173, TST+ Control Batch 1 Plate 3 ("14555.fcs_643665", "14561.fcs_725703", "14559.fcs_743061", "14557.fcs_719466")
#   RS102175, Persistent Neg Batch 1 Plate 3 ("14565.fcs_433686", "14567.fcs_499620", "14571.fcs_519123", "14569.fcs_556314")
#   RS102177, Persistent Neg Batch 1 Plate 3 ("14575.fcs_552750", "14581.fcs_620928", "14577.fcs_647790", "14579.fcs_669372")
#  Batch 2:
#   None

# Remove samples deemed unfit:
samples2Drop <- c("15712.fcs_751608", "15718.fcs_695046", "15716.fcs_613734", "15714.fcs_781110",
                  "14555.fcs_643665", "14561.fcs_725703", "14559.fcs_743061", "14557.fcs_719466",
                  "14565.fcs_433686", "14567.fcs_499620", "14571.fcs_519123", "14569.fcs_556314",
                  "14575.fcs_552750", "14581.fcs_620928", "14577.fcs_647790", "14579.fcs_669372")
gs1 <- subset(gs1, !(tmpRownames %in% samples2Drop))
gs2 <- subset(gs2, !(tmpRownames %in% samples2Drop))

# Finally, save the GatingSetList
prepare.gating.set.list.4.compass(gsList=list(gs1, gs2),
                                  outDir=gatingSetListOutDir)
