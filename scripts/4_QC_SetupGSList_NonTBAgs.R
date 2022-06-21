#!/usr/bin/env Rscript
library(CytoML)
library(here)
library(flowWorkspace)
library(tidyverse)
library(data.table)

projectDir <- here::here()
source(file.path(projectDir, "scripts/flowHelperFunctions.R")) # for boxplot.cell.counts and prepare.gating.set.list.4.compass

# Task: Read in FlowJo workspaces, perform QC, and prepare a GatingSetList for the COMPASS run.

flowJoXmlPath1 <- file.path(projectDir, "data/NonTBAgs/20180207_OMIP14_Batch1/20180209_NonTB_Antgns_OMIP14.xml")
flowJoXmlPath2 <- file.path(projectDir, "data/NonTBAgs/20180214_OMIP14_Batch2/20180216_NonTB_Antgns_OMIP14.xml")
batch1FcsDirPath <- file.path(projectDir, "data/NonTBAgs/20180207_OMIP14_Batch1/20180209_RSTR_NonTB_Atgns_Omip14_FCS")
batch2FcsDirPath <- file.path(projectDir, "data/NonTBAgs/20180214_OMIP14_Batch2/20180216_RSTR_NonTB_Atgns_Omip14_FCS")

file_map <- read.table(here::here("data/ImmPort_FCS_FileMapping.tsv"), sep = "\t", header = T,
                           colClasses = c("character", "character", "character", "character", "numeric"))

# Additional variables
qcOutDir <- file.path(projectDir, "out/QC/NonTBAgs")
patientStatusFilePath <- file.path(projectDir, "data/20170518_HiRisk_VisitA_Only_1.txt")
gatingSetListOutDir <- file.path(projectDir, "out/GatingSets/AllBatchesForCompass_NonTBAgs")

# First read in the flowJoXmlPaths with desired keywords
ws1 <- open_flowjo_xml(flowJoXmlPath1)
# getKeywords(ws1, "14521.fcs")
keywords2import=c("PATIENT ID", "Comp", "Antigen", "PLATE NAME", "TUBE NAME", "WELL ID")
sampleGroup <- "Samples"

cs1 <- cytoset()
quietly(pmap(file_map %>% dplyr::filter(grepl("^NonTBAgs", New_Folder_Path) & grepl("Batch1", New_Folder_Path) & !grepl("Compensation", Original_Name)) %>%
               dplyr::select(New_Folder_Path, Original_Name, flowJo_xml_sampleID),
             function(New_Folder_Path, Original_Name, flowJo_xml_sampleID) {
               cat(sprintf("Loading %s/%s for sampleID %s\n", New_Folder_Path, Original_Name, flowJo_xml_sampleID))
               cf_tmp <- load_cytoframe_from_fcs(here::here("data", New_Folder_Path, Original_Name))
               flowJo_xml_sample_name <- fj_ws_get_samples(ws1) %>% dplyr::filter(sampleID == flowJo_xml_sampleID) %>% dplyr::pull(name)
               # below, flowjo_to_gatingset will use the name stored in flowJo_xml_sample_name to match the cytoframe to the flowjo workspace entry
               cs_add_cytoframe(cs = cs1, sn = flowJo_xml_sample_name, cf = cf_tmp)  
             }))

gs1 <- flowjo_to_gatingset(ws1, name=sampleGroup, cytoset = cs1, keywords=keywords2import)

png(file.path(qcOutDir, "Batch1GatingTree.png"))
plot(gs1)
dev.off()

ws2 <- open_flowjo_xml(flowJoXmlPath2)
keywords2import=c("PATIENT ID", "Comp", "Antigen", "PLATE NAME", "TUBE NAME", "WELL ID")
sampleGroup <- "Samples"

cs2 <- cytoset()
quietly(pmap(file_map %>% dplyr::filter(grepl("^nonTBAgs", New_Folder_Path) & grepl("Batch2", New_Folder_Path) & !grepl("Compensation", Original_Name)) %>%
               dplyr::select(New_Folder_Path, Original_Name, flowJo_xml_sampleID),
             function(New_Folder_Path, Original_Name, flowJo_xml_sampleID) {
               cat(sprintf("Loading %s/%s for sampleID %s\n", New_Folder_Path, Original_Name, flowJo_xml_sampleID))
               cf_tmp <- load_cytoframe_from_fcs(here::here("data", New_Folder_Path, Original_Name))
               flowJo_xml_sample_name <- fj_ws_get_samples(ws2) %>% dplyr::filter(sampleID == flowJo_xml_sampleID) %>% dplyr::pull(name)
               # below, flowjo_to_gatingset will use the name stored in flowJo_xml_sample_name to match the cytoframe to the flowjo workspace entry
               cs_add_cytoframe(cs = cs2, sn = flowJo_xml_sample_name, cf = cf_tmp)  
             }))

gs2 <- flowjo_to_gatingset(ws2, name=sampleGroup, cytoset = cs2, keywords=keywords2import)

png(file.path(qcOutDir, "Batch2GatingTree.png"))
plot(gs2)
dev.off()

# The batches have the same Gating tree layout. However, remove the extra 8+/4+ gates from both batches
gs_pop_remove(gs1, "8+/4+")
gs_pop_remove(gs2, "8+/4+")

# Read in the patient status mappings and add to GatingSet metadata
patientStatusFile <- read.table(patientStatusFilePath, header = T, sep = "\t")[,c("RS_SUB_ACCESSION_NO", "STATUS_TST")]
colnames(patientStatusFile) <- c("PATIENT ID", "Status")
patientStatusFile$Status <- dplyr::recode(patientStatusFile$Status, "PERSISTENT NEG." = "RSTR", "TST+ CONTROL" = "LTBI")
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
    gs_pop_add(gs1, gs_pop_get_gate(gs1, nodeName), parent=parent)
    gs_pop_add(gs2, gs_pop_get_gate(gs2, nodeName), parent=parent)
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
cd4CountsBatch1 <- merge(gs_pop_get_count_fast(gs1)[which(gs_pop_get_count_fast(gs1)$Population == "/S/LV/L/3+/4+"),], pData(gs1), by.x = "name", by.y = "tmpRownames")
cd4CountsBatch1 <- cd4CountsBatch1[order(cd4CountsBatch1$Count),]
cd4CountsBatch2 <- merge(gs_pop_get_count_fast(gs2)[which(gs_pop_get_count_fast(gs2)$Population == "/S/LV/L/3+/4+"),], pData(gs2), by.x = "name", by.y = "tmpRownames")
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
