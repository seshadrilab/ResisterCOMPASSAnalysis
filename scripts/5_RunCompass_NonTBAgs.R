#!/usr/bin/env Rscript
library(here)
library(flowWorkspace)
library(COMPASS)
library(grid)

projectDir <- here::here()
source(file.path(projectDir, "scripts/flowHelperFunctions.R")) # for generic.compass.wrapper

# Task: Run Compass on the CD4, CD8, DNeg, and DPos nodes

# Additional variables
gsPath <- file.path(projectDir, "out/GatingSets/AllBatchesForCompass_NonTBAgs")
seed <- 20180221
mapMarkers <- list("IL2", "IL4", "IFNg", "TNF", "IL17a", "CD154", "CD107a")
individuals <- "PATIENT ID"
uniqueruns <- "Antigen"
grouping <- "Status"

# NOTE: You might want to run each of the four runs below in parallel to save time, as opposed to one after another

# Run COMPASS for CD4+
mapNodesCD4 <- c("4+/IL2+", "4+/IL4+", "4+/IFNg+", "4+/TNF+", "4+/IL17a+", "4+/CD154+", "4+/CD107a+")
outdir <- file.path(projectDir, "out/CompassOutput/NonTBAgs/CD4")
node <- "4+"
nodeMarkerMapCD4 <- mapMarkers
names(nodeMarkerMapCD4) <- mapNodesCD4
generic.compass.wrapper(path=gsPath,
                        seed=seed,
                        outdir=outdir,
                        cnode=node,
                        nodemarkermap=nodeMarkerMapCD4,
                        individuals=individuals,
                        grouping=grouping,
                        uniqueruns=uniqueruns,
                        countFilterThreshold=0)

# Run COMPASS for CD8+
mapNodesCD8 <- c("8+/IL2+", "8+/IL4+", "8+/IFNg+", "8+/TNF+", "8+/IL17a+", "8+/CD154+", "8+/CD107a+")
outdir <- file.path(projectDir, "out/CompassOutput/NonTBAgs/CD8")
node <- "8+"
nodeMarkerMapCD8 <- mapMarkers
names(nodeMarkerMapCD8) <- mapNodesCD8
generic.compass.wrapper(path=gsPath,
                        seed=seed,
                        outdir=outdir,
                        cnode=node,
                        nodemarkermap=nodeMarkerMapCD8,
                        individuals=individuals,
                        grouping=grouping,
                        uniqueruns=uniqueruns,
                        countFilterThreshold=0)

# Run COMPASS for DNeg
mapNodesDN <- c("DNeg/IL2+", "DNeg/IL4+", "DNeg/IFNg+", "DNeg/TNF+", "DNeg/IL17a+", "DNeg/CD154+", "DNeg/CD107a+")
outdir <- file.path(projectDir, "out/CompassOutput/NonTBAgs/DN")
node <- "DNeg"
nodeMarkerMapDN <- mapMarkers
names(nodeMarkerMapDN) <- mapNodesDN
generic.compass.wrapper(path=gsPath,
                        seed=seed,
                        outdir=outdir,
                        cnode=node,
                        nodemarkermap=nodeMarkerMapDN,
                        individuals=individuals,
                        grouping=grouping,
                        uniqueruns=uniqueruns,
                        countFilterThreshold=0)

# Run COMPASS for DPos
mapNodesDP <- c("DPos/IL2+", "DPos/IL4+", "DPos/IFNg+", "DPos/TNF+", "DPos/IL17a+", "DPos/CD154+", "DPos/CD107a+")
outdir <- file.path(projectDir, "out/CompassOutput/NonTBAgs/DP")
node <- "DPos"
nodeMarkerMapDP <- mapMarkers
names(nodeMarkerMapDP) <- mapNodesDP
generic.compass.wrapper(path=gsPath,
                        seed=seed,
                        outdir=outdir,
                        cnode=node,
                        nodemarkermap=nodeMarkerMapDP,
                        individuals=individuals,
                        grouping=grouping,
                        uniqueruns=uniqueruns,
                        countFilterThreshold=0)