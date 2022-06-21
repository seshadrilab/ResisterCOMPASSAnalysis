#!/usr/bin/env Rscript
library(here)
library(data.table)
library(grid)
library(ggplot2)
library(ggcyto)
library(plyr)
library(flowWorkspace)
library(extrafont)
library(cowplot)
library(svglite)
library(coin)
library(COMPASS)
library(grDevices)
library(stringr)
library(tidyr)
library(ggsignif)

projectDir <- here::here()
source(file.path(projectDir, "scripts/flowHelperFunctions.R"))

date <- 20180731
set.seed(date) # Not really necessary

# Additional variables
cytokine_annotation_colors <- c("black", "black", "black", "black", "black", "black", "black")
ResisterStatusColors <- c("RSTR" = "gray78", "LTBI" = "gray27")
row_ann_colors <- list(`Status` = ResisterStatusColors)
grouping <- "Status"
outdir <- file.path(projectDir, "out/PostCompassPlots/TBAgs")
gsListPath <- file.path(projectDir, "out/GatingSets/AllBatchesForCompass_TBAgs")
gs <- load_gs(gsListPath)

cd4CompassResultDir <- file.path(projectDir, "out/CompassOutput/TBAgs/CD4")
CD4PP1CompassResult <- readRDS(file.path(cd4CompassResultDir, "4+_Peptide Pool 1/COMPASSResult_4+_Peptide Pool 1.rds"))
CD4PP2CompassResult <- readRDS(file.path(cd4CompassResultDir, "4+_Peptide Pool 2/COMPASSResult_4+_Peptide Pool 2.rds"))
CD4TBLysCompassResult <- readRDS(file.path(cd4CompassResultDir, "4+_TB Lysate/COMPASSResult_4+_TB Lysate.rds"))
CD4SEBCompassResult <- readRDS(file.path(cd4CompassResultDir, "4+_SEB/COMPASSResult_4+_SEB.rds"))

#######################################

# Figure 3
# Peptide Pool 1 CD4 responses
# Generate Heatmap for Panel A
p <- print(plot(CD4PP1CompassResult,
                grouping, show_rownames = FALSE,
                fontsize=16, fontsize_row=16, fontsize_col=12,
                row_annotation_colors=row_ann_colors))

svg(filename=file.path(outdir, "Figure3A_PP1_CD4_Heatmap.svg"), width=10, height=11)
print(grid.draw(p))
dev.off()

svglite(file=file.path(outdir, "Figure3A_PP1_CD4_Heatmap_svglite.svg"), width=10, height=11)
print(grid.draw(p))
dev.off()

# Reload GatingSet to avoid weird segfault error in next function
gs <- load_gs(gsListPath)

# Generate Polyfunctionality Score Boxplot for Panel B
cd4pp1PFSPlot <- fs.plot(gsOrGsListOrPath=gs,
                         compassResultOrPath=CD4PP1CompassResult,
                         stratifyBy=c("Status"),
                         ylims=NULL,
                         showTitle=FALSE,
                         removeGridAndBg=TRUE,
                         themeBaseSize = 18,
                         polyfunctionality = T,
                         axestitle_fontsize = 18,
                         axestick_fontsize = 14,
                         pvalue_fontsize = 5,
                         font="Arial",
                         geom_jitter_width=0.1,
                         xaxis_title = NULL)
cd4pp1PFSPlot$plot + scale_fill_manual(values = ResisterStatusColors)


svglite(file=file.path(outdir, "Figure3B_PP1_CD4_Polyfunctionality_svglite.svg"), width=5, height=6)
print(plot_grid(cd4pp1PFSPlot$plot + scale_fill_manual(values = ResisterStatusColors) + theme(legend.position="none"), labels=c("B")))
dev.off()

# Panel C Boxplots
CD4PP1_compass.subset.comparisons.result <- compass.subset.comparisons(compassResultOrPath=CD4PP1CompassResult,
                                                                       gsOrGsListOrPath=gs,
                                                                       parentSubset="4+",
                                                                       antigenCol="Antigen",
                                                                       stimAntigen="Peptide Pool 1",
                                                                       controlAntigen="DMSO",
                                                                       stratifyBy="Status",
                                                                       stratifyByValueMinuend = "LTBI",
                                                                       stratifyByValueSubtrahend = "RSTR",
                                                                       sampleIDCol="PATIENT ID")
CD4PP1_cr_cats <- orderHeatmapColumnsByCytokinePresence(CD4PP1CompassResult, "IFNg",
                                                        cats_df_cytokine_order_override = c("CD154", "IL2", "TNF", "CD107a", "IL4", "IL17a", "IFNg"))$fit$categories
CD4PP1_compass_subset_boxplots <- boxplotsCompassSubsetBgCorrPropStratified(CD4PP1_compass.subset.comparisons.result,
                                                                            CD4PP1_cr_cats,
                                                                            stratifyBy="Status",
                                                                            sampleIDCol="PATIENT ID",
                                                                            parentSubset="4+",
                                                                            removeGridAndBg=T,
                                                                            ylim=c(0, 0.007),
                                                                            stratifyBy_colors=ResisterStatusColors,
                                                                            axestitle_fontsize = 18,
                                                                            axestick_fontsize = 12,
                                                                            pvalue_fontsize = 5,
                                                                            font="Arial",
                                                                            legend_fontsize = 16,
                                                                            legend_position = c(0.22,1.05))
CD4PP1_compass_subset_boxplots

svglite(file=file.path(outdir, "Figure3C_PP1_CD4_COMPASS_Subset_boxplots_svglite.svg"), width=7, height=4)
print(plot_grid(CD4PP1_compass_subset_boxplots, labels="C"))
dev.off()

# Reload GatingSet to avoid weird segfault error in next function
gs <- load_gs(gsListPath)

# Flowplots D
# CD154+ IFNg+ IL2+ TNF+ (Four-function IFNg+ subset)
PP1CD4_flowPlot_IFNgPos <- highlight.boolean.subset.flow.plot(gs=gs,
                                                              individualsCol="PATIENT ID",
                                                              individual=c("RS102376", "RS102293"),
                                                              conditioncol="Antigen",
                                                              conditioncol2="Status",
                                                              exp="Peptide Pool 1",
                                                              ctrl="DMSO",
                                                              parentsubset="4+",
                                                              boolsubset="!4+/CD107a+&4+/CD154+&4+/IFNg+&!4+/IL17a+&4+/IL2+&!4+/IL4+&4+/TNF+",
                                                              xaxis="CD154",
                                                              yaxis="IFNg",
                                                              facetorder=c("Peptide Pool 1", "DMSO"),
                                                              facetorder2=c("RSTR", "LTBI"),
                                                              geomTextX = 775,
                                                              geomTextY = 3050,
                                                              xlims=c(0,4000),
                                                              ylims=c(0,3250),
                                                              stripLegendGridAxesTitle=TRUE,
                                                              font="Arial",
                                                              axestitle_fontsize = 21,
                                                              percentage_fontsize = 7) +
  labs(title="IFNg+ CD154+ IL2+ TNF+") + theme(plot.title = element_text(size = 21, hjust = 0.5),
                                               strip.text = element_text(size = 18, margin = margin(0,0,0,0, "cm")),
                                               line = element_blank(),
                                               panel.border = element_rect(fill = NA, color = "black"),
                                               panel.spacing = unit(0,"line"))
PP1CD4_flowPlot_IFNgPos

# Reload GatingSet to avoid weird segfault error in next function
gs <- load_gs(gsListPath)

# CD154+, IL2+, TNF+ (i.e. IFNg-)
PP1CD4_flowPlot_IFNgNeg <- highlight.boolean.subset.flow.plot(gs=gs,
                                                              individualsCol="PATIENT ID",
                                                              individual=c("RS102376", "RS102293"),
                                                              conditioncol="Antigen",
                                                              conditioncol2="Status",
                                                              exp="Peptide Pool 1",
                                                              ctrl="DMSO",
                                                              parentsubset="4+",
                                                              boolsubset="!4+/CD107a+&4+/CD154+&!4+/IFNg+&!4+/IL17a+&4+/IL2+&!4+/IL4+&4+/TNF+",
                                                              xaxis="CD154",
                                                              yaxis="IFNg",
                                                              facetorder=c("Peptide Pool 1", "DMSO"),
                                                              facetorder2=c("RSTR", "LTBI"),
                                                              geomTextX = 775,
                                                              geomTextY = 3050,
                                                              xlims=c(0,4000),
                                                              ylims=c(0,3250),
                                                              stripLegendGridAxesTitle=TRUE,
                                                              font="Arial",
                                                              axestitle_fontsize = 21,
                                                              percentage_fontsize = 7) +
  labs(title="CD154+ IL2+ TNF+") + theme(plot.title = element_text(size = 21, hjust = 0.5),
                                         strip.text = element_text(size = 18, margin = margin(0,0,0,0, "cm")),
                                         line = element_blank(),
                                         panel.border = element_rect(fill = NA, color = "black"),
                                         panel.spacing = unit(0,"line"))
PP1CD4_flowPlot_IFNgNeg

svglite(file=file.path(outdir, "Figure3D_PP1_CD4_FlowPlots_svglite.svg"), width=11, height=5)
print(plot_grid(plotlist=list(as.ggplot(PP1CD4_flowPlot_IFNgPos), as.ggplot(PP1CD4_flowPlot_IFNgNeg)),
          nrow = 1, labels=c("", ""), rel_widths = c(1,1)))
dev.off()

###############################################################

# Figure 4
# Peptide Pool 2 CD4 responses

# Heatmap for Panel A
p <- print(plot(CD4PP2CompassResult,
                grouping, show_rownames = FALSE,
                fontsize=16, fontsize_row=16, fontsize_col=12,
                row_annotation_colors=row_ann_colors))

svg(filename=file.path(outdir, "Figure4A_PP2_CD4_Heatmap.svg"), width=10, height=11)
print(grid.draw(p))
dev.off()

svglite(file=file.path(outdir, "Figure4A_PP2_CD4_Heatmap_svglite.svg"), width=10, height=11)
print(grid.draw(p))
dev.off()

# Reload GatingSet to avoid weird segfault error in next function
gs <- load_gs(gsListPath)

# Generate Polyfunctionality Score Boxplot for Panel B
cd4PP2PFSPlot <- fs.plot(gsOrGsListOrPath=gs,
                         compassResultOrPath=CD4PP2CompassResult,
                         stratifyBy=c("Status"),
                         ylims=NULL,
                         showTitle=FALSE,
                         removeGridAndBg=TRUE,
                         themeBaseSize = 18,
                         polyfunctionality = T,
                         axestitle_fontsize = 18,
                         axestick_fontsize = 14,
                         pvalue_fontsize = 5,
                         font="Arial",
                         geom_jitter_width=0.1,
                         xaxis_title = NULL)
cd4PP2PFSPlot$plot

svglite(file=file.path(outdir, "Figure4B_PP2_CD4_Polyfunctionality_svglite.svg"), width=5, height=6)
print(plot_grid(cd4PP2PFSPlot$plot + scale_fill_manual(values = ResisterStatusColors) + theme(legend.position="none"), labels=c("B")))
dev.off()

# Reload GatingSet to avoid weird segfault error in next function
gs <- load_gs(gsListPath)

# Flowplots C
# CD154+ IFNg+ IL2+ TNF+ (Four-function IFNg+ subset)
PP2CD4_flowPlot_IFNgPos <- highlight.boolean.subset.flow.plot(gs=gs,
                                                              individualsCol="PATIENT ID",
                                                              individual=c("RS102376", "RS102293"),
                                                              conditioncol="Antigen",
                                                              conditioncol2="Status",
                                                              exp="Peptide Pool 2",
                                                              ctrl="DMSO",
                                                              parentsubset="4+",
                                                              boolsubset="!4+/CD107a+&4+/CD154+&4+/IFNg+&!4+/IL17a+&4+/IL2+&!4+/IL4+&4+/TNF+",
                                                              xaxis="CD154",
                                                              yaxis="IFNg",
                                                              facetorder=c("Peptide Pool 2", "DMSO"),
                                                              facetorder2=c("RSTR", "LTBI"),
                                                              geomTextX = 1100,
                                                              geomTextY = 2850,
                                                              xlims=c(0,4000),
                                                              ylims=c(0,3250),
                                                              stripLegendGridAxesTitle=TRUE,
                                                              font="Arial",
                                                              axestitle_fontsize = 21,
                                                              percentage_fontsize = 7) +
  labs(title="IFNg+ CD154+ IL2+ TNF+") + theme(plot.title = element_text(size = 21, hjust = 0.5),
                                               strip.text.x = element_text(size = 18, margin = margin(0,0,0,0, "cm")),
                                               strip.text.y = element_text(size = 16, margin = margin(0,0,0,0, "cm")),
                                               line = element_blank(),
                                               panel.border = element_rect(fill = NA, color = "black"),
                                               panel.spacing = unit(0,"line"))
PP2CD4_flowPlot_IFNgPos

svglite(file=file.path(outdir, "Figure4C_PP2_CD4_FlowPlots_svglite.svg"), width=4, height=4.5)
print(plot_grid(as.ggplot(PP2CD4_flowPlot_IFNgPos),
          ncol = 1, labels=c("")))
dev.off()

#####################################################

# Figure 4
# TB Lysate CD4 responses

# Heatmap for Panel D
p <- print(plot(CD4TBLysCompassResult,
                grouping, show_rownames = FALSE,
                fontsize=16, fontsize_row=16, fontsize_col=12,
                row_annotation_colors=row_ann_colors))
svg(filename=file.path(outdir, "Figure4D_TBLys_CD4_Heatmap.svg"), width=10, height=11)
print(grid.draw(p))
dev.off()

svglite(file=file.path(outdir, "Figure4D_TBLys_CD4_Heatmap_svglite.svg"), width=10, height=11)
print(grid.draw(p))
dev.off()

# Reload GatingSet to avoid weird segfault error in next function
gs <- load_gs(gsListPath)

# Generate Polyfunctionality Score Boxplot for Panel E
cd4TBLysPFSPlot <- fs.plot(gsOrGsListOrPath=gs,
                           compassResultOrPath=CD4TBLysCompassResult,
                           stratifyBy=c("Status"),
                           ylims=NULL,
                           showTitle=FALSE,
                           removeGridAndBg=TRUE,
                           themeBaseSize = 18,
                           polyfunctionality = T,
                           axestitle_fontsize = 18,
                           axestick_fontsize = 14,
                           pvalue_fontsize = 5,
                           font="Arial",
                           geom_jitter_width=0.1,
                           xaxis_title = NULL)
cd4TBLysPFSPlot$plot

svglite(file=file.path(outdir, "Figure4E_PP2_TBLys_Polyfunctionality_svglite.svg"), width=5, height=6)
print(plot_grid(cd4TBLysPFSPlot$plot + scale_fill_manual(values = ResisterStatusColors) + theme(legend.position="none"), labels=c("E")))
dev.off()

# Reload GatingSet to avoid weird segfault error in next function
gs <- load_gs(gsListPath)

# Flowplot F
# CD154+ IFNg+ IL2+ TNF+ (Four-function IFNg+ subset)
TBLysCD4_flowPlot_IFNgPos <- highlight.boolean.subset.flow.plot(gs=gs,
                                                                individualsCol="PATIENT ID",
                                                                individual=c("RS102376", "RS102293"),
                                                                conditioncol="Antigen",
                                                                conditioncol2="Status",
                                                                exp="TB Lysate",
                                                                ctrl="DMSO",
                                                                parentsubset="4+",
                                                                boolsubset="!4+/CD107a+&4+/CD154+&4+/IFNg+&!4+/IL17a+&4+/IL2+&!4+/IL4+&4+/TNF+",
                                                                xaxis="CD154",
                                                                yaxis="IFNg",
                                                                facetorder=c("TB Lysate", "DMSO"),
                                                                facetorder2=c("RSTR", "LTBI"),
                                                                geomTextX = 1100,
                                                                geomTextY = 2850,
                                                                xlims=c(0,4000),
                                                                ylims=c(0,3250),
                                                                stripLegendGridAxesTitle=TRUE,
                                                                font="Arial",
                                                                axestitle_fontsize = 21,
                                                                percentage_fontsize = 7) +
  labs(title="IFNg+ CD154+ IL2+ TNF+")+ theme(plot.title = element_text(size = 21, hjust = 0.5),
                                              strip.text.x = element_text(size = 18, margin = margin(0,0,0,0, "cm")),
                                              strip.text.y = element_text(size = 18, margin = margin(0,0,0,0, "cm")),
                                              line = element_blank(),
                                              panel.border = element_rect(fill = NA, color = "black"),
                                              panel.spacing = unit(0,"line"))
TBLysCD4_flowPlot_IFNgPos

svglite(file=file.path(outdir, "Figure4F_TBLys_CD4_FlowPlots_svglite.svg"), width=4, height=4.5)
print(plot_grid(as.ggplot(TBLysCD4_flowPlot_IFNgPos),
          ncol = 1, labels=c("")))
dev.off()

# Panel G and H
CD4TBLys_compass.subset.comparisons.result <- compass.subset.comparisons(compassResultOrPath=CD4TBLysCompassResult,
                                                                         gsOrGsListOrPath=gs,
                                                                         parentSubset="4+",
                                                                         antigenCol="Antigen",
                                                                         stimAntigen="TB Lysate",
                                                                         controlAntigen="DMSO",
                                                                         stratifyBy="Status",
                                                                         stratifyByValueMinuend = "LTBI",
                                                                         stratifyByValueSubtrahend = "RSTR",
                                                                         sampleIDCol="PATIENT ID")
CD4TBLys_cr_cats <- orderHeatmapColumnsByCytokinePresence(CD4TBLysCompassResult, "IFNg",
                                                          cats_df_cytokine_order_override = c("CD154", "IL2", "TNF", "CD107a", "IL4", "IL17a", "IFNg"))$fit$categories
# Boxplots for IFNg+ subsets
CD4TBLys_compass.subset.comparisons.result_IFNgPos <- CD4TBLys_compass.subset.comparisons.result
# The subsets that get plotted depends on names(CD4TBLys_compass.subset.comparisons.result$wilcox), so subset those...
CD4TBLys_compass.subset.comparisons.result_IFNgPos$wilcox <- CD4TBLys_compass.subset.comparisons.result_IFNgPos$wilcox[grep("!IFNg", names(CD4TBLys_compass.subset.comparisons.result_IFNgPos$wilcox), invert=T)]
CD4TBLys_compass_subset_boxplots_IFNgPos <- boxplotsCompassSubsetBgCorrPropStratified(CD4TBLys_compass.subset.comparisons.result_IFNgPos,
                                                                                      CD4TBLys_cr_cats,
                                                                                      stratifyBy="Status",
                                                                                      sampleIDCol="PATIENT ID",
                                                                                      parentSubset="4+",
                                                                                      removeGridAndBg=T,
                                                                                      stratifyBy_colors=ResisterStatusColors,
                                                                                      axestitle_fontsize = 18,
                                                                                      axestick_fontsize = 12,
                                                                                      pvalue_fontsize = 5,
                                                                                      font="Arial",
                                                                                      legend_fontsize = 16,
                                                                                      ylim=c(0, 0.031),
                                                                                      legend_position = c(0.16,1))
CD4TBLys_compass_subset_boxplots_IFNgPos

svglite(file=file.path(outdir, "Figure4G_TBLys_CD4_IFNgPos_COMPASS_Subset_boxplots_svglite.svg"), width=9, height=5)
print(plot_grid(CD4TBLys_compass_subset_boxplots_IFNgPos, labels="G"))
dev.off()

# Boxplots for IFNg- subsets
CD4TBLys_compass.subset.comparisons.result_IFNgNeg <- CD4TBLys_compass.subset.comparisons.result
# The subsets that get plotted depends on names(CD4TBLys_compass.subset.comparisons.result$wilcox), so subset those...
CD4TBLys_compass.subset.comparisons.result_IFNgNeg$wilcox <- CD4TBLys_compass.subset.comparisons.result_IFNgNeg$wilcox[grep("!IFNg", names(CD4TBLys_compass.subset.comparisons.result_IFNgNeg$wilcox))]
CD4TBLys_compass_subset_boxplots_IFNgNeg <- boxplotsCompassSubsetBgCorrPropStratified(CD4TBLys_compass.subset.comparisons.result_IFNgNeg,
                                                                                      CD4TBLys_cr_cats,
                                                                                      stratifyBy="Status",
                                                                                      sampleIDCol="PATIENT ID",
                                                                                      parentSubset="4+",
                                                                                      removeGridAndBg=T,
                                                                                      stratifyBy_colors=ResisterStatusColors,
                                                                                      axestitle_fontsize = 18,
                                                                                      axestick_fontsize = 12,
                                                                                      pvalue_fontsize = 5,
                                                                                      font="Arial",
                                                                                      legend_fontsize = 16,
                                                                                      ylim=c(0, 0.031),
                                                                                      legend_position = c(0.14,1))
CD4TBLys_compass_subset_boxplots_IFNgNeg
svglite(file=file.path(outdir, "Figure4H_TBLys_CD4_IFNgNeg_COMPASS_Subset_boxplots_svglite.svg"), width=11, height=5)
print(plot_grid(CD4TBLys_compass_subset_boxplots_IFNgNeg, labels="H"))
dev.off()

#####################################################

# Supplemental Figure 1
# SEB CD4 Responses

# Heatmap for Panel D
p <- print(plot(CD4SEBCompassResult,
                grouping, show_rownames = FALSE,
                fontsize=16, fontsize_row=16, fontsize_col=12,
                row_annotation_colors=row_ann_colors))
print(grid.draw(p))

svg(filename=file.path(outdir, "FigureS1D_SEB_CD4_Heatmap.svg"), width=10, height=11)
print(grid.draw(p))
dev.off()

svglite(file=file.path(outdir, "FigureS1D_SEB_CD4_Heatmap_svglite.svg"), width=10, height=11)
print(grid.draw(p))
dev.off()

# Reload GatingSet to avoid weird segfault error in next function
gs <- load_gs(gsListPath)

# Generate Polyfunctionality Score Boxplot for Panel E
cd4sebPFSPlot <- fs.plot(gsOrGsListOrPath=gs,
                         compassResultOrPath=CD4SEBCompassResult,
                         stratifyBy=c("Status"),
                         ylims=NULL,
                         showTitle=FALSE,
                         removeGridAndBg=TRUE,
                         themeBaseSize = 18,
                         polyfunctionality = T,
                         axestitle_fontsize = 18,
                         axestick_fontsize = 14,
                         pvalue_fontsize = 5,
                         font="Arial",
                         geom_jitter_width=0.1,
                         xaxis_title = NULL)
cd4sebPFSPlot$plot

svglite(file=file.path(outdir, "FigureS1E_SEB_CD4_Polyfunctionality_svglite.svg"), width=7, height=6)
plot_grid(cd4sebPFSPlot$plot + scale_fill_manual(values = ResisterStatusColors) + theme(legend.position="none"), labels=c("E"))
dev.off()

# Reload GatingSet to avoid weird segfault error in next function
gs <- load_gs(gsListPath)

# Panel F
# CD154+ IFNg+ IL2+ TNF+ (Four-function IFNg+ subset)
SEBCD4_flowPlot_IFNgPos <- highlight.boolean.subset.flow.plot(gs=gs,
                                                              individualsCol="PATIENT ID",
                                                              individual=c("RS102376", "RS102293"),
                                                              conditioncol="Antigen",
                                                              conditioncol2="Status",
                                                              exp="SEB",
                                                              ctrl="DMSO",
                                                              parentsubset="4+",
                                                              boolsubset="!4+/CD107a+&4+/CD154+&4+/IFNg+&!4+/IL17a+&4+/IL2+&!4+/IL4+&4+/TNF+",
                                                              xaxis="CD154",
                                                              yaxis="IFNg",
                                                              facetorder=c("SEB", "DMSO"),
                                                              facetorder2=c("RSTR", "LTBI"),
                                                              geomTextX = 1000,
                                                              geomTextY = 3000,
                                                              xlims=c(0,4000),
                                                              ylims=c(0,3250),
                                                              stripLegendGridAxesTitle=TRUE,
                                                              font="Arial",
                                                              axestitle_fontsize = 21,
                                                              percentage_fontsize = 7,
                                                              overlayDotSize = 0.2,
                                                              geom_hex_bins = 120) +
  labs(title="IFNg+ CD154+ IL2+ TNF+") + theme(plot.title = element_text(size = 21, hjust = 0.5),
                                               strip.text = element_text(size = 18, margin = margin(0,0,0,0, "cm")),
                                               line = element_blank(),
                                               panel.border = element_rect(fill = NA, color = "black"),
                                               panel.spacing = unit(0,"line"))
SEBCD4_flowPlot_IFNgPos

svglite(file=file.path(outdir, "FigureS1F_SEB_CD4_FlowPlots_IFNgPos_svglite.svg"), width=4, height=4.5)
plot_grid(as.ggplot(SEBCD4_flowPlot_IFNgPos), labels=c(""))
dev.off()

######################################################

# Alternative COMPASS heatmaps (these are like what are in the final paper)
# Requires installation of COMPASS via devtools::install_github("malisas/COMPASS", ref="plotColorsCytokineHeight")
#
# The parameters cytokine_annotation_colors and cytokine_height_multiplier are not yet part of the main COMPASS repository,
# so a regular COMPASS installation will not support these options

# Figure 3
# Peptide Pool 1 CD4 responses
# Generate Heatmap for Panel A
p <- print(plot(orderHeatmapColumnsByCytokinePresence(CD4PP1CompassResult, "IFNg",
                                                      cats_df_cytokine_order_override = c("CD154", "IL2", "TNF", "CD107a", "IL4", "IL17a", "IFNg")),
                grouping, show_rownames = FALSE,
                fontsize=16, fontsize_row=16, fontsize_col=19,
                row_annotation_colors=row_ann_colors,
                cytokine_annotation_colors=cytokine_annotation_colors,
                cytokine_height_multiplier=2.5, order_by_max_functionality = F,
                fixed_column_order = T))

svg(filename=file.path(outdir, "Figure3A_PP1_CD4_Heatmap_alt.svg"), width=10, height=11)
print(grid.draw(p))
dev.off()

svglite(file=file.path(outdir, "Figure3A_PP1_CD4_Heatmap_alt_svglite.svg"), width=10, height=11)
print(grid.draw(p))
dev.off()

# Figure 4
# Peptide Pool 2 CD4 responses
# Heatmap for Panel A
p <- print(plot(orderHeatmapColumnsByCytokinePresence(CD4PP2CompassResult, "IFNg",
                                                      cats_df_cytokine_order_override = c("CD154", "IL2", "TNF", "CD107a", "IL4", "IL17a", "IFNg")),
                grouping, show_rownames = FALSE,
                fontsize=16, fontsize_row=16, fontsize_col=19,
                row_annotation_colors=row_ann_colors,
                cytokine_annotation_colors=cytokine_annotation_colors,
                cytokine_height_multiplier=2.5, order_by_max_functionality = F,
                fixed_column_order = T))

svg(filename=file.path(outdir, "Figure4A_PP2_CD4_Heatmap_alt.svg"), width=10, height=11)
print(grid.draw(p))
dev.off()

svglite(file=file.path(outdir, "Figure4A_PP2_CD4_Heatmap_alt_svglite.svg"), width=10, height=11)
print(grid.draw(p))
dev.off()

# Figure 4
# TB Lysate CD4 responses
# Heatmap for Panel D
p <- print(plot(orderHeatmapColumnsByCytokinePresence(CD4TBLysCompassResult, "IFNg",
                                                      cats_df_cytokine_order_override = c("CD154", "IL2", "TNF", "CD107a", "IL4", "IL17a", "IFNg")),
                grouping, show_rownames = FALSE,
                fontsize=16, fontsize_row=16, fontsize_col=19,
                row_annotation_colors=row_ann_colors,
                cytokine_annotation_colors=cytokine_annotation_colors,
                cytokine_height_multiplier=2.5, order_by_max_functionality = F,
                fixed_column_order = T))
svg(filename=file.path(outdir, "Figure4D_TBLys_CD4_Heatmap_alt.svg"), width=10, height=11)
print(grid.draw(p))
dev.off()

svglite(file=file.path(outdir, "Figure4D_TBLys_CD4_Heatmap_alt_svglite.svg"), width=10, height=11)
print(grid.draw(p))
dev.off()

# Supplemental Figure 1
# SEB CD4 Responses
# Heatmap for Panel D
p <- print(plot(orderHeatmapColumnsByCytokinePresence(CD4SEBCompassResult, "IFNg",
                                                      cats_df_cytokine_order_override = c("CD154", "IL2", "TNF", "CD107a", "IL4", "IL17a", "IFNg")),
                grouping, show_rownames = FALSE,
                fontsize=16, fontsize_row=16, fontsize_col=19,
                row_annotation_colors=row_ann_colors,
                cytokine_annotation_colors=cytokine_annotation_colors,
                cytokine_height_multiplier=2.5, order_by_max_functionality = F,
                fixed_column_order = T))
print(grid.draw(p))

svg(filename=file.path(outdir, "FigureS1D_SEB_CD4_Heatmap_alt.svg"), width=10, height=11)
print(grid.draw(p))
dev.off()

svglite(file=file.path(outdir, "FigureS1D_SEB_CD4_Heatmap_alt_svglite.svg"), width=10, height=11)
print(grid.draw(p))
dev.off()
