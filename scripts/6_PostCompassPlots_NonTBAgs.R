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
outdir <- file.path(projectDir, "out/PostCompassPlots/NonTBAgs")
gsListPath <- file.path(projectDir, "out/GatingSets/AllBatchesForCompass_NonTBAgs")
gs <- load_gs(gsListPath)

cd8CompassResultDir <- file.path(projectDir, "out/CompassOutput/NonTBAgs/CD8")
CD8CEFCompassResult <- readRDS(file.path(cd8CompassResultDir, "8+_CEF/COMPASSResult_8+_CEF.rds"))

# CEF CD8 plots

# Heatmap Panel A
p <- print(plot(CD8CEFCompassResult,
                grouping, show_rownames = FALSE,
                fontsize=16, fontsize_row=16, fontsize_col=12,
                row_annotation_colors=row_ann_colors))

svg(filename=file.path(outdir, "FigureS1A_CEF_CD8_Heatmap.svg"), width=10, height=11)
grid.draw(p)
dev.off()

svglite(file=file.path(outdir, "FigureS1A_CEF_CD8_Heatmap_svglite.svg"), width=10, height=11)
grid.draw(p)
dev.off()

# Polyfunctionality Boxplot Panel B
cd8cefPFSPlot <- fs.plot(gsOrGsListOrPath=gs,
                         compassResultOrPath=CD8CEFCompassResult,
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
cd8cefPFSPlot$plot

svglite(file=file.path(outdir, "FigureS1B_CEF_CD8_Polyfunctionality_svglite.svg"), width=7, height=6)
plot_grid(cd8cefPFSPlot$plot + scale_fill_manual(values = ResisterStatusColors) + theme(legend.position="none"), labels=c("B"))
dev.off()

# Reload GatingSet to avoid weird segfault error in next function
gs <- load_gs(gsListPath)

# Flowplot IFNg+
# Panel C
# RS102393 is Persistent Neg, RS102415 is TST+ Control
CEFCD8_flowPlot_IFNgPos <- highlight.boolean.subset.flow.plot(gs=gs,
                                                              individualsCol="PATIENT ID",
                                                              individual=c("RS102393", "RS102415"),
                                                              conditioncol="Antigen",
                                                              conditioncol2="Status",
                                                              exp="CEF",
                                                              ctrl="DMSO",
                                                              parentsubset="8+",
                                                              boolsubset="8+/CD107a+&!8+/CD154+&8+/IFNg+&!8+/IL17a+&8+/IL2+&!8+/IL4+&8+/TNF+",
                                                              xaxis="TNF",
                                                              yaxis="IFNg",
                                                              facetorder=c("CEF", "DMSO"),
                                                              facetorder2=c("RSTR", "LTBI"),
                                                              geomTextX = 450,
                                                              geomTextY = 3175,
                                                              xlims=c(-500,3500),
                                                              ylims=c(-500,4000),
                                                              stripLegendGridAxesTitle=TRUE,
                                                              font="Arial",
                                                              axestitle_fontsize = 21,
                                                              percentage_fontsize = 7) +
  labs(title="IFNg+ CD107a+ IL2+ TNF+") + theme(plot.title = element_text(size = 21, hjust = 0.5),
                                                strip.text = element_text(size = 19, margin = margin(0,0,0,0, "cm")),
                                                line = element_blank(),
                                                panel.border = element_rect(fill = NA, color = "black"),
                                                panel.spacing = unit(0,"line"))
CEFCD8_flowPlot_IFNgPos

svglite(file=file.path(outdir, "FigureS1C_CEF_CD8_FlowPlots_IFNgPos_svglite.svg"), width=4, height=4)
plot_grid(as.ggplot(CEFCD8_flowPlot_IFNgPos), labels=c(""))
dev.off()


######################################################

# Alternative COMPASS heatmap (this is like what is in the final paper)
# Requires installation of COMPASS via devtools::install_github("malisas/COMPASS", ref="plotColorsCytokineHeight")
#
# The parameters cytokine_annotation_colors and cytokine_height_multiplier are not yet part of the main COMPASS repository,
# so a regular COMPASS installation will not support these options

# CEF CD8
# Heatmap Panel A
p <- print(plot(orderHeatmapColumnsByCytokinePresence(CD8CEFCompassResult, "IFNg",
                                                      cats_df_cytokine_order_override = c("CD154", "IL2", "TNF", "CD107a", "IL4", "IL17a", "IFNg")),
                grouping, show_rownames = FALSE,
                fontsize=16, fontsize_row=16, fontsize_col=19,
                row_annotation_colors=row_ann_colors,
                cytokine_annotation_colors=cytokine_annotation_colors,
                cytokine_height_multiplier=2.5, order_by_max_functionality = F,
                fixed_column_order = T))

svg(filename=file.path(outdir, "FigureS1A_CEF_CD8_Heatmap_alt.svg"), width=10, height=11)
grid.draw(p)
dev.off()

svglite(file=file.path(outdir, "FigureS1A_CEF_CD8_Heatmap_alt_svglite.svg"), width=10, height=11)
grid.draw(p)
dev.off()
