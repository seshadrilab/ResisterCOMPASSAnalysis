# Collection of functions from the WIP https://github.com/seshadrilab/flowHelpers for use in this project

#' Plot cell counts for each sample (default CD3+ cells)
#'
#' Plot the number of cells that are CD3+ (or your specified population) for each sample.
#' Run this function once for each batch to be combined later. Can use to look for batch effects.
#'
#' @param flowJoXmlPath full path to FlowJo Xml file for import (flowJoXmlPath or gatingSetPath or gatingSet is required)
#' @param gatingSetPath path to saved GatingSet directory for import (flowJoXmlPath or gatingSetPath or gatingSet is required)
#' @param gatingSet a GatingSet object (flowJoXmlPath or gatingSetPath or gatingSet is required)]
#' @param batch the batch name, required if data is passed in as a gatingSet
#' @param stratifyByLevel1 Required keyword on which to stratify boxplot data (usually "PATIENT ID")
#' @param fcsPath (Optional) directory where FCS files reside, if not in the same directory as flowJoXmlPath
#' @param keywords2import (Optional) keywords to import from flowJo into GatingSet metadata, required for flowJoXmlPath
#' @param outdir (Optional) directory to save boxplot in. If not given, returns the boxplot instead.
#' @param sampleGroup (Optional) Specify this if you are providing a flowJoXmlPath and the flowJo sample group is not 3
#' @param subpopulation (Optional) node name to use for counts, if not 3+
#' @param stratifyByLevel2 (Optional) additional keyword on which to stratify boxplot data
#' @param threshold (Optional) where to draw the threshold cutoff line for number of cells, default 25,000
#' @param yUpperExpand (Optional) This will expand the yaxis to have this upper limit, if applicable.
#' @return Boxplot of cell counts stratified by stratifyByLevel1, unless outdir is specified.
#' @export boxplot.cell.counts
#' @keywords QC counts
#' @usage boxplot.cell.counts <- function(flowJoXmlPath=NULL, gatingSetPath=NULL,
#'                     fcsPath=if (!.isnull(flowJoXmlPath)) { dirname(flowJoXmlPath) },
#'                     outdir=NULL, sampleGroup=3, subpopulation="3+",
#'                     keywords2import=NULL, stratifyByLevel1, stratifyByLevel2=NULL)
#' @examples
#' \dontrun{
#' boxplot.cell.counts(flowJoXmlPath="/home/Batch1Directory/Batch1FlowJo.xml",
#'                     keywords2import=c("PATIENT ID", "Barcode"),
#'                     stratifyByLevel1="PATIENT ID",
#'                     stratifyByLevel2="Barcode")
#'                     }
boxplot.cell.counts <- function(flowJoXmlPath=NULL,
                                gatingSetPath=NULL,
                                gatingSet=NULL,
                                fcsPath=if (!.isnull(flowJoXmlPath)) { dirname(flowJoXmlPath) },
                                outdir=NULL,
                                sampleGroup=3,
                                subpopulation="3+",
                                keywords2import=NULL,
                                stratifyByLevel1,
                                stratifyByLevel2=NULL,
                                batch=NULL,
                                threshold=25000,
                                yUpperExpand=NULL
) {
  # Check that required arguments are provided
  if (is.null(flowJoXmlPath) & is.null(gatingSetPath) & is.null(gatingSet)) {
    stop("flowJoXmlPath or gatingSetPathor gatingSet parameter must be provided.")
  }
  if (!is.null(flowJoXmlPath)) {
    if (is.null(keywords2import)) {
      stop("keywords2import parameter must be provided.")
    }
    if (!(stratifyByLevel1 %in% keywords2import)) {
      stop("stratifyByLevel1 must exist within keywords2import.")
    }
    if (!(is.null(stratifyByLevel2))) {
      if (!(stratifyByLevel2 %in% keywords2import)) {
        stop("stratifyByLevel2 must exist within keywords2import.")
      }
    }
  } else if (!is.null(gatingSet)) {
    if (is.null(batch)) {
      stop("batch parameter must be provided.")
    }
  }
  if (is.null(stratifyByLevel1)) {
    stop("stratifyByLevel1 parameter must be provided.")
  }
  
  gs <- gatingSet
  if (!is.null(flowJoXmlPath)) {
    # Read in the workspace
    cat(paste(c("Opening ", flowJoXmlPath, "\n"), collapse=""))
    ws <- flowWorkspace::openWorkspace(flowJoXmlPath)
    
    # Read in the sample fcs files as a GatingSet
    # The keywords option, a character vector, specifies the keywords to be extracted as pData of GatingSet
    gs <- flowWorkspace::parseWorkspace(ws, name=sampleGroup, path=fcsPath, keywords=keywords2import)
  } else if (!is.null(gatingSetPath)) {
    gs <- flowWorkspace::load_gs(gatingSetPath)
  }
  # Make the name the rownames
  flowWorkspace::pData(gs)[,"name"] <- rownames(pData(gs))
  
  # Obtain the subpopulation cell counts
  # Counts indicate flowCore recomputed counts, not FlowJo amounts
  popStats <- flowWorkspace::getPopStats(gs, subpopulations = subpopulation)
  
  # Merge the pData and popStats data tables together.
  annotatedCounts <- merge(flowWorkspace::pData(gs), popStats, by="name")
  # Add a column labeling points which have less than threshold cells
  annotatedCounts$pointLabels <- ifelse(annotatedCounts$Count < threshold, annotatedCounts$name, as.numeric(NA))
  
  if (!is.null(flowJoXmlPath)) {
    closeWorkspace(ws)
  }
  
  batchName <- if (!is.null(batch)) {
    batch
  } else if (!is.null(flowJoXmlPath)) {
    tools::file_path_sans_ext(basename(flowJoXmlPath))
  } else {
    tools::file_path_sans_ext(basename(gatingSetPath))
  }
  yExpandedLims <- if(is.null(yUpperExpand)) { 0 } else { c(0, yUpperExpand)}
  countsBoxplot <- {
    if (is.null(stratifyByLevel2)) {
      plotTitle <- paste(c(subpopulation, " flowCore Counts\nfor ", batchName,
                           "\nGrouped by ", stratifyByLevel1), collapse="")
      plotCounts <- ggplot2::ggplot(annotatedCounts, ggplot2::aes(x=factor(get(stratifyByLevel1)), y=Count)) +
        ggplot2::geom_boxplot() +
        ggplot2::stat_summary(fun.y=mean, geom="point", shape=8, size=4) +
        ggplot2::geom_point(shape=16) +
        ggplot2::labs(title=plotTitle) + ggplot2::xlab(stratifyByLevel1) +
        ggplot2::geom_hline(yintercept=threshold) +
        ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5, size=30), axis.text=ggplot2::element_text(size=16),
                       axis.title=ggplot2::element_text(size=22,face="bold"), legend.position="none", strip.text.x = ggplot2::element_text(size = 15)) +
        ggplot2::geom_text(ggplot2::aes(label=pointLabels), na.rm = TRUE) +
        ggplot2::expand_limits(y = yExpandedLims)
      plotCounts
    } else {
      plotTitle <- paste(c(subpopulation, " flowCore Counts\nfor ", batchName,
                           "\nGrouped by ", stratifyByLevel1, " and ", stratifyByLevel2), collapse="")
      plotCounts <- ggplot2::ggplot(annotatedCounts, ggplot2::aes(x=factor(get(stratifyByLevel2)), y=Count)) +
        ggplot2::geom_boxplot() +
        ggplot2::stat_summary(fun.y=mean, geom="point", shape=8, size=4) +
        ggplot2::geom_point(shape=16) +
        ggplot2::labs(title=plotTitle) + ggplot2::xlab(stratifyByLevel1) +
        ggplot2::geom_hline(yintercept=threshold) +
        ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5, size=30), axis.text=ggplot2::element_text(size=16),
                       axis.title=ggplot2::element_text(size=22,face="bold"), legend.position="none", strip.text.x = ggplot2::element_text(size = 15)) +
        ggplot2::facet_grid(as.formula(paste(c("~ ", "`", stratifyByLevel1, "`"), collapse="")), scales="free_x") +
        ggplot2::geom_text(ggplot2::aes(label = pointLabels), na.rm = TRUE) +
        ggplot2::expand_limits(y = yExpandedLims)
      plotCounts
    }
  }
  
  if (is.null(outdir)) {
    countsBoxplot
  } else {
    subpopulationFmtd <- gsub("/", "_", subpopulation)
    pngName <- { if (is.null(stratifyByLevel2)) {
      paste(c("QC_Boxplot_", subpopulationFmtd, "_Counts_By_", stratifyByLevel1, "_", c(batchName, ".png")), collapse="")
    } else {
      paste(c("QC_Boxplot_", subpopulationFmtd, "_Counts_By_", stratifyByLevel1, "_and_", stratifyByLevel2, "_", batchName, ".png"), collapse="")
    }}
    cat(paste(c("Saving png to ", file.path(outdir, pngName), "\n"), collapse=""))
    png(file.path(outdir, pngName), width=1500, height=900)
    print(countsBoxplot)
    dev.off()
    
    # Write annotatedCounts to file
    countsFileName <- paste(c("QC_Annotated_Counts_", subpopulationFmtd, "_", batchName, ".txt"), collapse="")
    write.table(annotatedCounts, file=file.path(outdir, countsFileName), sep="\t")
  }
}

### PrepareGatingSetList4Compass.R #######################################################

#' Prepare flow data for COMPASS
#'
#' This function accepts GatingSets and/or reads in multiple FlowJo workspace xml files and their associated FCS files.
#' It then combines them all into a single GatingSetList, extracting the user-provided keywords
#' and checking that the batches are combine-able along the way.
#' GatingSets should already have keywords.
#' The final GatingSetList is saved to a user-provided output directory by default. You can then use it directly with COMPASS or modify it further.
#'
#' Note 1: All the batches need to have the same gating tree, so this function will drop unshared
#' nodes and channels from the batches to be merged. TODO: make option to turn off default
#' Note 2: Marker names and channel names must all be the same across the batches.
#'
#' @param xmlFiles character vector (xmlFiles AND/OR gsDirs AND/OR gsList is required)
#' @param gsDirs character vector (xmlFiles AND/OR gsDirs AND/OR gsList is required)
#' @param gsList list of GatingSet objects (xmlFiles AND/OR gsDirs AND/OR gsList is required)
#' @param outDir directory in which to save the final GatingSetList
#' @param fcsFiles (Optional) character vector. Only needed if fcs and xml files exist in different directories.
#' @param sampleGroups (Optional) numeric vector. Specify this if the flowJo sample group is not 3 for all batches.
#' @param keywords2import (Optional) character vector. List of keywords to import from FlowJo workspace and into pData.
#' @param keyword4samples2exclude (Optional) keyword used to identify samples for exclusion
#' @param samples2exclude (Optional) character vector. When making GatingSetList, excludes samples whose keyword4samples2exclude column matches patterns in this vector (using `grepl`).
#' @param returnGSList (Optional) Return the GatingSetList instead of saving it. Default FALSE
#' @return Nothing, or GatingSetList if requested
#' @keywords GatingSetList COMPASS
#' @export
#' @examples
#' \dontrun{
#' prepare.gating.set.list.4.compass(xmlFiles=c("/home/Batch1Files/Batch1FlowJoWorkspace.xml",
#'                                              "/home/Batch2Files/Batch2FlowJoWorkspace.xml",
#'                                              "/home/Batch3Files/Batch3FlowJoWorkspace.xml",
#'                                              "/home/Batch4Files/Batch4FlowJoWorkspace.xml"),
#'                                   outDir="/home/GatingSetListOutDir",
#'                                   keywords2import=c("Barcode", "Antigen", "PATIENT ID"),
#'                                   keyword4samples2exclude="Barcode",
#'                                   samples2exclude=c('1234567890', '1234567891', '1234567892', '1234567893', '1234567894'))
#'                                   }
#' TODO: Write tests. - When keyword4samples2exclude and gsList and/or gsDirs are provided. - When keyword4samples2exclude and xmlFiles are provided
prepare.gating.set.list.4.compass <- function(xmlFiles=NULL,
                                              fcsFiles=NULL,
                                              gsDirs=NULL,
                                              gsList=list(),
                                              outDir=NULL,
                                              sampleGroups=NULL,
                                              keywords2import=c(),
                                              keyword4samples2exclude=NULL,
                                              samples2exclude=NULL,
                                              returnGSList=FALSE) {
  if(is.null(xmlFiles) && is.null(gsDirs) && length(gsList) == 0) {stop("One or more of xmlFiles or gsDirs or gsList must be provided")}
  if(!returnGSList && is.null(outDir)) {stop("outDir must be provided")}
  if(!is.null(outDir) && length(list.files(outDir)) > 0) {
    stop("Please delete outDir")
  }
  
  # Load in all the GatingSets from gsDirs
  gsList <- gsList # not sure if this is necessary
  gsListLen <- length(gsList)
  if (!is.null(gsDirs)){
    for (i in seq_along(gsDirs)) {
      gsList[[gsListLen + i]] <- load_gs(gsDirs[i])
    }
  }
  
  # Read in all the FlowJo workspaces and their fcs files.
  if (is.null(sampleGroups)) {
    message("sampleGroups not found, defaulting to \"All Samples\"")
    sampleGroups <- rep("All Samples", length(xmlFiles))
  }
  wsList <- list()
  gsListLen <- length(gsList)
  if (!is.null(xmlFiles)) {
    for (i in 1:length(xmlFiles)) {
      wsList[[i]] <- flowWorkspace::openWorkspace(xmlFiles[i])
      if (!is.null(fcsFiles)) {
        # TODO: use parseWorkspace subset argument instead of subset.GatingSet   subset=grepl(paste(samples2exclude, collapse="|"), factor(get(keyword4samples2exclude)), invert=TRUE)
        gsList[[gsListLen + i]] <- flowWorkspace::parseWorkspace(wsList[[i]], name=sampleGroups[i], path=fcsFiles[i],
                                                                 keywords=unique(append(keywords2import, keyword4samples2exclude)))
        if (!is.null(keyword4samples2exclude) & !is.null(samples2exclude)) {
          gsList[[gsListLen + i]] <- subset.GatingSet(gsList[[gsListLen + i]], !grepl(paste(samples2exclude, collapse="|"), factor(get(keyword4samples2exclude))))
        }
      } else {
        gsList[[gsListLen + i]] <- parseWorkspace(wsList[[i]], name=sampleGroups[i],
                                                  keywords=unique(append(keywords2import, keyword4samples2exclude)))
        if (!is.null(keyword4samples2exclude) & !is.null(samples2exclude)) {
          gsList[[gsListLen + i]] <- subset.GatingSet(gsList[[gsListLen + i]], !grepl(paste(samples2exclude, collapse="|"), factor(get(keyword4samples2exclude))))
        }
      }
    }
  }
  
  # Now all the GatingSets are in gsList.
  # Next drop redundant nodes and channels in the GatingSets to help make them merge-able.
  gsGroups <- flowWorkspace::groupByTree(gsList)
  nodes2Remove <- flowWorkspace::checkRedundantNodes(gsGroups)
  if (!(length(nodes2Remove) == 1 & length(nodes2Remove[[1]]) == 0)) {
    paste(as.character(Sys.time()), "WARNING: Removing nodes:\n")
    paste(nodes2Remove)
    paste("\n")
  }
  flowWorkspace::dropRedundantNodes(gsGroups, nodes2Remove) # original GatingSets in gsList are modified via external pointers
  
  # And then drop any redundant channels in the GatingSets.
  # This specifically removes channels from a GatingSet if there are no nodes which have
  # the given channel as an associated gating channel/marker.
  # TODO: I am modifying these in place and overwriting the old gsList contents. Not sure if that's kosher.
  for (i in seq_along(gsList)) {
    gsList[[i]] <- flowWorkspace::dropRedundantChannels(gsList[[i]])
  }
  
  # There is no surefire way to obtain pairs of matched marker names and channel names,
  # but some sanity checking is possible:
  # Obtain the common marker and channel pairings:
  commonMarkerChannelMappings <- Reduce(merge, lapply(gsList, function(x) {
    # parameters of the first flowFrame in x, the GatingSet
    pData(parameters(flowWorkspace::getData(x[[1]])))[,1:2] }))
  paste("Common marker and channel pairings across all GatingSets:")
  print(commonMarkerChannelMappings)
  
  # The GatingSets should now be ready to combine into one GatingSetList
  gsList4COMPASS <- flowWorkspace::GatingSetList(gsList)
  # Subset the GatingSetList based on keyword4samples2exclude, if applicable
  # It might be more correct to do this when the GatingSets are first read in (prior to the marker/gating checking above), but this way is quicker
  if (!is.null(keyword4samples2exclude) & !is.null(samples2exclude)) {
    if (keyword4samples2exclude %in% colnames(pData(gsList4COMPASS))) {
      gsList4COMPASS <- subset(gsList4COMPASS, !grepl(paste(samples2exclude, collapse="|"), factor(get(keyword4samples2exclude))))
    } else {
      message(paste(keyword4samples2exclude, " not in colnames(pData(gsList4COMPASS))"))
    }
  }
  
  if (returnGSList) {
    gsList4COMPASS
  } else {
    # Save the new GatingSetList to disk
    # Delete outDir and its subdirectories
    paste(c(as.character(Sys.time()), " Overwriting outDir: ", outDir, "\n"), collapse="")
    unlink(outDir, recursive=TRUE)
    flowWorkspace::save_gslist(gsList4COMPASS, path=outDir)
    # Close all the workspaces, if applicable
    if (length(wsList) > 0) {
      for (i in 1:length(wsList)) {
        flowWorkspace::closeWorkspace(wsList[[i]])
      }
    }
  }
}

#' Run COMPASS once
#'
#' This function runs COMPASS once, saving output to disk.
#' Intended for use by the generic.compass.wrapper function below.
#' @param gs The GatingSet of GatingSetList
#' @param cnode Node on which to run COMPASS
#' @param nodemarkermap List mapping nodes to marker names
#' @param outdir Directory in which to save output, e.g. heatmaps
#' @param individuals pData column containing individual identifiers (rows of heatmap)
#' @param seed (Optional) Number to set seed to [default NULL]
#' @param grouping (Optional) pData columns on which to group rows in heatmap, as a character vector [default NULL]
#' @param uniqueruns (Optional) pData column identifying unique runs. Use if you need multiple runs. [default NULL]
#' @param lineplotxvar (Optional) pData column which defines groups along x-axis in FS-score line plot, e.g. Time [default NULL]
#' @param iter (Optional) Number of COMPASS iterations to perform on each repitition (8 repetitions total) [default 40,000]
#' @param lineplotgroupby (Optional) This should be specified if lineplotxvar is given. pData column which defines which values to connect in the line plot (usually something like "PTID")
#' @return Nothing
#' @keywords COMPASS
#' @import grid
#' @export
run.compass.once <- function(gs,
                             cnode,
                             individuals,
                             nodemarkermap,
                             iter=40000,
                             lineplotxvar=NULL,
                             run=NULL,
                             outdir,
                             uniqueruns=NULL,
                             grouping=NULL,
                             lineplotgroupby=NULL,
                             countFilterThreshold=0,
                             seed=NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Create a COMPASSContainer from the GatingSetList.
  CC <- COMPASS::COMPASSContainerFromGatingSet(gs, node=cnode, individual_id=individuals,
                                               mp=nodemarkermap, countFilterThreshold=countFilterThreshold)
  
  # Run COMPASS for this unique run
  fit <- COMPASS::COMPASS( CC,
                           treatment=trt == "Treatment",
                           control=trt == "Control",
                           iterations=iter
  )
  
  FS <- COMPASS::FunctionalityScore(fit)
  PFS <- COMPASS::PolyfunctionalityScore(fit)
  
  # Initialize subdirectory to save output for this run
  fileSuffix <- if (is.null(run)) { cnode } else { paste(cnode, run, sep="_") }
  subDir <- fileSuffix
  subDirPath <- file.path(outdir, subDir)
  dir.create(subDirPath, showWarnings = FALSE)
  setwd(subDirPath)
  
  # Save the COMPASS run as an RDS file for safekeeping
  saveRDS(fit, paste(c("COMPASSResult_", fileSuffix, ".rds"), collapse=""))
  
  # Initialize output text files for FS and PFS
  columnsFormat <- paste("CellSubset", uniqueruns, sep="\t")
  fsFile <- paste(paste("FS", fileSuffix, sep="_"), ".txt", sep="")
  pfsFile <- paste("P", fsFile, sep="")
  write(paste(columnsFormat, paste(names(FS), collapse="\t"), sep="\t"), file=fsFile, append=TRUE)
  write(paste(columnsFormat, paste(names(PFS), collapse="\t"), sep="\t"), file=pfsFile, append=TRUE)
  
  # Write results to file, one for each statistic type
  write(paste(paste(cnode, run, sep="\t"), paste(FS, collapse="\t"), sep="\t"), file=fsFile, append=TRUE)
  write(paste(paste(cnode, run, sep="\t"), paste(PFS, collapse="\t"), sep="\t"), file=pfsFile, append=TRUE)
  
  plotTitleSuffix <- paste(c(",\n", run, ", ", cnode, " Cells"), collapse="")
  cytokine_annotation_colors <- c("black", "black", "black", "black", "black", "black", "black")
  
  # Plot a heatmap of the mean probability of response, to visualize differences
  # in expression for each category
  png(filename=paste(c("HeatmapMeanProbResponse", "_", cnode, run, ".png"), collapse=""),
      width=800, height=650)
  try(grid::grid.draw(print(plot(fit, grouping, show_rownames = TRUE,
                                 main = paste("Heatmap of Mean Probability of Response", plotTitleSuffix, sep=""),
                                 fontsize=14, fontsize_row=13, fontsize_col=11,
                                 cytokine_annotation_colors=cytokine_annotation_colors))))
  dev.off()
  
  # Log scale of previous data, smaller changes show up better
  png(filename=paste(c("HeatmapLogPostResponse", "_", cnode, run, ".png"), collapse=""),
      width=800, height=650)
  try(grid::grid.draw(print(plot(fit, grouping, show_rownames = TRUE,
                                 measure=COMPASS::PosteriorLogDiff(fit), threshold=0,
                                 main = paste("Heatmap of Log Posterior Differences in Response", plotTitleSuffix, sep=""),
                                 fontsize=14, fontsize_row=13, fontsize_col=11,
                                 cytokine_annotation_colors=cytokine_annotation_colors))))
  dev.off()
  
  # If applicable, create line plot of Functionality Score vs. lineplotxvar
  try(if (!is.null(lineplotxvar)) {
    # First format FS into a data.table
    # The order of rows in metadata is not necessarily the same as those of FSdf. Perform a table merge
    FSdf <- data.frame(names(FS), FS, row.names=NULL)
    names(FSdf) <- c(individuals, "FunctionalityScore")
    metasub <- pData(gsListForCOMPASSsub)
    fsplotdf <- merge(x=metasub[metasub[,uniqueruns]==run,], y=FSdf, by.x=individuals, by.y=individuals)
    # Draw the line plot
    lineplot <- ggplot2::ggplot(data=fsplotdf, aes(x=factor(get(lineplotxvar)), y=FunctionalityScore, group=factor(get(lineplotgroupby)))) +
      ggplot2::geom_point() + ggplot2::geom_line() +
      ggplot2::labs(title=paste(c("Functionality Score vs. ", lineplotxvar, " Line Plot,\n", cnode, " ", run), collapse=""),
                    x=lineplotxvar) +
      ggplot2::theme_set(theme_gray(base_size = 25))
    ggplot2::ggsave(filename=paste(c("LinePlot_FSv", lineplotxvar, "_", cnode, run, ".png"), collapse=""),
                    plot=lineplot,
                    width=6.66, height=7.85)
  })
  
  # Call the garbage collector to free up memory
  gc()
}


#' COMPASS Wrapper
#'
#' This function runs COMPASS once for each of the unique values defined by the uniqueruns argument (if provided)
#' @param path The path to the folder in which the GatingSetList or GatingSet is saved
#' @param cnode Node on which to run COMPASS
#' @param nodemarkermap List mapping nodes to marker names
#' @param outdir Directory in which to save output, e.g. heatmaps
#' @param individuals pData column containing individual identifiers (rows of heatmap)
#' @param seed (Optional) Number to set seed to [default NULL]
#' @param grouping (Optional) pData columns on which to group rows in heatmap, as a character vector [default NULL]
#' @param uniqueruns (Optional) pData column identifying unique runs. Use if you need multiple runs. [default NULL]
#' @param lineplotxvar (Optional) pData column which defines groups along x-axis in FS-score line plot, e.g. Time [default NULL]
#' @param iter (Optional) Number of COMPASS iterations to perform on each repitition (8 repetitions total) [default 40,000]
#' @param lineplotgroupby (Optional) This should be specified if lineplotxvar is given. pData column which defines which values to connect in the line plot (usually something like "PTID")
#' @return Nothing
#' @keywords COMPASS
#' @import grid
#' @export
#' @examples
#' \dontrun{
#' generic.compass.wrapper(path="/home/path/to/GatingSetList",
#'                         seed=1,
#'                         outdir="/path/to/OutDirectory",
#'                         cnode="8+",
#'                         nodemarkermap=list("8+/154+" = "CD154",
#'                                            "8+/IFNg+" = "IFNg",
#'                                            "8+/IL4+" = "IL4",
#'                                            "8+/TNFa+" = "TNFa",
#'                                            "8+/IL22+" = "IL22",
#'                                            "8+/IL17+" = "IL17a",
#'                                            "8+/IL2+" = "IL2"),
#'                         individuals="PATIENT ID",
#'                         uniqueruns="Peptide")
#'                         }
generic.compass.wrapper <- function(path=NULL,
                                    seed=NULL,
                                    outdir=NULL,
                                    cnode=NULL,
                                    nodemarkermap=NULL,
                                    individuals=NULL,
                                    grouping=NULL,
                                    uniqueruns=NULL,
                                    lineplotxvar=NULL,
                                    iter=40000,
                                    lineplotgroupby=NULL,
                                    countFilterThreshold=0) {
  # Assumptions: pData trt column contains "Treatment" and "Control" labels
  # TODO: run in parallel
  # TODO: test for single run from command line
  
  # Check that required arguments are provided
  if (is.null(path)) {
    stop("Path parameter must be provided.")
  }
  if (is.null(outdir)) {
    stop("Outdir parameter must be provided.")
  }
  if (is.null(nodemarkermap)) {
    stop("Nodemarkermap parameter must be provided.")
  }
  if (is.null(individuals)) {
    stop("Individuals parameter must be provided.")
  }
  if (is.null(cnode)) {
    stop("cnode parameter must be provided.")
  }
  if (!is.null(lineplotxvar) & is.null(lineplotgroupby)) {
    stop("Please specify lineplotgroupby parameter, which must be provided if lineplotxvar is provided.")
  }
  
  cat(paste(as.character(Sys.time()), "Loading GatingSetList\n"))
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Load the saved GatingSetList or GatingSet:
  loadGSListOrGS <- function (path) {
    out <- try(flowWorkspace::load_gslist(path))
    if (class(out) == "try-error") {
      cat("Caught an error during flowWorkspace::load_gslist, trying flowWorkspace::load_gs.\n")
      out <- flowWorkspace::load_gs(path)
    }
    out
  }
  
  gsListForCOMPASS <- loadGSListOrGS(path)
  meta <- flowWorkspace::pData(gsListForCOMPASS)
  
  # Run COMPASS once or multiple times depending on whether uniqueruns is given
  if (is.null(uniqueruns)) {
    cat(paste(c(as.character(Sys.time()), " Running COMPASS for", cnode, "cells", "...\n"), collapse=" "))
    try(run.compass.once(gs=gsListForCOMPASS,
                         cnode=cnode,
                         individuals=individuals,
                         nodemarkermap=nodemarkermap,
                         iter=iter,
                         lineplotxvar=lineplotxvar,
                         run=NULL,
                         outdir=outdir,
                         uniqueruns=NULL,
                         grouping=grouping,
                         countFilterThreshold=countFilterThreshold))
  } else {
    # Get a list of values identifying unique runs from the "uniqueruns" column, minus the control value
    uniqueRunsList <- unique(meta[meta[,"trt"]!="Control",][,uniqueruns])
    # Obtain the value in the "uniqueruns" column that corresponds to the control value
    controlval <- unique(meta[meta[,"trt"]=="Control",][,uniqueruns])
    
    for (run in uniqueRunsList) {
      cat(paste(c(as.character(Sys.time()), " Running COMPASS for", cnode, "cells ", run, "...\n"), collapse=" "))
      
      # Subset the GatingSetList for the desired unique run
      rowSubsetBooleans <- meta[,uniqueruns] %in% c(run, controlval)
      rownamesAsChar <- as.character(rownames(meta[rowSubsetBooleans,]))
      gsListForCOMPASSsub <- gsListForCOMPASS[rownamesAsChar]
      
      try(run.compass.once(gs=gsListForCOMPASSsub,
                           cnode=cnode,
                           individuals=individuals,
                           nodemarkermap=nodemarkermap,
                           iter=iter,
                           lineplotxvar=lineplotxvar,
                           run=run,
                           outdir=outdir,
                           uniqueruns=uniqueruns,
                           grouping=grouping,
                           lineplotgroupby=lineplotgroupby,
                           countFilterThreshold=countFilterThreshold))
    }
  }
  cat(paste(as.character(Sys.time()), " All COMPASS runs Done\n"))
}

#' Order COMPASSresult categories columns for cytokine of interest
#' 
#' This function rearranges the mean_gamma and categories matrices so that subsets containing the cytokine of interest appear on the right.
#' The categories legend is rearranged such that the cytokine of interest is on the top row (so, top right cells will be shaded in).
#' The other cytokines are also rearranged in a way to be as aesthetically plasing as possible.
#' Or, specify your own cytokine order for the legend with cats_df_cytokine_order_override (doesn't override the cytokine-of-interest though).
#' 
#' @param cr a COMPASSresult
#' @param cytokine The cytokine of interest
#' @param cats_df_cytokine_order_override Use this if you want to override the sorting of the categories matrix cytokine. Don't include "Counts".
#' @param threshold
#' @param minimum_dof
#' @param maximum_dof
#' @export
#' @import data.table
#' @examples 
#' \dontrun{
#' new_cr <- orderHeatmapColumnsByCytokinePresence(myCompassResult, "IFNg")
#' # Good for the fixed_column_order = T option from modified plot function
#' p <- print(plot(new_cr, order_by_max_functionality = F,
#'                 fixed_column_order = T))
#' grid.draw(p)
#' }
orderHeatmapColumnsByCytokinePresence <- function(cr, cytokine, cats_df_cytokine_order_override=NULL, threshold=0.01, minimum_dof=1, maximum_dof=Inf) {
  # Return a version of the COMPASSResult where:
  # - cr$fit$mean_gamma columns are ordered so that cytokine-of-interest-positive subsets are listed last
  #     - Don't need to deal with the fully-cytokine-negative column, the plotting function filters it out somewhere
  #     - cr$fit$categories rows are re-ordered to be the same as the new mean_gamma columns order
  #     - Order categories-df columns by functionality after cytokine-presence sorting (so, keep as much of original functionality-sorted order as possible)
  # - cr$fit$categories cytokine columns are ordered by popularity (i.e. # of subsets each cytokine appears in)
  
  # First assign cr$fit$categories rownames to make re-ordering easier later
  rownames(cr$fit$categories) <- colnames(cr$fit$mean_gamma)
  
  # Make the column with the cytokine of interest show up second-to-last (before the Counts column) in the cr$fit$categories matrix columns so that it is on the top of the heatmap cytokine annotations section
  # Also order the other cytokines by the number of subsets they appear in
  if(is.null(cats_df_cytokine_order_override)) {
    # Calculate the order based on the filtered subsets
    #
    # like https://github.com/malisas/COMPASS/blob/b11fe4ddd6c4bfa3518dd812ba40c745dd5c72eb/R/plotMeanGamma.R#L231
    ind <- 1:ncol(cr$fit$mean_gamma)
    dof <- cr$fit$categories[, "Counts"]
    dof_ind <- which(dof >= minimum_dof & dof <= maximum_dof)
    ind <- intersect(ind, dof_ind)
    # much like https://github.com/malisas/COMPASS/blob/b11fe4ddd6c4bfa3518dd812ba40c745dd5c72eb/R/plotMeanGamma.R#L235
    m <- apply(cr$fit$mean_gamma, 2, function(x) {
      mean(x, na.rm = TRUE)
    })
    keep <- m > threshold
    ind <- intersect(ind, which(keep))
    
    cats_tmp <- cr$fit$categories[ind,]
    col_order_part1 <- setdiff(colnames(cats_tmp)[order(colSums(cats_tmp), decreasing=T)],
                               c(cytokine, "Counts"))
    cr$fit$categories <- cr$fit$categories[,c(col_order_part1, c(cytokine, "Counts"))]
  } else {
    cr$fit$categories <- cr$fit$categories[,c(cats_df_cytokine_order_override, "Counts")]
  }
  
  # First order the mean_gamma columns by 1) number of positive cytokines,
  # and then 2) based on the order of col_order_part1
  cats_tmp_2 <- as.data.table(cr$fit$categories, keep.rownames = T)[,1:(ncol(cr$fit$categories)-1)]
  setorderv(cats_tmp_2, colnames(cats_tmp_2)[-1])
  subsets_order_2 <- match(colnames(cr$fit$mean_gamma), rev(cats_tmp_2$rn)) # the subsets hierarchically ordered by presence in cr$fit$categories columns
  mean_gamma_cols_preordered <- colnames(cr$fit$mean_gamma)[order(-nchar(gsub("[^!]", "", colnames(cr$fit$mean_gamma))), subsets_order_2)]
  # Now order the mean_gamma columns based on presence of the cytokine
  no_cytokine_cols <- grep(paste0("!", cytokine), mean_gamma_cols_preordered,
                           fixed = T, value = T)
  cytokine_cols <- setdiff(mean_gamma_cols_preordered, no_cytokine_cols)
  new_order <- c(no_cytokine_cols, cytokine_cols)
  # Assign new order
  cr$fit$mean_gamma <- cr$fit$mean_gamma[,new_order]
  cr$fit$categories <- cr$fit$categories[new_order,]
  
  cr
}

#' Plot of Functionality Score vs some other axis
#' 
#' TODO: Make another function for paired values (i.e. line plot)
#' @param compassResultOrPath the COMPASSResult object or path to a COMPASSResult RDS object on disk
#' @param gsOrGsListOrPath (optional) Either a GatingSet, a GatingSetList, or path to one of these objects on disk. You can specify this if you want to use its metadata instead of the metadata in the COMPASSResult object
#' @param plotTextExtra Extra text (e.g. "CD4+, Peptide Pool 1") to place in plot title and filename. Commas and spaces, etc will be removed for filename.
#' @param outdir If given, plot will be saved to this directory instead of being returned in a list along with plot data
#' @param stratifyBy Character vector of metadata columns to stratify plot by
#' @param themeBaseSize
#' @param removeGridAndBg
#' @param showTitle
#' @param showSignificanceBracket
#' @param polyfunctionality Display Polyfunctionality Score instead of Functionality Score
#' @param pvalue_fontsize
#' @param axestitle_fontsize
#' @param axestick_fontsize
#' @param font
#' @param geom_jitter_width
#' @param xaxis_title
#' @export
#' @import coin
#' @import COMPASS
#' @import flowWorkspace
#' @import ggplot2
#' @import plyr
#' @examples 
#' \dontrun{
#' fsPlotList <- fs.plot(gsOrGsListOrPath="path/to/GatingSet/Folder",
#'        compassResultOrPath="path/to/CompassResultRDSFile.rds",
#'        stratifyBy="DiseaseState",
#'        plotTextExtra="",
#'        plotWilcox=TRUE)
#' }
fs.plot <- function(compassResultOrPath,
                    stratifyBy,
                    ylims=NULL,
                    outdir=NULL,
                    plotTextExtra="",
                    plotWilcox=FALSE,
                    themeBaseSize=18,
                    removeGridAndBg=FALSE,
                    showTitle=TRUE,
                    showSignificanceBracket=TRUE,
                    gsOrGsListOrPath=NULL,
                    polyfunctionality=F,
                    pvalue_fontsize=NULL,
                    axestitle_fontsize=NULL,
                    axestick_fontsize=NULL,
                    font=NULL,
                    geom_jitter_width=0.15,
                    xaxis_title=stratifyBy) {
  gs <- if (!is.null(gsOrGsListOrPath)) {
    if(class(gsOrGsListOrPath) == "GatingSet" || class(gsOrGsListOrPath) == "GatingSetList") {
      gsOrGsListOrPath
    } else {
      try(if(!(class(gsOrGsListOrPath) == "character")) stop("gsOrGsListOrPath must be either a GatingSet, a GatingSetList, or the path to the folder containing one of these objects on disk"))
      # Load the saved GatingSetList or GatingSet:
      loadGSListOrGS <- function (gsOrGsListOrPath) {
        out <- try(flowWorkspace::load_gslist(gsOrGsListOrPath))
        if (class(out) == "try-error") {
          cat("Caught an error during flowWorkspace::load_gslist, trying flowWorkspace::load_gs.\n")
          out <- flowWorkspace::load_gs(gsOrGsListOrPath)
        }
        out
      }
      loadGSListOrGS(gsOrGsListOrPath)
    }}
  
  cr <- if(class(compassResultOrPath) == "COMPASSResult") {
    compassResultOrPath
  } else {
    try(if(!(class(compassResultOrPath) == "character")) stop("compassResultOrPath must be a COMPASSResult object or path to a COMPASSResult rds file on disk"))
    # Load the saved COMPASSResult
    readRDS(compassResultOrPath)
  }
  fsTable <- if(polyfunctionality) {as.data.table(PolyfunctionalityScore(cr), keep.rownames = TRUE)} else {as.data.table(FunctionalityScore(cr), keep.rownames = TRUE)}
  individualIdentifier <- cr$data$individual_id # pData(gs)/cr$data$meta and the fsTable should have this column
  colnames(fsTable) <- c(individualIdentifier, "Score")
  
  meta <- if (!is.null(gsOrGsListOrPath)) { pData(gs) } else { cr$data$meta  }
  
  pData4Plot <- as.data.table(meta[,c(individualIdentifier, stratifyBy)])
  setkeyv(pData4Plot, individualIdentifier)
  pData4Plot <- unique(pData4Plot)
  fsTable <- merge(fsTable, pData4Plot, by=individualIdentifier)
  
  xaxis <- stratifyBy
  yaxis <- "Score"
  groupName <- individualIdentifier
  ylimits <- if(is.null(ylims)) { c(min(fsTable$Score) - 0.004, max(fsTable$Score) + 0.004) } else { ylims }
  p <- ggplot2::ggplot(data=fsTable, ggplot2::aes_string(x=xaxis, y=yaxis))
  p <- p + ggplot2::geom_boxplot(inherit.aes=FALSE, ggplot2::aes_string(x=xaxis, y=yaxis, fill = xaxis), colour = "black", outlier.shape = NA)
  p <- p +
    ggplot2::geom_jitter(width=geom_jitter_width) +
    ggplot2::labs(x=xaxis_title,
                  y=if(polyfunctionality) {"Polyfunctionality Score"} else {"Functionality Score"}) +
    ggplot2::theme_set(ggplot2::theme_gray(base_size = themeBaseSize)) +
    ggplot2::coord_cartesian(ylim=ylimits) +
    ggplot2::theme(axis.text=element_text(colour="black"))
  if(!is.null(axestick_fontsize)) {
    p <- p + ggplot2::theme(axis.text=element_text(size=axestick_fontsize))
  }
  if(!is.null(axestitle_fontsize)) {
    p <- p + ggplot2::theme(axis.title=element_text(size=axestitle_fontsize))
  }
  if(!is.null(font)) {
    p <- p + ggplot2::theme(text = element_text(family=font))
  }
  
  if(showTitle) {
    title <- sprintf("%s Score Boxplot\nBy %s\n%s",
                     if(polyfunctionality) {"Polyfunctionality"} else {"Functionality"},
                     stratifyBy, plotTextExtra)
    p <- p + ggplot2::labs(title=title)
  }
  if(removeGridAndBg) {
    p <- p + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                            panel.background = ggplot2::element_blank(), axis.line = element_line(colour = "black"))
  }
  
  # Wilcox rank sum test between groups
  testResult <- if(length(stratifyBy) == 1 && length(levels(as.factor(fsTable[,get(stratifyBy)]))) == 2) {
    fsTable[,stratifyBy] <- as.factor(fsTable[,get(stratifyBy)])
    tr <- coin::wilcox_test(Score ~ get(stratifyBy), data=fsTable)
    # p_value_text <- if(coin::pvalue(tr) < 0.001) {"p<0.001"} else {paste0("p=", signif(coin::pvalue(tr), digits=3))}
    p_value_text <- if(coin::pvalue(tr) < 0.001) {"p<0.001"} else {paste0("p=", format(round(coin::pvalue(tr), 3), nsmall = 3))}
    if(plotWilcox) {
      p <- p + ggplot2::geom_text(label = p_value_text,
                                  x = 1.5,
                                  y = max(fsTable$Score) + 0.003,
                                  colour="black",
                                  parse=FALSE,
                                  size=5)
    }
    
    if(showSignificanceBracket) {
      y_Bracket <- max(fsTable$Score)*1.04
      if(!is.null(pvalue_fontsize)) {
        if(!is.null(font)) {
          p <- p + ggsignif::geom_signif(annotation=p_value_text,
                                         y_position=y_Bracket, xmin=0.85, xmax=2.15, 
                                         tip_length = c(0.01, 0.01),
                                         textsize=pvalue_fontsize,
                                         family=font)
        } else {
          p <- p + ggsignif::geom_signif(annotation=p_value_text,
                                         y_position=y_Bracket, xmin=0.85, xmax=2.15, 
                                         tip_length = c(0.01, 0.01),
                                         textsize=pvalue_fontsize)
        }
      } else {
        if(!is.null(font)) {
          p <- p + ggsignif::geom_signif(annotation=p_value_text,
                                         y_position=y_Bracket, xmin=0.85, xmax=2.15, 
                                         tip_length = c(0.01, 0.01),
                                         family=font)
        } else {
          p <- p + ggsignif::geom_signif(annotation=p_value_text,
                                         y_position=y_Bracket, xmin=0.85, xmax=2.15, 
                                         tip_length = c(0.01, 0.01))
        }
      }
    }
    tr
  }
  
  if(!is.null(outdir)) {
    filePrefix <- gsub(" ", "_", gsub("[`!@#$%^&*(),?]", "", plotTextExtra)) # replace spaces with underscores and other symbols with an empty string
    filename <- sprintf("%s_Compass_%s_Plot.png", filePrefix, if(polyfunctionality) {"PFS"} else {"FS"}) 
    ggplot2::ggsave(filename=file.path(outdir, filename),
                    plot=p,
                    width=6.66, height=7.85)
  } else {
    list(plot = p, data = fsTable, test = testResult)
  }
}

#' Overlay and Highlight Polyfunctional Cell Subsets on a Flow plot
#'
#' Defines a function to Overlay and Highlight Polyfunctional Cell Subsets on a Flow plot
#' This function assumes there are two columns, conditioncol and conditioncol2, upon which you are stratifying the plots
#'
#' @param path path to directory holding GatingSetList or GatingSet
#' @param gsOrGsList GatingSet or GatingSetList object
#' @param individualsCol column which defines individual
#' @param individual value of individual(s) in individualsCol whose data you want to plot
#' @param conditioncol name of the column that defines the main experimental condition, e.g. Antigen
#' @param exp experimental value in conditioncol, e.g. ESAT-6
#' @param ctrl control value in conditioncol, e.g. DMSO
#' @param conditioncol2 second condition on which to stratify Flow plots, e.g. "PATIENT ID"
#' @param parentsubset unique name of parent node to use for plots
#' @param boolsubset the full boolean subset to be used by booleanfilter()
#' @param xaxis a marker name to plot on the x-axis
#' @param yaxis a marker name to plot on the y-axis
#' @param width width, in inches, of the plot
#' @param outdir (Optional) saves image in output directory, if given
#' @param facetorder (Optional) the levels of conditioncol (e.g. Antigen) in the order you want displayed
#' @param facetorder2 (Optional) the levels of conditioncol2 (e.g. Time) in the order you want displayed
#' @param overlayDotSize
#' @param themeBaseSize
#' @param xlims
#' @param ylims
#' @param axestitle_fontsize
#' @param font
#' @param percentage_fontsize
#' @param geom_hex_bins
#' @return Flow plot, unless outdir is specified
#' @import data.table
#' @import flowWorkspace
#' @import grDevices
#' @import svglite
#' @export
#' @keywords Flow Plot Polyfunctional Subset
#' @examples
#' \dontrun{
#' highlight.boolean.subset.flow.plot(path="/home/path/to/GatingSetListAllBatches",
#'                                    individualsCol="PTID",
#'                                    individual=12345678,
#'                                    conditioncol="Antigen",
#'                                    exp="ESAT-6",
#'                                    ctrl="DMSO",
#'                                    conditioncol2="Time",
#'                                    parentsubset="8+",
#'                                    boolsubset="8+/TNFa+&!8+/IFNg+&!8+/IL2+&!8+/IL4+",
#'                                    xaxis="TNFa",
#'                                    yaxis="IFNg",
#'                                    facetorder=c("DMSO", "ESAT-6"))
#'                                    }
highlight.boolean.subset.flow.plot <- function(path,
                                               gsOrGsList=NULL,
                                               outdir=NULL,
                                               individualsCol,
                                               individual,
                                               conditioncol,
                                               exp,
                                               ctrl,
                                               conditioncol2=".",
                                               parentsubset,
                                               boolsubset,
                                               xaxis,
                                               yaxis,
                                               facetorder=NULL,
                                               facetorder2=NULL,
                                               geomTextY=5,
                                               geomTextX=200,
                                               pngORsvg="png",
                                               width=NULL,
                                               overlayDotSize=0.4,
                                               themeBaseSize=18,
                                               stripLegendGridAxesTitle=FALSE,
                                               xlims=NULL,
                                               ylims=NULL,
                                               axestitle_fontsize=NULL,
                                               font=NULL,
                                               percentage_fontsize=NULL,
                                               geom_hex_bins=120
                                               
) {
  # TODO: check all required parameters exist
  #library(flowWorkspace) # flowWorkspace::add doesn't seem to work w/o this line
  
  gs <- if(!is.null(gsOrGsList)) {
    gsOrGsList
  } else {
    # Load the saved GatingSetList or GatingSet:
    loadGSListOrGS <- function (path) {
      out <- try(flowWorkspace::load_gslist(path))
      if (class(out) == "try-error") {
        cat("Caught an error during flowWorkspace::load_gslist, trying flowWorkspace::load_gs.\n")
        out <- flowWorkspace::load_gs(path)
      }
      out
    }
    loadGSListOrGS(path)
  }
  
  metaSub <- flowWorkspace::pData(gs)[intersect(which(flowWorkspace::pData(gs)[,individualsCol] %in% individual),
                                                union(which(flowWorkspace::pData(gs)[conditioncol] == exp),
                                                      which(flowWorkspace::pData(gs)[conditioncol] == ctrl))),]
  gsSub <- gs[rownames(metaSub)]
  
  boolsubsetName <- gsub("&", "and", gsub("!", "not_", gsub("/", ":", boolsubset)))
  addBooleanGate(gs=gs, booleanSubset=boolsubset, parentGate=parentsubset, overrideGate=FALSE, booleanGateName=boolsubsetName)
  boolsubsetPopStats <- flowWorkspace::getPopStats(gsSub, flowJo=FALSE, subpopulations=c(boolsubsetName))
  
  # When we plot the parent subset, we don't want to plot the events which are in the overlayed boolean subset as well as in the parent subset twice.
  # So for the purposes of plotting only, define a new parent subset which has the overlayed plots removed from it.
  parentsubset_ForPlot_name <- sprintf("%s_no%s_ForPlot", parentsubset, boolsubsetName)
  addBooleanGate(gs=gs, booleanSubset=sprintf("%s&!%s", parentsubset, boolsubsetName), parentGate=parentsubset, overrideGate=FALSE, booleanGateName=parentsubset_ForPlot_name)
  
  # gsSubMetaData <- flowWorkspace::pData(gsSub)[,2:length(colnames(flowWorkspace::pData(gsSub)))]
  gsSubMetaData <- flowWorkspace::pData(gsSub)
  gsSubMetaData <- cbind(gsSubMetaData, rownames(gsSubMetaData))
  colnames(gsSubMetaData)[length(colnames(gsSubMetaData))] <- "row.names"
  gsSubMetaDataCols <- if (conditioncol2 == ".") { c("row.names", conditioncol) } else
  {c("row.names", conditioncol, conditioncol2) }
  
  # print(gsSubMetaDataCols)
  # print(head(boolsubsetPopStats[, c("name", "Population", "Count", "ParentCount")]))
  # print(head(gsSubMetaData[, gsSubMetaDataCols]))
  
  boolsubsetPopStatsMerge <- merge(x=boolsubsetPopStats[, c("name", "Population", "Count", "ParentCount")], y=gsSubMetaData[, gsSubMetaDataCols], by.x="name", by.y="row.names")
  
  library(data.table)
  byCols <- c(conditioncol, conditioncol2)#, individualsCol)
  byCols <- unique(byCols[byCols != "."])
  boolsubsetPopStatsMergeCollapsed <- rbind(boolsubsetPopStatsMerge[, {
    # cols2makeUnique <- .BY #.BY[, c("Day", "Treatment")]
    # cols2makeUnique <- unlist(lapply(cols2makeUnique, function(x) {
    #   x <- unique(x[!is.na(x)])
    #   if(length(x) == 1) as.character(x)
    #   else if(length(x) == 0) NA_character_
    #   else "multiple"
    # }))
    cols2Sum <- .SD[, c("Count", "ParentCount")]
    cols2Sum <- colSums(cols2Sum)
    data.table(t(unlist(cols2Sum)))
    #cbind(data.table(t(unlist(cols2makeUnique))), data.table(t(unlist(cols2Sum))))
    #cbind(.BY, data.table(t(unlist(cols2Sum))))
    #cbind(data.table(t(unlist(cols2makeUnique))), data.table(t(unlist(cols2Sum))))
  },
  by=byCols])
  
  boolsubsetPopStatsMergeCollapsed[, "Proportion"] <- boolsubsetPopStatsMergeCollapsed[, "Count"] / boolsubsetPopStatsMergeCollapsed[, "ParentCount"]
  # boolsubsetPopStatsMergeCollapsed[, "Percent"] <- sapply(formatC(base::round(with(boolsubsetPopStatsMergeCollapsed, Count / ParentCount) * 100, 3), 3, format="f"), function(x) paste(x, "%", sep=""), USE.NAMES=FALSE)
  boolsubsetPopStatsMergeCollapsed[, "Percent"] <- sapply(with(boolsubsetPopStatsMergeCollapsed, Count / ParentCount) * 100, function(x) { if(x < 0.005) {"0%"} else { paste0(format(round(x, 2), nsmall = 2), "%") }}, USE.NAMES = F)
  
  flowtitle <- if (conditioncol2 == ".") {
    paste(c(individualsCol, " ", paste(individual, collapse=", "), ", ", parentsubset, " cells\nResponse to ", exp), collapse="")
  } else {
    paste(c(individualsCol, " ", paste(individual, collapse=", "), ", ", parentsubset, " cells\nResponse to ", exp, " vs ", conditioncol2), collapse="")
  }
  
  # Simplify boolean subset for display
  subsetsmpl <- strsplit(boolsubset, split="&")[[1]]
  # What we're selecting positively FOR
  possubset <- subsetsmpl[grep("!", subsetsmpl, invert=TRUE)]
  possubsetFmtd <- paste("Pos: ", paste(lapply(possubset, function(x) {x[length(x)][[1]]}), collapse=", "), sep="")
  # What we're selecting AGAINST
  negsubset <- subsetsmpl[grep("!", subsetsmpl)]
  negsubsetFmted <- paste("Neg: ", paste(lapply(negsubset, function(x) {splt <- strsplit(x, "!")[[1]]; splt[length(splt)][[1]]}), collapse=", "), sep="")
  plottitle <- paste(c("Difference in cell subset proportions between ", exp, " and ", ctrl), collapse="")
  subtitle1 <- paste(c("Full Boolean Subset: \n       ", possubsetFmtd, "\n       ", negsubsetFmted), collapse="")
  
  # Order the levels of conditioncol and conditioncol2, if facet order is provided
  if (!is.null(facetorder)) {
    flowWorkspace::pData(gsSub)[,conditioncol] <- factor(flowWorkspace::pData(gsSub)[,conditioncol], levels=facetorder)
    if (conditioncol2 != "." && is.null(facetorder2)) {
      flowWorkspace::pData(gsSub)[,conditioncol2] <- factor(flowWorkspace::pData(gsSub)[,conditioncol2])
    }
    boolsubsetPopStatsMergeCollapsed <- as.data.frame(boolsubsetPopStatsMergeCollapsed)
    boolsubsetPopStatsMergeCollapsed[,conditioncol] <- factor(boolsubsetPopStatsMergeCollapsed[,conditioncol], levels=facetorder)
  }
  if (!is.null(facetorder2)) {
    flowWorkspace::pData(gsSub)[,conditioncol2] <- factor(flowWorkspace::pData(gsSub)[,conditioncol2], levels=facetorder2)
    if (conditioncol != "." && is.null(facetorder)) {
      flowWorkspace::pData(gsSub)[,conditioncol] <- factor(flowWorkspace::pData(gsSub)[,conditioncol])
    }
    boolsubsetPopStatsMergeCollapsed <- as.data.frame(boolsubsetPopStatsMergeCollapsed)
    boolsubsetPopStatsMergeCollapsed[,conditioncol2] <- factor(boolsubsetPopStatsMergeCollapsed[,conditioncol2], levels=facetorder2)
  }
  
  # Rename conditioncol2 column in case it contains characters like spaces, which mess up formulas
  colnames(pData(gsSub))[which(colnames(pData(gsSub)) == conditioncol2)] <- "conditioncol2"
  colnames(boolsubsetPopStatsMergeCollapsed)[which(colnames(boolsubsetPopStatsMergeCollapsed) == conditioncol2)] <- "conditioncol2"
  conditioncol2 <- "conditioncol2"
  
  flowplot <- ggcyto::ggcyto(gsSub, ggplot2::aes_string(x=xaxis, y=yaxis, alpha=0.5), subset=parentsubset_ForPlot_name) +
    ggplot2::geom_hex(bins = geom_hex_bins) +
    ggcyto::labs_cyto("marker") +
    ggplot2::facet_grid(stats::as.formula(paste(conditioncol, "~", conditioncol2))) +
    ggcyto::geom_overlay(boolsubsetName, col="red", size=overlayDotSize, alpha=1) +
    ggplot2::geom_text(data=boolsubsetPopStatsMergeCollapsed, ggplot2::aes_string(x=get("geomTextX"), y=get("geomTextY"), label="Percent"),
                       colour="black", parse=FALSE, inherit.aes=FALSE,
                       size=if(is.null(percentage_fontsize)) { max(1, themeBaseSize-13) } else { percentage_fontsize }) +
    ggplot2::scale_alpha(guide = 'none') +
    ggplot2::theme_set(ggplot2::theme_gray(base_size = themeBaseSize))
  # ggplot2::theme(plot.title=ggplot2::element_text(vjust=-0.8, hjust=0.5, size=19),
  #                plot.subtitle=ggplot2::element_text(size=12),
  #                axis.text=ggplot2::element_text(size=14),
  #                axis.title=ggplot2::element_text(size=18),
  #                strip.text=ggplot2::element_text(size=16),
  #                legend.title=ggplot2::element_text(size=15),
  #                legend.text=ggplot2::element_text(size=12))
  if(!is.null(xlims) && is.null(ylims)) {
    flowplot <- flowplot + ggplot2::coord_cartesian(xlim = xlims) 
  } else if(!is.null(ylims) && is.null(xlims)) {
    flowplot <- flowplot +  ggplot2::coord_cartesian(ylim = ylims) 
  } else if(!is.null(xlims) && !is.null(ylims)) {
    flowplot <- flowplot +  ggplot2::coord_cartesian(xlim = xlims, ylim = ylims) 
  }
  if(stripLegendGridAxesTitle) {
    flowplot <- flowplot + 
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"),
                     legend.position="none", axis.text=ggplot2::element_blank(), axis.ticks=ggplot2::element_blank()) +
      ggplot2::labs(title="")
    # TODO: subset label seems to be plotted by default, not sure how to override it so title is not allocated space
  } else {
    flowplot <- flowplot + 
      ggplot2::labs(title=flowtitle, subtitle=subtitle1)
  }
  if(!is.null(axestitle_fontsize)) {
    flowplot <- flowplot + ggplot2::theme(axis.title=element_text(size=axestitle_fontsize))
  }
  if(!is.null(font)) {
    flowplot <- flowplot + ggplot2::theme(text = element_text(family=font))
  }
  
  width <- if (is.null(width)) { if (conditioncol2 == ".") { 5 } else { 9 } } else { width }
  
  if (!is.null(outdir)) {
    # Simplify subset for file name
    subsetsmpl <- strsplit(boolsubset, split="&")[[1]]
    # What we're selecting positively FOR
    possubset <- subsetsmpl[grep("!", subsetsmpl, invert=TRUE)]
    # Rewrite as one string
    possubset <- paste(lapply(possubset, function(x) {splt <- strsplit(x, "/")[[1]]; splt[length(splt)][[1]]}), collapse="")
    ext <- pngORsvg # if(pngORsvg == "png") {"png" } else { "svg" }
    ggplot2::ggsave(filename=paste(c("FlowPlot_", individualsCol, "_", individual, "_", parentsubset, "_", exp, "_", possubset, ".", ext), collapse=""),
                    plot=flowplot, path=outdir, device=if(pngORsvg == "png") {"png" } else { grDevices::svg() }, width=width, height=8, units="in")
  } else {
    flowplot
  }
}

#' COMPASS Subset Comparisons
#'
#' Obtain background-corrected proportions for COMPASS subsets. Perform Wilcoxon rank sum test comparing mean for one group to mean for another group.
#' Defaults to only subsets above threshold 0.01.
#' Saves a csv with the above data to the outdir, if provided.
#' 
#' It is assumed that any given GatingSetList will have the same gates across all GatingSets within.
#'
#' @param compassResultOrPath 
#' @param gsOrGsListOrPath 
#' @param parentSubset 
#' @param threshold (optional) Filters subsets where the average mean_gamma is greater than the threshold
#' @param antigenCol 
#' @param stimAntigen 
#' @param sampleIDCol
#' @param controlAntigen 
#' @param stratifyBy column of pData to stratify data by
#' @param outdir (optional)
#' @param stratifyByValueMinuend (optional) stratifyBy value to be used as selector for minuend (A in A - B)
#' @param stratifyByValueSubtrahend (optional) stratifyBy value to be used as selector for subtrahend (B in A - B)
#' @import flowWorkspace
#' @import stringr
#' @import tidyr
#' @import coin
#' @import data.table
#' @export
#' @examples 
#' \dontrun{
#' ifngResults <- compass.subset.comparisons(compassResultOrPath=myCompassResult,
#'        gsOrGsListOrPath=gs,
#'        parentSubset="4+",
#'        antigenCol="Antigen",
#'        stimAntigen="Peptide Pool 1",
#'        controlAntigen="DMSO",
#'        stratifyBy="Status",
#'        sampleIDCol="PATIENT ID")
#' }
compass.subset.comparisons <- function(compassResultOrPath,
                                       gsOrGsListOrPath,
                                       parentSubset,
                                       threshold=0.01,
                                       antigenCol,
                                       stimAntigen,
                                       controlAntigen,
                                       stratifyBy,
                                       outdir=NULL,
                                       stratifyByValueMinuend=NULL,
                                       stratifyByValueSubtrahend=NULL,
                                       sampleIDCol) {
  compassResult <- if(class(compassResultOrPath) == "character") {
    readRDS(compassResultOrPath)
  } else if(class(compassResultOrPath) == "COMPASSResult") {
    compassResultOrPath
  }
  gs <- if(class(gsOrGsListOrPath) == "GatingSet" || class(gsOrGsListOrPath) == "GatingSetList") {
    gsOrGsListOrPath
  } else {
    try(if(!(class(gsOrGsListOrPath) == "character")) stop("gsOrGsListOrPath must be either a GatingSet, a GatingSetList, or the path to the folder containing one of these objects on disk"))
    # Load the saved GatingSetList or GatingSet:
    loadGsOrGsList <- function (gsOrGsListOrPath) {
      out <- try(flowWorkspace::load_gs(gsOrGsListOrPath))
      if (class(out) == "try-error") {
        cat("Caught an error during flowWorkspace::load_gs, trying flowWorkspace::load_gslist.\n")
        out <- flowWorkspace::load_gslist(gsOrGsListOrPath)
      }
      out
    }
    loadGsOrGsList(gsOrGsListOrPath)
  }
  
  # We're only interested in the subsets where the average mean_gamma is greater than the threshold (default in heatmap and this function is 0.01)
  # Reapply the filter here...
  m <- apply(compassResult$fit$mean_gamma, 2, function(x) {
    mean(x, na.rm = TRUE)
  })
  keep <- m > threshold
  compassSubsetsFiltered <- names(which(keep))
  # Get rid of the subset with 0 positive markers (equal number of ! as there are markers)
  numMarkers <- ncol(compassResult$fit$categories) - 1
  compassSubsetsFiltered <- compassSubsetsFiltered[lengths(regmatches(compassSubsetsFiltered, gregexpr("!", compassSubsetsFiltered))) != numMarkers]
  # Make markers alphabetical
  makeBoolSubsetAlphabetical <- function(subset) {
    and_split <- stringr::str_split(subset, "&")[[1]]
    rawMarkerNames <- unlist(lapply(and_split, function(str) { tmpvec <- stringr::str_split(str, "!")[[1]]; tmpvec[[length(tmpvec)]] }))
    negMarkers <- rawMarkerNames[grep("!", and_split)]
    rawMarkerNamesSorted <- sort(rawMarkerNames)
    # This next line is very specific to the way the gates are named in the RSTR data (with "+"). Unfortunately prevents the function from general use.
    and_split_sorted <- unlist(lapply(rawMarkerNamesSorted, function(str) { if(str %in% negMarkers) { paste0("!", str, "+") } else { paste0(str, "+") } }))
    paste(and_split_sorted, collapse="&")
  }
  compassSubsetsFilteredAlpha <- unlist(lapply(compassSubsetsFiltered, makeBoolSubsetAlphabetical))
  
  # Add all the COMPASS subsets of interest to the GatingSet/GatingSetList
  newNodeAdded <- FALSE
  for(boolSubset in compassSubsetsFilteredAlpha) {
    call <- substitute(flowWorkspace::booleanFilter(v), list(v = as.symbol(boolSubset)))
    g <- eval(call)
    boolSubsetNodeName <- paste0(parentSubset, ":", boolSubset)
    if(boolSubsetNodeName %in% flowWorkspace::getNodes(gs[[1]], path="auto")) {
      # If the gate already exists, keep it
      message(paste0("Gate ", boolSubsetNodeName, " already exists. Keeping old gate..."))
    } else {
      flowWorkspace::add(gs, g, parent = parentSubset, name=boolSubsetNodeName)
      newNodeAdded <- TRUE
    }
  }
  if(newNodeAdded) {
    flowWorkspace::recompute(gs)
  } else {
    message("No new nodes added, not recomputing gates.")
  }
  
  compassPopStats <- lapply(paste0(parentSubset, ":", compassSubsetsFilteredAlpha), function(boolSubsetNodeName) {
    d <- flowWorkspace::getPopStats(gs, subpopulation=boolSubsetNodeName)
    d[,boolSubsetNodeName] <- d$Count / d$ParentCount
    d[, c("name", boolSubsetNodeName), with=FALSE]
  })
  mergedCompassPopStats <- Reduce(function(...) merge(...), compassPopStats)
  flowWorkspace::pData(gs)$name <- rownames(flowWorkspace::pData(gs))
  compassPopStatsMeta <- merge(flowWorkspace::pData(gs), mergedCompassPopStats, by="name")
  
  compassPopStatsMeta_ctrl <- compassPopStatsMeta[compassPopStatsMeta[,antigenCol] == controlAntigen,]
  compassPopStatsMeta_stim <- compassPopStatsMeta[compassPopStatsMeta[,antigenCol] == stimAntigen,]
  # Subtract Control proportions from Stim proportions
  compassPopStatsMetaBgCorr <- merge(compassPopStatsMeta_ctrl, compassPopStatsMeta_stim, by=sampleIDCol, suffixes = c(".Ctrl", ".Stim"))
  for(boolSubsetNodeName in paste0(parentSubset, ":", compassSubsetsFilteredAlpha)) {
    stim_colname <- paste0(boolSubsetNodeName, ".Stim")
    ctrl_colname <- paste0(boolSubsetNodeName, ".Ctrl")
    compassPopStatsMetaBgCorr[, paste0(boolSubsetNodeName, ".BgCorr")] <- unlist(purrr::map2(compassPopStatsMetaBgCorr[,stim_colname], compassPopStatsMetaBgCorr[,ctrl_colname],
                                                                                             function(stimProp, ctrlProp) {
                                                                                               max(stimProp - ctrlProp, 0)
                                                                                             }))
  }
  
  stopifnot(compassPopStatsMetaBgCorr[,paste0(stratifyBy, ".Ctrl")] == compassPopStatsMetaBgCorr[,paste0(stratifyBy, ".Stim")])
  colnames(compassPopStatsMetaBgCorr)[which(colnames(compassPopStatsMetaBgCorr) == paste0(stratifyBy, ".Ctrl"))] <- stratifyBy
  
  # Make the stratifyBy column a factor, sometimes needed
  compassPopStatsMetaBgCorr[, stratifyBy] <- as.factor(compassPopStatsMetaBgCorr[, stratifyBy])
  tests <- sapply(paste0(parentSubset, ":", compassSubsetsFilteredAlpha, ".BgCorr"), function(colname) {
    coin::wilcox_test(data = compassPopStatsMetaBgCorr, get(colname) ~ get(stratifyBy))
  })
  minuend <- if(!is.null(stratifyByValueMinuend)) { stratifyByValueMinuend } else { unique(compassPopStatsMetaBgCorr[,stratifyBy])[[1]] }
  subtrahend <- if(!is.null(stratifyByValueSubtrahend)) { stratifyByValueSubtrahend } else { unique(compassPopStatsMetaBgCorr[,stratifyBy])[[2]] }
  minuendMeanMag <- sapply(paste0(parentSubset, ":", compassSubsetsFilteredAlpha, ".BgCorr"), function(colname) {
    mean(compassPopStatsMetaBgCorr[which(compassPopStatsMetaBgCorr$Status == minuend), colname])
  })
  subtrahendMeanMag <- sapply(paste0(parentSubset, ":", compassSubsetsFilteredAlpha, ".BgCorr"), function(colname) {
    mean(compassPopStatsMetaBgCorr[which(compassPopStatsMetaBgCorr$Status == subtrahend), colname])
  })
  deltaMeanMagnitudes <- minuendMeanMag - subtrahendMeanMag
  
  bgCorrProportionsTestStratified <- data.table::rbindlist(lapply(compassSubsetsFilteredAlpha, function(subset) {
    and_split <- stringr::str_split(subset, "&")[[1]]
    rawMarkerNames <- unlist(lapply(and_split, function(str) { tmpvec <- stringr::str_split(str, "!")[[1]]; tmpvec[[length(tmpvec)]] }))
    subsetRepresentation <- rep("+", length(and_split))
    subsetRepresentation[grep("!", and_split)] <- "-"
    subsetRepresentation <- data.table::as.data.table(t(subsetRepresentation))
    colnames(subsetRepresentation) <- rawMarkerNames
    subsetRepresentation
  }))
  bgCorrProportionsTestStratified <- cbind(bgCorrProportionsTestStratified, p.adjust(sapply(tests, coin::pvalue), method="bonferroni"),
                                           minuendMeanMag, subtrahendMeanMag, deltaMeanMagnitudes)
  colnames(bgCorrProportionsTestStratified) <- c(unlist(lapply(stringr::str_split(compassSubsetsFilteredAlpha[[1]], "&")[[1]], function(str) {
    tmpvec <- stringr::str_split(str, "!")[[1]];
    tmp <- tmpvec[[length(tmpvec)]];
    substr(tmp, 1, nchar(tmp)-1)
  })), "p-value (adj)", paste0("Mean ", minuend), paste0("Mean ", subtrahend), "Diff Mean Bg-Corr Prop")
  bgCorrProportionsTestStratified <- bgCorrProportionsTestStratified[order(bgCorrProportionsTestStratified$`p-value (adj)`),]
  
  # Save output to outdir if given
  if(!is.null(outdir)) {
    file_prefix <- paste0(gsub(" ", "", stimAntigen), "_", parentSubset, "_")
    # Save the wilcox tests as rds
    saveRDS(tests, file=file.path(outdir, paste0(file_prefix, "WilcoxTestsBy", stratifyBy, ".rds")))
    # Save the pvalue table as csv
    write.csv(bgCorrProportionsTestStratified, file = file.path(outdir, paste0(file_prefix, "bgCorrProportionsTestBy", stratifyBy, ".csv")))
    # Save the table of bg-corrected proportions as csv
    write.csv(compassPopStatsMetaBgCorr, file = file.path(outdir, paste0(file_prefix, "bgCorrProportions.csv")))
  }
  
  # Return a list with 3 items: 1) the wilcox test objects, 2) the p-value table, 3) the table of bg-corrected proportions
  list("wilcox" = tests,
       "pValueTable" = bgCorrProportionsTestStratified,
       "bgCorrProportions" = compassPopStatsMetaBgCorr)
}

#' Create Stratified Boxplots of the Background Corrected Proportions for COMPASS subsets
#' 
#' Subsets will be plotted in the order in which they appear in the cr_cats matrix.
#' Only subsets which appear in compass.subset.comparisons.result will be plotted.
#' Intended to be shown with the categories legend from the corresponding COMPASS heatmap displayed below.
#' 
#' This function is hacked together.
#'
#' @param compass.subset.comparisons.result The output from a run of COMPASSHelpers::compass.subset.comparisons(). Contains the p-values and background-corrected proportions for each subset.
#' @param cr_cats a COMPASSResult categories matrix (cr$fit$categories) with rows defining the order in which you want the boxplots to appear. Can be obtained by a call to COMPASSHelpers::orderHeatmapColumnsByCytokinePresence(), for example.
#' @param stratifyBy what metadata column to stratify the boxplots by, e.g. "Status"
#' @param sampleIDCol e.g. "PATIENT ID"
#' @param parentSubset e.g. "4+"
#' @param themeBaseSize
#' @param removeGridAndBg
#' @param showTitle
#' @param showSignificanceBracket
#' @param pvalue_fontsize
#' @param axestitle_fontsize
#' @param axestick_fontsize
#' @param font
#' @param geom_jitter_width
#' @param xaxis_title
#' @param show_points
#' @param padj_method
#' @param ylim e.g. c(0,0.01)
#' @param p_alpha The alpha level to use for choosing which subsets to plot p-values for
#' @param stratifyBy_colors e.g. c("RSTR" = "white", "LTBI" = "black")
#' @param box_fill_alpha
#' @param legend_fontsize
#' @export
#' @import data.table
#' @import flowWorkspace
#' @import ggplot2
#' @import ggsignif
#' @import tidyr
#' @examples 
#' \dontrun{
#' compass.subset.comparisons.result <- compass.subset.comparisons(compassResultOrPath=CD4PP1CompassResult,
#'                                                              gsOrGsListOrPath=gs,
#'                                                              parentSubset="4+",
#'                                                              antigenCol="Antigen",
#'                                                              stimAntigen="Peptide Pool 1",
#'                                                              controlAntigen="DMSO",
#'                                                              stratifyBy="Status",
#'                                                              stratifyByValueMinuend = "LTBI",
#'                                                              stratifyByValueSubtrahend = "RSTR",
#'                                                              sampleIDCol="PATIENT ID")
#' cr_cats <- orderHeatmapColumnsByCytokinePresence(CD4PP1CompassResult, "IFNg",
#'                                       cats_df_cytokine_order_override = c("CD154", "IL2", "TNF", "CD107a", "IL4", "IL17a", "IFNg"))$fit$categories
#' boxplotsCompassSubsetBgCorrPropStratified(compass.subset.comparisons.result,
#'                                           cr_cats,
#'                                           stratifyBy="Status",
#'                                           sampleIDCol,
#'                                           parentSubset,
#'                                           removeGridAndBg=T,
#'                                           ylim=c(0, 0.01),
#'                                           stratifyBy_colors=c("RSTR" = "white", "LTBI" = "black"))
#' }
boxplotsCompassSubsetBgCorrPropStratified <- function(compass.subset.comparisons.result,
                                                      cr_cats,
                                                      stratifyBy,
                                                      sampleIDCol,
                                                      parentSubset,
                                                      themeBaseSize=18,
                                                      removeGridAndBg=FALSE,
                                                      showSignificanceBracket=TRUE,
                                                      pvalue_fontsize=NULL,
                                                      axestitle_fontsize=NULL,
                                                      axestick_fontsize=NULL,
                                                      font=NULL,
                                                      geom_jitter_width=0.05,
                                                      xaxis_title=stratifyBy,
                                                      dot_size=0.75,
                                                      legend_position = c(1,1),
                                                      show_points=F,
                                                      padj_method="bonferroni",
                                                      ylim=NULL,
                                                      p_alpha=0.05,
                                                      stratifyBy_colors=NULL,
                                                      box_fill_alpha=1,
                                                      legend_fontsize=NULL) {
  
  # Put the cr_cats columns in the order they appear in compass.subset.comparisons.result$pValueTable columns
  # This is to recreate the subset names using the same cytokine order as those in the compass.subset.comparisons.result object
  cr_cats_mod <- ifelse(cr_cats[,colnames(compass.subset.comparisons.result$pValueTable)[1:ncol(cr_cats) - 1]]==0, "!", "")
  # The subsets in plotting order:
  # The "+" addition makes this RSTR study specific...
  subsetsInOrder <- apply(cr_cats_mod, 1, function(x) {paste0(parentSubset, ":", paste(paste0(x, paste(colnames(cr_cats_mod), "+", sep="")), collapse="&"), ".BgCorr")})
  subsetsInOrder2Plot <- subsetsInOrder[which(subsetsInOrder %in% names(compass.subset.comparisons.result$wilcox))]
  data <- compass.subset.comparisons.result$bgCorrProportions[, c(sampleIDCol, stratifyBy, subsetsInOrder2Plot)]
  data_long <- gather(data, Subset, BgCorrProp, subsetsInOrder2Plot)
  # Make the data_long Subset column a factor so that it preserves the order desired in subsetsInOrder2Plot
  data_long$Subset <- factor(data_long$Subset, levels=subsetsInOrder2Plot)
  
  p <- ggplot(data_long, aes(x=Status, y=BgCorrProp, fill=Status)) +
    geom_boxplot(outlier.shape=NA, position = position_dodge2(preserve = "total"))
  if(show_points) {
    p <- p + geom_jitter(width = geom_jitter_width, size=dot_size)
  }
  p <- p +
    facet_grid(. ~ Subset) + #, switch="x") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank()) +
    theme(legend.justification = c(1,1), legend.position = legend_position) +
    labs(y="Background corrected proportion") +
    guides(fill=guide_legend(title=NULL)) +
    ggplot2::theme(axis.text=element_text(colour="black"))
  if(!is.null(stratifyBy_colors)) {
    p <- p +
      scale_alpha_manual(values=names(stratifyBy_colors)) +
      scale_fill_manual(values=alpha(stratifyBy_colors, box_fill_alpha))
  }
  if(!is.null(axestick_fontsize)) {
    p <- p + ggplot2::theme(axis.text=element_text(size=axestick_fontsize))
  }
  if(!is.null(axestitle_fontsize)) {
    p <- p + ggplot2::theme(axis.title=element_text(size=axestitle_fontsize))
  }
  if(!is.null(legend_fontsize)) {
    p <- p + ggplot2::theme(legend.text=element_text(size=legend_fontsize))
  }
  if(!is.null(font)) {
    p <- p + ggplot2::theme(text = element_text(family=font))
  }
  if(removeGridAndBg) {
    p <- p + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                            panel.background = ggplot2::element_blank(), axis.line = element_line(colour = "black"))
  }
  if(!is.null(ylim)) {
    p <- p +
      ylim(ylim)
  }
  if(showSignificanceBracket) {
    subsets_p_adj <- p.adjust(sapply(as.vector(subsetsInOrder2Plot), function(subset) { coin::pvalue(compass.subset.comparisons.result$wilcox[[subset]]) }), method=padj_method)
    for(subset in subsetsInOrder2Plot) {
      if(subsets_p_adj[[subset]] < p_alpha) {
        # The following code is necessary to get the brackets to show but don't do anything: Status=data_long[[stratifyBy]][[1]], BgCorrProp=0, group=1), aes(group=group)
        # TODO: deal with outliers....i.e. just instead of ymax do 75% range thing quantile
        #p_text <- if(subsets_p_adj[[subset]] < 0.001) {"p<0.001"} else {paste0("p=", signif(subsets_p_adj[[subset]], digits=3))}
        # Show 3 decimal places, or "p<0.001" if less than 0.001
        p_text <- if(subsets_p_adj[[subset]] < 0.001) {"p<0.001"} else {paste0("p=", format(round(subsets_p_adj[[subset]], 3), nsmall = 3))}
        find_whisker_max <- function(x) { quantile(x, probs = c(0.75)) + IQR(x) }
        whisker_max <- max(aggregate(data_long[which(data_long$Subset==subset),"BgCorrProp"], by = list(data_long[which(data_long$Subset==subset), stratifyBy]), FUN = find_whisker_max)$x)
        y_pos <- whisker_max + if(!is.null(ylim)) {ylim[[2]]/30} else {max(data_long$BgCorrProp)/30}
        if(!is.null(font)) {
          if(!is.null(pvalue_fontsize)) {
            p <- p +
              geom_signif(data=data.frame(Subset=c(subset),
                                          Status=data_long[[stratifyBy]][[1]], BgCorrProp=0, group=1),
                          aes(group=group, family=font),
                          annotations=p_text,
                          y_position=y_pos,
                          xmin=1, xmax=2,
                          tip_length = c(0.005, 0.005),
                          textsize=pvalue_fontsize)
          } else {
            p <- p +
              geom_signif(data=data.frame(Subset=c(subset),
                                          Status=data_long[[stratifyBy]][[1]], BgCorrProp=0, group=1),
                          aes(group=group, family=font),
                          annotations=p_text,
                          y_position=y_pos,
                          xmin=1, xmax=2,
                          tip_length = c(0.005, 0.005))
          }
        } else {
          if(!is.null(pvalue_fontsize)) {
            p <- p +
              geom_signif(data=data.frame(Subset=c(subset),
                                          Status=data_long[[stratifyBy]][[1]], BgCorrProp=0, group=1),
                          aes(group=group),
                          annotations=p_text,
                          y_position=y_pos,
                          xmin=1, xmax=2,
                          tip_length = c(0.005, 0.005),
                          textsize=pvalue_fontsize)
          } else {
            p <- p +
              geom_signif(data=data.frame(Subset=c(subset),
                                          Status=data_long[[stratifyBy]][[1]], BgCorrProp=0, group=1),
                          aes(group=group),
                          annotations=p_text,
                          y_position=y_pos,
                          xmin=1, xmax=2,
                          tip_length = c(0.005, 0.005))
          }
        }
      }
    }
  }
  p
}

#' Adds a boolean gate to the GatingSet or GatingSetList
#' 
#' @param gs Either a GatingSet, a GatingSetList, or path to one of these objects on disk
#' @param booleanSubset The booleanSubset (a combination of existing gates) in string format, e.g. "8+/GMM+&!8+/GAMMADELTA"
#' @param parentGate The gate under which the booleansubset is to be added
#' @param booleanGateName optional. What to call the new gate
#' @import flowWorkspace
#' @export addBooleanGate
addBooleanGate <- function(gs,
                           booleanSubset,
                           parentGate,
                           overrideGate=FALSE,
                           booleanGateName=NULL) {
  try(if(missing(gs) || missing(booleanSubset) || missing(parentGate)) stop("Required arguments missing.") )
  # Check if booleanSubset already exists under parentGate
  booleanSubsetName <- if(is.null(booleanGateName)) {
    gsub("/", ":", paste0(parentGate, ":", booleanSubset))
  } else {
    booleanGateName
  }
  message("Checking for gate ", booleanSubsetName, "\n in the first GatingHierarchy")
  if(booleanSubsetName %in% getNodes(gs[[1]], path="auto")) {
    # If the gate already exists, decide what to do with it.
    if(overrideGate) {
      message(paste0("Gate ", booleanSubsetName, " already exists. Deleting old gate and adding new gate..."))
      call <- substitute(flowWorkspace::booleanFilter(v), list(v = as.symbol(booleanSubset)))
      g <- eval(call)
      flowWorkspace::Rm(booleanSubsetName, gs)
      flowWorkspace::add(gs, g, parent = parentGate, name=booleanSubsetName)
      flowWorkspace::recompute(gs, booleanSubsetName)
    } else {
      message(paste0("Gate ", booleanSubsetName, " already exists. Keeping old gate..."))
    }
  } else {
    # If the gate doesn't exist, add it.
    message(paste0("Adding new gate ", booleanSubsetName, " to GatingSet/GatingSetList..."))
    call <- substitute(flowWorkspace::booleanFilter(v), list(v = as.symbol(booleanSubset)))
    g <- eval(call)
    flowWorkspace::add(gs, g, parent = parentGate, name=booleanSubsetName)
    flowWorkspace::recompute(gs, booleanSubsetName)
  }
}