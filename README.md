# COMPASS analysis of M. tuberculosis Resisters and Non-Resisters from Uganda

This repository contains code to perform COMPASS analysis on Intracellular Cytokine Staining data from the Ugandan Resister and Non-Resister blood samples.  
See the accompanying paper published in Nature Medicine in May 2019, titled [IFN-γ-independent immune markers of Mycobacterium tuberculosis exposure](https://www.nature.com/articles/s41591-019-0441-3)

## Running the analysis

1. First clone or download this repository. The code depends on the current directory structure.  
2. Next, install the R packages listed under the Dependencies section.  
3. Download the [Resister cohort flow data from ImmPort study accession SDY1385](https://www.immport.org/shared/study/SDY1385), place it into the `data` subfolder, and unzip it. You will want to download the folder `Browse Shared Data > SDY1385 > ResultFiles > Flow_cytometry_result (375 files)`. After running the `0_Copy_and_Rename_FCS_Files.R` script in step 4, the `data` subfolder will have the following structure:  

```
├── data  
    ├── 20170518_HiRisk_VisitA_Only_1.txt (This file maps each sample to disease status)  
    ├── ImmPort_FCS_FileMapping.tsv (This file is used to move and rename downloaded FCS files)
    ├── flowJoWorkspace_sampleID_FileMapping.tsv (Maps FCS files to FlowJo workspace sampleID. Needed due to modification of FCS keywords by Immport.)
    ├── Flow_cytometry_result (Folder containing FCS files downloaded from Immport)
    ├── NonTBAgs  
    |   ├── 20180207_OMIP14_Batch1  
    |   |   ├── 20180209_NonTB_Antgns_OMIP14.xml  
    |   |   └── 20180209_RSTR_NonTB_Atgns_Omip14_FCS (this folder will contain FCS files)  
    |   └── 20180214_OMIP14_Batch2  
    |       ├── 20180216_NonTB_Antgns_OMIP14.xml  
    |       └── 20180216_RSTR_NonTB_Atgns_Omip14_FCS (this folder will contain FCS files)  
    └── TBAgs  
        ├── 20170605_RSTR_OMIP14_ICS_Batch1  
        |   ├── 20170607_RSTR_ICS_Batch1.xml  
        |   └── 20170607_RSTR_OMIP14_ICS (this folder contains FCS files)  
        └── 20170612_RSTR_OMIP14_ICS_Batch2  
            ├── 20170614_RSTR_ICS_Batch2.xml  
            └── 20170614_RSTR_OMIP14_ICS_Batch2 (this folder contains FCS files)  
```

4. Create an R project in this top-level folder. Open it in RStudio and run the following scripts (see `scripts` subfolder) in order:

```
0_Copy_and_Rename_FCS_Files.R  
1_QC_SetupGSList_TBAgs.R  
2_RunCompass_TBAgs.R  
3_PostCompassPlots_TBAgs.R  
4_QC_SetupGSList_NonTBAgs.R  
5_RunCompass_NonTBAgs.R  
6_PostCompassPlots_NonTBAgs.R  
```

The output will get placed into the `out` subfolder.

## Dependencies

This work depends on the following R packages:

```
BH # required by flowWorkspace
RcppArmadillo
coin
cowplot
data.table
extrafont
ggplot2
ggsignif # for significance bars
here # for path management
plyr
stringr
survival # required by coin
svglite # for saving plots
tidyr
tidyverse
```

You can use bioconductor or `devtools::install_github()` to install the relevant flow cytometry packages:

```
COMPASS
CytoML
flowCore
flowWorkspace
ggcyto
ncdfFlow
```

[openCyto](https://bioconductor.org/packages/release/bioc/html/openCyto.html) should install flowWorkspace, flowCore, and ncdfFlow.  
[ggcyto](https://bioconductor.org/packages/release/bioc/html/ggcyto.html) for bivariate flow dot plots.  
[COMPASS](https://bioconductor.org/packages/release/bioc/html/COMPASS.html)  
[CytoML](https://bioconductor.org/packages/release/bioc/html/CytoML.html)

Quick install:
```
# Note: on ubuntu 20.04, you may need to install libfontconfig1-dev first
install.packages(c("BH", "RcppArmadillo", "coin", "cowplot", "data.table", "extrafont",
  "ggplot2", "ggsignif", "here", "plyr", "stringr", "survival",
  "svglite", "tidyr", "tidyverse))
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("openCyto")
BiocManager::install("CytoML")
BiocManager::install("ggcyto")
BiocManager::install("COMPASS")
```
