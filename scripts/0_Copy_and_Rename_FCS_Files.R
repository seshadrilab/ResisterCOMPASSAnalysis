# This script copies the downloaded ImmPort FCS files to the appropriate folders for data analysis and renames them to the original FCS file names. This is necessary in order to perform the rest of the downstream analyses.
# Before you run this script, download the flow cytometry result files for study SDY1385 from ImmPort (https://www.immport.org/home) and place the downloaded `Flow_cytometry_result` folder into the `ResisterCOMPASSAnalysis > data` subfolder (unzipped).
# As of 20190426, these FCS files are on ImmPort under `Browse Shared Data > SDY1385 > ResultFiles > Flow_cytometry_result (375 files)` https://browser.immport.org/browser?path=SDY1385%2FResultFiles.

library(here)
projectDir <- here::here()

file_map <- read.table(file.path(projectDir, "data/ImmPort_FCS_FileMapping.tsv"), sep = "\t", header = T, colClasses = "character")
# We don't need the compensation files for the data analysis, so exclude those entries.
# (If you want to look at them anyway, you can comment out the next line. Make sure you download the `Flow_cytometry_compensation_or_control` folder from ImmPort and copy the contents into the unzipped `Flow_cytometry_result` folder so this simple script can find them)
file_map <- file_map[grep("Compensation", file_map$Original_Name, invert = T),]

ImmPort_dl_path <- file.path(projectDir, "data/Flow_cytometry_result")
subfolder_paths <- unique(file_map$New_Folder_Path)
for(subfolder_path in subfolder_paths) {
  subfolder_path_abs <- file.path(projectDir, "data", subfolder_path)
  message(sprintf("Copying and renaming files to %s", subfolder_path_abs))
  if(!dir.exists(subfolder_path_abs)) {
    dir.create(subfolder_path_abs, recursive = T)
  }
  curRows <- subset(file_map, New_Folder_Path == subfolder_path)
  file.copy(file.path(ImmPort_dl_path, curRows$ImmPort_Name), subfolder_path_abs)
  file.rename(file.path(subfolder_path_abs, curRows$ImmPort_Name), file.path(subfolder_path_abs, curRows$Original_Name))
}
