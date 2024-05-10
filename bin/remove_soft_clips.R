#!/usr/bin/env Rscript

library(tidyverse)

#Since I am only looking for soft clips at end and begining of read, this should work. 
# Change the parent directory to the place where your taylor output bam files are
# parent_directory <- "/Volumes/lizso_backup_drive/GSU_CL_VOCBeta_Swift/April2024/swift_out/"

args <- commandArgs(trailingOnly=TRUE)
parent_directory <- args[1]

P01sam <- read_tsv(paste0(parent_directory,"to_clip.sam"), col_names=FALSE, comment="@")
n_soft <- str_count(P01sam$X6, "S")
S_start <- as.data.frame(str_split_fixed(P01sam$X6, "S", 2))
S_start$start_clip <- ifelse(S_start$V2 != "", S_start$V1, 0)
S_end <- data.frame(ends=str_sub(P01sam$X6, start=(str_length(P01sam$X6) - 2), end=str_length(P01sam$X6)))
#types <- S_end %>% group_by(ends) %>% summarize(n=n())
S_end$stop_clip <- ifelse(str_detect(S_end$ends, "S") == TRUE, str_sub(S_end$ends, 1,2), 0)
#examined, only leading char are X or =, so just gsub
S_end$stop_clip <- gsub("X", "", S_end$stop_clip)
S_end$stop_clip <- gsub("=", "", S_end$stop_clip)
S_end$stop_clip <- gsub("D", "", S_end$stop_clip)
S_end$stop_clip <- gsub("I", "", S_end$stop_clip)
S_end$stop_clip <- gsub("M", "", S_end$stop_clip)
S_end$stop_clip <- gsub("N", "", S_end$stop_clip)
S_end$stop_clip <- gsub("P", "", S_end$stop_clip)

S_start$start_clip <- as.numeric(S_start$start_clip)
S_end$stop_clip <- as.numeric(S_end$stop_clip)

S_end[is.na(S_end$stop_clip),]

A <- str_sub(P01sam$X10, S_start$start_clip + 1, str_length(P01sam$X10))
B <- str_sub(A, 1, str_length(A) - S_end$stop_clip)

C <- str_sub(P01sam$X11, S_start$start_clip + 1, str_length(P01sam$X11))
D <- str_sub(C, 1, str_length(C) - S_end$stop_clip)
#str_length(C)
#str_length(D)

new_sam <- P01sam
new_sam$X10 <- B
new_sam$X11 <- D
new_sam$X6 <- gsub("S", "H", new_sam$X6)

new_sam <- new_sam[,1:11]

write.table(new_sam, paste0(parent_directory,"hard_clipped.sam"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

