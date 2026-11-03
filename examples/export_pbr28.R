# R base only; the source RData is fetched separately from the pinned kinfitr URL.
args <- commandArgs(trailingOnly=TRUE)
load(args[1])
dir.create(args[2],recursive=TRUE,showWarnings=FALSE)
summary <- data.frame()
for (i in seq_len(nrow(pbr28))) {
  name <- as.character(pbr28$PET[i])
  folder <- file.path(args[2],name)
  dir.create(folder,showWarnings=FALSE)
  write.csv(as.data.frame(pbr28$tacs[[i]]),file.path(folder,'tacs.csv'),row.names=FALSE)
  write.csv(as.data.frame(pbr28$input[[i]]),file.path(folder,'input.csv'),row.names=FALSE)
  summary <- rbind(summary,data.frame(scan=name,subject=as.character(pbr28$Subjname[i]),
                                    visit=pbr28$PETNo[i],genotype=as.character(pbr28$Genotype[i])))
}
write.csv(summary,file.path(args[2],'scans.csv'),row.names=FALSE)
