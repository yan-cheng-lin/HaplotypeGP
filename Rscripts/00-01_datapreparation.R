#### Primary Setting ####
### Working Directory ###
setwd("/home/ge47yov/GP/Analysis")

### Libraries ###
library("tidyverse")

### Functions ###

##### Genotypic Data #####
### Load Data ###
GenoRDATA <- "Raw_Data/GenotypicData/geno_map_F2.GC.S0.DH.BL.merged.RData"
PhenoRDATA <- "Raw_Data/PhenotypicData/pheno_Klimafit.MAZE_BLUEs.acrossEnv_EIN.ROG_2017.2018.txt"
load(GenoRDATA)
pheno <- read.table(PhenoRDATA, header = T)

id.geno <- rownames(geno)
id.pheno <- pheno[["Genotype"]]

## DH: DH lines; S0: Original landrace; S1 and S2: GC lines; EL: elite parents
table(substr(id.geno, start = 1,stop = 5)) 
table(substr(id.pheno, start = 1,stop = 5)) 

### Save Map ###
save(map, file = "Geno/map.RData")

#### Split and Transformate the GenoMatrix ####
# Split data genotypic data into three parts, DH, GC and REST.
# Each dataset has two format (coding):
# .geno coded as 0,1,2
# .geno.*.gamete coded as 0,1 for gametes *_1 and *_2

## DH Lines ##
keep <- intersect(grep(id.geno, pattern = "DH", value = T), id.pheno)
geno.DH <- geno[ keep, ]

geno.DH.gamete <- matrix( nrow = nrow(geno.DH) * 2,
                          ncol = ncol(geno.DH), 
                          dimnames = list( paste(rep(rownames(geno.DH), each = 2), rep(c(1,2), times = nrow(geno.DH)), sep = "_"),
                                           colnames(geno.DH)))
geno.DH.gamete[ seq(1, nrow(geno.DH.gamete), by = 2) , ] <- (geno.DH / 2)
geno.DH.gamete[ seq(2, nrow(geno.DH.gamete), by = 2) , ] <- (geno.DH / 2)

save(list = c("geno.DH", "map"), file = "Geno/geno_map_DH.RData")
save(list = c("geno.DH.gamete", "map"), file = "Geno/geno_map_DH_gamete.RData")

## GC Lines ##
keep <- intersect(grep(id.geno, pattern = "S1|S2", value = T),
                  c(paste(id.pheno, "_1", sep = ""), paste(id.pheno, "_2", sep = "")))
geno.GC.gamete <- (geno[ keep, ] / 2)

geno.GC <- geno.GC.gamete[ seq(1, nrow(geno.GC.gamete), by = 2) , ] + geno.GC.gamete[ seq(2, nrow(geno.GC.gamete), by = 2) , ]
rownames(geno.GC) <- substr( rownames(geno.GC.gamete[ seq(1, nrow(geno.GC.gamete), by = 2) , ]), start = 1, stop = 9)
save(list = c("geno.GC.gamete", "map"), file = "Geno/geno_map_GC_gamete.RData")
save(list = c("geno.GC", "map"), file = "Geno/geno_map_GC.RData")

## Rest Lines ##
geno.s0.gamete <- (geno[grep(id.geno, pattern = "S0"),] / 2)
geno.inbred <- geno[grep(id.geno, pattern = "EL|F2"),]

geno.s0 <- geno.s0.gamete[ seq(1, nrow(geno.s0.gamete), by = 2) , ] + geno.s0.gamete[ seq(2, nrow(geno.s0.gamete), by = 2) , ]
rownames(geno.s0) <- substr( rownames(geno.s0.gamete[ seq(1, nrow(geno.s0.gamete), by = 2) , ]), start = 1, stop = 9)

geno.inbred.gamete <- matrix( nrow = nrow(geno.inbred) * 2,
                              ncol = ncol(geno.inbred), 
                              dimnames = list( paste(rep(rownames(geno.inbred), each = 2), rep(c(1,2), times = nrow(geno.inbred)), sep = "_"),
                                               colnames(geno.inbred)))
geno.inbred.gamete[ seq(1, nrow(geno.inbred.gamete), by = 2) , ] <- (geno.inbred / 2)
geno.inbred.gamete[ seq(2, nrow(geno.inbred.gamete), by = 2) , ] <- (geno.inbred / 2)

geno.REST <- rbind(geno.s0, geno.inbred)
geno.REST.gamete <- rbind(geno.s0.gamete, geno.inbred.gamete)

save(list = c("geno.REST","map"),file = "Geno/geno_map_F2.S0.EL.RData")
save(list = c("geno.REST.gamete","map"), file = "Geno/geno_map_F2.S0.EL_gamete.RData")

### Overall Geno Matrix ###
load("Geno/geno_map_GC.RData")
load("Geno/geno_map_DH.RData")

geno <- rbind(geno.DH, geno.GC)
key <- apply(geno, 2, function(x){length(unique(x))}) > 1
table(key)

#### Heritability ####
dh.h2 <- read.table("heritability_2017and2018_HQDataonly.csv", sep = ";", header = T)
gc.h2 <- read.table("LP_heritability_2017and2018.csv", sep = ";", header = T)

dh.h2$Pop <- paste("DH", dh.h2$Pop, sep = "_")
H2 <- rbind(dh.h2[ , colnames(gc.h2) ], gc.h2)

H2$Idx <- gsub(H2$Pop, pattern = "S2", replacement = "GC")

save(H2, file = "Heritability.RData")


#### Write out data ####
### Pheno ###
traits <- c( "EV_V6", "PH_V6", "PH_final", "DtSILK", "RL_R6")
pheno$Genotype <- gsub(pheno$Genotype, pattern = "S1_PE", replacement = "GC_PE") %>%
  gsub(pattern = "S1_KE", replacement = "GC_KE")

pheno_out <- 
pheno[, c("Genotype", traits)] %>%
  subset(grepl(Genotype, pattern = "DH_KE|DH_PE|GC_KE|GC_PE"))
write.table(pheno_out, "/home/ge47yov/GP/Analysis/Publication/Data/PhenotypicData.txt",
            row.names = F, quote = F, sep = "\t")

### Heritability ###
rownames(H2) <- H2$Idx
write.table(H2[ c("DH_KE", "DH_PE", "GC_KE", "GC_PE"), traits],
            "/home/ge47yov/GP/Analysis/Publication/Data/Heritability.txt",
            quote = F, sep = "\t")

### Geno ###
load("Geno/geno_map_GC_gamete.RData")
load("Geno/geno_map_DH_gamete.RData")
load("Geno/geno_map_F2.S0.EL_gamete.RData")

geno.sel.gamete <- rbind(
  geno.DH.gamete[grep(rownames(geno.DH.gamete), pattern = "DH_KE|DH_PE"), ],
  geno.GC.gamete[grep(rownames(geno.GC.gamete), pattern = "S1_KE|S1_PE"), ],
  geno.REST.gamete[c("F2_1", "F2_2"), ]
)

rownames(geno.sel.gamete) <- gsub(rownames(geno.sel.gamete),
                                  pattern = "S1", replacement = "GC")
rownames(geno.sel.gamete) <- gsub(rownames(geno.sel.gamete),
                                  pattern = "F2", replacement = "FV2")

save(list = c("geno.sel.gamete", "map"),
     file = "/home/ge47yov/GP/Analysis/Publication/Data/GenotypicData.RData")

