#### Load data, libraries and functions ####
### Data ###
load("../Data/GenotypicData.RData")

### Library ###
library("magrittr")
library("stringr")

### Functions ###
MAF_Gamete <- function(x){
  af0 <- sum(x)/(length(x))
  if(af0 > 0.5 ){
    maf <- 1-af0
  } else {
    maf <- af0
  }
  return(maf)
}

#### Settings ####
# Samples for haplotype construction #
id.geno <- grep(rownames(geno.sel.gamete), pattern = "DH_KE", value = T) %>%
  substr(start = 1, stop = 9) %>%
  unique()
HVMethod <- "GAB" # Block construction methods: GAB, GAM or SPI


#### Data preparation  ####
minMAF <- 0.01
if (HVMethod == "GAB"){
  HVList <- "GABRIELblocks"
} else if (HVMethod == "GAM") {
  HVList <- "4GAMblocks"
} else if (HVMethod == "SPI") {
  HVList <- "SPINEblocks"
} else {
  print("Methods Wrong")
}


### MAF filtering for SNP (population based) ###
MAFfilter <- c()
for(pop in c("DH_KE", "DH_PE", "GC_KE", "GC_PE")){
  
  ids <- grep(id.geno, pattern = pop, value = T)
  if(length(ids) == 0){
  } else{
    gamete <- paste0(rep(ids, each = 2),
                     rep(c("_1", "_2"), times = length(ids)))
    MAFfilter <- union(MAFfilter, which(apply(geno.sel.gamete[gamete,], 2, MAF_Gamete) > minMAF))
  }
  
}

id.gamete <- paste0(rep(id.geno, each = 2),
                    rep(c("_1", "_2"), times = length(id.geno)))

geno <- geno.sel.gamete[id.gamete, sort(MAFfilter)]

samples <- data.frame(FamilyID = 0,
                      SampleID = rownames(geno))
mr.geno <- colnames(geno)


### Write out lists of samples and markers for HaploView construction ###
write.table(samples, file = "../Data/Sample.list",
            sep = "\t", row.names = F, quote = F, col.names = F)
write.table(mr.geno, file = "../Data/Marker.list",
            sep = "\t", row.names = F, quote = F, col.names = F)

#### Run HaploView via Run_HaploView.sh ####
system(paste("sh 01-02_Run_HaploView.sh", HVMethod))

#### Load HaploView Results ####
### Marker Info ###
info.temp <- list()
for ( chr in 1:10){
  f <- paste("../Data/geno_forHaploView.chr-", chr, ".info", sep = "")
  info.temp[[chr]] <- read.table(f, sep = "\t", col.names = c("Marker", "POS"))
}

info <- do.call(rbind, info.temp)
info$MRIDX <- paste(str_pad(rep(1:10, times = sapply(info.temp, nrow)), width = 2, side = "left", pad = "0"),
                    str_pad(unlist(sapply(info.temp, function(x){ rownames(x) })), width = 6, side = "left", pad = "0"),
                    sep = "_")

### Block Info ###


dat_bc <- read.table(paste0("../Data/HVBlock.", HVList, ".list"),
                     col.names = c("CHR", "BLOCK", "NMR"))
dat_bc$MRIDX <- paste(str_pad(dat_bc$CHR, width = 2, side = "left", pad = "0"),
                      str_pad(dat_bc$NMR, width = 6, side = "left", pad = "0"), sep = "_")

dat_bc$BCIDX <- paste(str_pad(dat_bc$CHR, width = 2, side = "left", pad = "0"),
                      str_pad(dat_bc$BLOCK, width = 5, side = "left", pad = "0"), sep = "_")

dat <- merge(dat_bc, info) 

Info <- dat[, c("Marker", "CHR", "POS", "BCIDX")]


#### Make Haplotype genotype matrix ####
geno.gamete <- geno

## Transform to haplotype matrix ##
geno.bc.list <- list()
geno.bc.list <- tapply(dat$Marker, 
                       dat$BCIDX, 
                       function(marker){
                         haplotypes <- apply(geno.gamete[ , marker], 1, paste, collapse = "" )
                         tags <- unique(haplotypes)
                         mat <- matrix(nrow = length(haplotypes), ncol = length(tags),
                                       dimnames = list(names(haplotypes), tags))
                         for ( i in tags){
                           mat[, i] <- as.numeric( haplotypes == i )
                         }
                         colnames(mat) <- str_pad(1:ncol(mat), width = 3, pad = "0")
                         return(mat)
                       })

geno.bc.mat <- do.call(cbind, geno.bc.list)

bc.name <- gsub(names(geno.bc.list), pattern = "_", replacement = "_wind0")
colnames(geno.bc.mat) <- paste("Win",
                               (rep(bc.name, times = sapply(geno.bc.list, ncol))),
                               colnames(geno.bc.mat),
                               sep = "_")
## Modify col, rownames ##
haplo_binary_mat.gamete <- geno.bc.mat # Gamete haplotype genotype matrix
HB <- geno.bc.mat[ seq(1, nrow(geno.bc.mat), 2) , ] +
  geno.bc.mat[ seq(2, nrow(geno.bc.mat), 2) , ]
rownames(HB) <- gsub(rownames(HB), pattern = "_1$", replacement = "")


#### Output Data ####
save(list = c("HB", "haplo_binary_mat.gamete", "Info", "id.geno"),
     file = "../Data/HaploView_HMat.RData")


# HB: Haplotype matrix for genomic relationship matrix (GRM) calculation.
# Naming system of haplotype allele: Win_01_wind000002_003 indicates third allele of haplotype block 2 on chromosome 1.
# haplo_binary_mat.gamete: Similar to "HB", but in encoded by gametes
# Info: Information of maker name (Marker), physical position (CHR, POS), haplotype block (BCIDX)
# id.geno: Genotypes in this dataset 


