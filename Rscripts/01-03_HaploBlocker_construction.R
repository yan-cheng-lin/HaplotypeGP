#### Load data, libraries and functions ####
### Data ###
load("../Data/GenotypicData.RData")

### Library ###
library('HaploBlocker')
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
id.geno <- grep(rownames(geno.sel.gamete), pattern = "GC_KE", value = T) %>%
  substr(start = 1, stop = 9) %>%
  unique()
win <- 20 # Window size setting, can be a vecter, etc. c(5,10,20,50), for multi-window size mode 
tc <- 0.90 # Target coverage setting (0-1), use NULL for not using 
Minsubgroup <- 1 # Min Subgroup setting
MCMB <- 5000 # MCMB setting


#### Data preparation  ####
minMAF <- 0.01
MultiWindowMode <- length(win) > 1 # If window size is a vector, use multi-window size mode


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

pmap <- map[colnames(geno), ]


### Sub-group List ###
subgroup <- list()
for (pop in c("DH_KE", "DH_PE", "GC_KE", "GC_PE")){
  id.pop <- grep(id.gamete, pattern = pop)
  
  if(length(id.pop) > 0){
    subgroup[[pop]] <- id.pop
  } else {}
}

# Put FV2 into second subgroup (should be GC) #
if("FV2" %in% id.geno){
  subgroup[[2]] <- grep(id.gamete, pattern = "FV2")
} else {
  
}

# Change parameter setting, if there is only one subgroup #
if(length(subgroup) == 0){ 
  subgroup <- NULL
  Minsubgroup <- NULL
} else {
  tc <- NULL
}

### Modify format for HaploBlocker ###
geno.chr <- list()
for ( chr in unique(pmap[["chr"]])){
  mr.chr <- rownames(pmap)[pmap$chr == chr]
  gamete.pop <- geno[, mr.chr] 
  geno.chr[[chr]] <- t( gamete.pop[ rownames(gamete.pop) , ] ) # Trans genotype matrix
} 


#### Run block_calculation ####
bc <- list()
for ( chr in unique(pmap$chr) ){
  
  bpmap <- pmap[rownames(geno.chr[[chr]]), "pos"]
  
  bc[[chr]] <- block_calculation(geno.chr[[chr]],
                                 bp_map = bpmap,
                                 window_size = win,
                                 multi_window_mode = MultiWindowMode,
                                 target_coverage = tc,
                                 overlap_remove = FALSE,
                                 subgroups = subgroup,
                                 min_per_subgroup = Minsubgroup,
                                 min_majorblock = MCMB,
                                 verbose = FALSE,
                                 merging_error = 0)
  
  print(paste(chr,"Finished"))
  
}

#### Make Haplotype Matrix and Marker Lists ####
HB <- list()
MarkerList <- list()
for ( chr in 1:10){
  
  ## Make HB Matrix ##
  mat.hb.raw <- block_matrix_construction(bc[[chr]])
  mat.hb <- t(mat.hb.raw[ , seq(1,ncol(mat.hb.raw),by=2)] +  mat.hb.raw[ , seq(2,ncol(mat.hb.raw),by=2)] )
  rownames(mat.hb) <- id.geno
  colnames(mat.hb) <- paste(chr, "_", colnames(mat.hb), sep="")
  HB[[chr]] <- mat.hb
  
  ## List for BCIDX v.s. Markers ##
  mrlist <- sapply(bc[[chr]], function(x){names(x[[7]]$snp)})
  
  bc.nMR <- sapply(mrlist, length)
  nHB <- length(mrlist)
  
  MarkerList[[chr]] <- data.frame(Markers = unlist(mrlist),
                                  BCIDX = paste(str_pad(chr, 2, side = c("left"), pad = "0"),
                                                "HB",
                                                str_pad(rep(1:nHB, times = bc.nMR), 4, side = c("left"), pad = "0"),
                                                sep = "_"))
  
} 
HB <- do.call(cbind, HB)
MarkerList <- do.call(rbind, MarkerList)


#### Information of HaploBlocks ####
nblock <- sapply(bc, length)

hbidx <- paste(str_pad(rep(1:length(nblock), times = nblock), 2, side = c("left"), pad = "0"),
               "HB",
               str_pad(unlist(sapply(nblock, function(x){1:x})), 4, side = c("left"), pad = "0"),
               sep = "_")

bc.info <- data.frame(row.names = hbidx)

for ( t in c("snp", "window", "bp")){
  bc.INFO <- lapply(bc, blocklist_startend, type = t)
  
  bc.INFO <- as.data.frame(do.call(rbind, bc.INFO))
  
  bc.info[[t]] <- bc.INFO$end - bc.INFO$start
  
}

colnames(HB) <- hbidx


#### Output Data ####
save(list = c("HB", "bc", "bc.info", "MarkerList", "id.geno"),
     file = "../Data/HaploBlocker_HMat.RData")

# HB: Haplotype matrix for genomic relationship matrix (GRM) calculation.
# Naming system of haplotype allele: 01_HB_0002 indicates haplotype 2 on chromosome 1.
# bc: Original results of block_calculation from R/HaploBlocker
# bc.info: Length information of the haplotypes
# MarkerList: List of markers in each haplotype
# id.geno: Genotypes in this dataset 

