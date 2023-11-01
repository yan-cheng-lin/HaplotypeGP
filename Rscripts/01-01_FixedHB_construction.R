#### Load data, libraries and functions ####
### Data ###
load("../Data/GenotypicData.RData")

### Library ###
library(zoo)
library(stringr)


### Functions ###
# function for calculating minor allele frequency using gamete genotype matrix
MAF_Gamete <- function(x){
  af0 <- sum(x)/(length(x))
  if(af0 > 0.5 ){
    maf <- 1-af0
  } else {
    maf <- af0
  }
  return(maf)
}

# function for recoding geno into "haplotype alleles"
recode_hapAlleles_f <- function(markers, geno){
  geno_wind <- geno[, markers]
  geno_wind_unique <- unique(geno_wind)
  haps <- apply(geno_wind, 1, paste0, collapse ="")
  return(haps)
}

# function for writing a Markerlist containing the names of markers within sliding windows along each chromosome
FixedLengthMarkerlist_f <- function(map, nSNPs, steps) {
  Markerlist <- NULL
  for (Chr in unique(map$chr)){	
    SNPs_chr <- rownames(map)[map$chr == Chr]
    pos_chr <- map$pos[map$chr == Chr]
    names(pos_chr) <- SNPs_chr
    Markerlist <- rbind(Markerlist, rollapply(SNPs_chr, width = nSNPs, FUN = function(x){x}, by = steps))
  }
  return(Markerlist)
}

# function for generating Infolist of the corresponding blocks
makeInfolist_f <- function(x, map) {
  map_temp <- map[which(rownames(map) %in% x),]
  Chr <- unique(map_temp$chr)
  if(length(Chr) != 1){
    print(x)
    print("block goes across chromosomes!")
  }
  pos <- map_temp$pos
  names(pos) <- rownames(map_temp)
  info <- c(Chr, pos[x[1]], pos[x[length(x)]], pos[x[length(x)]] - pos[x[1]] + 1, length(x))
  names(info) <- c("Chr", "Pos_Start_bp", "Pos_End_bp", "Size_bp", "n_Markers")
  return(info)
}

# function for recoding vector of haplotype alleles (1,2,3,4,5,...) into binary format with one column per haplotype allele
# for windows with only 2 alleles, remove one of them (the second one would provide no additional information)

### Don't remove the alternative allele in this cases ###
# only keep alleles passing a certain MAF threshold
binary_coding_f <- function(x, MAF) {
  haps <- unique(x)
  if (length(haps) > 1){
    binary <- matrix(x, nrow = length(x), ncol = length(haps))
    binary <- t(apply(binary, 1, helper_binary_coding_f, haps = haps))
    colnames(binary) <- haps
    # if(ncol(binary) == 2){binary <- matrix(binary[,1], ncol = 1)} #
    maf <- apply(binary, 2, MAF_f)
    binary <- matrix(binary[ ,which(maf > MAF)], ncol = length(which(maf > MAF)))
    return(binary)
  }
}
helper_binary_coding_f <- function(y, haps){
  y[y != haps] <- 0
  y[y == haps] <- 1	# 2 -> 1
  return(y)
} ### Setting here is change into gamete type, so carry one a allele would be cod in 1

# function for MAF calculation
MAF_f <- function(x) {
  y <- sum(length(which(x==0))*2+length(which(x==1)))/(2*(length(x)-length(which(is.na(x)))))
  y <- ifelse(y > 0.5,return(1-y),return(y))
  return(y)
}





#### Settings ####
# Samples for haplotype construction #
id.geno <- c("DH_KE0001", "DH_KE0002", "DH_KE0003",
             "DH_PE0001", "DH_PE0002", "DH_PE0003") 
win <- 20 # Window size (nSNP)
minCount <- 0 # Frequency filtering for haplotype allele, use 0 to keep all haplotype allele.



#### Data preparation  ####
nSNPs <- win
steps <- win
minMAF <- 0.01


# Minor allele frequency filtering for SNP (population based) #
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

#### Run Haplotype construction ####
HBlist <- list()
WinMarker <- list()
Info <- list()
HBGlist <- list()
for (CHR in unique(pmap$chr)){
  
  ## Subset each CHR ##
  map.chr <- pmap[which(pmap$chr == CHR), ]
  geno.chr <- geno[ , rownames(map.chr)]

  ### define windows ###
  # write Markerlist according to the specified windows
  Marker_List <- FixedLengthMarkerlist_f(map = map.chr, nSNPs = nSNPs, steps = steps)
  rownames(Marker_List) <- paste("wind", str_pad(seq(1, nrow(Marker_List), 1), width = 6, pad = 0 ), sep="")

  # define Info_List for these windows (position, size, ...)
  Info_List <- t(apply(Marker_List, 1, makeInfolist_f, map = map.chr))

  # define haplotypes per window
  haplo <- apply(Marker_List, 1, recode_hapAlleles_f, geno = geno.chr)
  colnames(haplo) <- rownames(Marker_List)
  rownames(haplo) <- rownames(geno.chr)

  # code binary
  # set MAF threshold for binary coding according to minCount (minimum number of counts of a haplotype allele)
  thresh_maf <- minCount / nrow(haplo) 
  haplo_binary_list <- apply(haplo, 2, binary_coding_f, MAF = thresh_maf)
  
  # Convert to matrix #
  haplo_binary_mat.gamete <- do.call(cbind, haplo_binary_list)
  class(haplo_binary_mat.gamete) <- "numeric"
  rownames(haplo_binary_mat.gamete) <- rownames(haplo)
  
  ## Convert gamete matrix into diploid matrix ##
  haplo_binary_mat <- haplo_binary_mat.gamete[ seq(1, nrow(haplo_binary_mat.gamete), by = 2) , ] +
    haplo_binary_mat.gamete[ seq(2, nrow(haplo_binary_mat.gamete), by = 2) , ]
  rownames(haplo_binary_mat) <- id.geno
  
  ## Rename Windows: Haplotype ID: "wind_CHROMOSOENUMBER_WINDOWNUMBER_HAPLOTYPENUMBERWITHINWINDOW"
  hbnumber <- c()
  for ( i in sapply(haplo_binary_list, ncol) ){
    hbnumber <- c(hbnumber, 1:i)
  }
  
  colnames(haplo_binary_mat) <- paste("Win", 
                                      str_pad(CHR, width = 2, side = "left", pad = 0),
                                      rep(names(haplo_binary_list), times = sapply(haplo_binary_list, ncol)),
                                      str_pad(hbnumber, width = 3, side = "left", pad = 0),
                                      sep = "_")
  
  colnames(haplo_binary_mat.gamete) <- colnames(haplo_binary_mat) # Also save the data from haplo_binary_mat.gamete
  
  rownames(Marker_List) <- paste("Win",
                                 str_pad(CHR, width = 2, side = "left", pad = 0),
                                 rownames(Info_List),
                                 sep = "_")
  
  rownames(Info_List) <- paste("Win",
                               str_pad(CHR, width = 2, side = "left", pad = 0),
                               rownames(Info_List),
                               sep = "_")
  
  ### Collect data ###
  HBlist[[CHR]] <- haplo_binary_mat
  WinMarker[[CHR]] <- Marker_List
  Info[[CHR]] <- Info_List
  HBGlist[[CHR]] <- haplo_binary_mat.gamete
  
  print(paste("Chromosome", CHR, "finished"))
}

HB <- do.call(cbind, HBlist)
WinMarker <- do.call(rbind, WinMarker)
Info <- do.call(rbind, Info)
haplo_binary_mat.gamete <- do.call(cbind, HBGlist)

class(HB) <- "numeric"

#### Output ####

save(list = c("HB", "haplo_binary_mat.gamete", "WinMarker", "Info", "id.geno", ),
     file = "../Data/FixedHB_HMat.RData")

# HB: Haplotype matrix for genomic relationship matrix (GRM) calculation.
# Naming system of haplotype allele: Win_01_wind000002_003 indicates third allele of haplotype block 2 on chromosome 1.
# haplo_binary_mat.gamete: Similar to "HB", but in encoded by gametes
# WinMarker: Markers in the haplotype blocks (windows)
# Info: Length information of the haplotype blocks (windows)
# id.geno: Genotypes in this dataset 
