#### Library ####
library("HaploBlocker")
library("magrittr")
library("dplyr")

#### Settings ####
HBRdata <- "../Data/HaploBlocker_HMat.RData" # RData for the haplotype matrix from HaploBlocker
type <- "snp" # Length of the haplotype (snp, cM or gene)
s <- 0.5 # Scaling factor for the weight

#### Functions ####
VanRaden_kinship_f_Weight <- function(W,
                                      BlockInfo,
                                      MRList,
                                      Length,
                                      Scale,
                                      n.rm = 2) {
  #### Make Weighting Vector ####
  rownames(MRList) <- rownames(MRList$Markers)
  
  if ( Length == "snp" ){
    w.vector <- BlockInfo$snp
    
  } else if ( Length == "cM" ) {
    MRList$pos_cM <- map_chr_phy_gen[ as.character(MRList$Markers) , "pos_cM" ]
    w.vector <- tapply(MRList$pos_cM, MRList$BCIDX, function(x){ return(max(x, na.rm = T) - min(x, na.rm = T))})
    print(paste("Number of ", sum(is.infinite(w.vector))))
    w.vector[is.infinite(w.vector)] <- 0
    
  } else if ( Length == "gene" ) {
    Region <- as.data.frame(do.call(rbind, lapply(bc, blocklist_startend, type = "bp")))
    Region$HBIDX <- rownames(BlockInfo)
    Region$chr <- substr(Region$HBIDX, start = 1, stop = 2)
    gene_pos$middle_bp <- (gene_pos$start_bp + gene_pos$end_bp)/2
    
    w.vector <- apply(Region, 1, function(x){
      key1 = gene_pos$chr == as.numeric(x["chr"])
      key2 = gene_pos$middle_bp >= as.numeric(x["start"])
      key3 = gene_pos$middle_bp <= as.numeric(x["end"])
      return(sum(key1 & key2 & key3))
      
    })
    
  }
  
  
  
  # number of individuals
  N = nrow(W)
  
  # remove all markers with allele counts <= n.rm
  M_rm <- which(apply(W, 2, sum) < n.rm)
  if(length(M_rm) > 0){
    print(paste(length(M_rm), "markers removed due to allele counts <", n.rm))
    W <- W[ , -M_rm] 
  }
  
  # number of markers
  M = ncol(W)
  
  # minor allele frequencies * 2
  maf = colMeans(W)
  
  P = matrix(rep(maf, N),
             ncol = M,
             byrow = TRUE)
  
  # calculate Z
  Z <- W - P
  
  # Weight Matrix #
  WeightMat <- diag(w.vector^Scale)
  
  # calculate kin
  return(2*(Z %*% WeightMat %*% t(Z)) / sum( w.vector^Scale * maf * (2 - maf)))
}

RecodeMAF <- function(W){
  Mat <- W-1
  key <- apply(Mat, 2, sum) < 0
  vec <- (key + 0)*2-1
  trans_mat <- matrix(rep(vec, each = nrow(Mat)), byrow = F, nrow = nrow(Mat))
  Mat_new <- trans_mat * Mat + 1
  return(Mat_new)
}


### Geno Info ###
load(HBRdata)
load("../Data/GenePosition_B73v4.RData")
load("../Data/LinkageMap.RData")



#### Weighting GRM ####
HB <- RecodeMAF(HB)
Kin.Gv <- VanRaden_kinship_f_Weight(W = HB,
                                    BlockInfo = bc.info,
                                    MRList = MarkerList,
                                    Length = type,
                                    Scale = s,
                                    n.rm = 2)

# Kin.Gv is the weighted genomic relationship matrix (GRM).
# To run cross-validation with weighted GRM, 
# please replace function "VanRaden_kinship_f" with "VanRaden_kinship_f_Weight" in 02-01_CrossValidation_GBLUP.R

