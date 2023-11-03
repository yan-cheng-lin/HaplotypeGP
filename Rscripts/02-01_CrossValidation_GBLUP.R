#### Libraries and Functions ####
### Library ###
library("HaploBlocker")
library("magrittr")
library("dplyr")
library("asreml")

### Functions ###
## Construction of van Raden kinship ##
# input: genotype matrix coded as 0, 1, 2 (rows = individuals and columns = markers)
# n.rm = 2, remove haplotype alleles, which total number in the population < 2 #
VanRaden_kinship_f <- function(W, n.rm = 2) {
  
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
  
  # calculate kin
  return(2*(Z %*% t(Z)) / sum(maf * (2 - maf)))
  
}


# Making G a positive semidefinite matrix by bending
# Modify the negative and very small positive eigenvalues (<0.0001) of the relationship matrix 
# Jorjani, H., Klei, L. and Emanuelson, U., 2003. A simple method for weighted bending of genetic (co) variance matrices. Journal of dairy science, 86(2), pp.677-679.
possemiD <- function(matrix, epsilon = 0.0001){
  eigenval <- eigen(matrix)
  eigenval$values[eigenval$values < epsilon] <- epsilon
  n <- length(eigenval$values)
  Eigen <- matrix(nrow = n, ncol = n, data = 0)
  diag(Eigen) <- eigenval$values
  kmat2 <- (eigenval$vectors %*% Eigen) %*% t(eigenval$vectors)
  rownames(kmat2) <- rownames(matrix)
  colnames(kmat2) <- colnames(matrix)
  return(kmat2)
}

#### Settings ####
HBRData <- "../Data/HaploView_HMat.RData" # Path to the RData of haplotype matrix (*_HMat.RData), with object "HB"
TrainSet <- "DH_PE" # Training Set: DH_KE, DH_PE, GC_KE and GC_PE
ValidSet <- "DH_KE" # Validation Set: DH_KE, DH_PE, GC_KE and GC_PE

### Load Data ###
load(HBRData)
load("../Data/CV_SampleList.RData")
dat.pheno.raw <- read.table("../Data/PhenotypicData.txt",header = T)
traits <- c("EV_V6", "PH_V6", "PH_final", "DtSILK", "RL_R6")
TS <- unlist(strsplit(TrainSet, split = "_"))
PS <- unlist(strsplit(ValidSet, split = "_"))

idx <- sum(which(TS == PS)) # Scenarion1 = 3; Scenarion2 = 2; Scenarion3 = 1

#### Cross Validation ####
if (idx == 3){
  
  ### Scenario 1 ###
  
  ## Five-fold 10 times CrossValidation ##
  CV.mat <- matrix(nrow = 50, ncol = length(traits),
                   dimnames = list(1:50, traits))
  for (i in 1:50){
    for(trait in traits){
      ### Sampling ###
      id.select <- ID.Sel[[TrainSet]][[trait]][[i]][,1]
      id.valid <- with(ID.Sel[[TrainSet]][[trait]][[i]], ID[Set == "Valid"])
      
      ### Make GRM and G Inverse ###
      HB.sample <- HB[ id.select, ]
      Kin.Gv <- VanRaden_kinship_f(HB.sample)
      
      line <- id.select
      dat.pheno <- na.omit(dat.pheno.raw[ , c("Genotype", trait)])
      G <- Kin.Gv
      G <- possemiD(G)
      G.inverse <- solve(G)
      attr(G.inverse, "rowNames") <- line
      attr(G.inverse, "colNames") <- line
      attr(G.inverse, "INVERSE") <- TRUE
      
      ### Prepare Phenotypic Data ###
      y.raw <- dat.pheno[match( line , dat.pheno$Genotype ), trait]
      dat <- data.frame(y = y.raw, line = as.factor(line))
      valid <- dat$line %in% id.valid
      dat$y[ valid ] <- NA
      
      
      ### Build GBLUP by ASReml ###
      modelGBLUP <- tryCatch({
        
        asreml(fixed = y~1,
               random = ~vm(line, G.inverse),
               data = dat, trace = F)
        
      }, warning = function(w) {
        
        asreml(fixed = y~1,
               random = ~vm(line, G.inverse),
               data = dat, trace = F)
        
      }, error = function(e) {
        
        return(NA)
        
        print("error")
      })
      
      ### Prediction Accuracy ###
      if ( !is.list(modelGBLUP) ){
        
        CV.mat[ i , trait ] <- NA
        
      } else {
        
        BLUPs <- summary(modelGBLUP, coef = T)$coef.random[,"solution"] + modelGBLUP$coefficients$fixed[1]
        Cor <- cor(y.raw[valid],BLUPs[valid])
        
        CV.mat[ i , trait ] <- Cor
        
      }
    }
  }
  
  
} else if (idx == 2) {
  ### Scenario 2 ###
  
  direction <- paste0(TS[1], "to", PS[1])
  pop <- TS[2]
  
  ## 100 times CrossValidation ##
  CV.mat <- matrix(nrow = 100, ncol = length(traits),
                   dimnames = list(1:100, traits))
  Samples <- get(paste0(direction, ".ID.Sel"))[[pop]] # Samples for validation Ex: DHtoGC.ID.Sel[["KE"]]
  
  
  for (i in 1:100){
    for(trait in traits){
      id.select <- unlist(Samples[[trait]][i])
      id.DH.select <- grep(id.select, pattern = "DH", value = T)
      id.GC.select <- grep(id.select, pattern = "GC", value = T)
      
      if(direction == "DHtoGC"){
        id.valid <- id.GC.select
      } else if(direction == "GCtoDH") {
        id.valid <- id.DH.select
      } else {
        cat("Wrong validation set")
      }
      
      ### Make GRM and G Inverse ###
      HB.sample <- HB[ id.select, ]
      Kin.Gv <- VanRaden_kinship_f(HB.sample)
      
      line <- id.select
      dat.pheno <- na.omit(dat.pheno.raw[ , c("Genotype", trait)])
      G <- Kin.Gv
      G <- possemiD(G)
      G.inverse <- solve(G)
      attr(G.inverse, "rowNames") <- line
      attr(G.inverse, "colNames") <- line
      attr(G.inverse, "INVERSE") <- TRUE
      
      ### Prepare phenotypic data ###
      y.raw <- dat.pheno[match( line , dat.pheno$Genotype ), trait]
      dat <- data.frame(y = y.raw, line = as.factor(line))
      valid <- dat$line %in% id.valid
      dat$y[ valid ] <- NA
      
      ### Build GBLUP by ASReml ###
      modelGBLUP <- tryCatch({
        
        asreml(fixed = y~1,
               random = ~vm(line, G.inverse),
               data = dat, trace = F)
        
      }, warning = function(w) {
        
        asreml(fixed = y~1,
               random = ~vm(line, G.inverse),
               data = dat, trace = F)
        
      }, error = function(e) {
        
        return(NA)
        
        print("error")
      })
      
      ### Prediction Accuracy ###
      if ( !is.list(modelGBLUP) ){
        
        CV.mat[ i , trait ] <- NA
        
      } else {
        
        BLUPs <- summary(modelGBLUP, coef = T)$coef.random[,"solution"] + modelGBLUP$coefficients$fixed[1]
        Cor <- cor(y.raw[valid],BLUPs[valid])
        
        CV.mat[ i , trait ] <- Cor
        
      }
      
      
    }
  }
  

} else if (idx == 1) {
  
  ### Scenario 3 ###
  
  direction <- paste0(TS[2], "to", PS[2])
  pop <- TS[1]
  
  ## 100 times CrossValidation ##
  CV.mat <- matrix(nrow = 100, ncol = length(traits),
                   dimnames = list(1:100, traits))
  Samples <- get(paste0(direction, ".ID.Sel"))[[pop]] # Samples for validation Ex: DHtoGC.ID.Sel[["KE"]]
  
  
  for (i in 1:100){
    for(trait in traits){
      id.select <- unlist(Samples[[trait]][i])
      id.KE.select <- grep(id.select, pattern = "KE", value = T)
      id.PE.select <- grep(id.select, pattern = "PE", value = T)
      
      if(direction == "KEtoPE"){
        id.valid <- id.PE.select
      } else if(direction == "PEtoKE") {
        id.valid <- id.KE.select
      } else {
        cat("Wrong validation set")
      }
      
      ### Make GRM and G Inverse ###
      HB.sample <- HB[ id.select, ]
      Kin.Gv <- VanRaden_kinship_f(HB.sample)
      
      line <- id.select
      dat.pheno <- na.omit(dat.pheno.raw[ , c("Genotype", trait)])
      G <- Kin.Gv
      G <- possemiD(G)
      G.inverse <- solve(G)
      attr(G.inverse, "rowNames") <- line
      attr(G.inverse, "colNames") <- line
      attr(G.inverse, "INVERSE") <- TRUE
      
      ### Prepare phenotypic data ###
      y.raw <- dat.pheno[match( line , dat.pheno$Genotype ), trait]
      dat <- data.frame(y = y.raw, line = as.factor(line))
      valid <- dat$line %in% id.valid
      dat$y[ valid ] <- NA
      
      ### Build GBLUP by ASReml ###
      modelGBLUP <- tryCatch({
        
        asreml(fixed = y~1,
               random = ~vm(line, G.inverse),
               data = dat, trace = F)
        
      }, warning = function(w) {
        
        asreml(fixed = y~1,
               random = ~vm(line, G.inverse),
               data = dat, trace = F)
        
      }, error = function(e) {
        
        return(NA)
        
        print("error")
      })
      
      ### Prediction Accuracy ###
      if ( !is.list(modelGBLUP) ){
        
        CV.mat[ i , trait ] <- NA
        
      } else {
        
        BLUPs <- summary(modelGBLUP, coef = T)$coef.random[,"solution"] + modelGBLUP$coefficients$fixed[1]
        Cor <- cor(y.raw[valid],BLUPs[valid])
        
        CV.mat[ i , trait ] <- Cor
        
      }
      
      
    }
  }
  
} else {
  cat("Please check pop of the training set and validation set.")
}



