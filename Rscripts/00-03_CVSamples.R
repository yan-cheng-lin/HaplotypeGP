#### Loading ####
dat.pheno.raw <- read.table("../Data/PhenotypicData.txt", header = T, sep = "\t")
traits <- colnames(dat.pheno.raw)[-1]
ids <- dat.pheno.raw$Genotype

#### Scenario 1 ####
ID.Sel <- list()
for (pop in c("DH_KE", "GC_KE", "DH_PE", "GC_PE")){
  id.geno <- grep(ids, pattern = pop, value = T)
  
  for (trait in traits){
    dat.pheno <- na.omit(dat.pheno.raw[ , c("Genotype", trait)])
    lines <- intersect(id.geno, dat.pheno$Genotype)
    n <- length(lines)
    i <- 1
    for (batch in 1:10){
      set.seed(batch)
      RandomOrder <- sample(lines)
      foldsize <- round(n/5)
      
      for (fold in 1:5){
        
        if(fold == 5){
          val <- (foldsize*4+1):n
        } else {
          val <- (foldsize * (fold - 1) + 1):(foldsize * fold)
        }
        
        ID.Sel[[pop]][[trait]][[i]] <- data.frame(ID = lines,
                                                  Set = ifelse(lines %in% RandomOrder[val], "Valid", "Train"))
        
        i <- i+1
      }
      
    }
  }
  
}

#### Scenario 2 ####
DHtoGC.ID.Sel <- list()
GCtoDH.ID.Sel <- list()

### KE ###
for (trait in traits){
  dat.pheno <- na.omit(dat.pheno.raw[ , c("Genotype", trait)])
  lines <- grep(dat.pheno$Genotype, pattern = "KE", value = T)
  
  for (i in 1:100){
    set.seed(i)
    id.DH.select <- sample(x = grep(lines, pattern = "DH", value = T), 200)
    id.GC.select <- sample(x = grep(lines, pattern = "GC", value = T), 50)
    id.select <- c(id.DH.select, id.GC.select)
    DHtoGC.ID.Sel[["KE"]][[trait]][[i]] <- id.select
    
    set.seed(i)
    id.DH.select <- sample(x = grep(lines, pattern = "DH", value = T), 50)
    id.GC.select <- sample(x = grep(lines, pattern = "GC", value = T), 200)
    id.select <- c(id.DH.select, id.GC.select)
    GCtoDH.ID.Sel[["KE"]][[trait]][[i]] <- id.select
    
  }
}
### PE ###
for (trait in traits){
  dat.pheno <- na.omit(dat.pheno.raw[ , c("Genotype", trait)])
  lines <- grep(dat.pheno$Genotype, pattern = "PE", value = T)
  
  for (i in 1:100){
    set.seed(i)
    id.DH.select <- sample(x = grep(lines, pattern = "DH", value = T), 200)
    id.GC.select <- sample(x = grep(lines, pattern = "GC", value = T), 50)
    id.select <- c(id.DH.select, id.GC.select)
    DHtoGC.ID.Sel[["PE"]][[trait]][[i]] <- id.select
    
    set.seed(i)
    id.DH.select <- sample(x = grep(lines, pattern = "DH", value = T), 50)
    id.GC.select <- sample(x = grep(lines, pattern = "GC", value = T), 200)
    id.select <- c(id.DH.select, id.GC.select)
    GCtoDH.ID.Sel[["PE"]][[trait]][[i]] <- id.select
    
  }
}

#### Scenario 3 ####
KEtoPE.ID.Sel <- list()
PEtoKE.ID.Sel <- list()

### DH ###
for (trait in traits){
  dat.pheno <- na.omit(dat.pheno.raw[ , c("Genotype", trait)])
  lines <- grep(dat.pheno$Genotype, pattern = "DH", value = T)
  
  for (i in 1:100){
    set.seed(i)
    id.KE.select <- sample(x = grep(lines, pattern = "KE", value = T), 200)
    id.PE.select <- sample(x = grep(lines, pattern = "PE", value = T), 50)
    id.select <- c(id.KE.select, id.PE.select)
    KEtoPE.ID.Sel[["DH"]][[trait]][[i]] <- id.select
    
    set.seed(i)
    id.PE.select <- sample(x = grep(lines, pattern = "PE", value = T), 200)
    id.KE.select <- sample(x = grep(lines, pattern = "KE", value = T), 50)
    id.select <- c(id.PE.select, id.KE.select)
    PEtoKE.ID.Sel[["DH"]][[trait]][[i]] <- id.select
    
  }
}

### GC ###
for (trait in traits){
  dat.pheno <- na.omit(dat.pheno.raw[ , c("Genotype", trait)])
  lines <- grep(dat.pheno$Genotype, pattern = "GC", value = T)
  
  for (i in 1:100){
    set.seed(i)
    id.KE.select <- sample(x = grep(lines, pattern = "KE", value = T), 200)
    id.PE.select <- sample(x = grep(lines, pattern = "PE", value = T), 50)
    id.select <- c(id.KE.select, id.PE.select)
    KEtoPE.ID.Sel[["GC"]][[trait]][[i]] <- id.select
    
    set.seed(i)
    id.PE.select <- sample(x = grep(lines, pattern = "PE", value = T), 200)
    id.KE.select <- sample(x = grep(lines, pattern = "KE", value = T), 50)
    id.select <- c(id.PE.select, id.KE.select)
    PEtoKE.ID.Sel[["GC"]][[trait]][[i]] <- id.select
    
  }
}

#### Save CV Samples ####
save(list = c("ID.Sel", "DHtoGC.ID.Sel", "GCtoDH.ID.Sel", "KEtoPE.ID.Sel", "PEtoKE.ID.Sel"),
     file = "../Data/CV_SampleList.RData")
