ChunkFileCreator <- function(nSamples, ChunkSize) {
  numberSamples = nSamples
  totalChunks = round(numberSamples/ChunkSize)
  Skipnumber = 0
  fileNumber = 1
  ChunkSize = ChunkSize
  columnnames = colnames(fread("SNPsDF_recodeA.raw", nrow = 1, drop=c(1,3:6)))
  for (chunk in seq(1, totalChunks)) {
    
    
    if (chunk == 1) {
      fileChunk = fread("SNPsDF_recodeA.raw", nrow = ChunkSize - 1,
                        skip = Skipnumber, drop=c(1,3:6))
      saveRDS(fileChunk, file = paste0("part",fileNumber,".rds"))
    }
    
    else if ( Skipnumber + ChunkSize > numberSamples) {
      
      fileChunk = fread("SNPsDF_recodeA.raw", nrow = numberSamples - Skipnumber,
                        skip = Skipnumber, drop=c(1,3:6), col.names = columnnames)
      saveRDS(fileChunk, file = paste0("part",fileNumber,".rds"))
    }
    
    else {
      fileChunk = fread("SNPsDF_recodeA.raw", nrow = ChunkSize,
                        skip = Skipnumber, drop=c(1,3:6), col.names = columnnames)
      saveRDS(fileChunk, file = paste0("part",fileNumber,".rds"))
    }
    
    Skipnumber = Skipnumber + ChunkSize
    fileNumber = fileNumber + 1
  }
  
}



GWAS_survival <- function(clinical, PCs, i) {
  
  geneticChunk <- read.csv(paste0("GeneticChunksFiles/genetic_",i,".txt"), sep = "", header = TRUE)
  PCs <- PCs
  clinical <- clinical
  
  TABLE <- clinical %>%
    inner_join(PCs, by = "IID") %>%
    inner_join(geneticChunk, by = "IID")
  
  TABLE <- as.data.frame(TABLE)
  
  # We initialize the df where we will store the results of the analysis
  coefficients<-as.data.frame(matrix(ncol= 9))
  names(coefficients) <- c("SNP","Coeff", "se", "Pvalue", "Cox.zphPVal", "N", "ov.lik.ratio","logrank", "r2" )
  
  # Then we gon into the for loop to generate a cox model for each snp "z" from the "i" genetic chunk
  for (z in 25:ncol(TABLE)) {
    #print(colnames(TABLE)[z])
    snp <- TABLE[,c(z,1:24)]
    
    model.cox <- coxph(Surv(snp$time_Event_midpoint, snp$event_dyskinesia) ~
                         snp[,1] + snp$SEX + snp$AAO.std + snp$PC1 + snp$PC2 +snp$PC3 +snp$PC4 +snp$PC5
                       + duration.std + ldopaDose.std + motorBL.std, data = snp)
    kmz<- cox.zph(model.cox, transform = "km") # We test the PH assumption
    
    j = z - 24 # To access the coefficient table from the very beginning
    
    coefficients[j,1]<- paste(colnames(TABLE)[z])
    coefficients[j,2]<- summary(model.cox)$coefficients[1,1]
    coefficients[j,3]<- summary(model.cox)$coefficients[1,3]
    coefficients[j,4]<- summary(model.cox)$coefficients[1,5]
    coefficients[j,5]<- kmz$table[1,3] # We check the PH assumption
    coefficients[j,6]<- model.cox$n
    coefficients[j,7]<- summary(model.cox)$logtest[[1]]
    coefficients[j,8]<- summary(model.cox)$sctest[[1]]
    coefficients[j,9]<- summary(model.cox)$rsq[[1]] # nagelkerke r square
  }
  coefficients
}
