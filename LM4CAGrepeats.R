library(dplyr)
library(limma)
library(Seurat)
#####
## generate pseudobulk counts - function
#######
pseudo_bulk_donor_region <- function(obj){
  
  hd_donor <- names(which(table(obj@colData$Donor, obj@colData$Condition)[,2] != 0))
  donor_ind_oligo <- list()
  donor_ind_opc <- list()
  
  for(i in 1:length(hd_donor)){
    
    donor_ind_oligo[[i]] <- which(obj@colData$Donor == hd_donor[i] & obj@colData$Lineage == "Oligodendrocyte")
    donor_ind_opc[[i]] <- which(obj@colData$Donor == hd_donor[i] & obj@colData$Lineage == "OPC" )
    
  }
  
  names(donor_ind_oligo) <- hd_donor
  names(donor_ind_opc) <- hd_donor
  
  oligo_donor_sum <- lapply(donor_ind_oligo,function(x){
    
    apply(obj@assays@data$counts[,x],1,sum)
    
  })
  
  opc_donor_sum <- lapply(donor_ind_opc,function(x){
    
    apply(obj@assays@data$counts[,x],1,sum)
    
  })
  
  oligo_hd_expr <- do.call(cbind,oligo_donor_sum)
  opc_hd_expr <- do.call(cbind,opc_donor_sum)
  
  oligo_hd_seurat <- Seurat::CreateSeuratObject(counts = oligo_hd_expr)
  oligo_hd_seurat <- Seurat::NormalizeData(oligo_hd_seurat, normalization.method = "LogNormalize")
  
  opc_hd_seurat <- Seurat::CreateSeuratObject(counts = opc_hd_expr)
  opc_hd_seurat <- Seurat::NormalizeData(opc_hd_seurat, normalization.method = "LogNormalize")
  
  oli_expr_mat <- as.matrix(oligo_hd_seurat@assays$RNA@data) 
  opc_expr_mat <- as.matrix(opc_hd_seurat@assays$RNA@data)   
  
  name1 <- paste0("Oligo")
  name2 <- paste0("OPC")
  expr_mat <- list(name1 = oli_expr_mat, name2 = opc_expr_mat)
  names(expr_mat) <- c(name1,name2)
  
  return(expr_mat)
  
  
} 
######
### read object
##########
obj<-readRDS("oligo_opc_seurat_object.rds")
############
#Pseudo bulking for Oligodendrocyte and OPC 
#within each region using helper function
accumbens_hd_mat <- pseudo_bulk_donor_region(obj, "Accumbens")
cingulate_hd_mat <- pseudo_bulk_donor_region(obj, "Cingulate")
caudate_hd_mat <- pseudo_bulk_donor_region(obj, "Caudate")

pseudo_bulk_summary <- list("Oligo_Accumbens" = accumbens_hd_mat$Oligo_Accumbens,
                            "OPC_Accumbens" = accumbens_hd_mat$OPC_Accumbens,
                            "Oligo_Cingulate" = cingulate_hd_mat$Oligo_Cingulate,
                            "OPC_Cingulate" = cingulate_hd_mat$OPC_Cingulate,
                            "Oligo_Caudate" = caudate_hd_mat$Oligo_Caudate,
                            "OPC_Caudate" = caudate_hd_mat$OPC_Caudate)
pseudo_bulk_summary$Oligo_all<-cbind(accumbens_hd_mat$Oligo_Accumbens,cingulate_hd_mat$Oligo_Cingulate) %>% cbind(caudate_hd_mat$Oligo_Caudate)
pseudo_bulk_summary$OPC_all<-cbind(accumbens_hd_mat$OPC_Accumbens,cingulate_hd_mat$OPC_Cingulate) %>% cbind(caudate_hd_mat$OPC_Caudate)
#### Prep metadata ####
hd_donor <- which(table(obj@colData$Donor, obj@colData$Condition)[,2] != 0) %>% names()
hd_donor_table_ind <- which(table(obj@colData$Donor, obj@colData$Condition)[,2] != 0)
all_cag <-names(table(obj@colData$Genotype))
all_age <- names(table(obj@colData$age))
all_gender <- c("M","F")
all_batch <- names(table(obj@colData$Batch))
all_region <- names(table(obj@colData$region))
hd_donor_cag <- apply(table(obj@colData$Donor, obj@colData$Genotype)[hd_donor_table_ind,],1,function(x){
  
  all_cag[which(x != 0)]
  
}) %>% as.numeric()

hd_donor_age <- apply(table(obj@colData$Donor, obj@colData$age)[hd_donor_table_ind,],1, function(x){
  
  all_age[which(x != 0)]
  
}) %>% unlist() %>% as.numeric() 

hd_donor_gender <- apply(table(obj@colData$Donor, obj@colData$gender)[hd_donor_table_ind,],1,function(x){
  
  all_gender[which(x != 0)]
  
}) %>% factor()

hd_donor_batch <- apply(table(obj@colData$Donor, obj@colData$Batch)[hd_donor_table_ind,],1,function(x){
  
  all_batch[which(x != 0)]
  
}) %>% factor()

hd_donor_region <- apply(table(obj@colData$Donor, obj@colData$region)[hd_donor_table_ind,],1,function(x){
  
  all_region[which(x != 0)]
  
}) %>% factor()
#Model matrix 
mat <- model.matrix(~hd_donor_gender)
mat <- cbind(mat,hd_donor_age, hd_donor_batch, hd_donor_cag, hd_donor_region)

#CAG regression analysis
pseudo_bulk_summary1<-pseudo_bulk_summary[7:8]
pseudo_bulk_summary_fit <- lapply(pseudo_bulk_summary1, function(x){
  
  temp <- lmFit(x,mat)
  temp <- eBayes(temp)
  return(temp)
  
})

#Coefficient extraction
df_coeff <- lapply(pseudo_bulk_summary_fit1, function(x){
  
  all_gene <- rownames(pseudo_bulk_summary_fit$Oligo_all$coefficients)
  df <- data.frame("Gene" = all_gene, "Coefficient" = x$coefficients[,5], "p-val" = x$p.value[,5])
  rownames(df) <- 1:dim(df)[1]
  return(df)
  
})
