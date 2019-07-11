# --- Load the libraries of pheatmap and RColorBrewer---
  library(pheatmap)
  library(RColorBrewer)
  colorz = rev(brewer.pal(10,"RdBu"))
  # set variables as well as the locations of necessary files
  GeneUniverse <- readRDS(file="C:/Users/sarah/Dropbox/Kioussi-Brown/sol-gn-tib-manuscript/data/GeneUniverse.RDS")
  mouse_hash <- readRDS(file = "C:/Users/sarah/Dropbox/Kioussi-Brown/sol-gn-tib-manuscript/data/mouse-hash.rds")
  enrich_type <- "GO:BP"
  
# -- (7) Filter out genes associated with Metabolic Process (GO:0008152) via PANTHER
  # filter out genes in soleus data set 
  soleus_WO6_DEGenes <- deedgSOLWO6_flt$table[row.names(deedgSOLWO6_flt$table) %in% deedgSOLWO6_flt$sig_genes,]
  soleus_WO6_DEMetGenes <- soleus_WO6_DEGenes[row.names(soleus_WO6_DEGenes) %in% GeneUniverse, ]
  # filter out genes in tibialis data set
  tibialis_renamed_DEGenes <- deedgTA_flt_renamed$table[row.names(deedgTA_flt_renamed$table) %in% 
                                                          deedgTA_flt_renamed$sig_genes,]
  tibialis_renamed_DEMetGenes <- tibialis_renamed_DEGenes[row.names(tibialis_renamed_DEGenes) %in% GeneUniverse,]
  
# -- (8) GO term enrichment analysis--
  # load 03_catenrich.R function
  
  # run the following
  soleus_WO6_enrich <- catenrich(all.genes = GeneUniverse, sig.genes = row.names(soleus_WO6_DEMetGenes), 
                                 entrez_hash = mouse_hash, enrich_type = enrich_type)
  tibialis_renamed_enrich <- catenrich(all.genes = GeneUniverse, sig.genes = row.names(tibialis_renamed_DEMetGenes), 
                                       entrez_hash = mouse_hash, enrich_type = enrich_type)
  # identify shared GO terms
  SharedSOLvsTA <- Reduce(intersect, list(soleus_WO6_enrich$overrep$category, 
                                          tibialis_renamed_enrich$overrep$category))
  # load the genes associated with those shared GO terms
  soleus_sig_genes <- readRDS(file="C:/Users/sarah/Dropbox/Kioussi-Brown/sol-gn-tib-manuscript/data/soleus_overrep_genes.RDS")
  tibialis_sig_genes <- readRDS(file="C:/Users/sarah/Dropbox/Kioussi-Brown/sol-gn-tib-manuscript/data/tibialis_overrep_genes.RDS")
  # concatenate the list of genes
  soleus_sig_genes2 <- Reduce(union, soleus_sig_genes)
  tibialis_sig_genes2 <- Reduce(union, tibialis_sig_genes)
  # identify shared genes between the two data sets
  SharedSOLvsTA2 <- Reduce(intersect, list(soleus_sig_genes2, tibialis_sig_genes2))

# -- (9) identify the average log-fold expression combining soleus and tibialis anterior data sets --
  # create a table of the log-fold changes
  soleus_joined_sig_genes <- deedgSOLWO6_flt$table[row.names(deedgSOLWO6_flt$table) %in% SharedSOLvsTA2,]
  soleus_joined_sig_genes1 <- soleus_joined_sig_genes[order(row.names(soleus_joined_sig_genes)),]
  tibialis_joined_sig_genes <- deedgTA_flt_renamed$table[row.names(deedgTA_flt_renamed$table) %in% SharedSOLvsTA2,]
  tibialis_joined_sig_genes1 <- tibialis_joined_sig_genes[order(row.names(tibialis_joined_sig_genes)),]
  sharedFC <- cbind(Soleus_FC= soleus_joined_sig_genes1$logFC, Tibialis_FC= tibialis_joined_sig_genes1$logFC)
  row.names(sharedFC) <- row.names(soleus_joined_sig_genes1)
  sharedFC <- data.frame(sharedFC)
  # take the average of genes with same log-fold change sign between both data sets
  averages <- rowMeans(sharedFC)
  sharedFC <- cbind(sharedFC, Average=averages)
  # take the absolute value of the average to remove log-fold change less than 1
  removal <- abs(sharedFC$Average)
  sharedFC <- cbind(sharedFC, Absolute=removal)
  # compare log-fold change signs
  checkSign <- sign(sharedFC$Soleus_FC)*sign(sharedFC$Tibialis_FC)
  sharedFC <- cbind(sharedFC, Sign=checkSign)
  sharedFC <- sharedFC[with(sharedFC, order(sharedFC$Sign, sharedFC$Average, decreasing = TRUE)),]
  # remove opposing log-fold change signs
  sharedFC1 <- sharedFC[sharedFC$Sign == 1,]
  # remove average log-fold change less than the absolute value of one
  sharedFC1 <- sharedFC1[!sharedFC1$Absolute < 1,]
  
# -- (10) Generate heatmap of the shared genes and their log-fold changes
  # isolate the greater than 1 log-fold change of shared metabolic process between soleus and tibialis 
  sharedFC2 <- cbind(Soleus=sharedFC1$Soleus_FC, Tibialis=sharedFC1$Tibialis_FC)
  row.names(sharedFC2) <- row.names(sharedFC1)
  # creating Figure 3 heatmap of the log-fold changes
  breakList = seq(min(sharedFC2),max(sharedFC2), by=0.2)
  colorz = colorRampPalette(colorz)(length(breakList))
  pheatmap(sharedFC2, breaks=breakList, color = colorz, cluster_rows = FALSE, show_colnames = TRUE, 
           show_rownames = TRUE, treeheight_row = 0, treeheight_col = 0)

# -- (11) GO terms based on ENTREZ IDs and GO evidence
  # load 04_GO_genes.R function
  
  # run the following:
  sharedSOLvsTA3 <- cbind(logFC=sharedFC1$Average)
  row.names(sharedSOLvsTA3) <- row.names(sharedFC1)
  sharedSOLvsTA3 <- data.frame(sharedSOLvsTA3)
  sharedSOLvsTA_allgo <- GO_genes(gene_frame = sharedSOLvsTA3, double_hash = mouse_hash, 
                                  go_ids = SharedSOLvsTA)
  sharedSOLvsTA_allgo <- sharedSOLvsTA_allgo[sharedSOLvsTA_allgo$GENES %in% row.names(sharedSOLvsTA3),]

# -- (12) Chord plots of GO terms and corresponding genes
  # load 05_GO_assemble.R function
  
  # run the following:
  sharedSOLvsTA_goplot <- GO_assemble(gene_frame = sharedSOLvsTA3, double_hash = mouse_hash,
                                    all_go = sharedSOLvsTA_allgo, overrep_file = soleus_WO6_enrich,
                                    go_ids = SharedSOLvsTA)
  # set universal variables
  group_total <- sharedSOLvsTA_goplot$chord[ , colSums(sharedSOLvsTA_goplot$chord !=0)>0]
  group_logFC <- c("logFC")
  group_logFC2 <- group_total[ , (colnames(group_total) %in% group_logFC)]
  
  
# load 06_scr_chordplot2.R function and run the following:
#soleus_WO6_chord <- scr_chordplot2(gene_frame = sharedSOLvsTA, goplot_data = sharedSOLvsTA_goplot, file_name=output)
