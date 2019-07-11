# ---- Amide and Nitrogen Chord ----
  # --- Complete list of GO terms examined in this chord ---
  group_AmideNitro <- c("amide biosynthetic process", "cellular amide metabolic process", "positive regulation of cellular amide metabolic process", "nitrogen compound metabolic process", "organonitrogen compound biosynthetic process", "logFC")
  
  # --- Combining and Renaming Categories ---
  # -- Amide group --
  group_Amide2 <- c("amide biosynthetic process", "cellular amide metabolic process", "positive regulation of cellular amide metabolic process")
  groups_Amide <- group_total[ , (colnames(group_total) %in% group_Amide2)]
  groups_Amide3 <- as.numeric(as.logical(base::rowSums(groups_Amide, na.rm=TRUE)))
  
  # -- Nitrogen group --
  group_Nitro2 <- c("nitrogen compound metabolic process", "organonitrogen compound biosynthetic process")
  groups_Nitro <- group_total[ , (colnames(group_total) %in% group_Nitro2)]
  group_Nitro3 <- as.numeric(as.logical(base::rowSums(groups_Nitro, na.rm=TRUE)))
  
  # --- Combining groups and Creating Chord for Amide and Nitrogen ---
  groups_AmideNitro2 <- cbind(Amide = groups_Amide3, Nitrogen = group_Nitro3)
  groups_AmideNitro2 <- cbind(groups_AmideNitro2, logFC = group_logFC2)
  
  
  
  # load 06_scr_chordplot2.R function and run the following:
  output <- "BLANK.pdf"
  scr_chordplot2(gene_frame = groups_AmideNitro2, goplot_data = groups_AmideNitro2, file_name=output)
  
  pl <- GOChord(groups_AmideNitro2,limit= c(1,0),space=0.02,gene.order='logFC',gene.space=0.25,lfc.col=c('red','white','blue'),lfc.min=fc_min,lfc.max=fc_max,border.size=0, process.label = 10)
  pl <- pl + coord_fixed()
  print(pl)
  dev.off()

  
  
  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- Nucleic Acids Chord ----
  # --- Complete list of GO terms examined in this chord ---
  group_NucleicAcids <- c("nucleic acid metabolic process", "RNA processing", "ribosome biogenesis", "posttranscriptional regulation of gene expression", "ncRNA processing", "logFC")
  
  # --- Combining and Renaming Categories ---
  # -- DNA group --
  group_DNA2 <- c("nucleic acid metabolic process", "DNA replication")
  groups_DNA <- group_total[ , (colnames(group_total) %in% group_DNA2)]
  groups_DNA3 <- as.numeric(as.logical(base::rowSums(groups_DNA, na.rm=TRUE)))
  
  # -- RNA group --
  group_RNA2 <- c("RNA processing")
  groups_RNA <- group_total[ , (colnames(group_total) %in% group_RNA2)]
  
  # -- Ribosome group --
  group_Ribo2 <- c("ribosome biogenesis")
  groups_Ribo <- group_total[ , (colnames(group_total) %in% group_Ribo2)]
  
  # -- Post-transcription --
  group_PostTran2 <- c("posttranscriptional regulation of gene expression")
  groups_PostTran <- group_total[ , (colnames(group_total) %in% group_PostTran2)]
  
  # -- ncRNA group --
  group_ncRNA2 <- c("ncRNA processing")
  groups_ncRNA <- group_total[ , (colnames(group_total) %in% group_ncRNA2)]
  
  # --- Combining groups and Creating chord for DNA, RNA, Ribosome, Post-transcription, and ncRNA ---
  groups_NucleicAcids <- cbind(DNA=groups_DNA3, RNA= groups_RNA, Ribosome = groups_Ribo, PostTranscription= groups_PostTran, ncRNA = groups_ncRNA)
  groups_NucleicAcids <- cbind(groups_NucleicAcids, logFC = group_logFC2)
  pl <- GOChord(groups_NucleicAcids,limit= c(1,0),space=0.02,gene.order='logFC',gene.space=0.25,lfc.col=c('red','white','blue'),lfc.min=fc_min,lfc.max=fc_max,border.size=0, process.label = 10)
  pl <- pl + coord_fixed()
  print(pl)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- Cellular Activities Chord ----
# --- Complete list of GO terms examined in this chord ---
group_CellularOrg2 <- c("cellular component biogenesis", "organelle organization", "cellular component assembly", "cytoskeleton organization", "logFC")

# --- Combining and Renaming Categories ---
# -- Cellular Component Biogenesis and Assembly --
group_BandA <- c("cellular component biogenesis","cellular component assembly")
groups_BandA <- group_total[ , (colnames(group_total) %in% group_BandA)]
groups_BandA3 <- as.numeric(as.logical(base::rowSums(groups_BandA, na.rm=TRUE)))

# -- Organelle --
group_organelle <- c("organelle organization")
groups_organelle <- group_total[ , (colnames(group_total) %in% group_organelle)]

# -- Cytoskeleton --
group_cyto <- c("cytoskeleton organization")
groups_cyto <- group_total[ , (colnames(group_total) %in% group_cyto)]

# --- Combining groups and Creating chord for Cellular Activities ---
groups_CellularOrg2 <- cbind(BiogenesisAndAssembly= groups_BandA3, Organelle= groups_organelle, Cytoskeleton= groups_cyto)
groups_CellularOrg2 <- cbind(groups_CellularOrg2, logFC = group_logFC2)
pl <- GOChord(groups_CellularOrg2,limit= c(1,0),space=0.02,gene.order='logFC',gene.space=0.25,lfc.col=c('red','white','blue'),lfc.min=fc_min,lfc.max=fc_max,border.size=0, process.label = 10)
pl <- pl + coord_fixed()
print(pl)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# --- Protein Chord ---
# -- Complete list of GO terms examined in this chord --
group_Protein2 <- c("peptide metabolic process", "proteolysis involved in cellular protein catabolic process", "cellular protein catabolic process", "protein ubiquitination", "regulation of protein ubiquitination", "protein maturation", "translational initiation","creatine metabolic process", "logFC")

# --- Combining and Renaming Categories ---
# -- Peptide Metabolism --
group_PepMet <- c("peptide metabolic process")
groups_PepMet <- group_total[ , (colnames(group_total) %in% group_PepMet)]

# -- Proteolysis --
group_Prote <- c("proteolysis involved in cellular protein catabolic process", "cellular protein catabolic process")
groups_Prote <- group_total[ , (colnames(group_total) %in% group_Prote)]
groups_Prote2 <- as.numeric(as.logical(base::rowSums(groups_Prote, na.rm=TRUE)))

# -- Ubiquitination --
group_Ubiq <- c("protein ubiquitination", "regulation of protein ubiquitination")
groups_Ubiq <- group_total[ , (colnames(group_total) %in% group_Ubiq)]
groups_Ubiq2 <- as.numeric(as.logical(base::rowSums(groups_Ubiq, na.rm=TRUE)))

# -- Protein Maturation --
group_ProteinMat <- c("protein maturation")
groups_ProteinMat <- group_total[ , (colnames(group_total) %in% group_ProteinMat)]

# -- Translation --
group_Trans <- c("translational initiation")
groups_Trans <- group_total[ , (colnames(group_total) %in% group_Trans)]

# -- Creatine Metabolism --
group_CreatineMet <- c("creatine metabolic process")
groups_CreatineMet <- group_total[ , (colnames(group_total) %in% group_CreatineMet)]

# --- Combining groups and Creating chord for Protein ---
groups_Protein2 <- cbind(PeptideMetabolism= groups_PepMet, Proteolysis=groups_Prote2, Ubiquitination=groups_Ubiq2, ProteinMaturation=groups_ProteinMat, Translation= groups_Trans,CreatineMetabolism=groups_CreatineMet)
groups_Protein2 <- cbind(groups_Protein2, logFC = group_logFC2)
pl <- GOChord(groups_Protein2,limit= c(1,0),space=0.02,gene.order='logFC',gene.space=0.25,lfc.col=c('red','white','blue'),lfc.min=fc_min,lfc.max=fc_max,border.size=0, process.label = 10)
pl <- pl + coord_fixed()
print(pl)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# --- Signaling Chord ---
# -- Complete list of GO terms examined in this chord --
group_Signaling2 <- c("positive regulation of transferase activity","regulation of telomerase activity","positive regulation of canonical Wnt signaling pathway", "TOR signaling", "regulation of cell cycle G1/S phase transition", "steroid hormone mediated signaling pathway", "regulation of intracellular steroid hormone receptor signaling pathway", "logFC")

# --- Combining and Renaming Categories ---
# -- TOR Signaling --
group_TOR <- c("TOR signaling")
groups_TOR <- group_total[ , (colnames(group_total) %in% group_TOR)]

# -- Wnt Signaling --
group_WntSignaling <- c("positive regulation of canonical Wnt signaling pathway")
groups_Wnt <- group_total[ , (colnames(group_total) %in% group_WntSignaling)]

# -- Steroid Hormone Signaling --
group_Steroid <- c("steroid hormone mediated signaling pathway", "regulation of intracellular steroid hormone receptor signaling pathway")
groups_Steroid <- group_total[ , (colnames(group_total) %in% group_Steroid)]
groups_Steroid2 <- as.numeric(as.logical(base::rowSums(groups_Steroid, na.rm=TRUE)))

# -- Cell Cycle --
group_CellCycle <- c("regulation of cell cycle G1/S phase transition", "regulation of telomerase activity")
groups_CellCycle <- group_total[ , (colnames(group_total) %in% group_CellCycle)]
groups_CellCycle2 <- as.numeric(as.logical(base::rowSums(groups_CellCycle, na.rm=TRUE)))

# -- Transferase Activity --
group_Transferase <- c("positive regulation of transferase activity")
groups_Transferase <- group_total[ , (colnames(group_total) %in% group_Transferase)]


# --- Combining groups and Creating chord for Signaling ---
groups_Signaling2 <- cbind(TOR=groups_TOR, WntSignaling=groups_Wnt, SteroidHormone=groups_Steroid2, CellCycle=groups_CellCycle2)
groups_Signaling2 <- cbind(groups_Signaling2, logFC= group_logFC2)
pl <- GOChord(groups_Signaling2,limit= c(1,0),space=0.02,gene.order='logFC',gene.space=0.25,lfc.col=c('red','white','blue'),lfc.min=fc_min,lfc.max=fc_max,border.size=0, process.label = 10)
pl <- pl + coord_fixed()
print(pl)
