# --- genelist consist of gene symbol and log fold change ---
# --- folder of excel spreadsheets of the GO annotation from MGI ---
require(readxl)

# --- Change the CODE here ---
genelist = read.csv("C:/Users/sarah/Documents/RNASeq_analysis/2019_04_17/genelist.csv")
setwd("C:/Users/sarah/Documents/RNASeq_analysis/2019_04_17/MGI GO/")
files = list.files(".")
outputFile = "C:/Users/sarah/Downloads/AnnotatedGeneList4.csv"

# ---Leave this code (unless DeBugging)---
# --- Subsections of the three major GO categories ---
molecular_function = c("carbohydrate derivative binding", "cytoskeletal protein binding", 
                       "DNA binding", "enzyme regulator", "hydrolase", "ligase", "lipid binding", 
                       "oxidoreductase", "RNA binding", "signaling receptor activity", 
                       "signaling receptor binding", "transcription", "transferase", "transporter")
biological_function = c("carbohydrate derivative metabolism", "cell death", "cell differentiation", 
                        "cell population proliferation", "cellular component organization", 
                        "establishment of localization", "homeostatic process", 
                        "immune system process", "lipid metabolic process", "nucleic acid-templated transcription",
                        "protein metabolic process", "response to stimulus", "signaling", "system development")
cellular_component = c("cell projection", "cytoplasmic vesicle", "cytoskeleton", "cytosol", "endoplasmic reticulum", 
                       "endosome", "extracellular region", "Golgi apparatus", "mitochondrion", 
                       "non-membrane-bounded organelle", "nucleus", "organelle envelope", "organelle lumen", 
                       "plasma membrane", "synapse", "vacuole")

# --- Clear and reset the lists ---
pasted_molecular_function = c()
pasted_biological_function = c()
pasted_cellular_component = c()
pasted_references = c()


# --- Need to convert the gene list to a data frame ---
genelist = data.frame(genelist)

# --- For molecular function, biological process, and cellular component GO characteristics ---
# --- Extract and Paste terms together and create a column for MF, BP, and CC ---
for (i in 1: length(files)){
  data = read_excel(files[i], sheet = 1)
  
  # - only select rows with the appropriate ontology
  # - ignore any terms whose evidence comes from ND (no biological data available), 
  #   NAS (non-traceable author statement), or IEA (Inferred from electronic annotation)
  data <- data[data$Evidence != "ND" & data$Evidence != "NAS" & data$Evidence != "IEA",]
  
  list1 = c()
  list1_1 = c()
  list2 = c()
  list3 = c()
  
  # -- Molecular Function --
  for (func in molecular_function){
    if (func %in% data$Category){
      list1 = c(list1, func)
      list1_1 = c(list1_1, data$`Reference(s)`)
    }else{
      list1 = c(list1, NA)
    }
  }
  with_NA1 = unlist(list1)
  no_NA1 = with_NA1[!is.na(with_NA1)]
  pasted_molecular_function[i] = paste(no_NA1, collapse = ", ") 
  
  # -- Biological Process --
  for (func in biological_function){
    if (func %in% data$Category){
      list2 = c(list2, func)
      list1_1 = c(list1_1, data$`Reference(s)`)
    }else{
      list2 = c(list2, NA)
    }
  }
  with_NA2 = unlist(list2)
  no_NA2 = with_NA2[!is.na(with_NA2)]
  pasted_biological_function[i] = paste(no_NA2, collapse = ", ")  
  
  # -- Cellular Component --
  for (func in cellular_component){
    if (func %in% data$Category){
      list3 = c(list3, func)
      list1_1 = c(list1_1, data$`Reference(s)`)
    }else{
      list3 = c(list3, NA)
    }
  }
  with_NA3 = unlist(list3)
  no_NA3 = with_NA3[!is.na(with_NA3)]
  pasted_cellular_component[i] = paste(no_NA3, collapse = ", ") 
  
  with_NA4 = unlist(list1_1)
  no_NA4 = with_NA4[!is.na(with_NA4)]
  no_NA5 = Reduce(union, no_NA4)
  pasted_references[i] = paste(no_NA5, collapse = "; ")
  
}

genelist$MolecularFunction = pasted_molecular_function
genelist$BiologicalFunction = pasted_biological_function
genelist$CellularComponent = pasted_cellular_component
genelist$References = pasted_references

# --- Write to file
write.csv(genelist, file = outputFile)
