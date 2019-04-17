
# --- genelist consist of gene symbol and log fold change ---
# --- folder of excel spreadsheets of the GO annotation from MGI ---

require(readxl)

# --- Change the CODE here ---
genelist = read.csv("C:/Users/sarah/Downloads/genelist.csv")
setwd("C:/Users/sarah/Downloads/MGI GO/")
files = list.files(".")
# data = read_excel(files[1], sheet = 1)

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

pasted_molecular_function = c()
pasted_biological_function = c()
pasted_cellular_component = c()

# --- Need to convert the gene list to a data frame ---
genelist = data.frame(genelist)

# --- For molecular function, biological process, and cellular component GO characteristics
# --- Extract and Paste terms together 
# --- and create a column for MF, BP, and CC
for (i in 1: length(files)){
  data = read_excel(files[i], sheet = 1)
  
  list1 = c()
  list2 = c()
  list3 = c()
  
  for (func in molecular_function){
    if (func %in% data$Category){
      list1 = c(list1, func)
    }else{
      list1 = c(list1, NA)
    }
  }
  with_NA1 = unlist(list1)
  no_NA1 = with_NA1[!is.na(with_NA1)]
  pasted_molecular_function[i] = paste(no_NA1, collapse = ", ") 
  
  for (func in biological_function){
    if (func %in% data$Category){
      list2 = c(list2, func)
    }else{
      list2 = c(list2, NA)
    }
  }
  with_NA2 = unlist(list2)
  no_NA2 = with_NA2[!is.na(with_NA2)]
  pasted_biological_function[i] = paste(no_NA2, collapse = ", ")  
  
  for (func in cellular_component){
    if (func %in% data$Category){
      list3 = c(list3, func)
    }else{
      list3 = c(list3, NA)
    }
  }
  with_NA3 = unlist(list3)
  no_NA3 = with_NA3[!is.na(with_NA3)]
  pasted_cellular_component[i] = paste(no_NA3, collapse = ", ")  
}

genelist$MolecularFunction = pasted_molecular_function
genelist$BiologicalFunction = pasted_biological_function
genelist$CellularComponent = pasted_cellular_component


# --- Change CODE here ---
# --- Write to file
write.csv(genelist, file = "C:/Users/sarah/Downloads/AnnotatedGeneList2.csv")
