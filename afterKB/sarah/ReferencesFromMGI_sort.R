# --- genelist consist of gene symbol and log fold change ---
genelist = read.delim(file="C:/Users/sarah/Documents/RNASeq_analysis/2019_04_29/genelist2.txt",header=TRUE)
# --- folder of text files of the references rom MGI ---
setwd("C:/Users/sarah/Documents/RNASeq_analysis/2019_04_29/Gene_References_SharedGO_SolvsTA")
files = list.files(".")
outputFile = "C:/Users/sarah/Downloads/AnnotatedReferenceList.txt"

# --- Clear and reset the lists ---
list1 = c()
list2 = c()
list3 = c()

# --- For each reference per gene given by MGI ---
for (i in 1: length(files)){
  data = read.delim(file= files[i],header = FALSE, sep = "\t")
  data = data.frame(data)
  
  # list2: Take the first column from the references file with the Reference Number
  list2 = c(list2, as.character(data$V1))
  # list1: Take the gene name and make a list of gene name the same length as the references listed
  for (j in 1:length(data[,1])){
    list1 = c(list1, as.character(genelist$GENE[i]))
  }
  # list3: Take the second column from the references file with the Citation Information
  list3 = c(list3, as.character(data$V2))
  
}

# --- Need to convert the  to a data frame ---
GeneSymbol <- cbind(Genes = list1, RefNum = list2, Reference = list3)

# --- Write to file
write.table(GeneSymbol, file = outputFile)
