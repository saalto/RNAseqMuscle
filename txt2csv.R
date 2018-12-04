setwd(data_dir)
filelist = list.files(pattern = ".txt")
for (i in 1:length(filelist)){
  input<-filelist[i]
  output<-paste0(gsub("\\.txt$", "",input), ".csv")
  print(paste("Processing the file:", input))
  data = read.delim(input, header = TRUE)   
  setwd("C:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/soleus/CSV")
  write.table(data, file=output, sep=",", col.names=TRUE, row.names=FALSE)
  setwd(data_dir)
}
