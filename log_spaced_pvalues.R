x <- readRDS('de-deseq-results.rds')
x <- x[order(x$padj),]
y <- readRDS('de-edgeR-results.rds')

log_spaced <- log(x$padj)
df1 <- data.frame(log_spaced)
df1$condition <- "DESeq"

log_spaced <- log(y$PValue)
df2 <- data.frame(log_spaced)
df2$condition <- "edgeR"

df3 <- rbind(df1, df2)

library(reshape2)
m.df <- melt(df3, id.vars="condition")
sbox <- ggplot(m.df, aes(x=condition, y=value)) + geom_boxplot(aes(fill=variable)) + theme(legend.position = "none")
sbox + geom_hline(yintercept = log(0.05), linetype= "dashed", color = "red") + geom_hline(yintercept = log(0.01), linetype = "dashed", color = "green")

