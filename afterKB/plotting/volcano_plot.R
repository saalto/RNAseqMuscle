#' script for making Volcano plots for RNAseq data
#' @param gene_frame: data frame with fold change/p value info
#' @param fc_name: name of column corresponding to log fold change
#' @param p_name: name of column with p-values (adjusted)
#' @param p_levels: a sorted (smallest to largest) vector of three (only 3!) cutoffs for color coding
#' @param group_colors: a set of four colors (1 per p_level +1 for nonsignificant genes)
#' @param file_name: name of the file to write the result
#' @export
volcano_plot <- function(gene_frame,fc_name="logFC",p_name="FDR",p_levels=c(0.01,0.025,0.05),group_colors=c("#009E73","#CC79A7","#000000","#999999"),file_name){
  require(ggplot2)
  # append a column to gene_frame that computes -log(p)
  gene_frame[["neglogp"]] <- -1.0*log(gene_frame[[p_name]])
  # create a color column with three levels based on FDR
  # d marks nonsignificant genes
  # c marks genes >= level 1
  # b marks genes between levels 1 and 2
  # a marks genes between levels 2 and 3
  c_vec <- rep("d",dim(gene_frame)[1])
  c_vec[gene_frame[[p_name]] <= p_levels[3]] <- "c"
  c_vec[gene_frame[[p_name]] <= p_levels[2]] <- "b"
  c_vec[gene_frame[[p_name]] <= p_levels[1]] <- "a"
  # bind to the data frame
  gene_frame[["cutoff"]] <- c_vec
  # get overall x limits
  max_x <- max(abs(gene_frame[[fc_name]]))
  # get the biggest neglogp
  max_y <- max(gene_frame[["neglogp"]])
  pdf(file_name)
  pl <- ggplot(gene_frame, aes_string(x=fc_name,y="neglogp",colour="cutoff")) + geom_point(size=1.0,alpha=0.1) + xlim(-max_x,max_x) + ylim(0,1.01*max_y)
  # fix the colors
  pl <- pl + scale_colour_manual(labels=c(paste("p <",p_levels[1]),paste(p_levels[2],"> p >",p_levels[1]),paste(p_levels[3],"> p >",p_levels[2]),paste("p >",p_levels[3])),values=group_colors)
  # remove gray background
  pl <- pl + theme_bw()
  # horizontal line at FDR = 0.05
  #pl <- pl + geom_hline(yintercept=-1.0*log(0.05),color="black",linetype="dashed",size=1)
  # change x and y labels
  pl <- pl + labs(x = "log(fold change)") + labs(y = "-log(p)")
  # font sizes
  pl <- pl + theme(axis.title.x = element_text(size=18),axis.title.y = element_text(size=18))  # axis labels
  pl <- pl + theme(axis.text.x = element_text(size=14),axis.text.y = element_text(size=14)) # tick labels
  # legend stuff
  pl <- pl + theme(legend.title = element_text(size=16),legend.text=element_text(size=12),legend.position="top")
  # override alpha and marker size for legend
  pl <- pl + guides(colour = guide_legend(override.aes = list(alpha=1,size=6)))
  print(pl)
  dev.off()
}
