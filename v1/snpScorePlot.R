library(argparser)
library(data.table)
library(stringr)
library(omsCMplot)


options(stringsAsFactors = F)
p <- arg_parser("snp varscore plot")
p <- add_argument(
  p, '--input', 
  help = 'input table')
p <- add_argument(
  p, '--output', 
  help = 'output prefix')
argv <- parse_args(p)

var_table <- argv$input
output_prefix <- argv$output
var_table_df <- fread(var_table)


if (! 'varscore' %in% colnames(var_table_df)) {
  omsCMplot(var_table_df,plot.type="d",bin.size=1e6,
         col=c("darkgreen", "yellow", "red"),
         file="jpg", dpi=300, out.name = output_prefix)
  omsCMplot(var_table_df,plot.type="d",bin.size=1e6,
         col=c("darkgreen", "yellow", "red"),
         file="pdf", dpi=300, out.name = output_prefix) 
} else {
  var_table_df$Start <- var_table_df$Start + 1
  var_table_df$SNP <- paste(var_table_df$Chrom, var_table_df$Start, sep = ':')
  var_table_df$Chrom <-  str_replace(var_table_df$Chrom, 'chr', '')
  plot_data <- var_table_df[, c('SNP', 'Chrom', 'Start', 'varscore')]
  omsCMplot(plot_data,plot.type="m",LOG10=F,threshold=NULL,
             chr.den.col=NULL,file="jpg",memo="test",dpi=300,ylab = "Score",
             out.name = output_prefix, cex.axis = 0.8)
  omsCMplot(plot_data,plot.type="m",LOG10=F,threshold=NULL,
             chr.den.col=NULL,file="pdf",memo="test",dpi=300,ylab = "Score",
             out.name = output_prefix, cex.axis = 0.8)
}



