suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(omsCMplot))
suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(scales))
suppressMessages(library(RColorBrewer))
suppressMessages(library(argparser))
suppressMessages(library(gridExtra))
suppressMessages(library(omplotr))


options(stringsAsFactors = F)

p <- arg_parser("snp score plot")
p <- add_argument(
  p, '--input', 
  help = 'input table')
p <- add_argument(
  p, '--plot_type', 
  help = 'plot type.')
p <- add_argument(
  p, '--output', 
  help = 'output prefix')
p <- add_argument(
  p, '--title',
  help = 'plot title',
  default="")
p <- add_argument(
  p, '--chr_size',
  help = 'chr size file',
  default=NULL)
argv <- parse_args(p)

# var_table <- 'genome.window.w1000000.s500000.bed.var.score.txt'
# output_prefix <- 'test'
# plot_type <- 'var'

var_table <- argv$input
output_prefix <- argv$output
plot_type <- argv$plot_type
plot_title <- argv$title
chr.size <- argv$chr_size


wheat_cols <- c("#377EB8", "#4DAF4A", "#FF7F00")

snp_index_plot <- function(each_chr) {
  
  chrom_df_1 <- filter(m_snp_index_df, Chrom == each_chr & variable == 'mutant.FREQ')
  p1 <- ggplot(chrom_df_1, aes(start_m, value, color=variable)) +
    geom_point(alpha=0.5) + geom_smooth(color=line_col,method = "loess",alpha = 0.5,span = 0.2) +
    theme_onmath() +
    scale_x_continuous(limits = c(0, chrom_len_max),
                       breaks = seq(0,chrom_len_max,chrom_len_unit)) +
    scale_color_manual(values = point_cols) +
    guides(color=F) + xlab("") + ylab("") + ggtitle('mutant')
  
  chrom_df_2 <- filter(m_snp_index_df, Chrom == each_chr & variable == 'wild.FREQ')
  p2 <- ggplot(chrom_df_2, aes(start_m, value, color=variable)) +
    geom_point(alpha=0.5,colour=point_cols[2]) + geom_smooth(color=line_col,method = "loess",alpha = 0.5,span = 0.2) +
    theme_onmath() +
    scale_x_continuous(limits = c(0, chrom_len_max),
                       breaks = seq(0,chrom_len_max,chrom_len_unit)) +
    scale_color_manual(values = point_cols) +
    guides(color=F) + xlab("") + ylab("") + ggtitle('wild')
  
  chrom_df_3 <- filter(m_snp_index_df, Chrom == each_chr & variable == 'snp_score')
  p3 <- ggplot(chrom_df_3, aes(start_m, value, color=variable)) +
    geom_point(alpha=0.5,colour=point_cols[3]) + geom_smooth(color=line_col,method = "loess",alpha = 0.5,span = 0.2) +
    theme_onmath() +
    scale_x_continuous(limits = c(0, chrom_len_max),
                       breaks = seq(0,chrom_len_max,chrom_len_unit)) +
    scale_y_continuous(limits = c(-1, 1),
                       breaks = seq(-1,1,0.5)) +
    scale_color_manual(values = point_cols) +
    guides(color=F) + xlab(paste(each_chr,"(MB)")) + ylab("") + ggtitle('Diff')
  
  p <- grid.arrange(p1, p2, p3, nrow = 3,ncol = 1)
  
  
  plot_name <- file.path(output_prefix, each_chr)
  ggsave(filename = paste(plot_name, 'png', sep = '.'),
         plot=p,
         width = 8,
         height = 6)
  ggsave(filename = paste(plot_name, 'pdf', sep = '.'),
         plot=p,
         width = 12,
         height = 10)
}


var_table_df <- fread(var_table)

if (plot_type == 'density') {
  var_table_df$SNP <- paste(var_table_df$Chr, var_table_df$Pos, sep = ':')
  plot_data <- var_table_df[, c('SNP', 'Chr', 'Pos')]
  omsCMplot(plot_data,plot.type="d",bin.size=1e6,
         col=c("darkgreen", "yellow", "red"),
         file="jpg", dpi=300, out.name = output_prefix,
         plot.title=plot_title, chr.size=chr.size)
  omsCMplot(plot_data,plot.type="d",bin.size=1e6,
         col=c("darkgreen", "yellow", "red"),
         file="pdf", dpi=300, out.name = output_prefix,
         plot.title=plot_title, chr.size=chr.size) 
} else if (plot_type == 'snp_index') {
  
  if ( ! dir.exists(output_prefix)) {
    dir.create(output_prefix)
  }
  m_snp_index_df <- melt(var_table_df, id.vars = c('Chrom', 'Start', 'End'))
  m_snp_index_df$start_m <- m_snp_index_df$Start / 1000000
  chroms <- unique(m_snp_index_df$Chrom)
  set1_cols <- brewer.pal(9, 'Set1')
  point_cols <- c(set1_cols[5], set1_cols[2], set1_cols[3])
  line_col <- set1_cols[1]
  
  chrom_len_unit <- 10 ^ floor(log10(max(m_snp_index_df$start_m)))
  chrom_len_max <- ceiling(max(m_snp_index_df$start_m) / chrom_len_unit) * chrom_len_unit
  lapply(chroms, snp_index_plot)
} else {
  var_table_df$Start <- var_table_df$Start + 1
  var_table_df$SNP <- paste(var_table_df$Chrom, var_table_df$Start, sep = ':')
  plot_data <- var_table_df[, c('SNP', 'Chrom', 'Start', 'snp_score')]
  if (!(is.na(chr.size))) {
    chr.size.df <- read.delim(chr.size, header=F, col.names=c('Chrom', 'Start'))
    chr.size.df$SNP <- paste(chr.size.df$Chrom, chr.size.df$Start, sep = ':')
    chr.size.df$snp_score <- 1
    chr.size.df <- chr.size.df[, c('SNP', 'Chrom', 'Start', 'snp_score')]
    plot_data <- rbind(plot_data, chr.size.df)
    plot_data <- arrange(plot_data, Chrom, Start)
  }
  plot_data$Chrom <- str_remove(plot_data$Chrom, fixed('chr', ignore_case = T))
  omsCMplot(plot_data,plot.type="m",LOG10=F,threshold=NULL, col = wheat_cols,
             chr.den.col=NULL,file="jpg",memo="test",dpi=300,ylab = "Score",
             out.name = output_prefix, cex.axis = 0.8, plot.title=plot_title, )
  omsCMplot(plot_data,plot.type="m",LOG10=F,threshold=NULL, col = wheat_cols,
             chr.den.col=NULL,file="pdf",memo="test",dpi=300,ylab = "Score",
             out.name = output_prefix, cex.axis = 0.8, plot.title=plot_title)
}



