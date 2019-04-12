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

p <- arg_parser("Plot snp index")
p <- add_argument(
  p, '--snp_index_file',
  help = 'snp-index')
p <- add_argument(
  p, '--plot_dir',
  help = 'directory store plot data.')
argv <- parse_args(p)

# snp_index_file <- 'snp.all.genome_window.snp_index.score.txt'
# plot_dir <- 'snp_index'

snp_index_file <- argv$snp_index_file
plot_dir <- argv$plot_dir

if ( ! dir.exists(plot_dir)) {
  dir.create(plot_dir)
}


snp_index_df <- fread(snp_index_file)
add_col <- c('mutant_bulk', 'wild_bulk') 
m_snp_index_df <- melt(snp_index_df, id.vars = c('Chrom', 'Start', 'End', 'Label'))
m_snp_index_df$start_m <- m_snp_index_df$Start / 1000000


chroms <- unique(m_snp_index_df$Chrom)
set1_cols <- brewer.pal(9, 'Set1')
point_cols <- c(set1_cols[5], set1_cols[2], set1_cols[3])
line_col <- set1_cols[1]

chrom_len_unit <- 10 ^ floor(log10(max(m_snp_index_df$start_m)))
chrom_len_max <- ceiling(max(m_snp_index_df$start_m) / chrom_len_unit) * chrom_len_unit

snp_plot_fun <- function(each_chr) {

  chrom_df_1 <- filter(m_snp_index_df, Chrom == each_chr & variable == add_col[1])
  p1 <- ggplot(chrom_df_1, aes(start_m, value, color=variable)) +
    geom_point(alpha=0.5) + geom_smooth(color=line_col,method = "loess",alpha = 0.5,span = 0.2) +
    theme_onmath() +
    scale_x_continuous(limits = c(0, chrom_len_max),
      breaks = seq(0,chrom_len_max,chrom_len_unit)) +
    scale_color_manual(values = point_cols) +
    guides(color=F) + xlab("") + ylab("") + ggtitle(add_col[1])

  chrom_df_2 <- filter(m_snp_index_df, Chrom == each_chr & variable == add_col[2])
  p2 <- ggplot(chrom_df_2, aes(start_m, value, color=variable)) +
    geom_point(alpha=0.5,colour=point_cols[2]) + geom_smooth(color=line_col,method = "loess",alpha = 0.5,span = 0.2) +
    theme_onmath() +
    scale_x_continuous(limits = c(0, chrom_len_max),
                       breaks = seq(0,chrom_len_max,chrom_len_unit)) +
    scale_color_manual(values = point_cols) +
    guides(color=F) + xlab("") + ylab("") + ggtitle(add_col[2])

  chrom_df_3 <- filter(m_snp_index_df, Chrom == each_chr & variable == 'varscore')
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


  plot_name <- file.path(plot_dir, each_chr)
  ggsave(filename = paste(plot_name, 'png', sep = '.'),
    plot=p,
    width = 8,
    height = 6)
  ggsave(filename = paste(plot_name, 'pdf', sep = '.'),
    plot=p,
    width = 12,
    height = 10)
}

lapply(chroms, snp_plot_fun)
