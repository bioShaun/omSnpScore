suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(argparser))

options(stringsAsFactors = F)

p <- arg_parser("for var compare heatmap")
p <- add_argument(p, "--var_density_file", help = "var density by window")
p <- add_argument(p, "--out_prefix", help = "output plot prefix")
argv <- parse_args(p)


var_density_file <- argv$var_density_file
out_prefix <- argv$out_prefix

# var_density_file <- './TCE000007.K1-1_TCE000008.R1-1-depth_1-window_1.0M.csv'
# out_prefix <- 'test-plot'

var_density_df <- read.csv(var_density_file, check.names=F)
max_chr_len = ceiling(max(var_density_df$end)/(10^8))

var_cols = colnames(var_density_df)
melt_cols <- var_cols[1:3]
sample_names <- var_cols[4:length(var_cols)]
sample_num <- length(sample_names)

cor_plot_col <- colorRampPalette(brewer.pal(9, 'Reds'))(100)
m_var_density_df <- melt(var_density_df, id.vars = melt_cols, variable.name = 'sample_id', value.name = 'variant_count')

plot_label <- seq(1,3*sample_num -2,3)
plot_pos <- seq(sample_num)

names(plot_label) <- sample_names
names(plot_pos) <- sample_names
names(sample_names) <- plot_label

m_var_density_df$y_pos <- plot_pos[m_var_density_df$sample_id]
m_var_density_df$log_depth <- log2(m_var_density_df$variant_count)

p <- ggplot(m_var_density_df) + 
  geom_rect(aes(xmin=start, xmax=end, ymin=y_pos*3-3, 
                ymax=y_pos*3-1, fill = log_depth)) +
  facet_wrap(chrom~., nrow = 1, strip.position="bottom") +
  theme(strip.text.y = element_text(size=rel(.8), 
                                    face="bold",
                                    angle = 0),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, hjust =0, size=rel(0.75))) +
  scale_fill_gradientn(colours = cor_plot_col) +
  guides(fill=guide_colourbar(title='Log2(SNP) per 1Mb')) +
  xlab('') + ylab('') +
  coord_flip()+ scale_x_reverse() +
  scale_y_continuous(labels=sample_names, breaks = plot_label, position = 'right')

if (sample_num == 1) {
    p <- p + theme(axis.text.x =  element_blank())
}

chrom_list <- as.character(unique(m_var_density_df$chrom))
chrom_num <- length(chrom_list)
sample_num_scale = ifelse(sample_num < 2, 2, sample_num)
p_width = 2 * chrom_num * sample_num_scale / 7
if (p_width <= 50){
    ggsave(paste(out_prefix, 'png', sep='.'),
           plot = p, width = p_width, height = 12,
           dpi = 300, type = "cairo")
}
ggsave(paste(out_prefix, 'pdf', sep='.'), 
       plot = p, width = p_width, height = 12, limitsize=F)
