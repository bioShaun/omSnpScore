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

p <- add_argument(
    p, "--web", 
    help = "if plot for web server, do not output pdf file to save time", 
    flag=TRUE)  
argv <- parse_args(p)

# var_table <- 'genome.window.w1000000.s500000.bed.var.score.txt'
# output_prefix <- 'test'
# plot_type <- 'var'

var_table <- argv$input
output_prefix <- argv$output
plot_type <- argv$plot_type
plot_title <- argv$title
chr.size <- argv$chr_size
is_web <- argv$web

wheat_cols <- c("#377EB8", "#4DAF4A", "#FF7F00")

snp_index_plot <- function(each_chr) {
  
  chrom_df_1 <- filter(m_snp_index_df, Chrom == each_chr & variable == 'mutant')
  p1 <- ggplot(chrom_df_1, aes(start_m, value, color=variable)) +
    geom_point(alpha=0.5) + geom_smooth(color=line_col,method = "loess",alpha = 0.5,span = 0.2) +
    theme_onmath() +
    scale_x_continuous(limits = c(0, chrom_len_max),
                       breaks = seq(0,chrom_len_max,chrom_len_unit)) +
    scale_color_manual(values = point_cols) +
    guides(color=F) + xlab("") + ylab("") + ggtitle('mutant')
  
  chrom_df_2 <- filter(m_snp_index_df, Chrom == each_chr & variable == 'wild')
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
  if (! is_web ) {
  ggsave(filename = paste(plot_name, 'pdf', sep = '.'),
         plot=p,
         width = 12,
         height = 10)

  }
}


ed_plot <- function(plot.df, breaks, cutoff, labelpos) {
        par(mfrow=c(2,1))
        plot(1:nrow(plot.df), plot.df$fitted, type='l', 
            ylim=c(min(plot.df$fitted, na.rm=T),
            1.1*max(max(plot.df$fitted, na.rm=T),cutoff)), 
            ylab=substitute("ED"^p~ ~"(Loess fit)", list(p=4)), 
            xaxt='n', xaxs='i', xlab="Chromosome", cex=.6, cex.lab=.8, cex.axis=.8,
            main=plot_title)
        abline(v=(breaks[1:length(breaks)-1]+2), col="grey")
        abline(h=cutoff, col='red', lty=2)
        mtext(unique(plot.df$CHROM[!is.na(plot.df$CHROM)]), at = labelpos, side=1, cex=.5)
        plot(1:nrow(plot.df), plot.df$unfitted, pch=16, cex=.6, 
            col="#999999AA", ylim=c(min(plot.df$unfitted, na.rm=T),
            1.1*max(plot.df$unfitted, na.rm=T)), 
            ylab=substitute("ED"^p, list(p=power)), 
            xaxt='n', xaxs='i', xlab="Chromosome", 
            cex.lab=.8, cex.axis=.8)
        abline(v=(breaks[1:length(breaks)-1]+2), col="grey")
        mtext(unique(plot.df$CHROM[!is.na(plot.df$CHROM)]), at = labelpos, side=1, cex=.5)
    }


qtlseqr_data <- function(field) {
    var_table_df <- data.frame(var_table_df)
    var_table_df$SNP <- paste(var_table_df$CHROM, var_table_df$POS, sep = ':')    
    plot_data <- var_table_df[, c('SNP', 'CHROM', 'POS', field)]
    if (!(is.na(chr.size))) {
        chr.size.df <- read.delim(chr.size, header=F, col.names=c('CHROM', 'POS'))
        chr.size.df$SNP <- paste(chr.size.df$CHROM, chr.size.df$POS, sep = ':')
        chr.size.df[, field] <- NA
        chr.size.df <- chr.size.df[, c('SNP', 'CHROM', 'POS', field)]
        plot_data <- rbind(plot_data, chr.size.df)
        plot_data <- arrange(plot_data, CHROM, POS)
    }
    plot_data$CHROM <- str_remove(plot_data$CHROM, fixed('chr', ignore_case = T))
    plot_data <- filter(plot_data, CHROM != "Un")
    return(plot_data)
}


qtlseqr_plot <- function(field, ylab, suffix, fdrT=0.05) {

    if (field %in% colnames(var_table_df)) {
        plot_data <- qtlseqr_data(field)
        output_prefix = paste(output_prefix, suffix, sep='.')
        omsCMplot(plot_data,plot.type="m",LOG10=F,col = wheat_cols,
                    chr.den.col=NULL,file="jpg",memo="test",dpi=300,ylab = ylab,
                    out.name = output_prefix, cex.axis = 0.8, plot.title=plot_title,
                    amplify=F, signal.col=NULL,cex=0.6)
        if ( ! is_web ) {
            omsCMplot(plot_data,plot.type="m",LOG10=F,col = wheat_cols,
                        chr.den.col=NULL,file="pdf",memo="test",dpi=300,ylab = ylab, out.name = output_prefix, cex.axis = 0.8, plot.title=plot_title)
        }   
    }

}

qtlseqr_snp_index_plot <- function() {
    var_table_df <- filter(var_table_df, CHROM != "chrUn")
    var_table_df$POS <- var_table_df$POS / 1e6
    if ('chr1A' %in% var_table_df$CHROM) {

        var_table_df$chrom_genome <- str_extract(var_table_df$CHROM, "chr\\d")
        var_table_df$chrom_num <- str_remove(var_table_df$CHROM, "chr\\d")
        snp_index_df <- melt(var_table_df, id.vars = c("chrom_genome", "chrom_num", "POS"), 
                            measure.vars = c("tricubeDeltaSNP"))
        plot_width = 14
        plot_height = 6
    } else {
        snp_index_df <- var_table_df
        snp_index_df <- melt(var_table_df, id.vars = c("CHROM", "POS"), 
                            measure.vars = c("tricubeDeltaSNP"))
        chrom_num <- length(unique(snp_index_df$CHROM))
        plot_width <- ceiling(sqrt(chrom_num)) * 3
        plot_height = ceiling(sqrt(chrom_num)) * 1.5   
    }
    
    snp_index_col <- brewer.pal(3, 'Set1')
    names(snp_index_col) <- c('tricubeDeltaSNP', "CI_95", "CI_99")
    if ('chr1A' %in% var_table_df$CHROM) {
        ints_df <-
            dplyr::select(var_table_df, chrom_genome, chrom_num, POS, dplyr::matches("CI_")) %>% tidyr::gather(key = "Interval", value = "value",-chrom_genome, -chrom_num,-POS)
    } else {
        ints_df <-
            dplyr::select(var_table_df, CHROM, POS, dplyr::matches("CI_")) %>% tidyr::gather(key = "Interval", value = "value",-CHROM,-POS) 
    }

    p <- ggplot(snp_index_df, aes(POS, value, color=variable)) +
    geom_line() + 
    geom_line(data = ints_df, aes(POS, value, color=Interval)) +
    geom_line(data = ints_df, aes(POS, -value, color=Interval)) +
    theme_bw() + theme(strip.text.y = element_text(angle = 0)) +
    scale_color_manual(values = snp_index_col) +
    xlab("Genomic Position (Mb)") + ylab(expression(Delta * '(SNP-index)')) +
    guides(color=guide_legend(title = "")) + ggtitle(plot_title)

    if ('chr1A' %in% var_table_df$CHROM) {
      p <- p + facet_grid(chrom_num~chrom_genome, scales = "free_x")
    } else {
      p <- p + facet_wrap(.~CHROM, scales = "free_x")
    }
    output_prefix = paste(output_prefix, 'snpIndex.plot', sep='.')
    ggsave(filename = paste(output_prefix, 'jpg', sep = '.'), plot=p, width = plot_width, height = plot_height)
    if (! is_web) {
        ggsave(filename = paste(output_prefix, 'pdf', sep = '.'), plot=p, width = plot_width, height = plot_height)
    }


}

replace_zero <- function(num_count) {
  if (num_count == 0) {
    return(-Inf)
  } else {
    return(num_count)
  }
}

var_table_df <- fread(var_table)
if (!(is.na(chr.size))) {
    chr.size.df <- read.delim(chr.size, header=F, col.names=c('CHROM', 'Start'))
    if ('CHROM' %in% colnames(var_table_df)) {
        var_table_df$CHROM <- factor(var_table_df$CHROM, levels = chr.size.df$CHROM)
    }
}

if (plot_type == 'density') {
  var_table_df$SNP <- paste(var_table_df$CHROM, var_table_df$POS, sep = ':')
  plot_data <- var_table_df[, c('SNP', 'CHROM', 'POS')]
  plot_data <- arrange(plot_data, CHROM, POS)
  plot_data$CHROM <- str_remove(plot_data$CHROM, fixed('chr', ignore_case = T))
  plot_data <- filter(plot_data, CHROM != "Un")
  omsCMplot(plot_data,plot.type="d",bin.size=1e6,
         col=c("darkgreen", "yellow", "red"),
         file="jpg", dpi=300, out.name = output_prefix,
         plot.title=plot_title, chr.size=chr.size)
  if ( ! is_web ) {
    omsCMplot(plot_data,plot.type="d",bin.size=1e6,
            col=c("darkgreen", "yellow", "red"),
            file="pdf", dpi=300, out.name = output_prefix,
            plot.title=plot_title, chr.size=chr.size) 
  }

} else if (plot_type == 'snp_index_old') {
  
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
} else if (plot_type == 'Gprime'){
    qtlseqr_plot('Gprime', "G' value", 'Gprime.plot')
} else if (plot_type == 'snpIndex'){
    qtlseqr_snp_index_plot()
} else if (plot_type == 'ED-old'){
    if ('fitted' %in% colnames(var_table_df))  {
        var_table_df <- data.frame(var_table_df)
        var_table_df <- filter(var_table_df, CHROM != "chrUn")
        breaks <- NULL
        for (row in 1:(nrow(var_table_df)-1)){
            if (var_table_df$CHROM[row] != var_table_df$CHROM[row+1]) {
                breaks <- append(breaks, row+1)
            }
        }
        breaks <- append(breaks, nrow(var_table_df))
        lastbreak <- 1
        labelpos <- NULL
        plot.df <- data.frame(cbind(NA, NA, NA, NA, NA))
        colnames(plot.df) <- c('CHROM','POS','euc','fitted','unfitted')    
        for( b in breaks){
            plot.df <- rbind(plot.df, var_table_df[lastbreak:b,c('CHROM','POS','euc','fitted','unfitted')],c(NA, NA, NA, NA, NA))
            labelpos <- append(labelpos, ((b+1)-lastbreak)/2+lastbreak)
            lastbreak <- b
        }    
        plot.df$CHROM <- str_remove(plot.df$CHROM, fixed('chr', ignore_case = T))
        cutoff <- 3*(sd(var_table_df$fitted)+median(var_table_df$fitted))        
        output_prefix = paste(output_prefix, 'ED.plot', sep=".")
        save_general_plot(ed_plot(plot.df, breaks, cutoff, labelpos), 
        output_prefix, plot_type = 'png', width = 10)
        if (! is_web)  {
            save_general_plot(ed_plot(plot.df, breaks, cutoff, labelpos), 
            out_prefix, plot_type = 'pdf', width = 10)
        }
    }
} else if (plot_type == 'ED'){
    field <- c('fitted', 'unfitted')
    plot_data <- qtlseqr_data(field)
    ylab <- substitute("ED"^p, list(p=4))
    omsCMplot(plot_data,plot.type="m", multracks=TRUE, LOG10=F,col = wheat_cols,
                chr.den.col=NULL,file="jpg",memo="",dpi=300,ylab = ylab,
                out.name = output_prefix, cex.axis = 0.8, plot.title=plot_title,
                amplify=F, signal.col=NULL,cex=0.6)    
} else if (plot_type == 'var'){
  var_table_df$Start <- var_table_df$Start + 1
  var_table_df$SNP <- paste(var_table_df$CHROM, var_table_df$Start, sep = ':')
  plot_data <- var_table_df[, c('SNP', 'CHROM', 'Start', 'varBScore')]
  if (!(is.na(chr.size))) {
    chr.size.df <- read.delim(chr.size, header=F, col.names=c('CHROM', 'Start'))
    chr.size.df$SNP <- paste(chr.size.df$CHROM, chr.size.df$Start, sep = ':')
    chr.size.df$varBScore <- NA
    chr.size.df <- chr.size.df[, c('SNP', 'CHROM', 'Start', 'varBScore')]
    plot_data <- rbind(plot_data, chr.size.df)
    plot_data <- arrange(plot_data, CHROM, Start)
  }
  plot_data$CHROM <- str_remove(plot_data$CHROM, fixed('chr', ignore_case = T))
  plot_data <- filter(plot_data, CHROM != "Un")
  omsCMplot(plot_data,plot.type="m",LOG10=F,threshold=NULL, col = wheat_cols,
             chr.den.col=NULL,file="jpg",memo="test",dpi=300,ylab = "Score",
             out.name = output_prefix, cex.axis = 0.8, plot.title=plot_title, )
  if ( ! is_web ) {
    omsCMplot(plot_data,plot.type="m",LOG10=F,threshold=NULL, col = wheat_cols,
                chr.den.col=NULL,file="pdf",memo="test",dpi=300,ylab = "Score",
                out.name = output_prefix, cex.axis = 0.8, plot.title=plot_title)
  }

}  else if (plot_type == 'density-new') {
    var_table_df$chrom <- str_remove(var_table_df$chrom, fixed('chr', ignore_case = T))
    chr.size.df$CHROM <- str_remove(chr.size.df$CHROM, fixed('chr', ignore_case = T))
    var_table_df$chrom <- factor(var_table_df$chrom, levels = chr.size.df$CHROM)
    plot_data <- filter(var_table_df, chrom != "Un")

    chrom_list <- as.character(unique(plot_data$chrom))
    chrom_pos <- 2 *seq(length(chrom_list))
    plot_data$chrom_pos <- chrom_pos[plot_data$chrom]
    plot_data$variantCount2 <- sapply(plot_data$variantCount, replace_zero)

    mega_base = floor(log10(max(chr.size.df$Start)))
    mega_unit = str_pad('M', mega_base - 6 + 1,  pad = "0")
    max_chr_len = ceiling(max(chr.size.df$Start)/(10^mega_base))

    cor_plot_col <- colorRampPalette(brewer.pal(11, 'RdYlGn'))(100)

    p <- ggplot(plot_data) + 
        geom_rect(aes(xmin=start, xmax=end, ymin=0, 
                        ymax=1, fill = variantCount2)) +
        facet_grid(chrom~.) +
        theme(strip.text.y = element_text(size=rel(.8), 
                                            face="bold",
                                            angle = 0),
                axis.text.y = element_blank(),
                axis.ticks = element_blank(),
                panel.background = element_rect(fill = "white"),
                plot.title = element_text(hjust = 0.5, size=rel(.8))
                ) +
        scale_x_continuous(limits = c(0, max_chr_len * 10^mega_base), 
                            breaks = seq(0, max_chr_len * 10^mega_base, 10^mega_base),
                            labels = c(0, paste0(seq(1, max_chr_len), mega_unit, sep = ""))) +
        scale_fill_gradientn(colours = rev(cor_plot_col), na.value = "grey80") +
        guides(fill=guide_colourbar(title='Number of SNP')) +
        xlab('Chromosome Length') + ylab('') +
        ggtitle(plot_title)

    chrom_num <- length(chrom_list)
    p_height = 2 * chrom_num / 7
    ggsave(paste(output_prefix, 'jpg', sep='.'), 
    plot = p, width = 10, height = p_height,
    dpi = 300)
    if ( ! is_web ) {
        ggsave(paste(output_prefix, 'pdf', sep='.'), 
        plot = p, width = 10, height = p_height)    
    }

} else {
    print('unsuported plot type.')
}



