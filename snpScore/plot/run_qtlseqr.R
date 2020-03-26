suppressMessages(library(QTLseqr))
suppressMessages(library(argparser))
suppressMessages(library(stringr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(stringr))
suppressMessages(library(fANCOVA))
suppressMessages(library(omsCMplot))
suppressMessages(library(omplotr))

options(stringsAsFactors = F)

p <- arg_parser("run QTLseqr analysis")
p <- add_argument(
    p, '--input', 
    help = 'input table'
)
p <- add_argument(
    p, '--high_bulk', 
    help = 'The sample name of the High Bulk.'
)
p <- add_argument(
    p, '--low_bulk', 
    help = 'The sample name of the Low Bulk.'
)
p <- add_argument(
    p, '--out_dir',
    help = 'output directory.'
)
p <- add_argument(
    p, '--window', 
    help = 'snp window size for calculation.',
    default=1e7
)
p <- add_argument(
    p, '--ref_freq',
    help = 'filter out SNPs with a Reference Allele Frequency less than refAlleleFreq and greater than 1 - refAlleleFreq',
    default=0.3
)
p <- add_argument(
    p, '--min_sample_dp',
    help = 'The minimum read depth for a SNP in each bulk',
    default = 5
)
p <- add_argument(
    p, '--pop_stru',
    help = 'the population structure. Defaults to "F2" and assumes "RIL" otherwise.',
    default = 'RIL'
)
p <- add_argument(
    p, '--qtlseqr',
    help = 'Run qtlseqr analysis.',
    flag = TRUE
)
p <- add_argument(
    p, '--ed',
    help = 'Run snp ED analysis.',
    flag = TRUE
)
p <- add_argument(
    p, '--ed_plot',
    help = 'Run snp ED plot analysis, for test.',
    flag = TRUE
)
argv <- parse_args(p)

input_table <- argv$input
window <- argv$window
high_bulk <- str_split(argv$high_bulk, ',')[[1]]
low_bulk <- str_split(argv$low_bulk, ',')[[1]]
ref_freq <- argv$ref_freq
min_sample_dp <- argv$min_sample_dp
pop_stru <- argv$pop_stru
out_dir <- argv$out_dir
qtlseqr_flag <- argv$qtlseqr
ed_flag <- argv$ed
ed_plot_flag <- argv$ed_plot

plot_cols <- brewer.pal(9, 'Set1')
wheat_cols <- c(plot_cols[2], plot_cols[3], plot_cols[5])




df <- importFromTable(file=input_table, highBulk = high_bulk, lowBulk = low_bulk)
df_filt <- filterSNPs(SNPset=df, refAlleleFreq = ref_freq, minSampleDepth = min_sample_dp)
chrom_num <- length(unique(df_filt$CHROM))

table_name <- c()

if (qtlseqr_flag) {
    filter_stat = 1 - 2*ref_freq
    table_name <- c('qtlseqr', as.integer(window), pop_stru)
    outname <- paste('qtlSeqr', as.integer(window), ref_freq, min_sample_dp, pop_stru, sep='.')
    out_prefix <- file.path(out_dir, outname)
    df_filt <- runGprimeAnalysis(SNPset = df_filt, windowSize = window, outlierFilter = "deltaSNP")
    df_filt <- runQTLseqAnalysis(df_filt, windowSize = window, popStruc=pop_stru, bulkSize=50, filter=filter_stat)
    png(paste(out_prefix, 'deltaSNP.png', sep='.'), width=4 * chrom_num, height=6, res=300, unit='in')
    print(plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals = TRUE))
    dev.off()

    pdf(paste(out_prefix, 'deltaSNP.pdf', sep='.'), width=4 * chrom_num, height=6)
    print(plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals = TRUE))
    dev.off()

    png(paste(out_prefix, 'Gprime.png', sep='.'), width=4 * chrom_num, height=6, res=300, unit='in')
    print(plotQTLStats(SNPset = df_filt, var = "Gprime", plotIntervals = TRUE, q = 0.01))
    dev.off()

    pdf(paste(out_prefix, 'Gprime.pdf', sep='.'), width=4 * chrom_num, height=6)
    print(plotQTLStats(SNPset = df_filt, var = "Gprime", plotIntervals = TRUE, q = 0.01))
    dev.off()

    res <- try(getQTLTable(
        SNPset = df_filt,
        alpha = 0.1,
        export = TRUE,
        fileName = paste(out_prefix, "QTLseqr.csv", sep='.')
    ))
    if(inherits(res, "try-error"))
    {
       print("Can not find significant region using QTLseqr!")
    }
} 

if (ed_flag) {
    table_name <- c(table_name, 'ed')
    eu_power <- 4
    df_filt$euc<-sqrt(2 * (df_filt$SNPindex.LOW-df_filt$SNPindex.HIGH)^2)
    chrs <- unique(df_filt$CHROM)
    cat("Processing each chromosome:\nChr\tSnps\tSpan\n")
    for(chr in chrs) {
        e = df_filt$euc[df_filt$CHROM == chr]^eu_power
        p = df_filt$POS[df_filt$CHROM == chr]
        if(as.integer(length(p)) < 50){
            cat(chr, " has less than 50 snps (n=", length(p), ")\n", sep="")
            next
        }
        lo <- loess.as(p, e, degree = 1, family = 'symmetric', criterion = 'aicc')
        df_filt$fitted[df_filt$CHROM==chr] <- lo$fitted
        df_filt$unfitted[df_filt$CHROM==chr] <- lo$y
        usespan <- lo$pars$span
        cat(chr, length(p), round(usespan[1], digits=3), "\n", sep="\t")
    }    
    breaks <- NULL
    for (row in 1:(nrow(df_filt)-1)){
        if (df_filt$CHROM[row] != df_filt$CHROM[row+1]) {
            breaks <- append(breaks, row+1)
        }
    }
    breaks <- append(breaks, nrow(df_filt))
    lastbreak <- 1
    labelpos <- NULL
    plot.df <- data.frame(cbind(NA, NA, NA, NA, NA))
    colnames(plot.df) <- c('CHROM','POS','euc','fitted','unfitted')    
    for( b in breaks){
        plot.df <- rbind(plot.df, df_filt[lastbreak:b,c('CHROM','POS','euc','fitted','unfitted')],c(NA, NA, NA, NA, NA))
        labelpos <- append(labelpos, ((b+1)-lastbreak)/2+lastbreak)
        lastbreak <- b
    }
    plot.df$CHROM <- str_remove(plot.df$CHROM, fixed('chr', ignore_case = T))
    cutoff <- 3*(sd(df_filt$fitted)+median(df_filt$fitted))
    df_filt$dis2edcutoff <- df_filt$fitted - cutoff
    ed_plot <- function() {
        par(mfrow=c(2,1))
        plot(1:nrow(plot.df), plot.df$fitted, type='l', 
            ylim=c(min(plot.df$fitted, na.rm=T),
            1.1*max(plot.df$fitted, na.rm=T)), 
            ylab=substitute("ED"^p~ ~"(Loess fit)", list(p=eu_power)), 
            xaxt='n', xaxs='i', xlab="Chromosome", cex=.6, cex.lab=.8, cex.axis=.8)
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
    outname <- paste('ED', ref_freq, min_sample_dp, sep='.')
    out_prefix <- file.path(out_dir, outname)
    plot_width <- length(chrs) * 0.1 + 8
    save_general_plot(ed_plot(), out_prefix, plot_type = 'png', width = plot_width)
    save_general_plot(ed_plot(), out_prefix, plot_type = 'pdf', width = plot_width)
}
table_name <- c(table_name, ref_freq, min_sample_dp, 'csv')
table_name_flat <- paste(table_name, collapse = '.')
table_prefix <- file.path(out_dir, table_name_flat)
write.csv(df_filt, file = table_prefix, quote=F, row.names=F)
