suppressMessages(library(QTLseqr))
suppressMessages(library(argparser))
suppressMessages(library(stringr))

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
    p, '--out_prefix',
    help = 'output prefix'
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
argv <- parse_args(p)

input_table <- argv$input
window <- argv$window
high_bulk <- str_split(argv$high_bulk, ',')[[1]]
low_bulk <- str_split(argv$low_bulk, ',')[[1]]
ref_freq <- argv$ref_freq
min_sample_dp <- argv$min_sample_dp
pop_stru <- argv$pop_stru
out_prefix <- argv$out_prefix

df <- importFromTable(file=input_table, highBulk = high_bulk, lowBulk = low_bulk)
df_filt <- filterSNPs(SNPset=df, refAlleleFreq = ref_freq, minSampleDepth = min_sample_dp)
df_filt <- runGprimeAnalysis(SNPset = df_filt, windowSize = window, outlierFilter = "deltaSNP")
df_filt <- runQTLseqAnalysis(df_filt, windowSize = window, popStruc=pop_stru, bulkSize=50)

chrom_num <- length(unique(df_filt$CHROM))

png(paste(out_prefix, 'deltaSNP.png', sep='.'), width=4 * chrom_num, height=6, res=300, unit='in')
plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals = TRUE)
dev.off()

pdf(paste(out_prefix, 'deltaSNP.pdf', sep='.'), width=4 * chrom_num, height=6)
plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals = TRUE)
dev.off()

png(paste(out_prefix, 'Gprime.png', sep='.'), width=4 * chrom_num, height=6, res=300, unit='in')
plotQTLStats(SNPset = df_filt, var = "Gprime", plotIntervals = TRUE, q = 0.01)
dev.off()

pdf(paste(out_prefix, 'Gprime.pdf', sep='.'), width=4 * chrom_num, height=6)
plotQTLStats(SNPset = df_filt, var = "Gprime", plotIntervals = TRUE, q = 0.01)
dev.off()

getQTLTable(
    SNPset = df_filt,
    alpha = 0.01,
    export = TRUE,
    fileName = paste(out_prefix, "QTLseqr.csv", sep='.')
)