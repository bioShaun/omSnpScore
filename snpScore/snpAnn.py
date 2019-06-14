import sys
import pandas as pd
from io import StringIO
from snpScore import core
from loguru import logger
from pathlib import Path, PurePath
from pybedtools import BedTool

COL_HEADER_MAP = {'Chr': '#CHROM', 'Pos': 'POS', 'Alt': 'ALT'}


def make_gene_bed(gene_bed, genes, target_bed):
    gene_bed_df = pd.read_csv(gene_bed,
                              sep='\t',
                              header=None,
                              names=['chrom', 'start', 'end'],
                              index_col=3)
    gene_index = pd.Index(genes)
    intersect_genes = gene_index.intersection(gene_bed_df.index)
    missed_genes = gene_index.difference(gene_bed_df.index)
    if intersect_genes.empty:
        logger.error('None of input genes is in database.')
        sys.exit(1)
    if not missed_genes.empty:
        missed_genes = missed_genes.astype('str')
        logger.warning('Input genes {} not found.'.format(
            ','.join(missed_genes)))
    logger.info('Making target region bed from input genes...')
    target_bed_df = gene_bed_df.loc[intersect_genes]
    target_bed_df.to_csv(target_bed, sep='\t', index=False, header=False)


def make_pos_bed(position, target_bed):
    chrom, start_end = position.split(':')
    start, end = start_end.split('-')
    start = int(start)
    if start > 0:
        start = start - 1
    logger.info('Making target region bed from input position...')
    with open(target_bed, 'w') as bed_inf:
        bed_inf.write(f'{chrom}\t{start}\t{end}\n')


def find_target_region_snp(snp_bed, region_bed):
    snp_bed_obj = BedTool(str(snp_bed))
    region_bed_obj = BedTool(str(region_bed))
    region_snp_obj = snp_bed_obj.intersect(region_bed_obj)
    region_snp_str = StringIO(str(region_snp_obj))
    region_snp_df = pd.read_csv(region_snp_str,
                                sep='\t',
                                header=None,
                                names=['#CHROM', 'Start', 'POS'])
    region_snp_df.drop_duplicates(inplace=True)
    return region_snp_df.loc[:, ['#CHROM', 'POS']]


def merge_snpeff_anno(snp_anno_obj):
    return '|'.join(snp_anno_obj)


def extract_snpeff_anno(anno_line):
    anno_stats = []
    fileds = (1, 3, 6, 9, 10)
    gene_anno = anno_line.split(';')[0]
    anno_line_stats = gene_anno.split(",")
    for annStr in anno_line_stats:
        annDetailArray = annStr.split("|")
        filed_stats = []
        for filled_i in fileds:
            filed_stats.append(annDetailArray[filled_i])
        anno_stats.append(filed_stats)
    zip_anno_stats = list(map(merge_snpeff_anno, zip(*anno_stats)))
    return zip_anno_stats


def split_dataframe_rows(df, column_selectors, row_delimiter):
    # we need to keep track of the ordering of the columns
    def _split_list_to_rows(row, row_accumulator, column_selector,
                            row_delimiter):
        split_rows = {}
        max_split = 0
        for column_selector in column_selectors:
            split_row = row[column_selector].split(row_delimiter)
            split_rows[column_selector] = split_row
            if len(split_row) > max_split:
                max_split = len(split_row)

        for i in range(max_split):
            new_row = row.to_dict()
            for column_selector in column_selectors:
                try:
                    new_row[column_selector] = split_rows[column_selector].pop(
                        0)
                except IndexError:
                    new_row[column_selector] = ''
            row_accumulator.append(new_row)

    new_rows = []
    df.apply(_split_list_to_rows,
             axis=1,
             args=(new_rows, column_selectors, row_delimiter))
    new_df = pd.DataFrame(new_rows, columns=df.columns)
    return new_df


def snp_reads_inf(snp_freq_obj):
    ref_df = snp_freq_obj.grp_dep_df - snp_freq_obj.grp_alt_df
    alt_df = snp_freq_obj.grp_alt_df
    grp_reads_inf_list = list()
    groups = alt_df.columns
    for group_i in groups:
        grp_i_ref = ref_df.loc[:, group_i].astype('int').astype('str')
        grp_i_alt = alt_df.loc[:, group_i].astype('int').astype('str')
        grp_reads_inf = grp_i_ref.str.cat(grp_i_alt, sep=',')
        grp_reads_inf_list.append(grp_reads_inf)
    grp_reads_df = pd.concat(grp_reads_inf_list, axis=1)
    return grp_reads_df


def df2list(df):
    df_list = list()
    if not df.empty:
        df_list.append(list(df.columns))
        for index_i in df.index:
            df_list.append(list(df.loc[index_i]))
    return df_list


def check_df(df, item='SNP'):
    if df.empty:
        logger.warning('{} Not Found!'.format(item))
        sys.exit(1)


def snp_ann_pipe(gene_bed, snp_ann_dir, outdir, genes, position, vcf_table_files,
         group_labels, outfmt='table'):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    # step1 make selected gene/region bedfile
    target_region_bed = outdir / 'target.bed'
    if genes:
        gene_list = genes.split(',')
        make_gene_bed(gene_bed, gene_list, target_region_bed)
    elif position:
        make_pos_bed(position, target_region_bed)
    else:
        raise ValueError('one of genes and postion must specified.')

    # step2 intersect with snp ann bed
    snp_ann_dir = PurePath(snp_ann_dir)
    snp_bed_file = snp_ann_dir / 'snp.ann.table.bed'
    logger.info('Finding target region snp...')
    target_snp_df = find_target_region_snp(snp_bed_file, target_region_bed)
    check_df(target_snp_df, item='Target SNP')

    # step3 extract snp annotation
    logger.info('Loading snp ann db...')
    snp_ann_file = snp_ann_dir / 'snp.ann.table.pkl'
    snp_ann_df = pd.read_pickle(snp_ann_file)
    logger.info('Annotating target snp...')
    target_snp_ann_df = target_snp_df.merge(snp_ann_df)

    # step4 extract annotation table
    logger.info('Extracting snpEff annotationg...')
    snpeff_anno = list(target_snp_ann_df.INFO.map(extract_snpeff_anno))
    snpeff_anno_df = pd.DataFrame(snpeff_anno)
    snpeff_anno_df.columns = ['Feature', 'Gene', 'Transcript', 'Variant_DNA_Level', 'Variant_Protein_Level']
    target_snp_ann_df = pd.concat([target_snp_ann_df, snpeff_anno_df], axis=1)
    target_snp_ann_df.drop('INFO', axis=1, inplace=True)
    flat_target_snp_ann_df = split_dataframe_rows(
        target_snp_ann_df,
        column_selectors=['Feature', 'Gene', 'Transcript', 'Variant_DNA_Level', 'Variant_Protein_Level'],
        row_delimiter='|')

    # step5 filter annotation
    if genes:
        logger.info('Filtering snp annotation...')
        mask = flat_target_snp_ann_df.Gene.isin(gene_list)
        flat_target_snp_ann_df = flat_target_snp_ann_df[mask]
        check_df(flat_target_snp_ann_df, item='Target Gene SNP')
    target_snp_loc_df = flat_target_snp_ann_df.loc[:, ['#CHROM', 'POS', 'ALT'
                                                       ]].drop_duplicates()
    target_snp_index = target_snp_loc_df.set_index(['#CHROM', 'POS',
                                                    'ALT']).index

    # step6 add snp freq
    logger.info('Loading input sample snp freq data...')
    snp_freq_obj = core.SNPscore(vcf_table_files, snp_ann_file, group_labels, outdir)
    snp_freq_obj.load_stats()
    snp_freq_obj.snp_stats_df = snp_freq_obj.snp_stats_df.reindex(
        target_snp_index)
    snp_freq_obj.snp_stats_df.fillna(0, inplace=True)
    snp_dep_df = snp_freq_obj.snp_stats_df.loc[:, 'dep_count']
    reads_cov_snp = snp_dep_df.sum(1) > 0
    snp_freq_obj.snp_stats_df = snp_freq_obj.snp_stats_df[reads_cov_snp]
    if group_labels:
        snp_freq_obj.group_stats()
        freq_df = snp_reads_inf(snp_freq_obj)
    else:
        sample_names = [
            PurePath(file_i).stem for file_i in vcf_table_files.split(',')
        ]
        freq_df = snp_freq_obj.snp_stats_df.loc[:, sample_names]
        freq_df.columns = freq_df.columns.droplevel()
    freq_df.columns.name = ''
    freq_df = freq_df.reset_index()
    flat_target_snp_ann_freq_df = flat_target_snp_ann_df.merge(freq_df)
    check_df(flat_target_snp_ann_freq_df, item='Target SNP in samples')
    flat_target_snp_ann_freq_df.sort_values(['#CHROM', 'POS'], inplace=True)
    if outfmt == 'string':
        flat_target_snp_ann_freq_str = flat_target_snp_ann_freq_df.to_string(index=False)
        print(flat_target_snp_ann_freq_str)
    else:
        if not flat_target_snp_ann_freq_df.empty:
            target_region_snp_inf_file = outdir / 'target_region_snp_inf.txt'
            flat_target_snp_ann_freq_df.to_csv(target_region_snp_inf_file,
                                           sep='\t',
                                           index=False)
        else:
            logger.warning('No snp was covered by reads in target region.')


    logger.info('The End.')
