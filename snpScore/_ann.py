import sys
import delegator
import pandas as pd
from io import StringIO
from loguru import logger
from pathlib import Path, PurePath
from datetime import datetime
from pybedtools import BedTool
from ._load import snpAnnTable, snpAnnTableByChr
from ._score import snpAnnBox
from ._utils import extract_snpeff_anno, split_dataframe_rows

COL_HEADER_MAP = {'Chr': '#CHROM', 'Pos': 'POS', 'Alt': 'ALT'}


def make_gene_bed(gene_bed, genes, target_bed):
    gene_bed_df = pd.read_csv(gene_bed,
                              sep='\t',
                              header=None,
                              names=['chrom', 'start', 'end'],
                              index_col=3)
    tmp_target_bed = PurePath(target_bed).with_suffix('.tmp.bed')
    sorted_tmp_target_bed = tmp_target_bed.with_suffix('.sorted.bed')
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
    target_bed_df.to_csv(tmp_target_bed, sep='\t', index=False, header=False)
    sort_cmd = f'sort -k1,1 -k2,2n {tmp_target_bed} > {sorted_tmp_target_bed}'
    delegator.run(sort_cmd)
    merge_region_cmd = f'bedtools merge -i {sorted_tmp_target_bed} > {target_bed}'
    delegator.run(merge_region_cmd)


def make_pos_bed(position, target_bed):
    chrom, start_end = position.split(':')
    start, end = start_end.split('-')
    start = int(start)
    if start > 0:
        start = start - 1
    logger.info('Making target region bed from input position...')
    with open(target_bed, 'w') as bed_inf:
        bed_inf.write(f'{chrom}\t{start}\t{end}\n')


def df2list(df):
    df_list = list()
    if not df.empty:
        df_list.append(list(df.columns))
        for index_i in df.index:
            df_list.append(list(df.loc[index_i]))
    return df_list


def printdf(df):
    col_str = [str(col_i) for col_i in df.columns]
    print('\t'.join(col_str))
    for index_i in df.index:
        line_i = [str(col_i) for col_i in df.loc[index_i]]
        print('\t'.join(line_i))


def check_df(df, item='SNP'):
    if df.empty:
        logger.warning('{} Not Found!'.format(item))
        sys.exit(1)


def init_logger(logger_file):
    if logger_file:
        logger.add(logger_file)
        logger.remove(0)


def snp_ann_pipe(gene_bed,
                 snp_ann_pkl,
                 outdir,
                 genes,
                 position,
                 vcf_dir,
                 sample_list,
                 group_list,
                 min_depth,
                 outfmt='table',
                 logger_file=None,
                 chrom=None):
    init_logger(logger_file)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    # step1 make selected gene/region bedfile
    time_now_str = '-'.join(str(datetime.now()).split())
    target_region_bed = outdir / f'target.{time_now_str}.bed'
    if genes:
        gene_list = genes.split(',')
        make_gene_bed(gene_bed, gene_list, target_region_bed)
    elif position:
        make_pos_bed(position, target_region_bed)
    else:
        raise ValueError('one of genes and postion must specified.')

    # step2 annotate target region
    if chrom is None:
        snp_table_obj = snpAnnTable(out_dir=outdir,
                                    table_dirs=vcf_dir,
                                    samples=sample_list,
                                    sample_label=group_list,
                                    min_depth=min_depth,
                                    filter_dp_grp=group_list,
                                    merge_method='outer',
                                    save_table=False)
    else:
        snp_table_obj = snpAnnTableByChr(out_dir=outdir,
                                         table_dirs=vcf_dir,
                                         samples=sample_list,
                                         sample_label=group_list,
                                         min_depth=min_depth,
                                         filter_dp_grp=group_list,
                                         save_table=False,
                                         merge_method='outer',
                                         chrom=chrom)

    snp_inf_df = snp_table_obj.alt_freq_df
    selected_cols = [each for each in snp_inf_df.columns if 'FREQ' not in each]
    ad_cols = [each for each in snp_inf_df.columns if 'AD' in each]

    snp_ann_obj = snpAnnBox(alt_freq_df=snp_inf_df.loc[:, selected_cols],
                            grp_list=group_list,
                            method_list='',
                            outdir=outdir,
                            min_depth=0,
                            vcf_ann_file=snp_ann_pkl,
                            target_bed=target_region_bed,
                            chr_size='')

    # step3 filter annotation
    snpeff_ann_df = snp_ann_obj.snp_window_ann_df.copy()
    snpeff_anno = list(snpeff_ann_df.INFO.map(extract_snpeff_anno))
    snpeff_anno_df = pd.DataFrame(snpeff_anno)
    snpeff_anno_df.columns = [
        'Feature', 'Gene', 'Transcript', 'Variant_DNA_Level',
        'Variant_Protein_Level'
    ]
    flat_target_snp_ann_df = pd.concat([snpeff_ann_df, snpeff_anno_df], axis=1)
    flat_target_snp_ann_df.drop('INFO', axis=1, inplace=True)
    flat_target_snp_ann_df = split_dataframe_rows(flat_target_snp_ann_df,
                                                  column_selectors=[
                                                      'Feature', 'Gene',
                                                      'Transcript',
                                                      'Variant_DNA_Level',
                                                      'Variant_Protein_Level'
                                                  ],
                                                  row_delimiter='|')

    flat_target_snp_ann_df.drop(['Start', 'End'], axis=1, inplace=True)
    flat_target_snp_ann_df.fillna(0, inplace=True)
    order_cols = [
        each for each in flat_target_snp_ann_df.columns if each not in ad_cols
    ]
    order_cols.extend(ad_cols)
    flat_target_snp_ann_df = flat_target_snp_ann_df.loc[:, order_cols]
    if genes:
        logger.info('Filtering snp annotation...')
        mask = flat_target_snp_ann_df.Gene.isin(gene_list)
        flat_target_snp_ann_df = flat_target_snp_ann_df[mask]
        check_df(flat_target_snp_ann_df, item='Target SNP')

    flat_target_snp_ann_df.sort_values(['Chrom', 'Pos'], inplace=True)

    if outfmt == 'string':
        printdf(flat_target_snp_ann_df)
    elif outfmt == 'df':
        return flat_target_snp_ann_df
    else:
        if not flat_target_snp_ann_df.empty:
            target_region_snp_inf_file = outdir / 'target_region_snp_inf.txt'
            flat_target_snp_ann_df.to_csv(target_region_snp_inf_file,
                                          index=False)
        else:
            logger.warning('No snp was covered by reads in target region.')

    logger.info('The End.')
