import gzip
import attr
import delegator
import numpy as np
import pandas as pd
from loguru import logger
from pathlib import Path
from functools import reduce
from ._var import MUT_NAME, WILD_NAME
from ._var import VCF_SAMPLE_INDEX
from ._utils import async_batch_sh_jobs, check_app, table2vcf, extract_snpeff_annotation
from ._utils import SampleFileNotMatch, UnsupportedFormat
from ._utils import DuplicatedRecord
from ._utils import valid_grp


@attr.s
class tableFromVcf:
    vcf = attr.ib(converter=Path)
    out_dir = attr.ib(converter=Path)
    thread = attr.ib(default=1, converter=int)
    sample_prefix = attr.ib(default=None)
    out_format = attr.ib(default='pkl')

    @property
    def vcf_samples(self):
        if self.vcf.suffix == '.gz':
            vcf_inf = gzip.open(self.vcf, 'rt')
        else:
            vcf_inf = open(self.vcf)
        for eachline in vcf_inf:
            if eachline[:6] == '#CHROM':
                # WARNING: changing of VCF format could output wrong names
                return eachline.strip().split('\t')[VCF_SAMPLE_INDEX:]
        vcf_inf.close()

    @property
    def samples(self):
        table_samples = []
        for sp_i in self.vcf_samples:
            name = sp_i
            if self.sample_prefix:
                name = f'{self.sample_prefix}-{sp_i}'
            sp_i_pkl = self.out_dir / f'{name}.pkl'
            if sp_i_pkl.is_file():
                logger.warning(f'{sp_i_pkl} exsits, omit sample [{sp_i}].')
            else:
                table_samples.append(sp_i)
        return table_samples

    def _extract_from_vcf(self, sample_id):
        check_app('bcftools')
        check_app('table2pkl')
        name = sample_id
        if self.sample_prefix:
            name = f'{self.sample_prefix}-{sample_id}'
        cmd = (f'bcftools view --samples {sample_id} {self.vcf} | '
               f'table2pkl from_stdin --sample_id {sample_id} '
               f'--table-file-pkl {self.out_dir}/{name}.pkl')
        if self.out_format == 'csv':
            check_app('table2csv-mp')
            cmd = (f'bcftools view --samples {sample_id} {self.vcf} | '
                   f'table2csv-mp from_stdin --sample_id {sample_id} '
                   f'--csv-dir {self.out_dir}/{name}')
        return cmd

    @property
    def make_table(self):
        self.out_dir.mkdir(parents=True, exist_ok=True)
        logger.info('Extracting sample vcf start...')
        cmd_list = [self._extract_from_vcf(sp_i) for sp_i in self.samples]
        async_batch_sh_jobs(cmd_list, self.thread)
        logger.info('Extracting sample vcf done.')


@attr.s
class tableFromVcfMP(tableFromVcf):
    @property
    def samples(self):
        return self.vcf_samples

    def _extract_from_vcf(self, sample_id):
        check_app('bcftools')
        check_app('table2pkl')
        cmd = (f'bcftools view --samples {sample_id} {self.vcf} | '
               f'table2pkl-mp from_stdin --sample_id {sample_id} '
               f'--pkl-dir {self.out_dir}/{sample_id}')
        return cmd


@attr.s
class tableFromSelectTable(tableFromVcf):
    def __attrs_post_init__(self):
        self.st_df = pd.read_csv(self.vcf, sep='\t')
        self.pos_cols = list(self.st_df.columns[:3])

    @property
    def vcf_samples(self):
        return self.st_df.columns[3:]

    def check_table(self):
        if len(self.st_df[self.st_df.duplicated()]):
            raise DuplicatedRecord(
                'Duplicated records in snp table, please check!')

    def _extract_from_vcf(self, sample_id):
        self.check_table()
        check_app('table2pkl')
        sample_cols = self.pos_cols[:]
        sample_cols.append(sample_id)
        sample_df = self.st_df.loc[:, sample_cols]
        sample_df.columns = ['Chr', 'Pos', 'Alt', sample_id]
        sample_df.loc[:, sample_id] = [
            str(each).replace('|', ',') for each in sample_df.loc[:, sample_id]
        ]
        sample_table = self.out_dir / f'{sample_id}.table'
        sample_df.to_csv(sample_table, sep='\t', index=False)
        cmd = f'table2pkl from_file --table_file {sample_table}'
        return cmd


@attr.s
class snpTable:
    out_dir = attr.ib(converter=Path)
    table_dirs = attr.ib(factory=list)
    samples = attr.ib(factory=list)
    sample_label = attr.ib(factory=list)
    min_depth = attr.ib(default=5, converter=int)
    filter_dp_grp = attr.ib(default=[MUT_NAME, WILD_NAME])
    save_table = attr.ib(default=True)
    merge_method = attr.ib(default='inner')

    def __attrs_post_init__(self):
        self._ad_df = None
        self._grp_dep_df = None
        self._grp_alt_dep_df = None
        self._grp_ref_dep_df = None
        self._alt_freq_df = None
        self._qtlseqr_snp_table = self.out_dir / 'qtlseqr.csv'
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.valid_grp = valid_grp(self.sample_label)
        self.alt_freq_file = self.out_dir / 'snp.freq.csv'

    @property
    def snp_table_files(self):
        table_file_list = []
        for sp_i in self.samples:
            sp_dir = []
            for dir_i in self.table_dirs:
                sp_i_pkl = Path(dir_i) / f'{sp_i}.pkl'
                sp_i_h5 = Path(dir_i) / f'{sp_i}.h5'
                if sp_i_pkl.is_file():
                    sp_dir.append(dir_i)
                    table_file_list.append(sp_i_pkl)
                elif sp_i_h5.is_file():
                    sp_dir.append(dir_i)
                    table_file_list.append(sp_i_h5)
            if len(sp_dir) > 1:
                sp_dir_str = ', '.join(sp_dir)
                logger.error(f'{sp_i} in multiple directory: {sp_dir_str}')
            elif len(sp_dir) == 0:
                logger.error(f'{sp_i} not found.')
        if len(table_file_list) != len(self.samples):
            raise SampleFileNotMatch
        return table_file_list

    @property
    def ad_df(self):
        if self._ad_df is None:
            logger.info('Loading tables...')
            self.ad_dfs = []
            for table_i in self.snp_table_files:
                if table_i.suffix == '.pkl':
                    table_i_df = pd.read_pickle(table_i)
                elif table_i.suffix == '.h5':
                    table_i_df = pd.read_hdf(table_i)
                else:
                    raise UnsupportedFormat
                self.ad_dfs.append(table_i_df)
            logger.info('Concatinating tables...')
            self._ad_df = reduce(
                lambda x, y: pd.merge(
                    x, y, on=['Chr', 'Pos', 'Alt'], how=self.merge_method),
                self.ad_dfs)
            self._ad_df.fillna(0, inplace=True)
        return self._ad_df

    @property
    def grp_dep_df(self):
        if self._grp_dep_df is None:
            self.dep_df = self.ad_df.loc[:, 'dep_count'].copy()
            self.dep_df.columns = self.sample_label
            logger.info('Group depth reads...')
            self._grp_dep_df = self.dep_df.T.groupby(level=0).agg('sum').T
            self._grp_dep_df = self._grp_dep_df.loc[:, self.valid_grp]
        return self._grp_dep_df

    @property
    def grp_alt_dep_df(self):
        if self._grp_alt_dep_df is None:
            logger.info('Group alt reads...')
            self.alt_df = self.ad_df.loc[:, 'alt_count'].copy()
            self.alt_df.columns = self.sample_label
            self._grp_alt_dep_df = self.alt_df.T.groupby(level=0).agg('sum').T
            self._grp_alt_dep_df = self._grp_alt_dep_df.loc[:, self.valid_grp]
        return self._grp_alt_dep_df

    @property
    def grp_ref_dep_df(self):
        if self._grp_ref_dep_df is None:
            self._grp_ref_dep_df = self.grp_dep_df - self.grp_alt_dep_df
            self._grp_ref_dep_df = self._grp_ref_dep_df.astype('int')
        return self._grp_ref_dep_df

    @property
    def grp_ad_df(self):
        self.grp_ref_dep_df.columns = [
            f'{col_i}.REF.AD' for col_i in self.grp_ref_dep_df.columns
        ]
        self.grp_alt_dep_df.columns = [
            f'{col_i}.ALT.AD' for col_i in self.grp_alt_dep_df.columns
        ]
        grp_ad_df = self.grp_ref_dep_df.merge(self.grp_alt_dep_df,
                                              left_index=True,
                                              right_index=True)
        return grp_ad_df

    @property
    def alt_freq_df(self):
        if self._alt_freq_df is None:
            if self.alt_freq_file.is_file():
                self._alt_freq_df = pd.read_csv(self.alt_freq_file)
            else:
                dep_passed_snp = self.grp_dep_df.loc[:,
                                 self.filter_dp_grp].min(
                    1) >= self.min_depth
                self.passed_grp_dep_df = self.grp_dep_df[dep_passed_snp]
                self.passed_grp_alt_dep_df = self.grp_alt_dep_df[
                    dep_passed_snp]
                logger.info('Filtering allele depth...')
                self.passed_grp_dep_df.applymap(
                    lambda x: x if x >= self.min_depth else np.nan)
                logger.info('Calculating alt allele freq...')
                self._alt_freq_df = self.passed_grp_alt_dep_df / \
                                    self.passed_grp_dep_df
                self._alt_freq_df.columns = [
                    f'{col_i}.FREQ' for col_i in self._alt_freq_df.columns
                ]
                self._alt_freq_df = self._alt_freq_df.merge(self.grp_ad_df,
                                                            left_index=True,
                                                            right_index=True)
                self._alt_freq_df = self._alt_freq_df.reset_index()
                self._alt_freq_df.loc[:, 'Chr'] = self._alt_freq_df.Chr.astype(
                    'str')
                self._alt_freq_df.sort_values(['Chr', 'Pos'], inplace=True)
                if self.save_table:
                    self._alt_freq_df.to_csv(self.alt_freq_file, index=False)
        return self._alt_freq_df

    @property
    def qtlseqr_snp_table(self):
        if not self._qtlseqr_snp_table.is_file():
            ref_df = self.grp_ref_dep_df.copy()
            alt_df = self.grp_alt_dep_df.copy()
            ref_df.columns = [f'AD_REF.{sp_i}' for sp_i in ref_df.columns]
            alt_df.columns = [f'AD_ALT.{sp_i}' for sp_i in alt_df.columns]
            ref_df = ref_df.astype('int')
            qtlseqr_df = ref_df.merge(alt_df, on=['Chr', 'Pos', 'Alt'])
            qtlseqr_df.index.names = ['CHROM', 'POS', 'ALT']
            qtlseqr_df = qtlseqr_df[qtlseqr_df.sum(1) > 0]
            qtlseqr_df = qtlseqr_df.reset_index()
            qtlseqr_df.to_csv(self._qtlseqr_snp_table, index=False)
        return self._qtlseqr_snp_table


@attr.s
class snpTableMP(snpTable):
    chrom = attr.ib(default=None)

    @property
    def snp_table_files(self):
        table_file_list = []
        for sp_i in self.samples:
            sp_dir = []
            for dir_i in self.table_dirs:
                sp_i_pkl = Path(dir_i) / f'{sp_i}' / f'{self.chrom}.pkl'
                sp_i_h5 = Path(dir_i) / f'{sp_i}' / f'{self.chrom}.h5'
                sp_i_csv = Path(dir_i) / f'{sp_i}' / f'{self.chrom}.csv'
                if sp_i_pkl.is_file():
                    sp_dir.append(dir_i)
                    table_file_list.append(sp_i_pkl)
                elif sp_i_h5.is_file():
                    sp_dir.append(dir_i)
                    table_file_list.append(sp_i_h5)
                elif sp_i_csv.is_file():
                    sp_dir.append(dir_i)
                    table_file_list.append(sp_i_csv)
            if len(sp_dir) > 1:
                sp_dir_str = ', '.join(sp_dir)
                logger.error(f'{sp_i} in multiple directory: {sp_dir_str}')
            elif len(sp_dir) == 0:
                logger.error(f'{sp_i} not found.')
        if len(table_file_list) != len(self.samples):
            raise SampleFileNotMatch
        return table_file_list

    @property
    def ad_df(self):
        if self._ad_df is None:
            logger.info('Loading tables...')
            self.ad_dfs = []
            for table_i in self.snp_table_files:
                if table_i.suffix == '.pkl':
                    table_i_df = pd.read_pickle(table_i)
                elif table_i.suffix == '.h5':
                    table_i_df = pd.read_hdf(table_i)
                elif table_i.suffix == '.csv':
                    table_i_df = pd.read_csv(table_i)
                    sample_name = Path(table_i).parent.name
                    table_i_df.loc[:, 'sample_id'] = sample_name
                    # filter multi alt snp
                    table_i_df = table_i_df[~table_i_df.Alt.str.contains(",")]
                    table_i_df = table_i_df.set_index(
                        ['Chr', 'Pos', 'Ref', 'Alt',
                         'sample_id']).unstack('sample_id')
                else:
                    raise UnsupportedFormat
                self.ad_dfs.append(table_i_df)
            logger.info('Concatinating tables...')
            self._ad_df = reduce(
                lambda x, y: pd.merge(
                    x, y, on=['Chr', 'Pos', 'Ref', 'Alt'], how=self.merge_method),
                self.ad_dfs)
            self._ad_df.fillna(0, inplace=True)
            self._ad_df = self._ad_df.astype('int')
        return self._ad_df

    @property
    def grp_ref_dep_df(self):
        if self._grp_ref_dep_df is None:
            logger.info('Group ref reads...')
            self.ref_df = self.ad_df.loc[:, 'ref_count'].copy()
            self.ref_df.columns = self.sample_label
            self._grp_ref_dep_df = self.ref_df.T.groupby(level=0).agg('sum').T
            self._grp_ref_dep_df = self._grp_ref_dep_df.loc[:, self.valid_grp]
        return self._grp_ref_dep_df

    @property
    def grp_dep_df(self):
        if self._grp_dep_df is None:
            self._grp_dep_df = self.grp_ref_dep_df + self.grp_alt_dep_df
            self._grp_dep_df = self._grp_dep_df.astype('int')
        return self._grp_dep_df

@attr.s
class snpAnnTable(snpTable):
    @property
    def grp_dep_df(self):
        if self._grp_dep_df is None:
            self.dep_df = self.ad_df.loc[:, 'dep_count'].copy()
            self.dep_df.columns = self.sample_label
            logger.info('Group depth reads...')
            self._grp_dep_df = self.dep_df.T.groupby(level=0).agg('sum').T
        return self._grp_dep_df

    @property
    def grp_alt_dep_df(self):
        if self._grp_alt_dep_df is None:
            logger.info('Group alt reads...')
            self.alt_df = self.ad_df.loc[:, 'alt_count'].copy()
            self.alt_df.columns = self.sample_label
            self._grp_alt_dep_df = self.alt_df.T.groupby(level=0).agg('sum').T
        return self._grp_alt_dep_df

    @property
    def alt_freq_df(self):
        if self._alt_freq_df is None:
            if self.alt_freq_file.is_file():
                self._alt_freq_df = pd.read_csv(self.alt_freq_file)
            else:
                dep_passed_snp = self.grp_dep_df.loc[:,
                                 self.filter_dp_grp].max(
                    1) >= self.min_depth
                self.passed_grp_dep_df = self.grp_dep_df[dep_passed_snp]
                self.passed_grp_alt_dep_df = self.grp_alt_dep_df[
                    dep_passed_snp]
                logger.info('Filtering allele depth...')
                self.passed_grp_dep_df.applymap(
                    lambda x: x if x >= self.min_depth else np.nan)
                logger.info('Calculating alt allele freq...')
                self._alt_freq_df = self.passed_grp_alt_dep_df / \
                                    self.passed_grp_dep_df
                self._alt_freq_df.columns = [
                    f'{col_i}.FREQ' for col_i in self._alt_freq_df.columns
                ]
                self._alt_freq_df = self._alt_freq_df.merge(self.grp_ad_df,
                                                            left_index=True,
                                                            right_index=True)
                self._alt_freq_df = self._alt_freq_df.reset_index()
                self._alt_freq_df.loc[:, 'Chr'] = self._alt_freq_df.Chr.astype(
                    'str')
                self._alt_freq_df.sort_values(['Chr', 'Pos'], inplace=True)
                if self.save_table:
                    self._alt_freq_df.to_csv(self.alt_freq_file, index=False)
        return self._alt_freq_df


@attr.s
class snpAnnTableByChr(snpTableMP):
    @property
    def grp_alt_dep_df(self):
        if self._grp_alt_dep_df is None:
            logger.info('Group alt reads...')
            self.alt_df = self.ad_df.loc[:, 'alt_count'].copy()
            self.alt_df.columns = self.sample_label
            self._grp_alt_dep_df = self.alt_df.T.groupby(level=0).agg('sum').T
        return self._grp_alt_dep_df

    @property
    def grp_ref_dep_df(self):
        if self._grp_ref_dep_df is None:
            logger.info('Group ref reads...')
            self.ref_df = self.ad_df.loc[:, 'ref_count'].copy()
            self.ref_df.columns = self.sample_label
            self._grp_ref_dep_df = self.ref_df.T.groupby(level=0).agg('sum').T
        return self._grp_ref_dep_df

    @property
    def alt_freq_df(self):
        if self._alt_freq_df is None:
            if self.alt_freq_file.is_file():
                self._alt_freq_df = pd.read_csv(self.alt_freq_file)
            else:
                dep_passed_snp = self.grp_dep_df.loc[:,
                                 self.filter_dp_grp].max(
                    1) >= self.min_depth
                self.passed_grp_dep_df = self.grp_dep_df[dep_passed_snp]
                self.passed_grp_alt_dep_df = self.grp_alt_dep_df[
                    dep_passed_snp]
                logger.info('Filtering allele depth...')
                self.passed_grp_dep_df.applymap(
                    lambda x: x if x >= self.min_depth else np.nan)
                logger.info('Calculating alt allele freq...')
                self._alt_freq_df = self.passed_grp_alt_dep_df / \
                                    self.passed_grp_dep_df
                self._alt_freq_df.columns = [
                    f'{col_i}.FREQ' for col_i in self._alt_freq_df.columns
                ]
                self._alt_freq_df = self._alt_freq_df.merge(self.grp_ad_df,
                                                            left_index=True,
                                                            right_index=True)
                self._alt_freq_df = self._alt_freq_df.reset_index()
                self._alt_freq_df.loc[:, 'Chr'] = self._alt_freq_df.Chr.astype(
                    'str')
                self._alt_freq_df.sort_values(['Chr', 'Pos'], inplace=True)
                if self.save_table:
                    self._alt_freq_df.to_csv(self.alt_freq_file, index=False)
        return self._alt_freq_df
