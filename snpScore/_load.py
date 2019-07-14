import gzip
import attr
import numpy as np
import pandas as pd
from enum import Enum
from loguru import logger
from pathlib import Path
from functools import reduce
from ._utils import async_batch_sh_jobs, check_app
from ._utils import AmbigousSample


class SnpGroup(Enum):
    mut = 'mutant'
    wild = 'wild'
    mut_pa = 'mutant_parent'
    wild_pa = 'wild_parent'
    bg = 'background'


MUT_NAME = SnpGroup.mut.value
WILD_NAME = SnpGroup.wild.value
VCF_SAMPLE_INDEX = 9


@attr.s
class tableFromVcf:

    vcf = attr.ib(converter=Path)
    out_dir = attr.ib(converter=Path)
    thread = attr.ib(default=1)

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

    @property
    def samples(self):
        table_samples = []
        for sp_i in self.vcf_samples:
            sp_i_pkl = self.out_dir / f'{sp_i}.pkl'
            if sp_i_pkl.is_file():
                logger.warning(f'{sp_i_pkl} exsits, omit transfer.')
            else:
                table_samples.append(sp_i)
        return table_samples

    def _extract_from_vcf(self, sample_id):
        check_app('bcftools')
        check_app('table2pkl')
        cmd = (f'bcftools view --samples {sample_id} {self.vcf} | '
               f'table2pkl from_stdin --sample_id {sample_id} '
               f'--table-file-pkl {self.out_dir}/{sample_id}.pkl')
        return cmd

    @property
    def make_table(self):
        self.out_dir.mkdir(parents=True, exist_ok=True)
        logger.info('Extracting sample vcf start...')
        cmd_list = [self._extract_from_vcf(sp_i) for sp_i in self.samples]
        async_batch_sh_jobs(cmd_list, self.thread)
        logger.info('Extracting sample vcf done.')


@attr.s
class snpTable:

    out_dir = attr.ib(converter=Path)
    table_dirs = attr.ib(factory=list)
    samples = attr.ib(factory=list)
    sample_label = attr.ib(factory=list)
    min_depth = attr.ib(default=5)

    def __attrs_post_init__(self):
        self._ad_df = None
        self._grp_dep_df = None
        self._grp_alt_dep_df = None
        self._alt_freq_df = None
        self.alt_freq_file = self.out_dir / 'snp.freq.csv'
        self._qtlseqr_snp_table = self.out_dir / 'qtlseqr.csv'
        self.out_dir.mkdir(parents=True, exist_ok=True)

    @property
    def snp_table_files(self):
        table_file_list = []
        for sp_i in self.samples:
            sp_dir = []
            for dir_i in self.table_dirs:
                sp_i_pkl = Path(dir_i) / f'{sp_i}.pkl'
                if sp_i_pkl.is_file():
                    sp_dir.append(dir_i)
                    table_file_list.append(sp_i_pkl)
            if len(sp_dir) > 1:
                sp_dir_str = ', '.join(sp_dir)
                logger.error(f'{sp_i} in multiple directory: {sp_dir_str}')
        if len(table_file_list) != len(self.samples):
            raise AmbigousSample
        return table_file_list

    @property
    def ad_df(self):
        if self._ad_df is None:
            logger.info('Loading tables...')
            self.ad_dfs = [
                pd.read_pickle(table_i) for table_i in self.snp_table_files
            ]
            logger.info('Concatinating tables...')
            self._ad_df = reduce(
                lambda x, y: pd.merge(x, y, on=['Chr', 'Pos', 'Alt']),
                self.ad_dfs)
            self._ad_df.fillna(0, inplace=True)
        return self._ad_df

    @property
    def grp_dep_df(self):
        if self._grp_dep_df is None:
            logger.info('Group depth reads...')
            self.dep_df = self.ad_df.loc[:, 'dep_count'].copy()
            self.dep_df.columns = self.sample_label
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
    def grp_ref_dep_df(self):
        return self.grp_dep_df - self.grp_alt_dep_df

    @property
    def alt_freq_df(self):
        if self._alt_freq_df is None:
            if self.alt_freq_file.is_file():
                self._alt_freq_df = pd.read_csv(self.alt_freq_file)
            else:
                logger.info('Filtering allele depth...')
                mw_group = [MUT_NAME, WILD_NAME]
                dep_passed_snp = self.grp_dep_df.loc[:, mw_group].min(
                    1) >= self.min_depth
                self.passed_grp_dep_df = self.grp_dep_df[dep_passed_snp]
                self.passed_grp_alt_dep_df = self.grp_alt_dep_df[
                    dep_passed_snp]
                self.passed_grp_dep_df.applymap(lambda x: x if x >= self.
                                                min_depth else np.nan)
                logger.info('Calculating alt allele freq...')
                self._alt_freq_df = self.passed_grp_alt_dep_df / \
                    self.passed_grp_dep_df
                self._alt_freq_df = self._alt_freq_df.reset_index()
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
