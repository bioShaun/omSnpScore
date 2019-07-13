import gzip
import attr
from ._utils import async_batch_sh_jobs
from loguru import logger
from pathlib import Path

VCF_SAMPLE_INDEX = 9


@attr.s
class tableFromVcf:

    vcf = attr.ib(converter=Path)
    out_dir = attr.ib(converter=Path)
    thread = attr.ib(default=1)

    @property
    def samples(self):
        if self.vcf.suffix == '.gz':
            vcf_inf = gzip.open(self.vcf, 'rt')
        else:
            vcf_inf = open(self.vcf)
        for eachline in vcf_inf:
            if eachline[:6] == '#CHROM':
                # WARNING: changing of VCF format could output wrong names
                return eachline.strip().split('\t')[VCF_SAMPLE_INDEX:]

    def _extract_from_vcf(self, sample_id):
        cmd = (f'bcftools view --samples {sample_id} {self.vcf} | '
               f'table2pkl from_stdin --sample_id {sample_id} '
               f'--table-file-pkl {self.out_dir}/{sample_id}.pkl')
        print(cmd)
        return cmd

    @property
    def make_table(self):
        self.out_dir.mkdir(parents=True, exist_ok=True)
        logger.info('Extracting sample vcf start...')
        cmd_list = [self._extract_from_vcf(sp_i) for sp_i in self.samples]
        async_batch_sh_jobs(cmd_list, self.thread)
        logger.info('Extracting sample vcf done.')
