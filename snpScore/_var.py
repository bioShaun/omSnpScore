import pkg_resources
import numpy as np
from enum import Enum, IntEnum
from pathlib import PurePath, Path


class SnpGroup(Enum):
    mut = 'mutant'
    wild = 'wild'
    mut_pa = 'mutant_parent'
    wild_pa = 'wild_parent'
    bg = 'background'


class SnpGroupFreq(Enum):
    mut = 'mutant.FREQ'
    wild = 'wild.FREQ'
    mut_pa = 'mutant_parent.FREQ'
    wild_pa = 'wild_parent.FREQ'
    bg = 'background.FREQ'


class SnpRep(IntEnum):
    alt = 0
    ref = 1
    unkown = 2


class VarScoreParams(IntEnum):
    snp_number_step = 5
    snp_number_window = 20
    min_depth = 5
    qtlseqr_min_depth = 5


script_dir = Path(__file__).parent

SNP_SCORE_PLOT = script_dir / 'plot/snpScorePlot.R'
QTLSEQR_PLOT = script_dir / 'plot/run_qtlseqr.R'
QTLSEQR_PLOT_WEB = script_dir / 'plot/run_qtlseqr_web.R'
DATA_DIR = PurePath(pkg_resources.resource_filename('snpScore', 'data'))
CHR_SIZE = DATA_DIR / 'chr.size'
CHR_WINDOW = DATA_DIR / 'chr.1m.window.bed'

OFFSET = 1e-05
SNP_FREQ_BIAS = 0.1
ALT_FREQ = np.round(2 / 3 - SNP_FREQ_BIAS, 2)
REF_FREQ = np.round(1 / 3 + SNP_FREQ_BIAS, 2)
VCF_SAMPLE_INDEX = 9

GROUPS = tuple([grp_i.value for grp_i in SnpGroup.__members__.values()])
MUT_NAME = SnpGroup.mut.value
WILD_NAME = SnpGroup.wild.value

# OUTPUT TABLE

COLUMN_NAME_MAP = {
    'Chr': 'CHROM',
    'Chrom': 'CHROM',
    'Pos': 'POS',
    'Alt': 'ALT',
    'snp_score': 'varBScore',
    'mutant.FREQ': 'mutant.AF',
    'wild.FREQ': 'wild.AF',
    'LOW.FREQ': 'wild.AF',
    'DP.LOW': 'wild.DP',
    'SNPindex.LOW': 'wild.SNPindex',
    'HIGH.FREQ': 'mutant.AF',
    'DP.HIGH': 'mutant.DP',
    'SNPindex.HIGH': 'mutant.SNPindex',
    'REF_FRQ': 'REF_FRQ',
    'euc': 'ED'
}

SNP_BASIC_COL = [
    'CHROM', 'POS', 'ALT', 'wild.AF', 'mutant.AF', 'AFD(deltaSNP)', 'REF_FRQ',
    'wild.DP', 'mutant.DP'
]

# varscore
VAR_SCORE_OUT_COL = [
    'CHROM', 'Start', 'End', 'varBScore', 'wild.AF', 'mutant.AF',
    'AFD(deltaSNP)', 'REF_FRQ', 'wild.DP', 'mutant.DP', 'POS', 'REF', 'ALT',
    'Feature', 'Gene', 'Transcript'
]

## qtlseqr
SCIENTIFIC_NUMBER_COLS = [
    'pvalue', 'negLog10Pval', 'qvalue', 'ED', 'fitted', 'unfitted',
    'dis2edcutoff'
]

ED_SPECIFIC_COLS = ['ED', 'fitted', 'unfitted', 'dis2edcutoff']

QTLSEQR_SPECIFIC_COLS = [
    'nSNPs', 'tricubeDeltaSNP', 'Gprime', 'negLog10Pval', 'qvalue', 'minDP',
    'CI_95', 'CI_99'
]


# OUTPUT DIR
class VarScoreOutDirName(Enum):
    snp_density = '1.SNP_Density'
    var_score = '2.varBScore'
    ed = '3.ed'
    qtlseqr = '4.QTLseqr'
    circos = '5.circosPlot'


# result readme
DOC_DIR = script_dir / 'doc'


class VarScoreDocName(Enum):
    snp_density = DOC_DIR / 'SNP_Density.readme.txt'
    var_score = DOC_DIR / 'varBScore.readme.txt'
    ed = DOC_DIR / 'ed.readme.txt'
    qtlseqr = DOC_DIR / 'qtlseqr.readme.txt'
    circos = DOC_DIR / 'circos.readme.txt'
