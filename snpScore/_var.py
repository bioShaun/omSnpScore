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
