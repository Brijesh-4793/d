# -*- coding: utf-8 -*-
"""Step_by_step_simulation.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1SkPqAHh5bg2LZWD_GirJE83QzaASdYBr
"""

from google.colab import drive
drive.mount('/content/drive')

import sys
sys.path.insert(0,'/content/drive/MyDrive/DNA_Project/Model')

import config

import numpy as np

import Model

from Model import *

import numpy as np

from config import DEFAULT_PASSER, TM_NGS, TM_NNP

pip install reedsolo

from reedsolo import RSCodec

import sys
sys.path.insert(0,'/content/drive/MyDrive/DNA_Project/Encode')

import Helper_Functions
import DNAFountain

from Helper_Functions import preprocess, rs_decode, dna_to_int_array, load_dna
from DNAFountain import DNAFountain, Glass

import sys
sys.path.insert(0,'/content/drive/MyDrive/DNA_Project/Analysis')

import html_printer

from html_printer import html_table, color_bold,background_color_bold

import Analysis
from Analysis import inspect_distribution, examine_strand, save_simu_result, dna_chunk

import Fountain_analyzer
from Fountain_analyzer import error_profile

import matplotlib.pyplot as plt
import plotly.express as px
import plotly.io as pio

import warnings
warnings.simplefilter("ignore")

pip install -U scienceplots

import scienceplots

# Commented out IPython magic to ensure Python compatibility.
pio.templates.default = "plotly_white"
pio.renderers.default = 'notebook'
plt.style.use(['science','no-latex']) #
plt.rcParams['savefig.dpi'] = 150 
plt.rcParams['figure.dpi'] = 150
# %load_ext autoreload
# %autoreload 2

"""# Load Data"""

# load dna
in_dnas = load_dna('/content/drive/MyDrive/DNA_Project/files/lena.dna')
arg = DEFAULT_PASSER

"""# Step-by-step simulation

## Synthesis
"""

index = 1
SYN = Synthesizer(arg)
dnas_syn = SYN(in_dnas)
inspect_distribution(dnas_syn, show = True)
examine_strand(dnas_syn, index = index)

"""## Decay"""

DEC = Decayer(arg)
dnas_dec = DEC(dnas_syn)

inspect_distribution(dnas_dec, show = True)
examine_strand(dnas_dec, index = index)

"""## PCR"""

PCR = PCRer(N = 12, p = 0.8)
dnas_pcr = PCR(dnas_dec)

inspect_distribution(dnas_pcr, show = True)
examine_strand(dnas_pcr, index = index)

"""## Sampling"""

SAM = Sampler(0.0005)
dnas_sam = SAM(dnas_pcr)

inspect_distribution(dnas_sam, show = True)
examine_strand(dnas_sam, index = index)

"""## Sequencing

"""

SEQ = Sequencer(arg)
dnas_seq = SEQ(dnas_sam)

inspect_distribution(dnas_seq, show = True)
examine_strand(dnas_seq, index = index)