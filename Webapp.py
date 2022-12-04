import streamlit as st
import plotly.io as pio
import matplotlib.pyplot as plt
from Model.Model import * 
import plotly.express as px

from Analysis.Analysis import dna_chunk, plot_oligo_number_distribution, plot_error_distribution, save_simu_result
from Analysis.Fountain_analyzer import FT_Analyzer_Simplified
from Encode.DNAFountain import *
pio.templates.default = "plotly_white"
import Model.config as config


# ---------------------------- Choosing parameters ---------------------------- # 

def inspect(dnas, num_th = 1, inspect_index = 20):
    fig = plt.figure(figsize = (12,6))
    plt.subplot(1,2,1)
    plot_oligo_number_distribution(dnas)
    plt.subplot(1,2,2)
    plot_error_distribution(dnas,th = num_th)
    st.write(fig)
    dc = dna_chunk(dnas[inspect_index],'html')
    table = dc.plot_re_dnas()
    st.markdown(table, unsafe_allow_html = True)


# --------------------------- assigning parameters ------------------------ # 
arg = config.DEFAULT_PASSER
st.sidebar.subheader('Parameters of DNA data storage channel')
arg.syn_number = st.sidebar.slider('Syn number', min_value = 10, max_value = 60, value = 30)
arg.syn_sub_prob = st.sidebar.number_input('Syn Error rate', min_value = 0.0, max_value = 0.1, value = 0.01) / 3 # 3 kinds of substitutions
arg.syn_yield = st.sidebar.slider('Syn Yield', min_value = 0.98, max_value = 0.995, value = 0.99)

arg.pcrc = st.sidebar.slider('PCR cycle',min_value = 0, max_value =20,value = 12)
arg.pcrp = st.sidebar.number_input('PCR prob',min_value = 0.5, max_value = 1.0,value = 0.8)

arg.sam_ratio = st.sidebar.number_input('Sampling ratio',min_value = 0.0, max_value =1.0,value = 0.005)
arg.seq_depth = st.sidebar.slider('Seq Depth', min_value = 1, max_value = 100, value = 10)
seq_platform = st.sidebar.selectbox('Sequencing Platform',['Illumina Sequencing','Nanopore'])

index = st.sidebar.slider('inspect index', max_value = 600, value = 0)

st.sidebar.subheader('Parameters of Fountain code')
alpha = st.sidebar.slider('Alpha', min_value = 0.25, max_value = 1.0, value = 0.5)
rs = st.sidebar.slider('RS', min_value = 4, max_value = 10, value = 4)
num_th = int(rs / 2)

if seq_platform == 'Illumina Sequencing':
    arg.seq_TM = config.TM_NGS
else:
    arg.seq_TM = config.TM_NNP

# -------------------------- encoding ------------------------------------ #
st.header('DNA error simulation pipeline')
st.markdown('201901296 - Yash jain')
st.markdown('201901427 - Shubham Chudasama')
st.markdown('201901087 - Vivek Makawan')
st.markdown('201901473 - Brijesh Rathva')

st.markdown('You can assign parameters of DNA data storage channel and fountain code in the sidebar.\
    Play with different combinations of parameters to see how the noise structures and optimal encoding designs are influenced.') 

file_name = 'lena.jpg'
file_name, suffix = file_name.split('.')
file_name = 'files/' + file_name

in_file_name = file_name + '.' + suffix
in_dna_name = file_name + '.dna'
out_dna_name = file_name + '_simu.dna'
out_file_name = file_name + '_re.' + suffix

print(in_file_name)

st.header('Error simulation of the DNA data storage channel')

st.subheader('Implementation')
with open(in_dna_name) as f:
    dnas = f.readlines()
in_dnas = [dna.split('\n')[0] for dna in dnas]

'''Here we have coded error rate generation code in python to show error produced in stages of DNA storage stages like Synthesis  Decay PCR Sampling Sequencing. 
Then with that generated error rate we will find number of sequences with error w.r.t sequencing depth.
We have generated error in different levels of a DNA storage pipeline by linking DNA modules 
'''
'The simulation process now begins:'

st.subheader('Synthesis')
SYN = Synthesizer(arg)
dnas_syn = SYN(in_dnas)
inspect(dnas_syn,inspect_index = index)

st.subheader('Decay')
DEC = Decayer(arg)
dnas_dec = DEC(dnas_syn)
inspect(dnas_dec,inspect_index = index)

st.subheader('PCR')
PCR = PCRer(arg = arg)
dnas_pcr = PCR(dnas_dec)
inspect(dnas_pcr,inspect_index = index)

st.subheader('Sampling')
SAM = Sampler(arg = arg)
dnas_sam = SAM(dnas_pcr)
inspect(dnas_sam,inspect_index = index)

st.subheader('Sequencing')
SEQ = Sequencer(arg)
dnas_seq = SEQ(dnas_sam)
inspect(dnas_seq,inspect_index = index)
save_simu_result(dnas_seq,out_dna_name)
