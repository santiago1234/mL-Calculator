"""
Extract the YRI and IBS samples
"""
import sys

import pandas as pd

sys.path.append('../../../../proyectos/mxb-genomes/')
from mxbgenomes.utils import load_populations_info

popinfo = load_populations_info('../../../../proyectos/mxb-genomes/')
pops_to_keep = ['YRI', 'IBS']
popinfo = popinfo[popinfo['Subpopulation'].isin(pops_to_keep)]
popinfo.rename(columns={'Subpopulation': 'Population'}, inplace=True)
popinfo = popinfo[['Samplename', 'Population']]

# save data 
popinfo.to_csv('data/popinfo.csv', index=False)
# list of samples
popinfo.Samplename.to_csv('data/samples.txt', index=False, header=False)



