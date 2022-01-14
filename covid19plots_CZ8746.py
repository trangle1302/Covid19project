import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import os

SAMPLES_DEST = "/data/trang/covid19_data_CZ8746_annotation/segmentation"
# df = pd.read_csv('/data/hao/forTrang/cells_combined_nobigcell.csv') 
# df = pd.read_csv('/home/trangle/Desktop/Covid19project/celldata_fromHao/cells_combined_nobigcell.csv')
meta = pd.read_csv('/home/trangle/Desktop/Covid19project/Experiment_design/meta_ab_CZ8746.csv')
df = pd.read_csv(os.path.join(SAMPLES_DEST, "cells_combined_nobigcell.csv")) 
df['well_id'] = [image_id.split('/')[-1][:-2] for image_id in df.image_id]
NONINFECTED = ["A", "B", "C", "D"]
INFECTED = ["E", "F", "G", "H"]

antibodies = set(meta.Antibody)
fc_p_df= pd.DataFrame()
for area in ['nuclei', 'cytosol','cell']:
    colname_protein = 'protein-%s-mean' % area
    #colname_virus = 'virus-%s-integration' % area
    
    for ab in antibodies:
        ab_wells = set(meta[meta.Antibody==ab].well_id)
        tmp = df[df.well_id.isin(ab_wells)]
        well_meta= meta[meta.Antibody==ab]
        infected = tmp[tmp.Infected==1][colname_protein]
        ninfected = tmp[tmp.Infected==0][colname_protein]
        fc = np.mean(infected)/np.mean(ninfected)
        t2, p2 = stats.ttest_ind(infected,ninfected)
        p_info = {'well_id': ab_wells,
                  'Antibody': well_meta.Antibody.values[0], 
                  'Gene_first': well_meta.Gene.values[0].split(';')[0],
                  'Gene': well_meta.Gene.values[0], 
                  'ENSGID': well_meta.ENSGID.values[0],
                  'Region': area,
                  'n_infected': len(infected), 'n_noninfected':len(ninfected),
                  'avg_infected': np.mean(infected), 'avg_noninfected':np.mean(ninfected),
                  'fc': np.mean(fc), 'pval': p2}
        fc_p_df = fc_p_df.append(p_info, ignore_index=True)

fc_p_df.to_csv(os.path.join(SAMPLES_DEST, 'Foldchange_ab_meanintensity.csv'), index=False)


empty_wells=set(df.well_id) - set(meta.well_id)
well_ids = set(df.well_id) - empty_wells
fc_p_df= pd.DataFrame()
for area in ['nuclei', 'cytosol','cell']:
    colname_protein = 'protein-%s-mean' % area
    #colname_virus = 'virus-%s-integration' % area
    
    for well_id in well_ids:
        tmp = df[df.well_id == well_id]
        well_meta= meta[meta.well_id==well_id]
        infected = tmp[tmp.Infected==1][colname_protein]
        ninfected = tmp[tmp.Infected==0][colname_protein]
        fc = np.mean(infected)/np.mean(ninfected)
        t2, p2 = stats.ttest_ind(infected,ninfected)
        p_info = {'well_id': well_id,
                  'Antibody': well_meta.Antibody.values[0], 
                  'Gene_first': well_meta.Gene.values[0].split(';')[0],
                  'Gene': well_meta.Gene.values[0], 
                  'ENSGID': well_meta.ENSGID.values[0],
                  'Region': area,
                  'n_infected': len(infected), 'n_noninfected':len(ninfected),
                  'avg_infected': np.mean(infected), 'avg_noninfected':np.mean(ninfected),
                  'fc': np.mean(fc), 'pval': p2}
        fc_p_df = fc_p_df.append(p_info, ignore_index=True)
        
fc_p_df.to_csv(os.path.join(SAMPLES_DEST, 'Foldchange_well_meanintensity.csv'), index=False)