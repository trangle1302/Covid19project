import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import os
from statsmodels.sandbox.stats.multicomp import multipletests

EXPERIMENT = "CZ8780"
ACQUIRED_DATA_PATH = f"/data/trang/HPA_{EXPERIMENT}" #HPA_DV9903_Prescreen"#/DV9903_240323_preHPA_II__2023-03-24T12_10_32-Measurement_1b"
FORMATTED_DATA_PATH = f"{ACQUIRED_DATA_PATH}_formatted"
IMG_FOLDER = f"{ACQUIRED_DATA_PATH}_max_projection"
ANNOTATION_FOLDER = f"{ACQUIRED_DATA_PATH}_max_projection_annotation"
SAVE_DEST = "/home/trangle/Desktop/Covid19project"
df = pd.read_csv(f"{ANNOTATION_FOLDER}/segmentation/cells_combined_nobigcell.csv")
df['well_id'] = [image_id.rsplit('_',1)[0] for image_id in df.image_id]
#meta = pd.read_csv('/home/trangle/Desktop/Covid19project/Experiment_design/meta_ab.csv')
meta = pd.read_csv(f'/home/trangle/Desktop/Covid19project/Experiment_design/meta_ab_{EXPERIMENT}.csv')
meta["well_id"] = meta.plate + "_" + meta.well_id
meta["Gene"] = meta["gene_names"]
meta["ENSGID"] = meta["ensembl_ids"]

empty_wells=set(df.well_id) - set(meta.well_id)
well_ids = set(df.well_id) - empty_wells
fc_p_df= pd.DataFrame()
for area in ['nuclei', 'cytosol','cell']:
    colname_protein = 'protein-%s-mean' % area
    #colname_virus = 'virus-%s-integration' % area
    
    for well_id in well_ids:
        tmp = df[df.well_id == well_id]
        well_meta = meta[meta.well_id==well_id]
        infected = tmp[tmp.Infected==1][colname_protein]
        ninfected = tmp[tmp.Infected==0][colname_protein]
        fc = np.mean(infected)/np.mean(ninfected)
        t2, p2 = stats.ttest_ind(infected,ninfected)
        p_info = {'well_id': well_id,
                  'Antibody': well_meta.Antibody.values[0], 
                  'Gene_first': str(well_meta.Gene.values[0]).split(';')[0],
                  'Gene': well_meta.Gene.values[0], 
                  'ENSGID': well_meta.ENSGID.values[0],
                  'Region': area,
                  'n_infected': len(infected), 'n_noninfected':len(ninfected),
                  'avg_infected': np.mean(infected), 'avg_noninfected':np.mean(ninfected),
                  'fc': np.mean(fc), 'pval': p2}
        fc_p_df = fc_p_df.append(p_info, ignore_index=True)
print(fc_p_df)
fc_p_df["Log2FC"] = np.log2(fc_p_df.fc)
fc_p_df["pval_adjusted"] = multipletests(fc_p_df.pval, method='bonferroni')[1]

fc_p_df.to_csv(os.path.join(SAVE_DEST, f'Foldchange_meanintensity_{EXPERIMENT}.csv'), index=False)
#https://gist.github.com/dblyon/402a5fe1e1e211f34d5b230e0a4c93cc
"""
# Create Protein-Antibody-Well mapping
meta = pd.read_csv('/home/trangle/Desktop/Covid19project/Experiment_design/meta.csv')
total = pd.DataFrame()

for n in range(1,8):
    tmp = pd.read_csv('/home/trangle/Desktop/Covid19project/Experiment_design/Template_%s.csv'%n)
    tmp = pd.melt(tmp, id_vars=['Plate%s'%n])
    tmp['well_id'] = ['_'.join((str(n),row['Plate%s'%n]+row['variable'])) for i, row in tmp.iterrows()]
    tmp['Antibody'] = tmp['value']
    total = pd.concat([total, tmp])

total = total[~total.Antibody.isin(['0','rr','Ace2 1:10','Ace2 1:100'])]
total = total[~total.Antibody.isna()]

total['Antibody'] = [ab.strip() for ab in total['Antibody']]

meta = pd.merge(meta,total[['Antibody','well_id']], how='left',on='Antibody')
meta.to_csv('/home/trangle/Desktop/Covid19project/Experiment_design/meta_ab.csv',index=False)

# In R
library(ggplot2)
fold_changes <- c(rnorm(20000, 0, 2))
pvalues <- runif(n=20000, min=1e-50, max=.1)
dif <- data.frame(fc =fold_changes,pv =pvalues)
dif$thershold <- ifelse(dif$fc > 1 & dif$pv < 0.01, "red", 
                        ifelse(dif$fc < -1 & dif$pv < 0.01, -1, "blue"))
ggplot(data=dif, aes(x=fc, y=-log10(pv))) +
  geom_point( size=1 ,aes(color=as.factor(thershold))) +
  theme(legend.position = "none") +
  xlim(c(-10, 10)) + ylim(c(0, 15)) +
  xlab("log2 fold change") + ylab("-log10 p-value")  + theme_bw()+
  annotate("label", x =c(-8,5), y = 4.75, label = c("400","120"), col=c("red","steelblue"))+
  annotate("text", x =c(-8,5), y = 5, label = c("50 FC>4","8FC <-4"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),)
"""

  
  