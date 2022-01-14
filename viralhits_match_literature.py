# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 15:17:25 2020

@author: trang.le
"""

import pandas as pd

# Our intensity hits
intensityvar_df = pd.read_excel('./Downloads/Review_Analysis.xlsx', sheet_name=6, skiprows=2)
intensityvar_df=intensityvar_df[:178]

# Hits from Gordon et al 2020
gordon_ST2 = pd.read_excel('./Downloads/41586_2020_2286_MOESM6_ESM.xlsx', skiprows=1)
gordon_ST2['bait'] = [b.split(' ')[1].upper() for b in gordon_ST2.Bait]

# Hits from Stukalov et al., 2020
stukalov_STA = pd.read_excel('./Downloads/media-2.xlsx', sheet_name=1, skiprows=1)
stukalov_STA = stukalov_STA[stukalov_STA.organism == 'SARS-CoV-2']
stukalov_STA.loc[stukalov_STA.gene_name == 43891, 'gene_name'] = 'MTARC1'
stukalov_STA.gene_name = [g.replace('...', '') for g in stukalov_STA.gene_name]

intensityvar_df['Bait_match_Gordon'] = ''
intensityvar_df['Bait_match_Stukalov'] = ''
for i, row in intensityvar_df.iterrows():
    gene = row.Gene
    gene_gordon = gordon_ST2[gordon_ST2.PreyGene == gene]
    gene_stukalov = stukalov_STA[stukalov_STA.gene_name == gene]
    if len(gene_gordon) > 0:
        bait_match = ', '.join(list(gene_gordon.bait))
    else:
        bait_match = None
    intensityvar_df.loc[i,'Bait_match_Gordon'] = bait_match
    
    if len(gene_stukalov) > 0:
        bait_match = ', '.join(list(gene_stukalov.bait_name))
    else:
        bait_match = None
    intensityvar_df.loc[i,'Bait_match_Stukalov'] = bait_match
    
intensityvar_df.to_csv('./Downloads/Sheet6_matchedGordon_matchedStukalov.csv', index=False)

# Our spatial hits
spatialvar_df = pd.read_excel('./Downloads/Review_Analysis.xlsx', sheet_name=4, skiprows=2)
spatialvar_df = spatialvar_df[:84]

spatialvar_df['Bait_match_Gordon'] = ''
spatialvar_df['Bait_match_Stukalov'] = ''
for i, row in spatialvar_df.iterrows():
    gene = row.Gene
    gene_gordon = gordon_ST2[gordon_ST2.PreyGene == gene]
    gene_stukalov = stukalov_STA[stukalov_STA.gene_name == gene]
    if len(gene_gordon) > 0:
        bait_match = ', '.join(list(gene_gordon.bait))
    else:
        bait_match = None
    spatialvar_df.loc[i,'Bait_match_Gordon'] = bait_match
    
    if len(gene_stukalov) > 0:
        bait_match = ', '.join(list(gene_stukalov.bait_name))
    else:
        bait_match = None
    spatialvar_df.loc[i,'Bait_match_Stukalov'] = bait_match
    

spatialvar_df.to_csv('./Downloads/Sheet4_matchedGordon_matchedStukalov.csv', index=False)

# Our common hits = intersection of spatial and intensity
commonhits = pd.read_excel('./Downloads/Review_Analysis.xlsx', sheet_name=7, skiprows=1)
commonhits = commonhits[:50]

commonhits['Bait_match_Gordon'] = ''
commonhits['Bait_match_Stukalov'] = ''
for i, row in commonhits.iterrows():
    gene = row.Gene
    gene_gordon = gordon_ST2[gordon_ST2.PreyGene == gene]
    gene_stukalov = stukalov_STA[stukalov_STA.gene_name == gene]
    if len(gene_gordon) > 0:
        bait_match = ', '.join(list(gene_gordon.bait))
    else:
        bait_match = None
    commonhits.loc[i,'Bait_match_Gordon'] = bait_match
    
    if len(gene_stukalov) > 0:
        bait_match = ', '.join(list(gene_stukalov.bait_name))
    else:
        bait_match = None
    commonhits.loc[i,'Bait_match_Stukalov'] = bait_match
    
commonhits.to_csv('./Downloads/Sheet7_matchedGordon_matchedStukalov.csv', index=False)


mappings = pd.read_csv('X:/home/trangle/Desktop/Covid19project/Experiment_design/meta_ab.csv')
sheet5 = pd.read_excel('./Downloads/Review_Analysis.xlsx', sheet_name=5)
sheet5 = sheet5[0: 36]
sheet5 = pd.merge(sheet5, mappings[['well_id', 'Gene']], how='left', left_on='WellID', right_on='well_id')

sheet5['Bait_match_Gordon'] = ''
sheet5['Bait_match_Stukalov'] = ''
for i, row in sheet5.iterrows():
    gene = row.Gene.split(';')
    gene_gordon = pd.DataFrame()
    gene_stukalov = pd.DataFrame()
    for g in gene:
        gene_gordon = gene_gordon.append(gordon_ST2[gordon_ST2.PreyGene == g], ignore_index=True)
        gene_stukalov = gene_stukalov.append(stukalov_STA[stukalov_STA.gene_name == g], ignore_index=True)
        
    if len(gene_gordon) > 0:
        bait_match = ', '.join(list(gene_gordon.bait))
    else:
        bait_match = None
    sheet5.loc[i,'Bait_match_Gordon'] = bait_match
    
    if len(gene_stukalov) > 0:
        bait_match = ', '.join(list(gene_stukalov.bait_name))
    else:
        bait_match = None
    sheet5.loc[i,'Bait_match_Stukalov'] = bait_match

sheet5.to_csv('./Downloads/Sheet5_matchedGordon_matchedStukalov.csv', index=False)
  
## TODO: overlap stukalov and gordon

# Spatial hits 
hits_df = pd.read_excel('./Downloads/Review_Analysis.xlsx', sheet_name=9, skiprows=1, nrows=84)

# Intensity hits 
hits_df = pd.read_excel('./Downloads/Review_Analysis.xlsx', sheet_name=9, skiprows=88, nrows=178)

# Common hits
hits_df = pd.read_excel('./Downloads/Review_Analysis.xlsx', sheet_name=9, skiprows=269, nrows=50)

to_String = dict({"WASHC4":"KIAA1033",
                  "ELOB":"TCEB2",
                  "ELOC":"TCEB1",
                  "ATP5PD":"ATP5H",
                  "SPART": "SPG20"})



sf = pd.DataFrame(columns=['Source', 'Interaction', 'Target', 'Weight'])
for i, row in hits_df.iterrows():
    target = row.Gene
    target = target.split(',')[0]
    if target in to_String.keys():
        target = to_String[target]
    bait1 = str(row.Bait_match_Gordon).split(', ')
    bait2 = str(row.Bait_match_Stukalov).split(', ')
    baits = bait1 + bait2
    baits = [b.upper() for b in baits if b != 'nan']
    if len(baits) > 0:
        for bait in baits:
            interaction = {'Source': bait, 'Interaction': 'vp', 'Target': target, 'Weight': 1, 'Type': 'Common'}
            sf = sf.append(interaction, ignore_index=True)

sf = sf.drop_duplicates()
#sf.to_csv('X:/home/trangle/Desktop/Covid19project/viralhits_spatial.sf', index=False)
sf.to_csv('X:/home/trangle/Desktop/Covid19project/viralhits_intensity.sf', index=False)
sf.to_csv('X:/home/trangle/Desktop/Covid19project/viralhits_common.sf', index=False)


mito = ["PMPCB","PMPCA","PITRM1","GRPEL1","HSPA9","TIMM8B","GFER"]
TGN = ["GCC2","COG7","VTI1A","VAMP4","STX6"]
noteAttribute = pd.DataFrame()
noteAttribute['Node'] = list(set(sf.Source).union(set(sf.Target)))
noteAttribute['Type'] = ''
noteAttribute['Type'][noteAttribute.Node.isin(list(sf.Source))] = 'bait'
noteAttribute['Type'][noteAttribute.Node.isin(list(sf.Target))] = 'protein'
noteAttribute['HitType'] = ''
noteAttribute['HitType'][noteAttribute.Node.isin(list(sf[sf.Type=='Spatial'].Target))] = 'Spatial'
noteAttribute['HitType'][noteAttribute.Node.isin(list(sf[sf.Type=='Intensity'].Target))] = 'Intensity'
noteAttribute['HitType'][noteAttribute.Node.isin(list(sf[sf.Type=='Common'].Target))] = 'Common'

noteAttribute['PathwayEnriched'] = ''
noteAttribute['PathwayEnriched'][noteAttribute.Node.isin(mito)] = 'Mitochrondrial Protein Support'
noteAttribute['PathwayEnriched'][noteAttribute.Node.isin(TGN)] = 'Retrograde transport at the Trans-Golgi-Networ'

noteAttribute.to_csv('X:/home/trangle/Desktop/Covid19project/viralhits_noteAttribute.csv', index=False)

go = pd.read_csv('X:/home/trangle/Desktop/Covid19project/allproteins.bgo', skiprows=18,sep='\t')

#%% Updated hit list
import pandas as pd
df = pd.read_excel('./Downloads/Review_Analysis (1).xlsx', sheet_name=18, nrows=97)
#df = df[:96].append(df[df.Gene == 'AKR1B1'])
df["Spatial"] = [True if h.find('Spatial')!=-1 else False for h in df["Spatial/Intensity up/Intensity down"]]
df["Intensity_up"] = [True if h.find('up')!=-1 else False for h in df["Spatial/Intensity up/Intensity down"]]
df["Intensity_down"] = [True if h.find('down')!=-1 else False for h in df["Spatial/Intensity up/Intensity down"]]

# Hits from Gordon et al 2020
gordon_ST2 = pd.read_excel('./Downloads/41586_2020_2286_MOESM6_ESM.xlsx', skiprows=1)
gordon_ST2['bait'] = [b.split(' ')[1].upper() for b in gordon_ST2.Bait]

# Hits from Stukalov et al., 2020
stukalov_STA = pd.read_excel('./Downloads/media-2.xlsx', sheet_name=1, skiprows=1)
stukalov_STA = stukalov_STA[stukalov_STA.organism == 'SARS-CoV-2']
stukalov_STA.loc[stukalov_STA.gene_name == 43891, 'gene_name'] = 'MTARC1'
stukalov_STA.gene_name = [g.replace('...', '') for g in stukalov_STA.gene_name]

df['Bait_match_Gordon'] = ''
df['Bait_match_Stukalov'] = ''

to_String = dict({"WASHC4":"KIAA1033",
                  "ELOB":"TCEB2",
                  "ELOC":"TCEB1",
                  "ATP5PD":"ATP5H",
                  "SPART": "SPG20"})
    
for i, row in df.iterrows():
    genes = row.Gene.split(';')
    if genes[0] in to_String.keys():
        df.loc[i, 'Gene'] = to_String[genes[0]]
        genes += [to_String[genes[0]]]
        
    for gene in genes:
        gene_gordon = gordon_ST2[gordon_ST2.PreyGene == gene]
        gene_stukalov = stukalov_STA[stukalov_STA.gene_name == gene]
        bait_match = df.loc[i,'Bait_match_Gordon']
        if len(gene_gordon) > 0:
            bait_match = ', '.join(list(gene_gordon.bait) + [bait_match])
            if i in [25, 33, 57]:
                print(gene,bait_match)
        else:
            bait_match = bait_match
        df.loc[i,'Bait_match_Gordon'] = bait_match
        
        bait_match = df.loc[i,'Bait_match_Stukalov']
        if len(gene_stukalov) > 0:
            bait_match = ', '.join(list(gene_stukalov.bait_name) + [bait_match])
            if i in [25, 33, 57]:
                print(gene,bait_match)
        else:
            bait_match = bait_match
        df.loc[i,'Bait_match_Stukalov'] = bait_match

df.to_csv('X:/home/trangle/Desktop/Covid19project/viralhits_sheet17.csv', index=False)
        
#%%

sf = pd.DataFrame(columns=['Source', 'Interaction', 'Target', 'Weight'])

for i, row in df.iterrows():
    target = row.Gene
    target = target.split(',')[0]
    bait1 = str(row.Bait_match_Gordon).split(', ')
    bait2 = str(row.Bait_match_Stukalov).split(', ')
    baits = bait1 + bait2
    baits = set([b.upper() for b in baits if b not in ['None','nan', '']])
    if len(baits) > 0:
        for bait in baits:
            if row.Spatial and row.Intensity_up:
                interaction = {'Source': bait, 'Interaction': 'vp', 'Target': target, 'Weight': 1, 'Type': 'Intensity up;Spatial'}
                sf = sf.append(interaction, ignore_index=True)    
            elif row.Spatial:
                interaction = {'Source': bait, 'Interaction': 'vp', 'Target': target, 'Weight': 1, 'Type': 'Spatial'}
                sf = sf.append(interaction, ignore_index=True)
            elif row.Intensity_up:
                interaction = {'Source': bait, 'Interaction': 'vp', 'Target': target, 'Weight': 1, 'Type': 'Intensity up'}
                sf = sf.append(interaction, ignore_index=True)
            elif row.Intensity_down:
                interaction = {'Source': bait, 'Interaction': 'vp', 'Target': target, 'Weight': 1, 'Type': 'Intensity down'}
                sf = sf.append(interaction, ignore_index=True)
sf.to_csv('X:/home/trangle/Desktop/Covid19project/viralhits_sheet17.sf', index=False)



noteAttribute = pd.DataFrame()
noteAttribute['Node'] = list(set(sf.Source).union(set(df.Gene)))
noteAttribute['Node'] = [to_String[n] if n in to_String.keys() else n for n in noteAttribute.Node]
noteAttribute['Type'] = ''
noteAttribute['Type'][noteAttribute.Node.isin(list(sf.Source))] = 'bait'
noteAttribute['Type'][~noteAttribute.Node.isin(list(sf.Source))] = 'protein'
noteAttribute['HitType'] = ''
noteAttribute['HitType'][noteAttribute.Node.isin(list(sf.Source))] = 'bait'
noteAttribute['HitType'][noteAttribute.Node.isin(list(df.Gene[df.Spatial]))] = 'Spatial'
noteAttribute['HitType'][noteAttribute.Node.isin(list(df.Gene[df.Intensity_up]))] = 'Intensity up'
noteAttribute['HitType'][noteAttribute.Node.isin(list(df.Gene[df.Intensity_down]))] = 'Intensity down'
noteAttribute['HitType'][noteAttribute.Node.isin(list(set(df.Gene[df.Intensity_up]) & set(df.Gene[df.Spatial])))] = 'Intensity up;Spatial'
#noteAttribute['HitType'][~noteAttribute.Node.isin(list(set(sf.Source).union(set(sf.Target))))] = 'No baits'

noteAttribute.to_csv('X:/home/trangle/Desktop/Covid19project/viralhits_noteAttribute.csv', index=False)
