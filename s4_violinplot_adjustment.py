import os
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import requests
import scipy.stats
import seaborn as sns



SAMPLES_DEST = "/data/trang/covid19_data_CZ8746_annotation/segmentation"
url4labels = "https://raw.githubusercontent.com/haoxusci/imjoy-plugin-config/master/config/Covid19ImageAnnotator.imjoy.config.json"
NONINFECTED = ["A","B","C","D"]
INFECTED = ["E","F","G","H"]

def config_to_labels(imjoy_url):
    labels_covid_json = requests.get(imjoy_url)
    labels_dict = labels_covid_json.json()["labels"]["numerized_labels"]
    labels_dict = {str(v): k for k, v in labels_dict.items()}
    return labels_dict


def get_cells_filtered(well_cells):
    # From https://github.com/CellProfiling/Annotation_Cell_Segmentation/blob/covid-application/utils/violin_plot_adjustment-tmp_jay.ipynb
    well_NI_cells = well_cells[well_cells['Infected']==0]
    well_I_cells = well_cells[well_cells['Infected']==1]
    well_cells.loc[well_cells['Infected']==0, 'Infected'] = 'NI'
    well_cells.loc[well_cells['Infected']==1, 'Infected'] = 'I'
    print(well_cells)
    """# integ intensity
    well_cells_nuclei_integ = well_cells[['image_id', 'protein-nuclei-integration', 'Infected']].rename(
        columns={"protein-nuclei-integration": "Intensity"}).assign(**dict.fromkeys(['value_type'], 'integ(nuclei)'))
    tStat_nuclei_integ, pValue_nuclei_integ = scipy.stats.ttest_ind(
        well_cells_nuclei_integ[well_cells_nuclei_integ['Infected']=='NI']['Intensity'],
        well_cells_nuclei_integ[well_cells_nuclei_integ['Infected']=='I']['Intensity']
    )

    well_cells_cytosol_integ = well_cells[['image_id', 'protein-cytosol-integration', 'Infected']].rename(
        columns={"protein-cytosol-integration": "Intensity"}).assign(**dict.fromkeys(['value_type'], 'integ(cytosol)'))
    tStat_cytosol_integ, pValue_cytosol_integ = scipy.stats.ttest_ind(
        well_cells_cytosol_integ[well_cells_cytosol_integ['Infected']=='NI']['Intensity'],
        well_cells_cytosol_integ[well_cells_cytosol_integ['Infected']=='I']['Intensity']
    )

    well_cells_cell_integ = well_cells[['image_id', 'protein-cell-integration', 'Infected']].rename(
        columns={"protein-cell-integration": "Intensity"}).assign(**dict.fromkeys(['value_type'], 'integ(cell)'))
    tStat_cell_integ, pValue_cell_integ = scipy.stats.ttest_ind(
        well_cells_cell_integ[well_cells_cell_integ['Infected']=='NI']['Intensity'],
        well_cells_cell_integ[well_cells_cell_integ['Infected']=='I']['Intensity']
    )"""
    
    # mean intensity
    well_cells_nuclei_integ = well_cells[['image_id', 'protein-nuclei-mean', 'Infected']].rename(
        columns={"protein-nuclei-mean": "Intensity"}).assign(**dict.fromkeys(['value_type'], 'mean(nuclei)'))
    tStat_nuclei_integ, pValue_nuclei_integ = scipy.stats.ttest_ind(
        well_cells_nuclei_integ[well_cells_nuclei_integ['Infected']=='NI']['Intensity'],
        well_cells_nuclei_integ[well_cells_nuclei_integ['Infected']=='I']['Intensity']
    )

    well_cells_cytosol_integ = well_cells[['image_id', 'protein-cytosol-mean', 'Infected']].rename(
        columns={"protein-cytosol-mean": "Intensity"}).assign(**dict.fromkeys(['value_type'], 'mean(cytosol)'))
    tStat_cytosol_integ, pValue_cytosol_integ = scipy.stats.ttest_ind(
        well_cells_cytosol_integ[well_cells_cytosol_integ['Infected']=='NI']['Intensity'],
        well_cells_cytosol_integ[well_cells_cytosol_integ['Infected']=='I']['Intensity']
    )

    well_cells_cell_integ = well_cells[['image_id', 'protein-cell-mean', 'Infected']].rename(
        columns={"protein-cell-mean": "Intensity"}).assign(**dict.fromkeys(['value_type'], 'mean(cell)'))
    tStat_cell_integ, pValue_cell_integ = scipy.stats.ttest_ind(
        well_cells_cell_integ[well_cells_cell_integ['Infected']=='NI']['Intensity'],
        well_cells_cell_integ[well_cells_cell_integ['Infected']=='I']['Intensity']
    )
    
    """well_cells_cell_mean = well_cells[['image_id', 'protein-cell-mean', 'Infected']].rename(
        columns={"protein-cell-mean": "Intensity"}).assign(**dict.fromkeys(['value_type'], 'mean(cell)'))
    well_cells_cell_mean.Intensity = well_cells_cell_mean.Intensity * (well_cells['protein-cell-integration'].max()/well_cells_cell_mean.Intensity.max())
    tStat_cell_mean, pValue_cell_mean = scipy.stats.ttest_ind(
        well_cells_cell_mean[well_cells_cell_mean['Infected']=='NI']['Intensity'],
        well_cells_cell_mean[well_cells_cell_mean['Infected']=='I']['Intensity']
    )"""

    """return {
        'well_id': well_id,
        'pValue_nuclei_integ': pValue_nuclei_integ,
        'pValue_cytosol_integ': pValue_cytosol_integ,
        'pValue_cell_integ': pValue_cell_integ,
        'pValue_cell_mean': pValue_cell_mean
    }"""
    well_cells_filterted = pd.concat(
        [
            well_cells_nuclei_integ,
            well_cells_cytosol_integ,
            well_cells_cell_integ,
            #well_cells_cell_mean
        ]    
    )

    return well_cells_filterted, pValue_nuclei_integ, pValue_cytosol_integ, pValue_cell_integ #, pValue_cell_mean

def plot_violin(well_cells_filterted, pValue_nuclei_integ, pValue_cytosol_integ, pValue_cell_integ, save_dir, well_id):
    props = dict(boxstyle='round', facecolor='#c4d48c', alpha=.2)
    fig = plt.figure(constrained_layout=True, figsize=(8, 8))
    gs = fig.add_gridspec(1, 1)
    ax = fig.add_subplot(gs[0, 0])
    #well_cells_filterted, pValue_nuclei_integ, pValue_cytosol_integ, pValue_cell_integ = get_cells_filtered(df_cells, image_ids, well_id)
    #print(df_cells)
    sns.violinplot(
        x="value_type",
        y="Intensity",
        hue="Infected",
        data=well_cells_filterted,
        palette="Set2",
        split=True,
        #order=['integ(nuclei)', 'integ(cytosol)'],#, 'integ(cell)'],
        order=['mean(nuclei)', 'mean(cytosol)', 'mean(cell)'],
        hue_order=['NI', 'I'],
        orient='v',
        scale='count',
        inner='quartile',
        scale_hue=False,
        bw=.2
    )
    ax.tick_params(
        axis='both',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        left=True,
        labelleft=True,
        bottom=True,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=True)
    ax.tick_params(axis='both', which='major', labelsize=22)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.legend().set_visible(False)
    ax.set_xlabel("Value Type", fontsize=30)
    ax.set_ylabel("Intensity", fontsize=30)

    ax.text(
        0.3,
        0.6,
        "P value : " + str(round(pValue_nuclei_integ, 6)),
        transform=ax.transAxes,
        fontsize=15,
        horizontalalignment='center',
        verticalalignment='top',
        bbox=props
    )
    ax.text(
        0.68,
        0.75,
        "P value : " + str(round(pValue_cytosol_integ, 6)),
        transform=ax.transAxes,
        fontsize=15,
        horizontalalignment='center',
        verticalalignment='top',
        bbox=props
    )
    ax.text(
        0.63,
        0.7,
        "P value : " + str(round(pValue_cell_integ, 6)),
        transform=ax.transAxes,
        fontsize=15,
        horizontalalignment='center',
        verticalalignment='top',
        bbox=props
    )
    """
    ax.text(
        0.85,
        0.78,
        "P value : " + str(round(pValue_cell_mean, 6)),
        transform=ax.transAxes,
        fontsize=15,
        horizontalalignment='center',
        verticalalignment='top',
        bbox=props
    )"""
    ax.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
    ymin, ymax = ax.get_ylim()
    print(ymin, ymax)
    ymin = min(ymin, 0)
    plt.ylim(ymin, )
    plt.savefig(os.path.join(save_dir, well_id + '_mean.png'), dpi=100)
    plt.savefig(os.path.join(save_dir, well_id + '_mean.svg'), dpi=100)

if __name__ == "__main__":
    df_cells = pd.read_csv(os.path.join(SAMPLES_DEST, "cells_combined_nobigcell.csv")) 
    meta = pd.read_csv('/home/trangle/Desktop/Covid19project/Experiment_design/meta_ab_CZ8746.csv')
    violin_per_well = False
    violin_per_ab = True

    labels_dict = config_to_labels(url4labels)
    
    if violin_per_well:
        image_ids = list(set(df_cells['image_id']))
        well_ids = list(set(map(lambda item: '_'.join(item.split(os.sep)[-1].split('_')[:-1]), image_ids)))
        well_ids.sort()
        save_dir = os.path.join(os.path.dirname(SAMPLES_DEST), "violin")
        for well_id in well_ids:    
            well_image_ids = list(filter(lambda item: item if well_id == item[:-2] else None, image_ids))
            well_cells = df_cells[df_cells['image_id'].isin(well_image_ids)]
            well_cells_filterted, pValue_nuclei_integ, pValue_cytosol_integ, pValue_cell_integ = get_cells_filtered(well_cells)
            plot_violin(well_cells_filterted, pValue_nuclei_integ, pValue_cytosol_integ, pValue_cell_integ, save_dir, well_id)
    
    if violin_per_ab:
        df_cells['well_id'] = [image_id.split('/')[-1][:-2] for image_id in df_cells.image_id]
        df_cells['Ab'] = [meta[meta.well_id==well].Antibody.values[0] for well in df_cells.well_id]
        antibodies = list(set(meta.Antibody))
        save_dir = os.path.join(os.path.dirname(SAMPLES_DEST), "violin_ab")
        for ab in antibodies:
            well_cells = df_cells[df_cells.Ab == ab]
            well_cells_filterted, pValue_nuclei_integ, pValue_cytosol_integ, pValue_cell_integ = get_cells_filtered(well_cells)
            plot_violin(well_cells_filterted, pValue_nuclei_integ, pValue_cytosol_integ, pValue_cell_integ, save_dir, ab)
