import pandas as pd 
import json

label_idx = json.load("/data/hao/Covid-well-annotation/all_locations.json")
def idx_to_labels(idx_list):
    labels = []
    for idx in idx_list:
        if idx in ["!undefined", "!null"]:
            print("Not here idx is not a good value: ", idx)
            idx = 0
        for k, v in label_idx.items():
            try:
                if v == int(float(idx)):
                    labels += [k]
            except:
                labels += [idx]
    return list(set(labels))
meta = pd.read_csv('/home/trangle/Desktop/Covid19project/Experiment_design/meta_ab_CZ8780.csv')
meta['plate'] = meta.plate.str.replace('CZ8780_plate_II','22')
meta['plate'] = meta.plate.str.replace('CZ8780_plate_I','21')
meta['task_id'] = meta.plate + '_' + meta.well_id
task_in_progress = pd.read_csv('tasks-in-progress.csv')

task_in_progress['Infected_location'] = [idx_to_labels(eval(l)['Infected'].split('|')) for l in task_in_progress.labels]
task_in_progress['Not_Infected_location'] = [idx_to_labels(eval(l)['Non_infected'].split('|')) for l in task_in_progress.labels]

df = task_in_progress.merge(meta, on ='task_id', how='left')