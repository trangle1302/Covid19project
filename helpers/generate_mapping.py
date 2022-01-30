import json
import pandas as pd
from pandas.core.indexes import base

base_dir = "/home/trangle/Desktop/Covid19project"
row_mapping = {
    "r01": "A",
    "r02": "B",
    "r03": "C",
    "r04": "D",
    "r05": "E",
    "r06": "F",
    "r07": "G",
    "r08": "H",
}

if __name__ == "__main__":
    meta_ab = pd.read_csv(f"{base_dir}/Experiment_design/meta_ab.csv")
    meta_ab["platexx"] = [f.split('_')[0] for f in meta_ab.well_id]
    print(meta_ab[meta_ab.platexx == 8])
    for plate in range(1,9):
        if plate == 8:
            continue
        target = json.load(open(f"{base_dir}/covid_data_Scilifelab_portal/plate{plate}.ome.tif_filenames.json"))
        target_df = pd.DataFrame(target)
        target_df["well_id"] = ""
        target_df["sample_id"] = ""
        for i,row in target_df.iterrows():
            r = row_mapping[row.files[:3]]
            c = str(int(row.files[4:6]))
            fov = str(int(row.files[7:9]))
            plane = str(int(row.files[10:12]))
            row.well_id = "".join([str(plate),"_",r,c])
            row.sample_id = "".join([fov,"_",plane])
            #print(r,c,fov,plane)

        df = target_df.merge(meta_ab[["Antibody","Gene","well_id"]], on="well_id", how="left")
        df["new_names"] = df.Antibody + "_" + df.Gene + "_" + df.sample_id
        print(plate, df.new_names.value_counts())
        id_nan = df.new_names.isna()
        df["new_names"][id_nan] = df.well_id[id_nan] + '_' + df.sample_id
        target["files"] = list(df.new_names.values)
        with open(f"{base_dir}/covid_data_Scilifelab_portal/plate_{plate}.ome.tif_filenames.json", "w") as outfile:
            json.dump(target, outfile)
