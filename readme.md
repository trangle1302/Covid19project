Subcellualar spatialproteomic screen of SARS-CoV-2 infection
===============================

This repository contains the code used the analysis of spatial proteomics screen in the manuscript *"Subcellular mapping of the protein landscape of SARS-CoV-2 infected cells for target-centric drug repurposing"*.

## Steps to preprocess plates for annotation-tool:
Step 1. Convert file names and max-projection

```sh
python generate_covid_dataset-max-projection.py
```
Paths need to be updated: 
- `ACQUIRED_DATA_PATH`: Raw data acquisition path with *Images* folder (input folder)
- `FORMATTED_DATA_PATH`: Images with formatted names (intermediate folder)
- `IMG_FOLDER`: Images with max-projection (output folder)
- `ANNOTATION_FOLDER`: Folder for `annotation-tool` with Imjoy (output folder, only used when manual annotation is later needed)

Step 2. Segmentation with HPACellSegmentator and convert mask to geojson file
```sh
python generate_covid_dataset-max-projection_segmentation.py
python mask_to_json.py
```
Paths need to be updated: 
- `IMG_FOLDER`: Images with max-projection 
- `SEGMENTATION_FOLDER`: path to save segmentation masks
- `NUC_MODEL`: nuclei model path
- `CELL_MODEL`: cell model path

If signals are too low in a plate, they need to be magnified for the segmentation algirthm to work. Note: The increased intensity is only used to create masks, subsequent analysis rely on raw intensity.
```sh
python generate_covid_dataset-max-projection.py
python generate_covid_dataset-max-projection_segmentation.py
```

If the experiment needs to be uploaded to annotation tool later, change `SEGMENTATION_FOLDER` and run
```sh
python mask_to_json.py
```
to create `annotation.json` files for each FOV.

## Steps to analyse intensity hits:
Step 3. Covid and protein quantification, Violin plot
```sh
python s1_covid_quantify_well.py
python s2_combine_quantification_df.py
python s3_define_cellsize_and_virus_threshold.py
python s4_violinplot_adjustment.py
```
Paths need to be updated: 
- `IMG_FOLDER`: Images with max-projection
- `SAMPLE_DEST`: sample destination
- (s4_)`NONINFECTED`: Wells or well pattern with non-infected cells
- (s4_)`INFECTED`: Wells or well pattern with non-infected cells

Step 4. Fold change, volcano and circle plot
```sh
python covid19plots.py
```
